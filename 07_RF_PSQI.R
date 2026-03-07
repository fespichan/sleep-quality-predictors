############################################################
# 07_RF_PSQI.R
# Random Forest: SMOTE-balanced training, pooled metrics,
# variable importance with LR p-values, and PDP
#
# Project: Biological and Anthropometric Determinants of
#          Sleep Quality Across the Adult Life Course
# Journal: Sleep Medicine
#
# Requires: data_stacked (from 02_Multiple_Imputation.R)
#           final_table_logit.rds (from 06_LR_Pooled_PSQI.R)
#
# Output:
#   Tables:
#     Table_RF_Performance.xlsx
#     Table_RF_Variable_Importance.xlsx
#     Table_RF_HL.xlsx
#     Table_RF_CV_SMOTE.xlsx
#   Figures:
#     Fig2_RF_Combined_4Panel.tif  (Importance + PDP + ROC + CM)
#   Data:
#     auc_values_rf.rds
############################################################

library(dplyr)
library(ggplot2)
library(caret)
library(randomForest)
library(pROC)
library(pdp)
library(themis)
library(recipes)
library(ResourceSelection)
library(patchwork)
library(writexl)

n_imp <- max(data_stacked$.imp)

# =================================================================
# STEP 1 — PREPARE DATA
# =================================================================

# Rename psqi levels for caret compatibility
levels(data_stacked$psqi) <- c("Good_Sleep", "Poor_Sleep")

# Consistent train/test partition (same as LR script)
set.seed(123456)
first_imp <- data_stacked %>% filter(.imp == 1)
indices_train <- createDataPartition(first_imp$psqi, p = 0.80, list = FALSE)

# =================================================================
# STEP 2 — CONTAINERS
# =================================================================

model_list_rf    <- list()
test_list_rf     <- list()
importance_list  <- list()
confusion_list   <- list()
roc_df_all       <- data.frame()

auc_values_rf    <- numeric(n_imp)
sensitivity_vals <- numeric(n_imp)
specificity_vals <- numeric(n_imp)
accuracy_vals    <- numeric(n_imp)
kappa_vals       <- numeric(n_imp)
ppv_vals         <- numeric(n_imp)
npv_vals         <- numeric(n_imp)
bal_acc_vals     <- numeric(n_imp)
hl_pvalues_rf    <- numeric(n_imp)

# =================================================================
# STEP 3 — MAIN LOOP: SMOTE + RF + EVALUATION
# =================================================================

control_rf <- trainControl(
  method = "cv", number = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

for (i in 1:n_imp) {

  if (i %% 20 == 0) message(sprintf("Processing imputation %d of %d", i, n_imp))

  capa_i <- data_stacked %>% filter(.imp == i)
  train_i <- capa_i[indices_train, ] %>%
    dplyr::select(-.imp, -.id) %>% as.data.frame()
  test_i <- capa_i[-indices_train, ] %>%
    dplyr::select(-.imp, -.id) %>% as.data.frame()

  if (sum(is.na(train_i)) > 0) {
    warning(sprintf("Imputation %d has NAs — skipping", i)); next
  }

  # SMOTE-NC on training set only
  tryCatch({
    set.seed(123)
    smote_prep <- prep(
      recipe(psqi ~ ., data = train_i) %>%
        step_smotenc(psqi, over_ratio = 1, seed = 123, neighbors = 5),
      training = train_i, verbose = FALSE
    )
    train_i_balanced <- bake(smote_prep, new_data = NULL)

    if (i == 1) {
      message("Class distribution BEFORE SMOTE (imp 1):")
      print(table(train_i$psqi))
      message("Class distribution AFTER SMOTE (imp 1):")
      print(table(train_i_balanced$psqi))
    }
  }, error = function(e) {
    warning(sprintf("SMOTE failed for imputation %d — using unbalanced", i))
    train_i_balanced <<- train_i
  })

  # Train Random Forest
  set.seed(123)
  rf_model_i <- train(
    psqi ~ sex + age + bmic + bf + hb,
    data = train_i_balanced, method = "rf",
    trControl = control_rf, metric = "ROC",
    importance = TRUE, ntree = 500
  )

  # Predict on unbalanced test set
  pred_class_i <- predict(rf_model_i, newdata = test_i)
  pred_prob_i  <- predict(rf_model_i, newdata = test_i, type = "prob")

  # ROC and AUC
  roc_obj_i <- roc(test_i$psqi, pred_prob_i$Poor_Sleep, quiet = TRUE)
  auc_values_rf[i] <- as.numeric(roc_obj_i$auc)

  coords_opt <- coords(roc_obj_i, "best", best.method = "youden",
                        ret = c("sensitivity", "specificity"))
  sensitivity_vals[i] <- coords_opt$sensitivity
  specificity_vals[i] <- coords_opt$specificity

  # Confusion matrix metrics
  cm_i <- confusionMatrix(pred_class_i, test_i$psqi, positive = "Poor_Sleep")
  accuracy_vals[i] <- cm_i$overall["Accuracy"]
  kappa_vals[i]    <- cm_i$overall["Kappa"]
  ppv_vals[i]      <- cm_i$byClass["Pos Pred Value"]
  npv_vals[i]      <- cm_i$byClass["Neg Pred Value"]
  bal_acc_vals[i]  <- cm_i$byClass["Balanced Accuracy"]
  confusion_list[[i]] <- cm_i$table

  # Hosmer-Lemeshow
  outcome_i <- ifelse(test_i$psqi == "Poor_Sleep", 1, 0)
  tryCatch({
    hl_i <- hoslem.test(outcome_i, pred_prob_i$Poor_Sleep, g = 10)
    hl_pvalues_rf[i] <- hl_i$p.value
  }, error = function(e) { hl_pvalues_rf[i] <<- NA })

  # ROC curve data
  roc_df_all <- rbind(roc_df_all, data.frame(
    specificity = roc_obj_i$specificities,
    sensitivity = roc_obj_i$sensitivities,
    imputation  = i
  ))

  # Variable importance
  imp_i <- varImp(rf_model_i, scale = FALSE)$importance
  importance_list[[i]] <- data.frame(
    Variable   = rownames(imp_i),
    Importance = if ("Overall" %in% colnames(imp_i)) imp_i$Overall
                 else rowMeans(imp_i),
    stringsAsFactors = FALSE
  )

  model_list_rf[[i]] <- rf_model_i
  test_list_rf[[i]]  <- test_i
}

saveRDS(auc_values_rf, "auc_values_rf.rds")

# =================================================================
# STEP 4 — POOLED RESULTS
# =================================================================

# Performance metrics
performance_rf <- data.frame(
  Metric = c("AUC", "Sensitivity", "Specificity", "Accuracy",
             "Kappa", "PPV", "NPV", "Balanced Accuracy"),
  Mean = round(c(mean(auc_values_rf), mean(sensitivity_vals),
                 mean(specificity_vals), mean(accuracy_vals),
                 mean(kappa_vals), mean(ppv_vals, na.rm = TRUE),
                 mean(npv_vals, na.rm = TRUE), mean(bal_acc_vals)), 3),
  SD = round(c(sd(auc_values_rf), sd(sensitivity_vals),
               sd(specificity_vals), sd(accuracy_vals),
               sd(kappa_vals), sd(ppv_vals, na.rm = TRUE),
               sd(npv_vals, na.rm = TRUE), sd(bal_acc_vals)), 3)
)

message("\n=== RANDOM FOREST PERFORMANCE ===")
print(performance_rf, row.names = FALSE)
write_xlsx(performance_rf, "Table_RF_Performance.xlsx")

# Variable importance
importance_all <- do.call(rbind, importance_list)
importance_avg <- importance_all %>%
  group_by(Variable) %>%
  summarise(Mean_Importance = mean(Importance),
            SD_Importance = sd(Importance), .groups = "drop") %>%
  arrange(desc(Mean_Importance))

message("\n=== VARIABLE IMPORTANCE (Mean Decrease Gini) ===")
print(as.data.frame(importance_avg), row.names = FALSE)
write_xlsx(as.data.frame(importance_avg), "Table_RF_Variable_Importance.xlsx")

# Confusion matrix
confusion_avg_rf <- Reduce("+", confusion_list) / n_imp
message("\nAverage confusion matrix:")
print(round(confusion_avg_rf, 2))

# Hosmer-Lemeshow
hl_summary_rf <- data.frame(
  Statistic = c("Mean p-value", "SD p-value",
                "% adequate (p > 0.05)"),
  Value = c(round(mean(hl_pvalues_rf, na.rm = TRUE), 3),
            round(sd(hl_pvalues_rf, na.rm = TRUE), 3),
            round(100 * mean(hl_pvalues_rf > 0.05, na.rm = TRUE), 1))
)

message("\n=== H-L CALIBRATION ===")
print(hl_summary_rf, row.names = FALSE)
write_xlsx(hl_summary_rf, "Table_RF_HL.xlsx")

# =================================================================
# STEP 5 — PARTIAL DEPENDENCE PLOT (HEMOGLOBIN)
# =================================================================

message("\nComputing PDP for hemoglobin...")

compute_avg_pdp <- function(var_name, model_list, n_models = n_imp) {
  pdp_list <- list()
  for (j in 1:n_models) {
    pd_j <- pdp::partial(model_list[[j]], pred.var = var_name,
                         prob = TRUE, which.class = "Poor_Sleep")
    pdp_list[[j]] <- pd_j
  }
  pdp_all <- do.call(rbind, pdp_list)
  pdp_all %>%
    group_by(across(all_of(var_name))) %>%
    summarise(yhat_mean = mean(yhat), yhat_sd = sd(yhat),
              yhat_lo = quantile(yhat, 0.025),
              yhat_hi = quantile(yhat, 0.975), .groups = "drop")
}

pdp_hb <- compute_avg_pdp("hb", model_list_rf)

# =================================================================
# STEP 6 — FIGURE 2: COMBINED 4-PANEL
# =================================================================

# Load pooled LR p-values for importance plot annotations
final_table <- readRDS("final_table_logit.rds") %>%
  mutate(term = case_when(
    term == "bmicOverweight" ~ "Overweight",
    term == "bmicObesity"    ~ "Obesity",
    TRUE ~ term
  ))

logit_pvalues <- final_table %>%
  dplyr::select(term, p.value) %>%
  mutate(Logit_P = case_when(
    p.value < 0.001 ~ "p < 0.001",
    p.value < 0.01  ~ "p < 0.01",
    p.value < 0.05  ~ "p < 0.05",
    TRUE ~ "ns"
  )) %>%
  rename(Variable = term)

importance_with_p <- importance_avg %>%
  left_join(logit_pvalues, by = "Variable") %>%
  mutate(
    Logit_P = ifelse(is.na(Logit_P), "", Logit_P),
    Significance = case_when(
      Logit_P %in% c("p < 0.001") ~ "p < 0.001",
      Logit_P %in% c("p < 0.01")  ~ "p < 0.01",
      Logit_P %in% c("p < 0.05")  ~ "p < 0.05",
      TRUE ~ "ns"
    )
  )

sig_colors <- c("p < 0.001" = "#EEC900", "p < 0.01" = "#FF8C00",
                "p < 0.05" = "#2196F3", "ns" = "#95a5a6")

# Panel A: Variable Importance

importance_with_p <-importance_with_p %>%
  mutate(Variable = case_when(
    Variable == "bmicOverweight" ~ "Overweight",
    Variable == "bmicObesity"    ~ "Obesity",
    TRUE ~ Variable
  ))

panel_a <- ggplot(importance_with_p,
                  aes(x = Mean_Importance,
                      y = reorder(Variable, Mean_Importance),
                      fill = Significance)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = Mean_Importance - SD_Importance,
                     xmax = Mean_Importance + SD_Importance),
                 height = 0.3, linewidth = 0.4) +
  geom_text(aes(x = Mean_Importance + SD_Importance + 1.5,
                label = Logit_P),
            hjust = 0, size = 5, fontface = "italic") +
  scale_fill_manual(values = sig_colors,
                    name = "Logistic Regression\n(Pooled p-value)") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
  theme_classic(base_size = 16) +
  labs(title = "A. Variable Importance",
       subtitle = "RF: Mean Decrease Gini (± SD) | Labels: Pooled Logistic Regression p-values",
       x = "Mean Decrease Gini (± SD)", y = "") +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey40"),
    axis.text = element_text(color = "black", size = 14),
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14)
  )

# Panel B: PDP Hemoglobin
panel_b <- ggplot(pdp_hb, aes(x = hb, y = yhat_mean)) +
  geom_ribbon(aes(ymin = yhat_lo, ymax = yhat_hi),
              fill = "#2196F3", alpha = 0.2) +
  geom_line(color = "#2196F3", linewidth = 1.5) +
  theme_classic(base_size = 16) +
  labs(title = "B. Partial Dependence: Hemoglobin",
       x = "Hemoglobin (g/dL)", y = "Prob. Poor Sleep") +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14, color = "black")
  )

# Panel C: ROC Curve
panel_c <- ggplot(roc_df_all, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(aes(group = imputation), alpha = 0.8, color = "#00C5CD",
            linewidth = 0.4) +
  stat_summary(fun = mean, geom = "line", color = "#E74C3C",
               linewidth = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "#A020F0", linewidth = 1.5) +
  annotate("text", x = 0.65, y = 0.25,
           label = sprintf("Mean AUC = %.3f", mean(auc_values_rf)),
           size = 5.5, fontface = "bold") +
  theme_classic(base_size = 16) +
  labs(title = "C. ROC Curve", x = "1 - Specificity",
       y = "Sensitivity") +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14, color = "black")
  ) +
  coord_equal()

# Panel D: Confusion Matrix
conf_df_rf <- as.data.frame(confusion_avg_rf)
colnames(conf_df_rf) <- c("Reference", "Prediction", "Freq")
total_rf <- sum(conf_df_rf$Freq)
conf_df_rf$Label <- sprintf("n = %.1f\n(%.1f%%)",
                            conf_df_rf$Freq,
                            conf_df_rf$Freq / total_rf * 100)

panel_d <- ggplot(conf_df_rf,
                  aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(color = "black", linewidth = 1) +
  geom_text(aes(label = Label), color = "white", size = 5.5,
            fontface = "bold", lineheight = 0.85) +
  scale_fill_gradient(low = "#3498DB", high = "#E74C3C") +
  labs(title = "D. Confusion Matrix", x = "Actual", y = "Predicted") +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14, face = "bold"),
    legend.position = "none", panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Combine and save
fig_combined_rf <- (panel_a + panel_b) / (panel_c + panel_d)

ggsave("Fig2_RF_Combined_4Panel.tif", fig_combined_rf,
       device = "tiff", width = 26, height = 24, units = "cm",
       dpi = 1200, compression = "lzw", bg = "white")

# =================================================================
# STEP 7 — 10-FOLD CROSS-VALIDATION
# =================================================================

first_imp_rf <- data_stacked %>%
  filter(.imp == 1) %>%
  dplyr::select(-.imp, -.id) %>%
  mutate(across(where(is.ordered),
                ~factor(as.character(.x), ordered = FALSE))) %>%
  as.data.frame()

set.seed(123456)
cv_control_rf <- trainControl(
  method = "cv", number = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  sampling = "smote"
)

cv_model_rf <- train(
  psqi ~ sex + age + bmic + bf + hb,
  data = first_imp_rf, method = "rf",
  trControl = cv_control_rf, metric = "ROC",
  importance = TRUE, ntree = 500
)

best_idx <- which.max(cv_model_rf$results$ROC)
cv_summary_rf <- data.frame(
  Metric = c("CV_ROC", "CV_Sensitivity", "CV_Specificity", "Best_mtry"),
  Value = c(round(cv_model_rf$results$ROC[best_idx], 3),
            round(cv_model_rf$results$Sens[best_idx], 3),
            round(cv_model_rf$results$Spec[best_idx], 3),
            cv_model_rf$results$mtry[best_idx])
)

message("\n=== 10-FOLD CV WITH SMOTE ===")
print(cv_summary_rf, row.names = FALSE)
write_xlsx(cv_summary_rf, "Table_RF_CV_SMOTE.xlsx")

# =================================================================
# SUMMARY
# =================================================================

message("\n========================================")
message("RANDOM FOREST ANALYSIS COMPLETED")
message("========================================")
message(sprintf("Mean AUC: %.3f (SD = %.3f)",
                mean(auc_values_rf), sd(auc_values_rf)))
message(sprintf("H-L: mean p = %.3f (%s%% adequate)",
                mean(hl_pvalues_rf, na.rm = TRUE),
                round(100 * mean(hl_pvalues_rf > 0.05, na.rm = TRUE), 1)))
message("========================================")


