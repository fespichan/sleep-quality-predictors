############################################################
# 08_XGBoost_PSQI.R
# XGBoost: SMOTE within CV folds, hyperparameter tuning,
# pooled metrics across 100 imputed datasets
#
# Project: Biological and Anthropometric Determinants of
#          Sleep Quality Across the Adult Life Course
# Journal: Journal of Sleep Research
#
# Requires: data_stacked (from 02_Multiple_Imputation.R)
#           psqi levels: Good_Sleep / Poor_Sleep (from 07)
#
# Output:
#   Tables:
#     Table_XGB_Performance.xlsx
#     Table_XGB_Variable_Importance.xlsx
#     Table_XGB_HL.xlsx
#     Table_XGB_CV_Results.xlsx
#   Figures:
#     Fig_XGB_Combined_4Panel.tif
#   Data:
#     auc_values_xgb.rds
############################################################

library(dplyr)
library(ggplot2)
library(caret)
library(xgboost)
library(pROC)
library(themis)
library(recipes)
library(ResourceSelection)
library(patchwork)
library(writexl)

n_imp <- max(data_stacked$.imp)

# =================================================================
# STEP 1 — PREPARE DATA
# =================================================================

# Ensure psqi levels are caret-compatible
if (!all(levels(data_stacked$psqi) %in% c("Good_Sleep", "Poor_Sleep"))) {
  data_stacked$psqi <- factor(data_stacked$psqi,
                               levels = c("g.sleep", "p.sleep"))
  levels(data_stacked$psqi) <- c("Good_Sleep", "Poor_Sleep")
}

# Consistent train/test partition
set.seed(123456)
first_imp <- data_stacked %>% filter(.imp == 1)
indices_train <- createDataPartition(first_imp$psqi, p = 0.80, list = FALSE)

# =================================================================
# STEP 2 — CONTAINERS
# =================================================================

model_list_xgb  <- list()
test_list_xgb   <- list()
importance_list  <- list()
confusion_list   <- list()
roc_df_all       <- data.frame()

auc_values_xgb   <- numeric(n_imp)
sensitivity_vals  <- numeric(n_imp)
specificity_vals  <- numeric(n_imp)
accuracy_vals     <- numeric(n_imp)
kappa_vals        <- numeric(n_imp)
ppv_vals          <- numeric(n_imp)
npv_vals          <- numeric(n_imp)
bal_acc_vals      <- numeric(n_imp)
hl_pvalues_xgb   <- numeric(n_imp)

# =================================================================
# STEP 3 — MAIN LOOP: XGBoost + CV + SMOTE
# =================================================================

# SMOTE applied within each CV fold (not pre-applied)
control_xgb <- trainControl(
  method = "cv", number = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  sampling = "smote"
)

# Hyperparameter grid
xgb_grid <- expand.grid(
  nrounds          = c(100, 200),
  max_depth        = c(3, 6),
  eta              = c(0.1, 0.3),
  gamma            = 0,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample        = 0.8
)

for (i in 1:n_imp) {

  if (i %% 20 == 0) message(sprintf("Processing imputation %d of %d", i, n_imp))

  capa_i <- data_stacked %>% filter(.imp == i)
  train_i <- capa_i[indices_train, ] %>%
    dplyr::select(-.imp, -.id) %>%
    mutate(across(where(is.ordered),
                  ~factor(as.character(.x), ordered = FALSE))) %>%
    as.data.frame()
  test_i <- capa_i[-indices_train, ] %>%
    dplyr::select(-.imp, -.id) %>%
    mutate(across(where(is.ordered),
                  ~factor(as.character(.x), ordered = FALSE))) %>%
    as.data.frame()

  if (sum(is.na(train_i)) > 0) {
    warning(sprintf("Imputation %d has NAs — skipping", i)); next
  }

  # Train XGBoost with 10-fold CV and SMOTE within folds
  set.seed(123)
  xgb_model_i <- train(
    psqi ~ sex + age + bmic + bf + hb,
    data = train_i, method = "xgbTree",
    trControl = control_xgb, tuneGrid = xgb_grid,
    metric = "ROC", verbosity = 0
  )

  if (i == 1) {
    message("\nBest hyperparameters (imputation 1):")
    print(xgb_model_i$bestTune)
  }

  # Predict on unbalanced test set
  pred_class_i <- predict(xgb_model_i, newdata = test_i)
  pred_prob_i  <- predict(xgb_model_i, newdata = test_i, type = "prob")

  # ROC and AUC
  roc_obj_i <- roc(test_i$psqi, pred_prob_i$Poor_Sleep, quiet = TRUE)
  auc_values_xgb[i] <- as.numeric(roc_obj_i$auc)

  coords_opt <- coords(roc_obj_i, "best", best.method = "youden",
                        ret = c("sensitivity", "specificity"))
  sensitivity_vals[i] <- coords_opt$sensitivity
  specificity_vals[i] <- coords_opt$specificity

  # Confusion matrix
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
    hl_pvalues_xgb[i] <- hl_i$p.value
  }, error = function(e) { hl_pvalues_xgb[i] <<- NA })

  # ROC curve data
  roc_df_all <- rbind(roc_df_all, data.frame(
    specificity = roc_obj_i$specificities,
    sensitivity = roc_obj_i$sensitivities,
    imputation  = i
  ))

  # Variable importance
  imp_i <- varImp(xgb_model_i, scale = FALSE)$importance
  importance_list[[i]] <- data.frame(
    Variable   = rownames(imp_i),
    Importance = if ("Overall" %in% colnames(imp_i)) imp_i$Overall
                 else rowMeans(imp_i),
    stringsAsFactors = FALSE
  )

  model_list_xgb[[i]] <- xgb_model_i
  test_list_xgb[[i]]  <- test_i
}

saveRDS(auc_values_xgb, "auc_values_xgb.rds")

# =================================================================
# STEP 4 — POOLED RESULTS
# =================================================================

# Performance metrics
performance_xgb <- data.frame(
  Metric = c("AUC", "Sensitivity", "Specificity", "Accuracy",
             "Kappa", "PPV", "NPV", "Balanced Accuracy"),
  Mean = round(c(mean(auc_values_xgb), mean(sensitivity_vals),
                 mean(specificity_vals), mean(accuracy_vals),
                 mean(kappa_vals), mean(ppv_vals, na.rm = TRUE),
                 mean(npv_vals, na.rm = TRUE), mean(bal_acc_vals)), 3),
  SD = round(c(sd(auc_values_xgb), sd(sensitivity_vals),
               sd(specificity_vals), sd(accuracy_vals),
               sd(kappa_vals), sd(ppv_vals, na.rm = TRUE),
               sd(npv_vals, na.rm = TRUE), sd(bal_acc_vals)), 3)
)

message("\n=== XGBOOST PERFORMANCE ===")
print(performance_xgb, row.names = FALSE)
write_xlsx(performance_xgb, "Table_XGB_Performance.xlsx")

# Variable importance
importance_avg <- do.call(rbind, importance_list) %>%
  group_by(Variable) %>%
  summarise(Mean_Importance = mean(Importance),
            SD_Importance = sd(Importance), .groups = "drop") %>%
  arrange(desc(Mean_Importance))

write_xlsx(as.data.frame(importance_avg), "Table_XGB_Variable_Importance.xlsx")

# Confusion matrix
confusion_avg <- Reduce("+", confusion_list) / n_imp

# Hosmer-Lemeshow
hl_summary <- data.frame(
  Statistic = c("Mean p-value", "SD p-value",
                "% adequate (p > 0.05)"),
  Value = c(round(mean(hl_pvalues_xgb, na.rm = TRUE), 3),
            round(sd(hl_pvalues_xgb, na.rm = TRUE), 3),
            round(100 * mean(hl_pvalues_xgb > 0.05, na.rm = TRUE), 1))
)

message("\n=== H-L CALIBRATION ===")
print(hl_summary, row.names = FALSE)
write_xlsx(hl_summary, "Table_XGB_HL.xlsx")

# CV results (from first imputation)
cv_results <- model_list_xgb[[1]]$results
best_idx   <- which.max(cv_results$ROC)
cv_summary <- data.frame(
  Metric = c("CV_ROC", "CV_Sensitivity", "CV_Specificity",
             "Best_nrounds", "Best_max_depth", "Best_eta"),
  Value = c(round(cv_results$ROC[best_idx], 3),
            round(cv_results$Sens[best_idx], 3),
            round(cv_results$Spec[best_idx], 3),
            cv_results$nrounds[best_idx],
            cv_results$max_depth[best_idx],
            cv_results$eta[best_idx])
)

message("\n=== 10-FOLD CV ===")
print(cv_summary, row.names = FALSE)
write_xlsx(cv_summary, "Table_XGB_CV_Results.xlsx")

# =================================================================
# STEP 5 — FIGURE: COMBINED 4-PANEL
# =================================================================

# Panel A: Variable Importance
panel_a <- ggplot(importance_avg,
                  aes(x = Mean_Importance,
                      y = reorder(Variable, Mean_Importance))) +
  geom_col(fill = "#FF9800", width = 0.6, color = "black", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = Mean_Importance - SD_Importance,
                     xmax = Mean_Importance + SD_Importance),
                 height = 0.3, linewidth = 0.4) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_classic(base_size = 14) +
  labs(title = "A. Variable Importance", x = "Gain (± SD)", y = "") +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text = element_text(size = 12, color = "black"))

# Panel B: ROC Curve
panel_b <- ggplot(roc_df_all, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(aes(group = imputation), alpha = 0.6, color = "#FFB74D",
            linewidth = 0.3) +
  stat_summary(fun = mean, geom = "line", color = "#E65100",
               linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "#A020F0", linewidth = 1) +
  annotate("text", x = 0.65, y = 0.25,
           label = sprintf("Mean AUC = %.3f", mean(auc_values_xgb)),
           size = 4, fontface = "bold") +
  theme_classic(base_size = 14) +
  labs(title = "B. ROC Curve", x = "1 - Specificity",
       y = "Sensitivity") +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text = element_text(size = 12, color = "black")) +
  coord_equal()

# Panel C: Confusion Matrix
conf_df <- as.data.frame(confusion_avg)
colnames(conf_df) <- c("Reference", "Prediction", "Freq")
total <- sum(conf_df$Freq)
conf_df$Label <- sprintf("n = %.1f\n(%.1f%%)", conf_df$Freq,
                         conf_df$Freq / total * 100)

panel_c <- ggplot(conf_df, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(color = "white", linewidth = 2) +
  geom_text(aes(label = Label), color = "white", size = 5,
            fontface = "bold", lineheight = 0.85) +
  scale_fill_gradient(low = "#FFE0B2", high = "#E65100") +
  labs(title = "C. Confusion Matrix", x = "Actual", y = "Predicted") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text = element_text(size = 12, face = "bold"),
        legend.position = "none", panel.grid = element_blank())

# Panel D: Calibration
calib_all <- data.frame()
for (i in 1:n_imp) {
  preds_cal <- predict(model_list_xgb[[i]], newdata = test_list_xgb[[i]],
                       type = "prob")$Poor_Sleep
  outcome_cal <- ifelse(test_list_xgb[[i]]$psqi == "Poor_Sleep", 1, 0)
  calib_i <- data.frame(predicted = preds_cal, observed = outcome_cal)
  calib_i$decile <- factor(dplyr::ntile(calib_i$predicted, 10))
  calib_s <- calib_i %>%
    group_by(decile) %>%
    summarise(predicted_mean = mean(predicted),
              observed_mean = mean(observed), .groups = "drop")
  calib_all <- rbind(calib_all, calib_s)
}

calib_avg <- calib_all %>%
  group_by(decile) %>%
  summarise(predicted_mean = mean(predicted_mean),
            observed_mean = mean(observed_mean), .groups = "drop")

panel_d <- ggplot(calib_avg, aes(x = predicted_mean, y = observed_mean)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "#A020F0", linewidth = 1) +
  geom_point(size = 3, color = "#E65100") +
  geom_smooth(method = "loess", se = TRUE, color = "#FF9800",
              fill = "white", linewidth = 0.8) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  annotate("text", x = 0.25, y = 0.9,
           label = sprintf("H-L p = %.3f", mean(hl_pvalues_xgb, na.rm = TRUE)),
           size = 4, fontface = "bold") +
  theme_classic(base_size = 14) +
  labs(title = "D. Calibration Plot",
       x = "Predicted Probability", y = "Observed Probability") +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text = element_text(size = 12, color = "black"))

# Combine and save
fig_combined <- (panel_a + panel_b) / (panel_c + panel_d)

ggsave("Fig_XGB_Combined_4Panel.tif", fig_combined,
       device = "tiff", width = 24, height = 22, units = "cm",
       dpi = 1200, compression = "lzw", bg = "white")

# =================================================================
# SUMMARY
# =================================================================

message("\n========================================")
message("XGBOOST ANALYSIS COMPLETED")
message("========================================")
message(sprintf("Mean AUC: %.3f (SD = %.3f)",
                mean(auc_values_xgb), sd(auc_values_xgb)))
message(sprintf("H-L: mean p = %.3f (%s%% adequate)",
                mean(hl_pvalues_xgb, na.rm = TRUE),
                round(100 * mean(hl_pvalues_xgb > 0.05, na.rm = TRUE), 1)))
message("========================================")
