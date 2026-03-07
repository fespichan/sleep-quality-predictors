############################################################
# 06_LR_Pooled_PSQI.R
# Logistic Regression: Pooled results (Rubin's rules),
# predictive performance, and diagnostic plots
#
# Project: Biological and Anthropometric Determinants of
#          Sleep Quality Across the Adult Life Course
# Journal: Sleep Medicine
#
# Requires: data_stacked (from 02_Multiple_Imputation.R)
#
# Output:
#   Tables:
#     Table2_LR_Pooled_NoSMOTE.xlsx   (primary inference)
#     Table2_LR_Pooled_SMOTE.xlsx     (with SMOTE)
#     Table3_LR_Performance.xlsx      (classification metrics)
#     Table_LR_CV_Results.xlsx        (cross-validation)
#     Table_LR_HL_NoSMOTE.xlsx        (calibration without SMOTE)
#     Table_LR_HL_SMOTE.xlsx          (calibration with SMOTE)
#     Table_LR_VIF.xlsx               (multicollinearity)
#   Figures:
#     Fig1_LR_Combined_4Panel.tif     (ROC + Calibration + 2 CMs)
#     Fig_LR_Forest_Plot_Pooled.tif   (Forest plot)
############################################################

library(dplyr)
library(ggplot2)
library(mice)
library(caret)
library(pROC)
library(car)
library(themis)
library(recipes)
library(ResourceSelection)
library(cowplot)
library(writexl)

n_imp <- max(data_stacked$.imp)

# =================================================================
# STEP 1 — CONSISTENT TRAIN/TEST PARTITION
# =================================================================
# Same subjects in train/test across all 100 imputed datasets

set.seed(123456)
first_imp <- data_stacked %>% filter(.imp == 1)
indices_train <- createDataPartition(first_imp$psqi, p = 0.80, list = FALSE)

# =================================================================
# STEP 2 — POOLED OR WITHOUT SMOTE (PRIMARY INFERENCE)
# =================================================================
# Rubin's rules applied to full datasets (no train/test split)
# This is the inferential model reported in the manuscript

model_list_no_smote <- list()

for (i in 1:n_imp) {
  capa_i <- data_stacked %>% filter(.imp == i) %>%
    dplyr::select(-.imp, -.id) %>%
    mutate(across(where(is.ordered),
                  ~factor(as.character(.x), ordered = FALSE))) %>%
    as.data.frame()

  model_list_no_smote[[i]] <- glm(psqi ~ sex + age + bmic + bf + hb,
                                   data = capa_i, family = "binomial")
}

pooled_no_smote <- pool(as.mira(model_list_no_smote))
table_no_smote  <- summary(pooled_no_smote, exponentiate = TRUE,
                           conf.int = TRUE)

message("\n=== POOLED OR — WITHOUT SMOTE (PRIMARY) ===")
print(table_no_smote)
write_xlsx(as.data.frame(table_no_smote), "Table2_LR_Pooled_NoSMOTE.xlsx")

# =================================================================
# STEP 3 — TRAIN MODELS WITH SMOTE (PREDICTIVE)
# =================================================================

train_balanced_list <- list()
test_list  <- list()
model_list <- list()

for (i in 1:n_imp) {

  if (i %% 20 == 0) message(sprintf("Processing imputation %d of %d", i, n_imp))

  capa_i <- data_stacked %>% filter(.imp == i)

  # Verify completeness
  if (sum(is.na(capa_i)) > 0) {
    warning(sprintf("Imputation %d has NAs — skipping", i))
    next
  }

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

  # Apply SMOTE-NC to training set only
  tryCatch({
    smote_recipe <- recipe(psqi ~ ., data = train_i) %>%
      step_smotenc(psqi, over_ratio = 1, seed = 123, neighbors = 5)
    smote_prep <- prep(smote_recipe, training = train_i, verbose = FALSE)
    train_i_balanced <- bake(smote_prep, new_data = NULL)
  }, error = function(e) {
    warning(sprintf("SMOTE failed for imputation %d — using unbalanced data", i))
    train_i_balanced <<- train_i
  })

  # Train logistic regression
  model_i <- glm(psqi ~ sex + age + bmic + bf + hb,
                 data = train_i_balanced, family = "binomial")

  train_balanced_list[[i]] <- train_i_balanced
  test_list[[i]]  <- test_i
  model_list[[i]] <- model_i
}

# Pool SMOTE models (for comparison)
pooled_smote <- pool(as.mira(model_list))
table_smote  <- summary(pooled_smote, exponentiate = TRUE, conf.int = TRUE)

message("\n=== POOLED OR — WITH SMOTE ===")
print(table_smote)
write_xlsx(as.data.frame(table_smote), "Table2_LR_Pooled_SMOTE.xlsx")
saveRDS(as.data.frame(table_smote), "final_table_logit.rds")

# =================================================================
# STEP 4 — PERFORMANCE METRICS (AUC, SENSITIVITY, SPECIFICITY)
# =================================================================

auc_values_logit <- numeric(n_imp)
optimal_thresholds <- numeric(n_imp)

for (i in 1:n_imp) {
  preds_prob <- predict(model_list[[i]], newdata = test_list[[i]],
                        type = "response")
  roc_obj <- roc(test_list[[i]]$psqi, preds_prob, quiet = TRUE)
  auc_values_logit[i] <- as.numeric(roc_obj$auc)
  optimal_thresholds[i] <- coords(roc_obj, "best",
                                   best.method = "youden")$threshold
}

saveRDS(auc_values_logit, "auc_values_logit.rds")

# Detailed metrics using average optimal threshold
threshold_avg <- mean(optimal_thresholds)
detailed_metrics_list <- list()

for (i in 1:n_imp) {
  preds_prob  <- predict(model_list[[i]], newdata = test_list[[i]],
                         type = "response")
  preds_class <- factor(ifelse(preds_prob > threshold_avg,
                               "p.sleep", "g.sleep"),
                        levels = c("g.sleep", "p.sleep"))

  cm <- confusionMatrix(preds_class, test_list[[i]]$psqi,
                        positive = "p.sleep")

  detailed_metrics_list[[i]] <- data.frame(
    AUC              = auc_values_logit[i],
    Accuracy         = cm$overall["Accuracy"],
    Kappa            = cm$overall["Kappa"],
    Sensitivity      = cm$byClass["Sensitivity"],
    Specificity      = cm$byClass["Specificity"],
    PosPredValue     = cm$byClass["Pos Pred Value"],
    NegPredValue     = cm$byClass["Neg Pred Value"],
    BalancedAccuracy = cm$byClass["Balanced Accuracy"]
  )
}

all_metrics <- do.call(rbind, detailed_metrics_list)

metrics_summary <- data.frame(
  Metric = colnames(all_metrics),
  Mean   = round(colMeans(all_metrics, na.rm = TRUE), 3),
  SD     = round(apply(all_metrics, 2, sd, na.rm = TRUE), 3)
)

message("\n=== CLASSIFICATION METRICS ===")
print(metrics_summary, row.names = FALSE)
write_xlsx(metrics_summary, "Table3_LR_Performance.xlsx")

# =================================================================
# STEP 5 — VIF (MULTICOLLINEARITY ASSESSMENT)
# =================================================================

vif_list <- list()
for (i in 1:n_imp) {
  vif_i <- car::vif(model_list[[i]])
  vif_list[[i]] <- if (is.matrix(vif_i)) vif_i[, ncol(vif_i)] else vif_i
}

vif_df  <- as.data.frame(do.call(rbind, vif_list))
vif_avg <- round(colMeans(vif_df), 3)

message("\n=== AVERAGE VIF ===")
print(vif_avg)
write_xlsx(vif_df, "Table_LR_VIF.xlsx")

# =================================================================
# STEP 6 — HOSMER-LEMESHOW CALIBRATION
# =================================================================
# NOTE: Formal H-L calibration analysis (with and without SMOTE)
# is performed in 12_HL_Calibration.R using independent partitions.
# Results: Without SMOTE mean p = 0.378, 87% adequate
#          With SMOTE    mean p = 0.005, 2% adequate
# The hl_p_no_smote variable is set here for use in the
# calibration plot (Panel B of Figure 1).

hl_p_no_smote_for_plot <- 0.378  # from 12_HL_Calibration.R

# =================================================================
# STEP 7 — 10-FOLD CROSS-VALIDATION
# =================================================================

first_imp_cv <- data_stacked %>%
  filter(.imp == 1) %>%
  dplyr::select(-.imp, -.id) %>%
  mutate(across(where(is.ordered),
                ~factor(as.character(.x), ordered = FALSE))) %>%
  as.data.frame()

set.seed(123456)
cv_control <- trainControl(
  method = "cv", number = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  sampling = "smote"
)

cv_model <- train(
  psqi ~ sex + age + bmic + bf + hb,
  data = first_imp_cv, method = "glm", family = "binomial",
  trControl = cv_control, metric = "ROC"
)

cv_summary <- data.frame(
  Metric = c("CV_ROC", "CV_Sensitivity", "CV_Specificity"),
  Value  = c(cv_model$results$ROC,
             cv_model$results$Sens,
             cv_model$results$Spec)
)

message("\n=== 10-FOLD CV WITH SMOTE ===")
print(cv_summary, row.names = FALSE)
write_xlsx(cv_summary, "Table_LR_CV_Results.xlsx")

# =================================================================
# STEP 8 — FIGURE 1: COMBINED 4-PANEL PLOT
# =================================================================

# --- Panel A: ROC Curve ---
roc_df_all <- data.frame()
for (i in 1:n_imp) {
  preds_prob <- predict(model_list[[i]], newdata = test_list[[i]],
                        type = "response")
  roc_obj <- roc(test_list[[i]]$psqi, preds_prob, quiet = TRUE)
  roc_df_all <- rbind(roc_df_all, data.frame(
    specificity = roc_obj$specificities,
    sensitivity = roc_obj$sensitivities,
    imputation  = i
  ))
}

panel_a <- ggplot(roc_df_all, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(aes(group = imputation), alpha = 1, color = "#00C5CD",
            linewidth = 0.4) +
  stat_summary(fun = mean, geom = "line", color = "#E74C3C",
               linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "#A020F0", linewidth = 1) +
  annotate("text", x = 0.7, y = 0.3,
           label = sprintf("Mean AUC = %.3f", mean(auc_values_logit)),
           size = 4.5, fontface = "bold") +
  labs(title = "A. ROC Curve", x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.background = element_rect(fill = "white", color = "white")
  ) +
  coord_equal()

# --- Panel B: Calibration Plot ---
# Average calibration across imputations using ntile
calib_data <- data.frame()
for (i in 1:n_imp) {
  preds_prob <- predict(model_list[[i]], newdata = test_list[[i]],
                        type = "response")
  outcome_i <- ifelse(test_list[[i]]$psqi == "p.sleep", 1, 0)
  calib_i <- data.frame(predicted = preds_prob, observed = outcome_i)
  calib_i$decile <- factor(dplyr::ntile(calib_i$predicted, 10))
  calib_data <- rbind(calib_data, calib_i)
}

calib_summary <- calib_data %>%
  group_by(decile) %>%
  summarise(predicted_mean = mean(predicted),
            observed_mean  = mean(observed), .groups = "drop")

panel_b <- ggplot(calib_summary,
                  aes(x = predicted_mean, y = observed_mean)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "#A020F0", linewidth = 1) +
  geom_point(size = 3, color = "#E74C3C") +
  geom_smooth(method = "loess", se = TRUE, color = "#3498DB",
              fill = "white", linewidth = 0.8) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  annotate("text", x = 0.25, y = 0.9,
           label = sprintf("H-L, p = %.3f", hl_p_no_smote_for_plot),
           size = 4.5, fontface = "bold") +
  labs(title = "B. Calibration Plot",
       x = "Predicted Probability", y = "Observed Probability") +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.background = element_rect(fill = "white", color = "white")
  ) +
  coord_equal()

# --- Panels C & D: Confusion Matrices ---
confusion_list <- list()
for (i in 1:n_imp) {
  preds_prob  <- predict(model_list[[i]], newdata = test_list[[i]],
                         type = "response")
  preds_class <- factor(ifelse(preds_prob > threshold_avg,
                               "p.sleep", "g.sleep"),
                        levels = c("g.sleep", "p.sleep"))
  confusion_list[[i]] <- confusionMatrix(preds_class,
                                          test_list[[i]]$psqi)$table
}

confusion_avg <- Reduce("+", confusion_list) / n_imp
conf_df <- as.data.frame(confusion_avg)
colnames(conf_df) <- c("Reference", "Prediction", "Freq")

total_cases <- sum(conf_df$Freq)
conf_df$Percentage <- (conf_df$Freq / total_cases) * 100
conf_df$Label <- sprintf("n = %.1f\n(%.1f%%)", conf_df$Freq, conf_df$Percentage)

# Row percentages
conf_df_row <- conf_df %>%
  group_by(Reference) %>%
  mutate(Pct_Row = (Freq / sum(Freq)) * 100,
         Label_Row = sprintf("%.1f\n(%.1f%%)", Freq, Pct_Row)) %>%
  ungroup()

panel_c <- ggplot(conf_df, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(color = "white", linewidth = 2) +
  geom_text(aes(label = Label), color = "white", size = 5,
            fontface = "bold", lineheight = 0.85) +
  scale_fill_gradient(low = "#3498DB", high = "#E74C3C") +
  annotate("text", x = 0.65, y = 0.65, label = "TN",
           color = "yellow", size = 5, fontface = "bold") +
  annotate("text", x = 1.35, y = 0.65, label = "FP",
           color = "yellow", size = 5, fontface = "bold") +
  annotate("text", x = 0.65, y = 1.35, label = "FN",
           color = "yellow", size = 5, fontface = "bold") +
  annotate("text", x = 1.35, y = 1.35, label = "TP",
           color = "yellow", size = 5, fontface = "bold") +
  labs(title = "C. Confusion Matrix (Total %)",
       x = "Actual", y = "Predicted") +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 14, face = "bold"),
    legend.position = "none", panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.background = element_rect(fill = "white", color = "white")
  )

panel_d <- ggplot(conf_df_row,
                  aes(x = Reference, y = Prediction, fill = Pct_Row)) +
  geom_tile(color = "white", linewidth = 2) +
  geom_text(aes(label = Label_Row), color = "white", size = 5,
            fontface = "bold", lineheight = 0.9) +
  scale_fill_gradient(low = "#27AE60", high = "#C0392B") +
  labs(title = "D. Confusion Matrix (Row %)",
       x = "Actual", y = "Predicted") +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 14, face = "bold"),
    legend.position = "none", panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.background = element_rect(fill = "white", color = "white")
  )

# Combine and save
combined_plot <- plot_grid(panel_a, panel_b, panel_c, panel_d,
                           ncol = 2, nrow = 2, align = "hv", axis = "tblr")

ggsave("Fig1_LR_Combined_4Panel.tif", combined_plot,
       device = "tiff", width = 20, height = 20, units = "cm",
       dpi = 1200, compression = "lzw", bg = "white")

# =================================================================
# STEP 9 — FOREST PLOT (POOLED OR)
# =================================================================

forest_data <- table_smote %>%
  filter(term != "(Intercept)") %>%
  mutate(
    Significance = case_when(
      p.value < 0.001 ~ "p < 0.001",
      p.value < 0.01  ~ "p < 0.01",
      p.value < 0.05  ~ "p < 0.05",
      TRUE            ~ "ns"
    ),
    label = sprintf("%.2f (%.2f-%.2f)", estimate, `2.5 %`, `97.5 %`)
  )

forest_plot <- ggplot(forest_data,
                      aes(x = estimate, y = reorder(term, estimate))) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "#0000CD", linewidth = 1) +
  geom_errorbarh(aes(xmin = `2.5 %`, xmax = `97.5 %`),
                 height = 0.25, linewidth = 0.5) +
  geom_point(aes(color = Significance), size = 3.5) +
  geom_text(aes(x = `97.5 %`, label = label),
            hjust = -0.1, size = 3.5, fontface = "italic") +
  scale_color_manual(
    values = c("p < 0.001" = "#EEC900", "p < 0.01" = "#00CD00",
               "p < 0.05" = "#483D8B", "ns" = "#FF3030")
  ) +
  scale_x_log10() +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(10, 80, 10, 10)
  ) +
  labs(title = "Forest Plot: Pooled OR (Rubin's Rules, with SMOTE)",
       x = "Odds Ratio (95% CI, log scale)", y = "")

ggsave("Fig_LR_Forest_Plot_Pooled.tif", forest_plot,
       device = "tiff", width = 22, height = 14, units = "cm",
       dpi = 1200, compression = "lzw", bg = "white")

# =================================================================
# SUMMARY
# =================================================================

message("\n========================================")
message("LOGISTIC REGRESSION ANALYSIS COMPLETED")
message("========================================")
message(sprintf("Mean AUC: %.3f (SD = %.3f)", mean(auc_values_logit),
                sd(auc_values_logit)))
message(sprintf("Optimal threshold: %.3f", threshold_avg))
message("H-L calibration: see 12_HL_Calibration.R")
message("  Without SMOTE: mean p = 0.378, 87% adequate")
message("  With SMOTE:    mean p = 0.005, 2% adequate")
message("========================================")

