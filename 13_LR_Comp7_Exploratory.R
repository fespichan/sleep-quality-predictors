############################################################
# 13_LR_Comp7_Exploratory.R
# Exploratory logistic regression for PSQI Component 7
# (Daytime Dysfunction): pooled inference + predictive
# performance with SMOTE
#
# Project: Biological and Anthropometric Determinants of
#          Sleep Quality Across the Adult Life Course
# Journal: Sleep Medicine
#
# Requires: data_stacked (from 02_Multiple_Imputation.R)
#
# Output:
#   Tables:
#     Table_LR_Comp7_Pooled_NoSMOTE.xlsx  (primary inference)
#     Table_LR_Comp7_Pooled_SMOTE.xlsx
#     Table_LR_Comp7_Performance.xlsx
#     Table_LR_Comp7_HL.xlsx
#     Table_LR_Comp7_CV.xlsx
#     Table_LR_Comp7_SMOTE_Justification.xlsx
#   Figures:
#     Fig_S7_LR_Comp7_Combined_4Panel.tif
#   Data:
#     auc_values_lr_comp7.rds
#     final_table_lr_comp7.rds
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
library(patchwork)
library(writexl)


# Extract all imputations in long format include original data
data_stacked <- complete(imputed_data, action = "long", include = TRUE)

# Round PSQI components to integers and clamp to 0-3 scale
data_stacked <- data_stacked %>%
  mutate(across(all_of(likert_vars),
                ~as.integer(pmin(pmax(round(.x), 0), 3))))

message("\nPSQI components rounded to 0-3 scale:")
for (v in likert_vars) {
  message(sprintf("  %s: [%d, %d]", v,
                  min(data_stacked[[v]]), max(data_stacked[[v]])))
}

# Export single imputation (dataset 10) for exploratory use
sleep_data_imp <- complete(imputed_data, 10)
write_xlsx(sleep_data_imp, "sleep_data_imp_10.xlsx")


# COMPUTE DERIVED VARIABLES
# Applied to ALL 100 imputed datasets to ensure internal consistency.

# BMI from imputed weight and height
data_stacked <- data_stacked %>%
  mutate(bmi = weight / (height / 100)^2)

# PSQI global score (sum of 7 components)
data_stacked <- data_stacked %>%
  mutate(
    across(Comp.1:Comp.7, ~as.numeric(as.character(.x))),
    psqi.total = Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5 + Comp.6 + Comp.7
  )

# PSQI classification: good sleep (<=5) vs. poor sleep (>5)
data_stacked <- data_stacked %>%
  mutate(psqi = case_when(
    psqi.total <= 5 ~ "g.sleep",
    psqi.total > 5  ~ "p.sleep",
    TRUE ~ NA_character_
  ))

# BMI category (WHO criteria; underweight merged with Normal)
data_stacked <- data_stacked %>%
  mutate(bmic = case_when(
    bmi < 25.0                ~ "Normal",
    bmi >= 25.0 & bmi < 30.0 ~ "Overweight",
    bmi >= 30.0               ~ "Obesity",
    TRUE ~ NA_character_
  ))

data_stacked$bmic <- factor(data_stacked$bmic,
                            levels = c("Normal", "Overweight", "Obesity"))

# Body fat category (sex-specific, ACE/ACSM guidelines)
data_stacked <- data_stacked %>%
  mutate(bfc = case_when(
    sex == "female" & bf >= 5  & bf < 14 ~ "Essential.fat",
    sex == "female" & bf >= 14 & bf < 21 ~ "Athletes",
    sex == "female" & bf >= 21 & bf < 25 ~ "Fitness",
    sex == "female" & bf >= 25 & bf < 32 ~ "Average",
    sex == "female" & bf >= 32            ~ "Obese",
    sex == "male"   & bf >= 2  & bf < 6  ~ "Essential.fat",
    sex == "male"   & bf >= 6  & bf < 14 ~ "Athletes",
    sex == "male"   & bf >= 14 & bf < 18 ~ "Fitness",
    sex == "male"   & bf >= 18 & bf < 25 ~ "Average",
    sex == "male"   & bf >= 25            ~ "Obese",
    TRUE ~ NA_character_
  ))

# Anaemia classification (WHO thresholds by sex)
data_stacked <- data_stacked %>%
  mutate(hbc = case_when(
    sex == "female" & hb < 12.0  ~ "anemia",
    sex == "female" & hb >= 12.0 ~ "normal",
    sex == "male"   & hb < 13.0  ~ "anemia",
    sex == "male"   & hb >= 13.0 ~ "normal",
    TRUE ~ NA_character_
  ))

# Convert categorical variables to factors
data_stacked <- data_stacked %>%
  mutate(
    psqi = factor(psqi, levels = c("g.sleep", "p.sleep")),
    bmic = factor(bmic, levels = c("Normal", "Overweight", "Obesity")),
    bfc  = as.factor(bfc),
    hbc  = as.factor(hbc),
    sex  = as.factor(sex)
  )

# Reorder columns for logical structure
data_stacked <- data_stacked %>%
  dplyr::select(
    .imp, .id,
    sex, age, height, weight, bmi, wc, bf, mm, water, hb,
    Comp.1, Comp.2, Comp.3, Comp.4, Comp.5, Comp.6, Comp.7,
    psqi.total, psqi, bmic, bfc, hbc
  )
n_imp <- max(data_stacked$.imp)

# =================================================================
# STEP 1 — CREATE BINARY OUTCOME (COMP.7)
# =================================================================

data_stacked <- data_stacked %>%
  mutate(
    Comp7_bin = factor(
      ifelse(Comp.7 == 0, "No_dysfunction", "Has_dysfunction"),
      levels = c("No_dysfunction", "Has_dysfunction")
    )
  )

message("\n=== COMP.7 DISTRIBUTION (first imputation) ===")
dist_comp7 <- data_stacked %>%
  filter(.imp == 1) %>%
  count(Comp7_bin) %>%
  mutate(pct = round(n / sum(n) * 100, 1))
print(dist_comp7)

# Consistent train/test partition
set.seed(123456)
first_imp <- data_stacked %>% filter(.imp == 1)
indices_train <- createDataPartition(first_imp$Comp7_bin, p = 0.80, list = FALSE)

# =================================================================
# STEP 2 — POOLED OR WITHOUT SMOTE (PRIMARY INFERENCE)
# =================================================================

imputed_long <- data_stacked %>%
  dplyr::select(.imp, .id, sex, age, bmic, bf, hb, Comp7_bin) %>%
  mutate(across(where(is.ordered),
                ~factor(as.character(.x), ordered = FALSE)))

imputed_mids <- as.mids(imputed_long)

mira_no_smote <- with(imputed_mids,
                       glm(Comp7_bin ~ sex + age + bmic + bf + hb,
                           family = "binomial"))

pooled_no_smote <- pool(mira_no_smote)
table_no_smote  <- summary(pooled_no_smote, exponentiate = TRUE,
                           conf.int = TRUE)

message("\n=== POOLED OR — WITHOUT SMOTE ===")
print(table_no_smote)

saveRDS(table_no_smote, "final_table_lr_comp7_no_smote.rds")
write_xlsx(as.data.frame(table_no_smote), "Table_LR_Comp7_Pooled_NoSMOTE.xlsx")

# AUC without SMOTE (for SMOTE justification)
auc_no_smote  <- numeric(n_imp)
sens_no_smote <- numeric(n_imp)
spec_no_smote <- numeric(n_imp)

for (i in 1:n_imp) {
  capa_i <- data_stacked %>% filter(.imp == i)
  train_i <- capa_i[indices_train, ] %>% dplyr::select(-.imp, -.id) %>%
    mutate(across(where(is.ordered),
                  ~factor(as.character(.x), ordered = FALSE))) %>%
    as.data.frame()
  test_i <- capa_i[-indices_train, ] %>% dplyr::select(-.imp, -.id) %>%
    mutate(across(where(is.ordered),
                  ~factor(as.character(.x), ordered = FALSE))) %>%
    as.data.frame()

  m_i <- glm(Comp7_bin ~ sex + age + bmic + bf + hb,
             data = train_i, family = "binomial")
  p_i <- predict(m_i, newdata = test_i, type = "response")

  roc_i <- roc(test_i$Comp7_bin, p_i, quiet = TRUE)
  auc_no_smote[i] <- as.numeric(roc_i$auc)

  coords_i <- coords(roc_i, "best", best.method = "youden",
                      ret = c("sensitivity", "specificity"))
  sens_no_smote[i] <- coords_i$sensitivity
  spec_no_smote[i] <- coords_i$specificity
}

message(sprintf("\nWithout SMOTE: AUC = %.3f (SD = %.3f)",
                mean(auc_no_smote), sd(auc_no_smote)))

# =================================================================
# STEP 3 — TRAIN WITH SMOTE (PREDICTIVE)
# =================================================================

model_list_lr    <- list()
test_list_lr     <- list()
confusion_list   <- list()
roc_df_all       <- data.frame()

auc_values       <- numeric(n_imp)
sensitivity_vals <- numeric(n_imp)
specificity_vals <- numeric(n_imp)
accuracy_vals    <- numeric(n_imp)
kappa_vals       <- numeric(n_imp)
ppv_vals         <- numeric(n_imp)
npv_vals         <- numeric(n_imp)
bal_acc_vals     <- numeric(n_imp)
optimal_thresholds <- numeric(n_imp)
hl_pvalues       <- numeric(n_imp)

for (i in 1:n_imp) {

  if (i %% 20 == 0) message(sprintf("Processing imputation %d of %d", i, n_imp))

  capa_i <- data_stacked %>% filter(.imp == i)
  train_i <- capa_i[indices_train, ] %>% dplyr::select(-.imp, -.id) %>%
    mutate(across(where(is.ordered),
                  ~factor(as.character(.x), ordered = FALSE))) %>%
    as.data.frame()
  test_i <- capa_i[-indices_train, ] %>% dplyr::select(-.imp, -.id) %>%
    mutate(across(where(is.ordered),
                  ~factor(as.character(.x), ordered = FALSE))) %>%
    as.data.frame()

  # SMOTE on training set
  rec_i <- recipe(Comp7_bin ~ sex + age + bmic + bf + hb, data = train_i) %>%
    step_smotenc(Comp7_bin, over_ratio = 1, neighbors = 5)
  train_balanced <- prep(rec_i) %>% juice()

  # Train model
  model_i <- glm(Comp7_bin ~ sex + age + bmic + bf + hb,
                 data = train_balanced, family = "binomial")

  # Predict
  pred_prob_i <- predict(model_i, newdata = test_i, type = "response")
  roc_obj_i   <- roc(test_i$Comp7_bin, pred_prob_i, quiet = TRUE)
  auc_values[i] <- as.numeric(roc_obj_i$auc)

  coords_opt <- coords(roc_obj_i, "best", best.method = "youden",
                        ret = c("threshold", "sensitivity", "specificity"))
  optimal_thresholds[i] <- coords_opt$threshold
  sensitivity_vals[i]   <- coords_opt$sensitivity
  specificity_vals[i]   <- coords_opt$specificity

  # Confusion matrix
  threshold_avg <- mean(optimal_thresholds[1:i])
  pred_class_i  <- factor(
    ifelse(pred_prob_i > threshold_avg, "Has_dysfunction", "No_dysfunction"),
    levels = c("No_dysfunction", "Has_dysfunction")
  )
  cm_i <- confusionMatrix(pred_class_i, test_i$Comp7_bin,
                          positive = "Has_dysfunction")
  accuracy_vals[i] <- cm_i$overall["Accuracy"]
  kappa_vals[i]    <- cm_i$overall["Kappa"]
  ppv_vals[i]      <- cm_i$byClass["Pos Pred Value"]
  npv_vals[i]      <- cm_i$byClass["Neg Pred Value"]
  bal_acc_vals[i]  <- cm_i$byClass["Balanced Accuracy"]
  confusion_list[[i]] <- cm_i$table

  # Hosmer-Lemeshow
  outcome_i <- ifelse(test_i$Comp7_bin == "Has_dysfunction", 1, 0)
  tryCatch({
    hl_i <- hoslem.test(outcome_i, pred_prob_i, g = 10)
    hl_pvalues[i] <- hl_i$p.value
  }, error = function(e) { hl_pvalues[i] <<- NA })

  # ROC data
  roc_df_all <- rbind(roc_df_all, data.frame(
    specificity = roc_obj_i$specificities,
    sensitivity = roc_obj_i$sensitivities,
    imputation  = i
  ))

  model_list_lr[[i]] <- model_i
  test_list_lr[[i]]  <- test_i
}

saveRDS(auc_values, "auc_values_lr_comp7.rds")

# Pool SMOTE models
pooled_smote <- pool(as.mira(model_list_lr))
table_smote  <- summary(pooled_smote, exponentiate = TRUE, conf.int = TRUE)
Table_LR_Comp7_Pooled_SMOTE <- table_smote

message("\n=== POOLED OR — WITH SMOTE ===")
print(table_smote)

saveRDS(table_smote, "final_table_lr_comp7.rds")
write_xlsx(as.data.frame(table_smote), "Table_LR_Comp7_Pooled_SMOTE.xlsx")

# =================================================================
# STEP 4 — PERFORMANCE METRICS
# =================================================================

performance_lr <- data.frame(
  Metric = c("AUC", "Sensitivity", "Specificity", "Accuracy",
             "Kappa", "PPV", "NPV", "Balanced Accuracy"),
  Mean = round(c(mean(auc_values), mean(sensitivity_vals),
                 mean(specificity_vals), mean(accuracy_vals),
                 mean(kappa_vals), mean(ppv_vals, na.rm = TRUE),
                 mean(npv_vals, na.rm = TRUE), mean(bal_acc_vals)), 3),
  SD = round(c(sd(auc_values), sd(sensitivity_vals),
               sd(specificity_vals), sd(accuracy_vals),
               sd(kappa_vals), sd(ppv_vals, na.rm = TRUE),
               sd(npv_vals, na.rm = TRUE), sd(bal_acc_vals)), 3)
)

message("\n=== COMP.7 PERFORMANCE (WITH SMOTE) ===")
print(performance_lr, row.names = FALSE)
write_xlsx(performance_lr, "Table_LR_Comp7_Performance.xlsx")

confusion_avg <- Reduce("+", confusion_list) / n_imp

# Hosmer-Lemeshow summary
hl_summary <- data.frame(
  Statistic = c("Mean p-value", "SD p-value",
                "% adequate (p > 0.05)"),
  Value = c(round(mean(hl_pvalues, na.rm = TRUE), 3),
            round(sd(hl_pvalues, na.rm = TRUE), 3),
            round(100 * mean(hl_pvalues > 0.05, na.rm = TRUE), 1))
)
write_xlsx(hl_summary, "Table_LR_Comp7_HL.xlsx")

# =================================================================
# STEP 5 — SMOTE JUSTIFICATION
# =================================================================

wt_auc  <- wilcox.test(auc_no_smote, auc_values, paired = TRUE)
wt_sens <- wilcox.test(sens_no_smote, sensitivity_vals, paired = TRUE)

smote_table <- data.frame(
  Approach = c("Without SMOTE", "With SMOTE"),
  Mean_AUC = round(c(mean(auc_no_smote), mean(auc_values)), 3),
  SD_AUC   = round(c(sd(auc_no_smote), sd(auc_values)), 3),
  Mean_Sens = round(c(mean(sens_no_smote), mean(sensitivity_vals)), 3),
  Wilcoxon_AUC_p  = c(NA, round(wt_auc$p.value, 4)),
  Wilcoxon_Sens_p = c(NA, round(wt_sens$p.value, 4))
)

message("\n=== SMOTE JUSTIFICATION ===")
print(smote_table, row.names = FALSE)
write_xlsx(smote_table, "Table_LR_Comp7_SMOTE_Justification.xlsx")

# =================================================================
# STEP 6 — 10-FOLD CROSS-VALIDATION
# =================================================================

first_imp_cv <- data_stacked %>%
  filter(.imp == 1) %>%
  dplyr::select(-.imp, -.id) %>%
  mutate(across(where(is.ordered),
                ~factor(as.character(.x), ordered = FALSE))) %>%
  as.data.frame()

set.seed(123456)
cv_model <- train(
  Comp7_bin ~ sex + age + bmic + bf + hb,
  data = first_imp_cv, method = "glm", family = "binomial",
  trControl = trainControl(method = "cv", number = 10,
                           classProbs = TRUE,
                           summaryFunction = twoClassSummary,
                           sampling = "smote"),
  metric = "ROC"
)

cv_summary <- data.frame(
  Metric = c("CV_ROC", "CV_Sensitivity", "CV_Specificity"),
  Value  = round(c(cv_model$results$ROC,
                   cv_model$results$Sens,
                   cv_model$results$Spec), 3)
)

message("\n=== 10-FOLD CV ===")
print(cv_summary, row.names = FALSE)
write_xlsx(cv_summary, "Table_LR_Comp7_CV.xlsx")

# =================================================================
# STEP 7 — FIGURE S7: COMBINED 4-PANEL
# =================================================================
Table_LR_Comp7_Pooled_NO_SMOTE<-readRDS("final_table_lr_comp7_no_smote.rds")
Table_LR_Comp7_Pooled_NO_SMOTE <- Table_LR_Comp7_Pooled_NO_SMOTE %>%
  mutate(term = case_when(
    term == "bmicOverweight" ~ "Overweight",
    term == "bmicObesity"    ~ "Obesity",
    TRUE ~ term
  ))



# Panel A: Forest Plot (No SMOTE OR)
forest_data <- Table_LR_Comp7_Pooled_NO_SMOTE %>%
  filter(term != "(Intercept)") %>%
  mutate(
    Significance = case_when(
      p.value < 0.001 ~ "p < 0.001",
      p.value < 0.01  ~ "p < 0.01",
      p.value < 0.05  ~ "p < 0.05",
      TRUE ~ "ns"
    ),
    label = sprintf("%.2f (%.2f-%.2f)", estimate, `2.5 %`, `97.5 %`)
  )

panel_a <- ggplot(forest_data,
                  aes(x = estimate, y = reorder(term, estimate))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50",
             linewidth = 0.6) +
  geom_errorbarh(aes(xmin = `2.5 %`, xmax = `97.5 %`),
                 height = 0.25, linewidth = 0.6) +
  geom_point(aes(color = Significance), size = 3.5) +
  geom_text(aes(x = `97.5 %`, label = label),
            hjust = -0.1, size = 3.5, fontface = "italic") +
  scale_color_manual(
    values = c("p < 0.001" = "#EEC900", "p < 0.01" = "#C0FF3E",
               "p < 0.05" = "#2196F3", "ns" = "#EE3B3B")
  ) +
  scale_x_log10() +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 14) +
  labs(title = "A. Forest Plot: Pooled OR (Rubin's Rules, No SMOTE)",
       x = "Odds Ratio (95% CI, log scale)", y = "") +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(10, 80, 10, 10)
  )

# Panel B: ROC Curve
panel_b <- ggplot(roc_df_all, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(aes(group = imputation), alpha = 0.5, color = "#00C5CD",
            linewidth = 0.5) +
  stat_summary(fun = mean, geom = "line", color = "#E74C3C",
               linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "#A020F0", linewidth = 1) +
  annotate("text", x = 0.65, y = 0.25,
           label = sprintf("Mean AUC = %.3f", mean(auc_values)),
           size = 4, fontface = "bold") +
  labs(title = "B. ROC Curve (with SMOTE)",
       x = "1 - Specificity", y = "Sensitivity") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text = element_text(size = 12, color = "black")) +
  coord_equal()

# Panel C: Calibration Plot
calib_all <- data.frame()
for (i in 1:n_imp) {
  preds_cal <- predict(model_list_lr[[i]], newdata = test_list_lr[[i]],
                       type = "response")
  outcome_cal <- ifelse(test_list_lr[[i]]$Comp7_bin == "Has_dysfunction", 1, 0)
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

panel_c <- ggplot(calib_avg, aes(x = predicted_mean, y = observed_mean)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "#A020F0", linewidth = 1) +
  geom_point(size = 3, color = "#E74C3C") +
  geom_smooth(method = "loess", se = TRUE, color = "#E74C3C",
              fill = "white", linewidth = 0.8) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  annotate("text", x = 0.25, y = 0.9,
           label = sprintf("H-L p = %.3f", mean(hl_pvalues, na.rm = TRUE)),
           size = 4, fontface = "bold") +
  labs(title = "C. Calibration Plot",
       x = "Predicted Probability", y = "Observed Probability") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text = element_text(size = 12, color = "black"))

# Panel D: Confusion Matrix
conf_df <- as.data.frame(confusion_avg)
colnames(conf_df) <- c("Reference", "Prediction", "Freq")
total <- sum(conf_df$Freq)
conf_df$Label <- sprintf("n = %.1f\n(%.1f%%)", conf_df$Freq,
                         conf_df$Freq / total * 100)

panel_d <- ggplot(conf_df, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(color = "white", linewidth = 2) +
  geom_text(aes(label = Label), color = "white", size = 5,
            fontface = "bold", lineheight = 0.85) +
  scale_fill_gradient(low = "#3498DB", high = "#E74C3C") +
  labs(title = "D. Confusion Matrix", x = "Actual", y = "Predicted") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text = element_text(size = 12, face = "bold"),
        legend.position = "none", panel.grid = element_blank())

# Combine and save
fig_combined <- (panel_a + panel_b) / (panel_c + panel_d)

ggsave("Fig_S7_LR_Comp7_Combined_4Panel.tif", fig_combined,
       device = "tiff", width = 26, height = 22, units = "cm",
       dpi = 1200, compression = "lzw", bg = "white")

# =================================================================
# SUMMARY
# =================================================================

message("\n========================================")
message("COMP.7 (DAYTIME DYSFUNCTION) COMPLETED")
message("========================================")
message(sprintf("Without SMOTE: AUC = %.3f", mean(auc_no_smote)))
message(sprintf("With SMOTE:    AUC = %.3f", mean(auc_values)))
message(sprintf("Wilcoxon p = %.4f", wt_auc$p.value))
message(sprintf("Overweight OR (No SMOTE) = %.2f, p = %.4f",
                table_no_smote$estimate[table_no_smote$term == "bmicOverweight"],
                table_no_smote$p.value[table_no_smote$term == "bmicOverweight"]))
message("========================================")
