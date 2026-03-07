############################################################
# 11_SMOTE_Justification.R
# Formal evaluation of SMOTE's impact on AUC, sensitivity,
# and specificity using paired Wilcoxon tests
#
# Project: Biological and Anthropometric Determinants of
#          Sleep Quality Across the Adult Life Course
# Journal: Sleep Medicine
#
# Requires: data_stacked (from 02_Multiple_Imputation.R)
#
# Output: Table_S_SMOTE_Justification.xlsx
#         Fig_S3_SMOTE_Justification.tif  (3-panel: 2+1 layout)
############################################################

library(dplyr)
library(ggplot2)
library(pROC)
library(caret)
library(themis)
library(recipes)
library(gridExtra)
library(grid)
library(writexl)

n_imp <- max(data_stacked$.imp)
set.seed(123)

# =================================================================
# STEP 1 — COMPUTE METRICS WITH AND WITHOUT SMOTE
# =================================================================

auc_smote     <- numeric(n_imp)
auc_no_smote  <- numeric(n_imp)
sens_smote    <- numeric(n_imp)
sens_no_smote <- numeric(n_imp)
spec_smote    <- numeric(n_imp)
spec_no_smote <- numeric(n_imp)

for (i in 1:n_imp) {

  if (i %% 20 == 0) message(sprintf("Processing dataset %d of %d", i, n_imp))

  capa_i <- data_stacked %>% filter(.imp == i) %>%
    dplyr::select(-.imp, -.id) %>%
    mutate(across(where(is.ordered),
                  ~factor(as.character(.x), ordered = FALSE))) %>%
    as.data.frame()

  idx_train <- createDataPartition(capa_i$psqi, p = 0.8, list = FALSE)
  train_i <- capa_i[idx_train, ]
  test_i  <- capa_i[-idx_train, ]

  # --- Model WITHOUT SMOTE ---
  model_no <- glm(psqi ~ sex + age + bmic + bf + hb,
                  data = train_i, family = "binomial")
  pred_no  <- predict(model_no, newdata = test_i, type = "response")
  roc_no   <- roc(test_i$psqi, pred_no, quiet = TRUE,
                  levels = c("Good_Sleep", "Poor_Sleep"), direction = "<")
  auc_no_smote[i] <- as.numeric(auc(roc_no))

  coords_no     <- coords(roc_no, "best", best.method = "youden")
  pred_class_no <- factor(ifelse(pred_no >= coords_no$threshold,
                                 "Poor_Sleep", "Good_Sleep"),
                          levels = c("Good_Sleep", "Poor_Sleep"))
  cm_no <- confusionMatrix(pred_class_no, test_i$psqi, positive = "Poor_Sleep")
  sens_no_smote[i] <- cm_no$byClass["Sensitivity"]
  spec_no_smote[i] <- cm_no$byClass["Specificity"]

  # --- Model WITH SMOTE ---
  rec_smote   <- recipe(psqi ~ sex + age + bmic + bf + hb, data = train_i) %>%
    step_smotenc(psqi, over_ratio = 1) %>% prep()
  train_smote <- bake(rec_smote, new_data = NULL)

  model_sm <- glm(psqi ~ sex + age + bmic + bf + hb,
                  data = train_smote, family = "binomial")
  pred_sm  <- predict(model_sm, newdata = test_i, type = "response")
  roc_sm   <- roc(test_i$psqi, pred_sm, quiet = TRUE,
                  levels = c("Good_Sleep", "Poor_Sleep"), direction = "<")
  auc_smote[i] <- as.numeric(auc(roc_sm))

  coords_sm     <- coords(roc_sm, "best", best.method = "youden")
  pred_class_sm <- factor(ifelse(pred_sm >= coords_sm$threshold,
                                 "Poor_Sleep", "Good_Sleep"),
                          levels = c("Good_Sleep", "Poor_Sleep"))
  cm_sm <- confusionMatrix(pred_class_sm, test_i$psqi, positive = "Poor_Sleep")
  sens_smote[i] <- cm_sm$byClass["Sensitivity"]
  spec_smote[i] <- cm_sm$byClass["Specificity"]
}

# =================================================================
# STEP 2 — STATISTICAL TESTS
# =================================================================

wilcox_auc  <- wilcox.test(auc_smote, auc_no_smote, paired = TRUE)
wilcox_sens <- wilcox.test(sens_smote, sens_no_smote, paired = TRUE)
wilcox_spec <- wilcox.test(spec_smote, spec_no_smote, paired = TRUE)

message("\n=== SMOTE JUSTIFICATION ===")
message(sprintf("AUC:         Without = %.3f (±%.3f) | With = %.3f (±%.3f) | p = %.4f",
                mean(auc_no_smote), sd(auc_no_smote),
                mean(auc_smote), sd(auc_smote), wilcox_auc$p.value))
message(sprintf("Sensitivity: Without = %.3f (±%.3f) | With = %.3f (±%.3f) | p = %.4f",
                mean(sens_no_smote), sd(sens_no_smote),
                mean(sens_smote), sd(sens_smote), wilcox_sens$p.value))
message(sprintf("Specificity: Without = %.3f (±%.3f) | With = %.3f (±%.3f) | p = %.4f",
                mean(spec_no_smote), sd(spec_no_smote),
                mean(spec_smote), sd(spec_smote), wilcox_spec$p.value))

# =================================================================
# STEP 3 — EXPORT TABLE
# =================================================================

smote_table <- data.frame(
  Metric = c("AUC", "Sensitivity", "Specificity"),
  Without_SMOTE_Mean = round(c(mean(auc_no_smote), mean(sens_no_smote),
                               mean(spec_no_smote)), 3),
  Without_SMOTE_SD   = round(c(sd(auc_no_smote), sd(sens_no_smote),
                               sd(spec_no_smote)), 3),
  With_SMOTE_Mean    = round(c(mean(auc_smote), mean(sens_smote),
                               mean(spec_smote)), 3),
  With_SMOTE_SD      = round(c(sd(auc_smote), sd(sens_smote),
                               sd(spec_smote)), 3),
  Difference         = round(c(mean(auc_smote) - mean(auc_no_smote),
                               mean(sens_smote) - mean(sens_no_smote),
                               mean(spec_smote) - mean(spec_no_smote)), 3),
  Wilcoxon_p         = round(c(wilcox_auc$p.value, wilcox_sens$p.value,
                               wilcox_spec$p.value), 4)
)

write_xlsx(smote_table, "Table_S_SMOTE_Justification.xlsx")

# =================================================================
# STEP 4 — FIGURE S3: 3-PANEL (2 TOP + 1 BOTTOM CENTERED)
# =================================================================

make_panel <- function(values_no, values_yes, metric_name,
                       fill_no, fill_yes, p_val) {

  df <- data.frame(
    Value  = c(values_no, values_yes),
    Method = factor(rep(c("Without\nSMOTE", "With\nSMOTE"),
                        each = length(values_no)),
                    levels = c("Without\nSMOTE", "With\nSMOTE"))
  )

  ast <- ifelse(p_val < 0.001, "***",
         ifelse(p_val < 0.01, "**",
         ifelse(p_val < 0.05, "*", "ns")))

  y_max  <- max(df$Value, na.rm = TRUE)
  y_min  <- min(df$Value, na.rm = TRUE)
  rango  <- y_max - y_min
  y_cap  <- y_max + rango * 0.08
  y_bar  <- y_max + rango * 0.14
  y_text <- y_bar + rango * 0.08

  ggplot(df, aes(x = Method, y = Value, fill = Method)) +
    geom_boxplot(alpha = 0.75, outlier.shape = NA, width = 0.55,
                 color = "black", linewidth = 0.5) +
    geom_jitter(width = 0.12, alpha = 0.25, size = 1.5, color = "gray25") +
    scale_fill_manual(values = c("Without\nSMOTE" = fill_no,
                                 "With\nSMOTE" = fill_yes)) +
    # Capped bracket
    annotate("segment", x = 1, xend = 1, y = y_cap, yend = y_bar,
             linewidth = 0.6) +
    annotate("segment", x = 2, xend = 2, y = y_cap, yend = y_bar,
             linewidth = 0.6) +
    annotate("segment", x = 1, xend = 2, y = y_bar, yend = y_bar,
             linewidth = 0.6) +
    annotate("text", x = 1.5, y = y_text, label = ast,
             size = 8, fontface = "bold") +
    # Mean ± SD annotations
    annotate("text", x = 1, y = y_min - rango * 0.12,
             label = sprintf("%.3f ± %.3f", mean(values_no), sd(values_no)),
             size = 6.5, color = fill_no, fontface = "bold") +
    annotate("text", x = 2, y = y_min - rango * 0.12,
             label = sprintf("%.3f ± %.3f", mean(values_yes), sd(values_yes)),
             size = 6.5, color = fill_yes, fontface = "bold") +
    labs(title = metric_name, x = "", y = metric_name) +
    theme_minimal(base_size = 20) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      axis.text.x = element_text(size = 18, face = "bold", color = "black"),
      axis.text.y = element_text(size = 18, color = "black"),
      axis.title.y = element_text(face = "bold", size = 18),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = "black",
                                      linewidth = 0.8),
      legend.position = "none",
      plot.margin = ggplot2::margin(10, 15, 10, 15)
    ) +
    coord_cartesian(ylim = c(y_min - rango * 0.18, y_text + rango * 0.08))
}

# Create panels
p_auc  <- make_panel(auc_no_smote, auc_smote,
                     "A. AUC", "#1E90FF", "#FF4500", wilcox_auc$p.value)
p_sens <- make_panel(sens_no_smote, sens_smote,
                     "B. Sensitivity", "#1E90FF", "#FF4500", wilcox_sens$p.value)
p_spec <- make_panel(spec_no_smote, spec_smote,
                     "C. Specificity", "#1E90FF", "#FF4500", wilcox_spec$p.value)

# Layout: 2 top + 1 bottom centered
titulo <- textGrob("Effect of SMOTE on Classification Performance",
                   gp = gpar(fontsize = 24, fontface = "bold"))
subtitulo <- textGrob(
  sprintf("Logistic Regression | Wilcoxon paired tests | n = %d imputed datasets",
          n_imp),
  gp = gpar(fontsize = 20, col = "#556B2F")
)

row_top    <- arrangeGrob(p_auc, p_sens, ncol = 2)
row_bottom <- arrangeGrob(nullGrob(), p_spec, nullGrob(),
                          ncol = 3, widths = c(1, 2, 1))

fig_final <- arrangeGrob(titulo, subtitulo, row_top, row_bottom,
                         nrow = 4, heights = c(0.8, 0.4, 5, 5))

ggsave("Fig_S3_SMOTE_Justification.tif", fig_final,
       device = "tiff", width = 28, height = 28, units = "cm",
       dpi = 1200, compression = "lzw", bg = "white")

# =================================================================
# SUMMARY
# =================================================================

message("\n========================================")
message("SMOTE JUSTIFICATION COMPLETED")
message("========================================")
message("Files saved:")
message("  Table_S_SMOTE_Justification.xlsx")
message("  Fig_S3_SMOTE_Justification.tif")
message("========================================")

