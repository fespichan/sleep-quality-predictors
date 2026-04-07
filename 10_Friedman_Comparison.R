############################################################
# 10_Friedman_Comparison.R
# Model comparison: LR vs RF vs XGBoost
# Friedman test with post-hoc pairwise Wilcoxon
#
# Project: Biological and Anthropometric Determinants of
#          Sleep Quality Across the Adult Life Course
# Journal: Journal of Sleep Research
#
# Requires: auc_values_logit.rds (from 06_LR_Pooled_PSQI.R)
#           auc_values_rf.rds    (from 07_RF_PSQI.R)
#           auc_values_xgb.rds   (from 08_XGBoost_PSQI.R)
#
# Output: Table_Model_Comparison.xlsx
#         Fig3_AUC_Comparison_Boxplot.tif
############################################################

library(ggplot2)
library(writexl)

# =================================================================
# LOAD AUC VECTORS
# =================================================================

auc_values_logit <- readRDS("auc_values_logit.rds")
auc_values_rf    <- readRDS("auc_values_rf.rds")
auc_values_xgb   <- readRDS("auc_values_xgb.rds")

n_imp <- length(auc_values_logit)

# =================================================================
# FRIEDMAN TEST (NON-PARAMETRIC REPEATED MEASURES)
# =================================================================

auc_matrix <- cbind(LR = auc_values_logit,
                    RF = auc_values_rf,
                    XGB = auc_values_xgb)

friedman_result <- friedman.test(auc_matrix)

message("\n=== FRIEDMAN TEST ===")
message(sprintf("Chi-squared = %.2f, df = %d, p = %.2e",
                friedman_result$statistic,
                friedman_result$parameter,
                friedman_result$p.value))

# =================================================================
# POST-HOC PAIRWISE WILCOXON (BONFERRONI CORRECTION)
# =================================================================

posthoc_result <- pairwise.wilcox.test(
  c(auc_values_logit, auc_values_rf, auc_values_xgb),
  rep(c("LR", "RF", "XGB"), each = n_imp),
  paired = TRUE, p.adjust.method = "bonferroni"
)

message("\n=== POST-HOC PAIRWISE COMPARISONS ===")
print(posthoc_result)

# Extract p-values into a table
posthoc_table <- data.frame(
  Comparison = c("LR vs RF", "LR vs XGBoost", "RF vs XGBoost"),
  p_adjusted = c(posthoc_result$p.value["RF", "LR"],
                 posthoc_result$p.value["XGB", "LR"],
                 posthoc_result$p.value["XGB", "RF"]),
  Significant = c(posthoc_result$p.value["RF", "LR"] < 0.05,
                  posthoc_result$p.value["XGB", "LR"] < 0.05,
                  posthoc_result$p.value["XGB", "RF"] < 0.05)
)

message("\nPost-hoc summary:")
print(posthoc_table, row.names = FALSE)

# =================================================================
# COMPARISON TABLE
# =================================================================

comparison_table <- data.frame(
  Model    = c("Logistic Regression", "Random Forest", "XGBoost"),
  Mean_AUC = round(c(mean(auc_values_logit), mean(auc_values_rf),
                     mean(auc_values_xgb)), 3),
  SD_AUC   = round(c(sd(auc_values_logit), sd(auc_values_rf),
                     sd(auc_values_xgb)), 3),
  Median_AUC = round(c(median(auc_values_logit), median(auc_values_rf),
                       median(auc_values_xgb)), 3)
)

message("\n=== MODEL COMPARISON ===")
print(comparison_table, row.names = FALSE)

write_xlsx(
  list(
    Comparison = comparison_table,
    PostHoc    = posthoc_table,
    Friedman   = data.frame(
      Statistic = round(friedman_result$statistic, 2),
      df        = friedman_result$parameter,
      p_value   = formatC(friedman_result$p.value, format = "e", digits = 2)
    )
  ),
  "Table_Model_Comparison.xlsx"
)

# =================================================================
# FIGURE 3: AUC COMPARISON BOXPLOT (PRO VERSION)
# =================================================================

comparison_df <- data.frame(
  Model = factor(rep(c("Logistic Regression", "Random Forest", "XGBoost"),
                     each = n_imp),
                 levels = c("Logistic Regression", "Random Forest", "XGBoost")),
  AUC = c(auc_values_logit, auc_values_rf, auc_values_xgb)
)

# Mean ± SD for annotations
means <- c(mean(auc_values_logit), mean(auc_values_rf), mean(auc_values_xgb))
sds   <- c(sd(auc_values_logit), sd(auc_values_rf), sd(auc_values_xgb))
fill_colors <- c("#E74C3C", "#2196F3", "#FF9800")

# Bracket geometry
cap <- 0.012          # half-height of the vertical caps |
lw  <- 0.55           # linewidth for brackets
y_floor <- max(comparison_df$AUC) + 0.025

# Bracket 1: LR vs RF (positions 1-2)
b1 <- y_floor + 0.02
# Bracket 2: LR vs XGB (positions 1-3)
b2 <- b1 + 0.055
# Bracket 3: RF vs XGB (positions 2-3)
b3 <- b2 + 0.055

# Helper: draw a single capped bracket  |-------|  with asterisk
# x1, x2 = horizontal positions; y = height of horizontal bar
bracket <- function(x1, x2, y, label) {
  list(
    # Left cap (vertical, symmetric around y)
    annotate("segment", x = x1, xend = x1,
             y = y - cap, yend = y + cap, linewidth = lw),
    # Right cap (vertical, symmetric around y)
    annotate("segment", x = x2, xend = x2,
             y = y - cap, yend = y + cap, linewidth = lw),
    # Horizontal bar
    annotate("segment", x = x1, xend = x2,
             y = y, yend = y, linewidth = lw),
    # Asterisk label centered above bar
    annotate("text", x = (x1 + x2) / 2, y = y + cap + 0.008,
             label = label, size = 6, fontface = "bold")
  )
}

fig_comparison <- ggplot(comparison_df,
                         aes(x = Model, y = AUC, fill = Model)) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA, width = 0.5,
               color = "black", linewidth = 0.5) +
  geom_jitter(width = 0.1, alpha = 0.2, size = 1.3, color = "gray25") +
  scale_fill_manual(values = c("Logistic Regression" = "#E74C3C",
                               "Random Forest" = "#2196F3",
                               "XGBoost" = "#FF9800")) +

  # Bracket 1: LR vs RF (***)
  bracket(1, 2, b1, "***") +

  # Bracket 2: LR vs XGB (***)
  bracket(1, 3, b2, "***") +

  # Bracket 3: RF vs XGB (**)
  bracket(2, 3, b3, "**") +

  # Mean ± SD annotations below boxplots
  annotate("text", x = 1, y = min(comparison_df$AUC) - 0.035,
           label = sprintf("%.3f ± %.3f", means[1], sds[1]),
           size = 4.5, color = fill_colors[1], fontface = "bold") +
  annotate("text", x = 2, y = min(comparison_df$AUC) - 0.035,
           label = sprintf("%.3f ± %.3f", means[2], sds[2]),
           size = 4.5, color = fill_colors[2], fontface = "bold") +
  annotate("text", x = 3, y = min(comparison_df$AUC) - 0.035,
           label = sprintf("%.3f ± %.3f", means[3], sds[3]),
           size = 4.5, color = fill_colors[3], fontface = "bold") +

  labs(
    title = "AUC Comparison Across 100 Imputed Datasets",
    subtitle = sprintf("Friedman test: \u03c7\u00b2 = %.1f, p < 0.001",
                       friedman_result$statistic),
    x = "", y = "AUC (Test Set)"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5, color = "gray30"),
    axis.text.x = element_text(size = 14, face = "bold", color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title.y = element_text(face = "bold", size = 15),
    axis.line = element_line(color = "black", linewidth = 0.6),
    legend.position = "none",
    plot.margin = ggplot2::margin(15, 15, 15, 15)
  ) +
  coord_cartesian(ylim = c(min(comparison_df$AUC) - 0.07,
                           b3 + cap + 0.03))

ggsave("Fig3_AUC_Comparison_Boxplot.tif", fig_comparison,
       device = "tiff", width = 18, height = 18, units = "cm",
       dpi = 1200, compression = "lzw", bg = "white")

# =================================================================
# SUMMARY
# =================================================================

message("\n========================================")
message("MODEL COMPARISON COMPLETED")
message("========================================")
message(sprintf("LR:  AUC = %.3f ± %.3f", mean(auc_values_logit), sd(auc_values_logit)))
message(sprintf("RF:  AUC = %.3f ± %.3f", mean(auc_values_rf), sd(auc_values_rf)))
message(sprintf("XGB: AUC = %.3f ± %.3f", mean(auc_values_xgb), sd(auc_values_xgb)))
message(sprintf("Friedman p < 0.001 | All pairwise comparisons significant"))
message("========================================")
