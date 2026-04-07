############################################################
# 12_HL_Calibration.R
# Hosmer-Lemeshow calibration comparison:
# with SMOTE vs without SMOTE across 100 imputed datasets
#
# Project: Biological and Anthropometric Determinants of
#          Sleep Quality Across the Adult Life Course
# Journal: Journal of Sleep Research
#
# Requires: data_stacked (from 02_Multiple_Imputation.R)
#
# Output: Table_S_HL_Calibration.xlsx
#         Fig_S4_HL_Calibration.tif  (2-panel histogram)
############################################################

library(dplyr)
library(ggplot2)
library(caret)
library(themis)
library(recipes)
library(ResourceSelection)
library(gridExtra)
library(grid)
library(writexl)

n_imp <- max(data_stacked$.imp)
set.seed(123)

# =================================================================
# STEP 1 — COMPUTE H-L P-VALUES WITH AND WITHOUT SMOTE
# =================================================================

hl_p_no_smote <- numeric(n_imp)
hl_p_smote    <- numeric(n_imp)

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

  outcome_i <- as.numeric(test_i$psqi == "Poor_Sleep")

  # --- Without SMOTE ---
  model_no <- glm(psqi ~ sex + age + bmic + bf + hb,
                  data = train_i, family = "binomial")
  pred_no  <- predict(model_no, newdata = test_i, type = "response")

  tryCatch({
    hl_no <- hoslem.test(outcome_i, pred_no, g = 10)
    hl_p_no_smote[i] <- hl_no$p.value
  }, error = function(e) { hl_p_no_smote[i] <<- NA })

  # --- With SMOTE ---
  rec_smote   <- recipe(psqi ~ sex + age + bmic + bf + hb, data = train_i) %>%
    step_smotenc(psqi, over_ratio = 1) %>% prep()
  train_smote <- bake(rec_smote, new_data = NULL)

  model_sm <- glm(psqi ~ sex + age + bmic + bf + hb,
                  data = train_smote, family = "binomial")
  pred_sm  <- predict(model_sm, newdata = test_i, type = "response")

  tryCatch({
    hl_sm <- hoslem.test(outcome_i, pred_sm, g = 10)
    hl_p_smote[i] <- hl_sm$p.value
  }, error = function(e) { hl_p_smote[i] <<- NA })
}

# =================================================================
# STEP 2 — RESULTS
# =================================================================

hl_p_no_smote_clean <- hl_p_no_smote[!is.na(hl_p_no_smote)]
hl_p_smote_clean    <- hl_p_smote[!is.na(hl_p_smote)]

pct_adequate_no <- round(mean(hl_p_no_smote_clean > 0.05) * 100, 1)
pct_adequate_sm <- round(mean(hl_p_smote_clean > 0.05) * 100, 1)

message("\n=== HOSMER-LEMESHOW CALIBRATION ===")
message(sprintf("Without SMOTE: mean p = %.3f (SD = %.3f) | %s%% adequate",
                mean(hl_p_no_smote_clean), sd(hl_p_no_smote_clean),
                pct_adequate_no))
message(sprintf("With SMOTE:    mean p = %.3f (SD = %.3f) | %s%% adequate",
                mean(hl_p_smote_clean), sd(hl_p_smote_clean),
                pct_adequate_sm))

# =================================================================
# STEP 3 — EXPORT TABLE
# =================================================================

hl_table <- data.frame(
  Condition    = c("Without SMOTE", "With SMOTE"),
  Mean_p       = c(round(mean(hl_p_no_smote_clean), 3),
                   round(mean(hl_p_smote_clean), 3)),
  SD_p         = c(round(sd(hl_p_no_smote_clean), 3),
                   round(sd(hl_p_smote_clean), 3)),
  Median_p     = c(round(median(hl_p_no_smote_clean), 3),
                   round(median(hl_p_smote_clean), 3)),
  Pct_adequate = c(pct_adequate_no, pct_adequate_sm),
  N_valid      = c(length(hl_p_no_smote_clean), length(hl_p_smote_clean))
)

write_xlsx(hl_table, "Table_S_HL_Calibration.xlsx")

# =================================================================
# STEP 4 — FIGURE S4: 2-PANEL HISTOGRAM
# =================================================================

make_hl_panel <- function(p_values, title_text, fill_color, pct_ok) {

  df <- data.frame(p_value = p_values)

  ggplot(df, aes(x = p_value)) +
    geom_histogram(bins = 25, fill = fill_color, color = "white", alpha = 0.8) +
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "red",
               linewidth = 1) +
    annotate("text", x = 0.07, y = Inf, vjust = 2, hjust = 0,
             label = "p = 0.05", color = "red", size = 6.5,
             fontface = "italic") +
    annotate("label", x = 0.75, y = Inf, vjust = 1.5,
             label = sprintf("Mean p = %.3f\n%s%% adequate\n(p > 0.05)",
                             mean(p_values), pct_ok),
             size = 4.5, fontface = "bold", fill = "white",
             label.size = 0.5) +
    labs(title = title_text,
         x = "Hosmer-Lemeshow p-value", y = "Frequency") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    theme_minimal(base_size = 20) +
    theme(
      plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
      axis.text = element_text(size = 18, color = "black"),
      axis.title = element_text(face = "bold", size = 18),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = "black",
                                      linewidth = 0.8),
      plot.margin = ggplot2::margin(10, 15, 10, 15)
    )
}

p_no_smote <- make_hl_panel(hl_p_no_smote_clean,
                            "A. Without SMOTE", "#3498DB",
                            as.character(pct_adequate_no))

p_smote <- make_hl_panel(hl_p_smote_clean,
                         "B. With SMOTE", "#FF4500",
                         as.character(pct_adequate_sm))

titulo <- textGrob("Hosmer-Lemeshow Calibration: Effect of SMOTE",
                   gp = gpar(fontsize = 24, fontface = "bold"))

subtitulo <- textGrob(
  sprintf("Logistic Regression | n = %d imputed datasets | Dashed line = significance threshold",
          n_imp),
  gp = gpar(fontsize = 20, col = "gray40")
)

fig_s4 <- arrangeGrob(
  titulo, subtitulo,
  arrangeGrob(p_no_smote, p_smote, ncol = 2),
  nrow = 3, heights = c(0.8, 0.4, 5)
)

ggsave("Fig_S4_HL_Calibration.tif", fig_s4,
       device = "tiff", width = 28, height = 16, units = "cm",
       dpi = 1200, compression = "lzw", bg = "white")

# =================================================================
# SUMMARY
# =================================================================

message("\n========================================")
message("H-L CALIBRATION COMPARISON COMPLETED")
message("========================================")
message(sprintf("Without SMOTE: %s%% adequate (mean p = %.3f)",
                pct_adequate_no, mean(hl_p_no_smote_clean)))
message(sprintf("With SMOTE:    %s%% adequate (mean p = %.3f)",
                pct_adequate_sm, mean(hl_p_smote_clean)))
message("\nFiles saved:")
message("  Table_S_HL_Calibration.xlsx")
message("  Fig_S4_HL_Calibration.tif")
message("========================================")
