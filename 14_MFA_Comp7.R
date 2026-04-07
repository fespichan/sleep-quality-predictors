############################################################
# 14_MFA_Comp7.R
# Multiple Factor Analysis for PSQI Component 7
# (Daytime Dysfunction) — exploratory supplementary analysis
#
# Project: Biological and Anthropometric Determinants of
#          Sleep Quality Across the Adult Life Course
# Journal: Journal of Sleep Research
#
# Requires: data_stacked (from 02_Multiple_Imputation.R)
#
# Output:
#   Figures:
#     Fig_S8_MFA_Comp7_Individuals_bmic.tif (for manuscript)
#     Fig_MFA_Comp7_Circle_[bfc/bmic/hbc/sex].tif
#     Fig_MFA_Comp7_Individuals_[bfc/bmic/hbc/sex].tif
#   Tables:
#     Table_MFA_Comp7_Results_[bfc/bmic/hbc/sex].xlsx
#
# Note: Figure S8 in the manuscript uses only the bmic
#       configuration. All 4 configurations are computed
#       for completeness.
############################################################

library(FactoMineR)
library(factoextra)
library(ggplot2)
library(dplyr)
library(writexl)

# =================================================================
# STEP 1 — PREPARE DATA
# =================================================================

first_imp <- data_stacked %>%
  filter(.imp == 1) %>%
  dplyr::select(-.imp, -.id) %>%
  as.data.frame()

# Create dichotomous Comp.7 outcome
first_imp$Comp7_bin <- factor(
  ifelse(first_imp$Comp.7 == 0, "No_dysfunction", "Has_dysfunction")
)

# =================================================================
# STEP 2 — DEFINE VARIABLE GROUPS
# =================================================================

vars_physio <- c("age", "height", "weight", "bmi", "wc",
                 "bf", "mm", "water", "hb")

vars_psqi <- c("Comp.1", "Comp.2", "Comp.3", "Comp.4",
               "Comp.5", "Comp.6", "Comp.7", "psqi.total")

bbf_combinations <- list(
  "bfc"  = c("Comp7_bin", "bfc"),
  "bmic" = c("Comp7_bin", "bmic"),
  "hbc"  = c("Comp7_bin", "hbc"),
  "sex"  = c("Comp7_bin", "sex")
)

labels <- c(
  "bfc"  = "Comp.7 x Body Fat Category",
  "bmic" = "Comp.7 x BMI Category",
  "hbc"  = "Comp.7 x Hemoglobin Category",
  "sex"  = "Comp.7 x Sex"
)

palette_map <- list(
  "bfc"  = c("#7d3c98", "#f1c40f", "#FC4E07", "#28b463",
             "#3498db", "#00EEEE", "#D02090"),
  "bmic" = c("#2196F3", "#FF9800", "#E74C3C", "#FFD700", "#00FF7F"),
  "hbc"  = c("#E74C3C", "#2196F3", "#D15FEE", "#00EE00"),
  "sex"  = c("#E74C3C", "#2196F3", "#FFD700", "#BCEE68")
)

# =================================================================
# STEP 3 — HELPER FUNCTION
# =================================================================

prepare_mfa_data <- function(data, bbf_vars) {
  df <- data %>%
    dplyr::select(all_of(c(bbf_vars, vars_physio, vars_psqi)))
  for (v in bbf_vars) df[[v]] <- as.factor(df[[v]])
  numeric_cols <- df %>% dplyr::select(all_of(c(vars_physio, vars_psqi)))
  df$mean <- rowMeans(numeric_cols, na.rm = TRUE)
  df$std  <- apply(numeric_cols, 1, sd, na.rm = TRUE)
  return(df)
}

# =================================================================
# STEP 4 — RUN ALL 4 MFA CONFIGURATIONS
# =================================================================

mfa_results <- list()

for (combo in names(bbf_combinations)) {

  bbf_vars <- bbf_combinations[[combo]]
  message(sprintf("\n=== Running MFA: %s ===", labels[combo]))

  mfa_data <- prepare_mfa_data(first_imp, bbf_vars)

  res.mfa <- MFA(
    mfa_data,
    group = c(2, 9, 8, 2),
    type = c("n", "s", "s", "s"),
    name.group = c("Biological & Behavioral", "Physiological",
                   "Sleep PSQI", "Illustrative"),
    num.group.sup = c(4), graph = FALSE
  )

  mfa_results[[combo]] <- res.mfa

  # Eigenvalues
  message("Eigenvalues (top 5):")
  print(head(get_eigenvalue(res.mfa), 5))

  # Save correlation circle
  p_circle <- fviz_mfa_var(
    res.mfa, "quanti.var", palette = "d3",
    col.var.sup = "purple", repel = TRUE,
    arrowsize = 1.5, labelsize = 5
  ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(size = 14, color = "black"),
      panel.background = element_rect(fill = "white", color = "black")
    ) +
    labs(title = paste("Correlation Circle -", labels[combo]))

  ggsave(paste0("Fig_MFA_Comp7_Circle_", combo, ".tif"),
         p_circle, device = "tiff",
         width = 16, height = 16, units = "cm",
         dpi = 1200, compression = "lzw", bg = "white")

  # Save individuals plot
  p_indiv <- fviz_ellipses(
    res.mfa, 1:2, geom = "point", repel = FALSE,
    pointsize = 3, labelsize = 1,
    palette = palette_map[[combo]]
  ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(size = 14, color = "black"),
      panel.background = element_rect(fill = "white", color = "black"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold")
    ) +
    labs(title = paste("Individuals -", labels[combo]))

  ggsave(paste0("Fig_MFA_Comp7_Individuals_", combo, ".tif"),
         p_indiv, device = "tiff",
         width = 16, height = 16, units = "cm",
         dpi = 1200, compression = "lzw", bg = "white")

  # Save Figure S8 (bmic configuration only)
  if (combo == "bmic") {
    ggsave("Fig_S8_MFA_Comp7_Individuals_bmic.tif",
           p_indiv, device = "tiff",
           width = 18, height = 16, units = "cm",
           dpi = 1200, compression = "lzw", bg = "white")
    message("  Saved: Fig_S8_MFA_Comp7_Individuals_bmic.tif (for manuscript)")
  }

  message(sprintf("  Saved: %s (circle + individuals)", combo))
}

# =================================================================
# STEP 5 — EXPORT NUMERICAL RESULTS
# =================================================================

message("\n--- Exporting numerical results ---")

for (combo in names(mfa_results)) {

  res <- mfa_results[[combo]]

  eig <- as.data.frame(get_eigenvalue(res))
  eig$Dimension <- rownames(eig)

  qv <- get_mfa_var(res, "quanti.var")
  coords <- as.data.frame(round(qv$coord, 4))
  coords$Variable <- rownames(coords)

  contrib <- as.data.frame(round(qv$contrib, 4))
  contrib$Variable <- rownames(contrib)

  grp <- get_mfa_var(res, "group")
  grp_coord <- as.data.frame(round(grp$coord, 4))
  grp_coord$Group <- rownames(grp_coord)

  write_xlsx(
    list(
      Eigenvalues = eig,
      Var_Coordinates = coords,
      Var_Contributions = contrib,
      Group_Coordinates = grp_coord
    ),
    path = paste0("Table_MFA_Comp7_Results_", combo, ".xlsx")
  )

  message(sprintf("  Saved: Table_MFA_Comp7_Results_%s.xlsx", combo))
}

# =================================================================
# SUMMARY
# =================================================================

message("\n========================================")
message("MFA COMP.7 ANALYSIS COMPLETED")
message("========================================")
message("Figures: 8 files (circle + individuals x 4 configs)")
message("         + Fig_S8_MFA_Comp7_Individuals_bmic.tif")
message("Tables:  4 Excel files")
message("========================================")
