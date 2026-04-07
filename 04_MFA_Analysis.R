############################################################
# 04_MFA_Analysis.R
# Multiple Factor Analysis (MFA) for variable selection
# 4 configurations: PSQI × bfc, bmic, hbc, sex
#
# Project: Biological and Anthropometric Determinants of
#          Sleep Quality Across the Adult Life Course
# Journal: Journal of Sleep Research
#
# Requires: data_stacked (from 02_Multiple_Imputation.R)
#
# Output:
#   Individual plots (per configuration × 6 = 24 figures):
#     Fig_MFA_Circle_[bfc/bmic/hbc/sex].tif
#     Fig_MFA_Individuals_[bfc/bmic/hbc/sex].tif
#     Fig_MFA_Contrib_Dim1_[bfc/bmic/hbc/sex].tif
#     Fig_MFA_Contrib_Dim2_[bfc/bmic/hbc/sex].tif
#     Fig_MFA_Screeplot_[bfc/bmic/hbc/sex].tif
#     Fig_MFA_Groups_[bfc/bmic/hbc/sex].tif
#   Combined panels (for manuscript):
#     Fig4_MFA_Circles_Combined.tif
#     Fig5_MFA_Individuals_Combined.tif
#   Tables (per configuration = 4 Excel files):
#     Table_MFA_Results_[bfc/bmic/hbc/sex].xlsx
#
# Note: MFA was performed on the first imputed dataset as
#       an exploratory technique. No standardized pooling
#       framework exists for MFA with multiply imputed data
#       (van Buuren, 2018).
############################################################

library(FactoMineR)
library(factoextra)
library(ggplot2)
library(dplyr)
library(writexl)
library(patchwork)

# =================================================================
# STEP 1 — EXTRACT FIRST IMPUTATION
# =================================================================

first_imp <- data_stacked %>%
  filter(.imp == 1) %>%
  dplyr::select(-.imp, -.id) %>%
  as.data.frame()

first_imp <- first_imp %>%
  mutate(
    psqi = factor(psqi, levels = c("g.sleep", "p.sleep")),
    bmic = factor(bmic, levels = c("Normal", "Overweight", "Obesity")),
    bfc  = as.factor(bfc),
    hbc  = as.factor(hbc),
    sex  = as.factor(sex)
  )

# =================================================================
# STEP 2 — DEFINE VARIABLE GROUPS
# =================================================================

# Physiological variables (Group 2) — 9 quantitative
vars_physio <- c("age", "height", "weight", "bmi", "wc",
                 "bf", "mm", "water", "hb")

# PSQI components (Group 3) — 8 quantitative
vars_psqi <- c("Comp.1", "Comp.2", "Comp.3", "Comp.4",
               "Comp.5", "Comp.6", "Comp.7", "psqi.total")

# Categorical variable pairs (Group 1) — 4 configurations
bbf_combinations <- list(
  "bfc"  = c("psqi", "bfc"),
  "bmic" = c("psqi", "bmic"),
  "hbc"  = c("psqi", "hbc"),
  "sex"  = c("psqi", "sex")
)

panel_titles <- c(
  "bfc"  = "PSQI x Body Fat Category",
  "bmic" = "PSQI x BMI Category",
  "hbc"  = "PSQI x Hemoglobin Category",
  "sex"  = "PSQI x Sex"
)

# Color palettes for individual plots
palette_map <- list(
  "bfc"  = c("#7d3c98", "#f1c40f", "#FC4E07", "#28b463",
             "#3498db", "#76EEC6", "#C0FF3E"),
  "bmic" = c("#2196F3", "#FF9800", "#E74C3C", "#E066FF", "#C0FF3E"),
  "hbc"  = c("#E74C3C", "#2196F3", "#00EE00", "#FFC125"),
  "sex"  = c("#E74C3C", "#2196F3", "#FFD700", "#B3EE3A")
)

# =================================================================
# STEP 3 — HELPER FUNCTIONS
# =================================================================

prepare_mfa_data <- function(data, bbf_vars) {
  # Select columns in order: BBF | Physiological | PSQI
  df <- data %>%
    dplyr::select(all_of(c(bbf_vars, vars_physio, vars_psqi)))

  for (v in bbf_vars) df[[v]] <- as.factor(df[[v]])

  # Supplementary illustrative variables
  numeric_cols <- df %>% dplyr::select(all_of(c(vars_physio, vars_psqi)))
  df$mean <- rowMeans(numeric_cols, na.rm = TRUE)
  df$std  <- apply(numeric_cols, 1, sd, na.rm = TRUE)

  return(df)
}

run_mfa <- function(data, name_label, bbf_vars) {

  message(sprintf("\n======== Running MFA: %s ========", name_label))

  mfa_data <- prepare_mfa_data(data, bbf_vars)

  # Group structure:
  #   Group 1 (BBF):           2 categorical variables  (cols 1-2)
  #   Group 2 (Physiological): 9 quantitative variables (cols 3-11)
  #   Group 3 (Sleep PSQI):    8 quantitative variables (cols 12-19)
  #   Group 4 (Illustrative):  2 supplementary variables (cols 20-21)

  res.mfa <- MFA(
    mfa_data,
    group         = c(2, 9, 8, 2),
    type          = c("n", "s", "s", "s"),
    name.group    = c("Biological & Behavioral", "Physiological",
                      "Sleep PSQI", "Illustrative"),
    num.group.sup = c(4),
    graph         = FALSE
  )

  # Report key results
  eig.val <- get_eigenvalue(res.mfa)
  message("\nEigenvalues (top 5):")
  print(head(eig.val, 5))

  quanti.var <- get_mfa_var(res.mfa, "quanti.var")
  message("\nVariable coordinates (Dim 1-2):")
  print(round(quanti.var$coord[, 1:2], 3))

  return(res.mfa)
}

# =================================================================
# STEP 4 — GENERATE INDIVIDUAL PLOTS
# =================================================================

generate_mfa_plots <- function(res.mfa, name_label, bbf_vars) {

  hab_var <- bbf_vars[2]
  colors_to_use <- palette_map[[hab_var]]

  # Screeplot
  p_scree <- fviz_screeplot(res.mfa) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 13),
      axis.text  = element_text(size = 12, color = "black"),
      panel.background = element_rect(fill = "white", color = "black")
    ) +
    labs(title = paste("Screeplot -", name_label))

  # Correlation circle
  p_circle <- fviz_mfa_var(
    res.mfa, "quanti.var",
    palette = "d3", col.var.sup = "purple",
    repel = TRUE, arrowsize = 1.5,
    labelsize = 5, axes.linetype = "solid"
  ) +
    theme_minimal(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 16),
      axis.text  = element_text(size = 14, color = "black"),
      panel.background = element_rect(fill = "white", color = "black")
    ) +
    labs(title = paste("Correlation Circle -", name_label))

  # Individuals with ellipses
  p_individuals <- fviz_ellipses(
    res.mfa, 1:2,
    geom = "point", repel = FALSE,
    pointsize = 3, labelsize = 1,
    palette = colors_to_use
  ) +
    theme_minimal(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 16),
      axis.text  = element_text(size = 14, color = "black"),
      panel.background = element_rect(fill = "white", color = "black"),
      legend.text  = element_text(size = 14),
      legend.title = element_text(size = 14, face = "bold")
    ) +
    labs(title = paste("Individuals -", name_label))

  # Contributions to Dim 1 and Dim 2
  p_contrib1 <- fviz_contrib(res.mfa, choice = "quanti.var",
                             axes = 1, top = 20, palette = "jco") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.background = element_rect(fill = "white", color = "black")
    ) +
    labs(title = paste("Contribution to Dim 1 -", name_label))

  p_contrib2 <- fviz_contrib(res.mfa, choice = "quanti.var",
                             axes = 2, top = 20, palette = "jco") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.background = element_rect(fill = "white", color = "black")
    ) +
    labs(title = paste("Contribution to Dim 2 -", name_label))

  # Group representation
  p_groups <- fviz_mfa_var(res.mfa, "group",
                           labelsize = 5, repel = FALSE) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      panel.background = element_rect(fill = "white", color = "black")
    ) +
    labs(title = paste("Groups -", name_label))

  return(list(
    screeplot = p_scree, circle = p_circle,
    individuals = p_individuals,
    contrib_dim1 = p_contrib1, contrib_dim2 = p_contrib2,
    groups = p_groups
  ))
}

# =================================================================
# STEP 5 — RUN ALL 4 MFA CONFIGURATIONS
# =================================================================

mfa_results <- list()
mfa_plots   <- list()

for (combo in names(bbf_combinations)) {
  bbf_vars <- bbf_combinations[[combo]]
  label    <- panel_titles[[combo]]

  mfa_results[[combo]] <- run_mfa(first_imp, label, bbf_vars)
  mfa_plots[[combo]]   <- generate_mfa_plots(mfa_results[[combo]],
                                              label, bbf_vars)
}

# =================================================================
# STEP 6 — SAVE INDIVIDUAL FIGURES
# =================================================================

message("\n--- Saving individual figures ---")

for (combo in names(mfa_plots)) {

  label_short <- gsub("psqi_", "", combo)

  ggsave(paste0("Fig_MFA_Circle_", label_short, ".tif"),
         mfa_plots[[combo]]$circle, device = "tiff",
         width = 16, height = 16, units = "cm",
         dpi = 1200, compression = "lzw", bg = "white")

  ggsave(paste0("Fig_MFA_Individuals_", label_short, ".tif"),
         mfa_plots[[combo]]$individuals, device = "tiff",
         width = 16, height = 16, units = "cm",
         dpi = 1200, compression = "lzw", bg = "white")

  ggsave(paste0("Fig_MFA_Contrib_Dim1_", label_short, ".tif"),
         mfa_plots[[combo]]$contrib_dim1, device = "tiff",
         width = 18, height = 12, units = "cm",
         dpi = 1200, compression = "lzw", bg = "white")

  ggsave(paste0("Fig_MFA_Contrib_Dim2_", label_short, ".tif"),
         mfa_plots[[combo]]$contrib_dim2, device = "tiff",
         width = 18, height = 12, units = "cm",
         dpi = 1200, compression = "lzw", bg = "white")

  ggsave(paste0("Fig_MFA_Screeplot_", label_short, ".tif"),
         mfa_plots[[combo]]$screeplot, device = "tiff",
         width = 14, height = 10, units = "cm",
         dpi = 1200, compression = "lzw", bg = "white")

  ggsave(paste0("Fig_MFA_Groups_", label_short, ".tif"),
         mfa_plots[[combo]]$groups, device = "tiff",
         width = 14, height = 14, units = "cm",
         dpi = 1200, compression = "lzw", bg = "white")

  message(sprintf("  Saved: %s (6 plots)", combo))
}

# =================================================================
# STEP 7 — COMBINED 2×2 PANELS (FOR MANUSCRIPT)
# =================================================================

message("\n--- Creating combined panels ---")

# Rebuild circle and individual plots with panel labels (A-D)
panel_labels <- c("A", "B", "C", "D")
circle_panels <- list()
indiv_panels  <- list()

for (idx in seq_along(names(bbf_combinations))) {

  combo    <- names(bbf_combinations)[idx]
  bbf_vars <- bbf_combinations[[combo]]
  mfa_data <- prepare_mfa_data(first_imp, bbf_vars)

  res.mfa <- MFA(
    mfa_data,
    group = c(2, 9, 8, 2),
    type = c("n", "s", "s", "s"),
    name.group = c("Biological & Behavioral", "Physiological",
                   "Sleep PSQI", "Illustrative"),
    num.group.sup = c(4), graph = FALSE
  )

  # Correlation circle panel
  circle_panels[[idx]] <- fviz_mfa_var(
    res.mfa, "quanti.var",
    palette = "d3", col.var.sup = "purple",
    repel = TRUE, arrowsize = 1.5,
    labelsize = 5, axes.linetype = "solid"
  ) +
    theme_minimal(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 15),
      axis.text  = element_text(size = 13, color = "black"),
      panel.background = element_rect(fill = "white", color = "black"),
      legend.position = "none"
    ) +
    labs(title = paste0(panel_labels[idx], ". ", panel_titles[combo]))

  # Individual plot panel
  indiv_panels[[idx]] <- fviz_ellipses(
    res.mfa, 1:2,
    geom = "point", repel = FALSE,
    pointsize = 1.5, labelsize = 1,
    palette = palette_map[[combo]]
  ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 11),
      axis.text  = element_text(size = 9, color = "black"),
      panel.background = element_rect(fill = "white", color = "black"),
      legend.text  = element_text(size = 9),
      legend.title = element_text(size = 9, face = "bold"),
      legend.key.size = unit(0.4, "cm"),
      legend.position = "bottom"
    ) +
    labs(title = paste0(panel_labels[idx], ". ", panel_titles[combo]))
}

# Save combined correlation circles (Figure 4)
fig_circles <- (circle_panels[[1]] + circle_panels[[2]]) /
  (circle_panels[[3]] + circle_panels[[4]])

ggsave("Fig4_MFA_Circles_Combined.tif",
       fig_circles, device = "tiff",
       width = 28, height = 28, units = "cm",
       dpi = 1200, compression = "lzw", bg = "white")

# Save combined individual plots (Figure 5)
fig_individuals <- (indiv_panels[[1]] + indiv_panels[[2]]) /
  (indiv_panels[[3]] + indiv_panels[[4]])

ggsave("Fig5_MFA_Individuals_Combined.tif",
       fig_individuals, device = "tiff",
       width = 28, height = 28, units = "cm",
       dpi = 1200, compression = "lzw", bg = "white")

message("  Saved: Fig4_MFA_Circles_Combined.tif")
message("  Saved: Fig5_MFA_Individuals_Combined.tif")

# =================================================================
# STEP 8 — EXPORT NUMERICAL RESULTS
# =================================================================

message("\n--- Exporting numerical results ---")

for (combo in names(mfa_results)) {

  res <- mfa_results[[combo]]

  eig    <- as.data.frame(get_eigenvalue(res))
  eig$Dimension <- rownames(eig)

  qv     <- get_mfa_var(res, "quanti.var")
  coords <- as.data.frame(round(qv$coord, 4))
  coords$Variable <- rownames(coords)

  cos2 <- as.data.frame(round(qv$cos2, 4))
  cos2$Variable <- rownames(cos2)

  contrib <- as.data.frame(round(qv$contrib, 4))
  contrib$Variable <- rownames(contrib)

  grp       <- get_mfa_var(res, "group")
  grp_coord <- as.data.frame(round(grp$coord, 4))
  grp_coord$Group <- rownames(grp_coord)

  grp_contrib <- as.data.frame(round(grp$contrib, 4))
  grp_contrib$Group <- rownames(grp_contrib)

  write_xlsx(
    list(
      Eigenvalues       = eig,
      Var_Coordinates   = coords,
      Var_Cos2          = cos2,
      Var_Contributions = contrib,
      Group_Coordinates = grp_coord,
      Group_Contributions = grp_contrib
    ),
    path = paste0("Table_MFA_Results_", combo, ".xlsx")
  )

  message(sprintf("  Saved: Table_MFA_Results_%s.xlsx", combo))
}

# =================================================================
# SUMMARY
# =================================================================

message("\n========================================")
message("MFA ANALYSIS COMPLETED")
message("========================================")
message("\nFigures (individual): 24 files (6 per configuration)")
message("Figures (combined):   Fig4_MFA_Circles_Combined.tif")
message("                      Fig5_MFA_Individuals_Combined.tif")
message("Tables:               4 Excel files with eigenvalues,")
message("                      coordinates, cos2, and contributions")
message("========================================")
