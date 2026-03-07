############################################################
# 02_Multiple_Imputation.R
# Multiple Imputation by Chained Equations (MICE)
# with random forest algorithm (m=100, maxit=50)
#
# Project: Biological and Anthropometric Determinants of
#          Sleep Quality Across the Adult Life Course
# Journal: Sleep Medicine
#
# Input:  sleep_data_clean_for_imputation.csv
#         (output from 01_Clinical_Audit.R)
#
# Output: data_stacked.xlsx (100 imputed datasets, long format)
#         sleep_data_imp_10.xlsx (single imputation, dataset 10)
#         Fig_S1_Missing_Pattern.tif
#         Fig_S2a_TracePlot_Physiological.tif
#         Fig_S2b_TracePlot_PSQI_Components.tif
#         Fig_S2c_DensityPlot_Physiological.tif
#         Fig_S2d_DensityPlot_PSQI_Components.tif
############################################################

# =================================================================
# LOAD PACKAGES
# =================================================================

library(mice)
library(ggmice)
library(ggplot2)
library(writexl)
library(lattice)
library(dplyr)

# =================================================================
# STEP 1 — LOAD AND EXPLORE DATA
# =================================================================

sleep_raw_data <- read.csv("sleep_data_clean_for_imputation.csv")

str(sleep_raw_data)
summary(sleep_raw_data)

# ----- Missing data pattern (Figure S1) -----

plot_pattern <- plot_pattern(sleep_raw_data) +
  theme_minimal() +
  theme(
    text = element_text(family = "sans"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,
                               colour = "black", size = 10),
    axis.text.y = element_text(colour = "black", size = 10),
    axis.title = element_text(face = "bold", size = 11),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("Fig_S1_Missing_Pattern.tif",
       plot = plot_pattern, device = "tiff",
       width = 20, height = 15, units = "cm",
       dpi = 1200, compression = "lzw")

# =================================================================
# STEP 2 — CONFIGURE MICE
# =================================================================

# Define variable groups
likert_vars <- c("Comp.1", "Comp.2", "Comp.3",
                 "Comp.4", "Comp.5", "Comp.6", "Comp.7")

physio_vars <- c("hb", "bf", "mm", "water", "height", "weight", "wc", "age")

# Set variable types
sleep_raw_data$sex <- as.factor(sleep_raw_data$sex)
sleep_raw_data[likert_vars] <- lapply(sleep_raw_data[likert_vars], as.integer)

# ----- Initialize MICE to extract default settings -----

init <- mice(sleep_raw_data, maxit = 0)
meth <- init$method
pred <- init$predictorMatrix

# ----- Imputation methods -----
# RF for continuous physiological variables (captures non-linear relationships)
# PMM for PSQI components (preserves discrete 0-3 scale with floor effects)
meth[physio_vars]  <- "rf"
meth[likert_vars]  <- "pmm"
meth["sex"]        <- ""    # Sex is complete; do not impute

# ----- Predictor matrix (block structure) -----
# Two independent blocks to avoid cross-contamination:
#   Block 1: Physiological variables predict each other + sex
#   Block 2: PSQI components predict each other + sex

pred[, ] <- 0

# Sex predicts all variables, but is not predicted by any
pred[, "sex"]  <- 1
pred["sex", ]  <- 0

# Physiological block: predict each other, not PSQI
pred[physio_vars, physio_vars]  <- 1
pred[physio_vars, likert_vars]  <- 0

# PSQI block: predict each other, not physiological
pred[likert_vars, likert_vars]  <- 1
pred[likert_vars, physio_vars]  <- 0

# Diagonal must be zero (no self-prediction)
diag(pred) <- 0

# ----- Pre-imputation verification -----

message("\n--- IMPUTATION CONFIGURATION ---")
message("Methods assigned:")
print(meth[meth != ""])

message("\nPredictor matrix verification:")
message(sprintf("  Physiological vars predict each other: %s",
                all(pred[physio_vars, physio_vars][row(pred[physio_vars, physio_vars]) !=
                    col(pred[physio_vars, physio_vars])] == 1)))
message(sprintf("  PSQI vars predict each other: %s",
                all(pred[likert_vars, likert_vars][row(pred[likert_vars, likert_vars]) !=
                    col(pred[likert_vars, likert_vars])] == 1)))
message(sprintf("  PSQI does NOT predict physiological: %s",
                all(pred[physio_vars, likert_vars] == 0)))
message(sprintf("  Physiological does NOT predict PSQI: %s",
                all(pred[likert_vars, physio_vars] == 0)))
message(sprintf("  Diagonal is zero: %s", all(diag(pred) == 0)))

# Missing values per variable
na_count <- sapply(sleep_raw_data, function(x) sum(is.na(x)))
message("\nMissing values per variable:")
print(na_count[na_count > 0])

# =================================================================
# STEP 3 — RUN IMPUTATION
# =================================================================

set.seed(123456)

imputed_data <- mice(
  sleep_raw_data,
  method          = meth,
  predictorMatrix = pred,
  m               = 100,
  maxit           = 50,
  printFlag       = FALSE
)

message("\nImputation completed: 100 datasets x 50 iterations")

# Summary of imputed variables
missing_summary <- data.frame(
  Variable    = names(na_count[na_count > 0]),
  N_missing   = na_count[na_count > 0],
  Pct_missing = round(100 * na_count[na_count > 0] / nrow(sleep_raw_data), 1)
)
print(missing_summary, row.names = FALSE)

# =================================================================
# STEP 4 — DIAGNOSTIC PLOTS
# =================================================================

# Variable groups for plotting
vars_physio <- names(imputed_data$nmis)[3:9]
vars_psqi   <- names(imputed_data$nmis)[10:16]

# ----- Trace plots: convergence assessment (Figure S2a) -----

tiff("Fig_S2a_TracePlot_Physiological.tif",
     width = 17, height = 20, units = "cm",
     res = 1200, compression = "lzw")
p <- plot(imputed_data, vars_physio, layout = c(2, 7))
print(p)
dev.off()

# ----- Trace plots: PSQI components (Figure S2b) -----

tiff("Fig_S2b_TracePlot_PSQI_Components.tif",
     width = 17, height = 20, units = "cm",
     res = 1200, compression = "lzw")
p <- plot(imputed_data, vars_psqi, layout = c(2, 7))
print(p)
dev.off()

# ----- Density plots: distributional fidelity (Figure S2c) -----

tiff("Fig_S2c_DensityPlot_Physiological.tif",
     width = 20, height = 12, units = "cm",
     res = 1200, compression = "lzw")
densityplot(imputed_data,
            ~ height + weight + wc + bf + mm + water + hb,
            layout = c(4, 2), thicker = 3)
dev.off()

# ----- Density plots: PSQI components (Figure S2d) -----

tiff("Fig_S2d_DensityPlot_PSQI_Components.tif",
     width = 24, height = 12, units = "cm",
     res = 1200, compression = "lzw")
densityplot(imputed_data,
            ~ Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5 + Comp.6 + Comp.7,
            layout = c(4, 2), thicker = 3)
dev.off()

message("\nDiagnostic plots saved (trace + density)")

# =================================================================
# STEP 5 — EXTRACT AND POST-PROCESS IMPUTED DATA
# =================================================================

# Extract all imputations in long format
data_stacked <- complete(imputed_data, action = "long", include = FALSE)

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

# =================================================================
# STEP 6 — COMPUTE DERIVED VARIABLES
# =================================================================
# Applied to ALL 100 imputed datasets to ensure internal consistency.

# 6.1 BMI from imputed weight and height
data_stacked <- data_stacked %>%
  mutate(bmi = weight / (height / 100)^2)

# 6.2 PSQI global score (sum of 7 components)
data_stacked <- data_stacked %>%
  mutate(
    across(Comp.1:Comp.7, ~as.numeric(as.character(.x))),
    psqi.total = Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5 + Comp.6 + Comp.7
  )

# 6.3 PSQI classification: good sleep (<=5) vs. poor sleep (>5)
data_stacked <- data_stacked %>%
  mutate(psqi = case_when(
    psqi.total <= 5 ~ "g.sleep",
    psqi.total > 5  ~ "p.sleep",
    TRUE ~ NA_character_
  ))

# 6.4 BMI category (WHO criteria; underweight merged with Normal)
data_stacked <- data_stacked %>%
  mutate(bmic = case_when(
    bmi < 25.0                ~ "Normal",
    bmi >= 25.0 & bmi < 30.0 ~ "Overweight",
    bmi >= 30.0               ~ "Obesity",
    TRUE ~ NA_character_
  ))

data_stacked$bmic <- factor(data_stacked$bmic,
                            levels = c("Normal", "Overweight", "Obesity"))

# 6.5 Body fat category (sex-specific, ACE/ACSM guidelines)
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

# 6.6 Anaemia classification (WHO thresholds by sex)
data_stacked <- data_stacked %>%
  mutate(hbc = case_when(
    sex == "female" & hb < 12.0  ~ "anemia",
    sex == "female" & hb >= 12.0 ~ "normal",
    sex == "male"   & hb < 13.0  ~ "anemia",
    sex == "male"   & hb >= 13.0 ~ "normal",
    TRUE ~ NA_character_
  ))

# 6.7 Convert categorical variables to factors
data_stacked <- data_stacked %>%
  mutate(
    psqi = factor(psqi, levels = c("g.sleep", "p.sleep")),
    bmic = factor(bmic, levels = c("Normal", "Overweight", "Obesity")),
    bfc  = as.factor(bfc),
    hbc  = as.factor(hbc),
    sex  = as.factor(sex)
  )

# 6.8 Reorder columns for logical structure
data_stacked <- data_stacked %>%
  dplyr::select(
    .imp, .id,
    sex, age, height, weight, bmi, wc, bf, mm, water, hb,
    Comp.1, Comp.2, Comp.3, Comp.4, Comp.5, Comp.6, Comp.7,
    psqi.total, psqi, bmic, bfc, hbc
  )

# =================================================================
# STEP 7 — POST-IMPUTATION VERIFICATION
# =================================================================

message("\n========================================")
message("POST-IMPUTATION VERIFICATION")
message("========================================")

message(sprintf("Total rows: %d", nrow(data_stacked)))
message(sprintf("Imputations: %d", length(unique(data_stacked$.imp))))
message(sprintf("Subjects per imputation: %d",
                nrow(data_stacked) / length(unique(data_stacked$.imp))))

# Check for NAs in derived variables
na_check <- data_stacked %>%
  summarise(across(c(bmi, psqi.total, psqi, bmic, bfc, hbc),
                   ~sum(is.na(.x))))
message("\nMissing values in derived variables:")
print(as.data.frame(na_check))

# PSQI distribution across 100 imputations
psqi_dist <- data_stacked %>%
  group_by(.imp, psqi) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(psqi) %>%
  summarise(mean_n = mean(n), sd_n = sd(n), .groups = "drop")

message("\nPSQI distribution across 100 imputations:")
print(as.data.frame(psqi_dist))

# BMI category distribution across 100 imputations
bmic_dist <- data_stacked %>%
  group_by(.imp, bmic) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(bmic) %>%
  summarise(mean_n = mean(n), sd_n = sd(n), .groups = "drop")

message("\nBMI category distribution across 100 imputations:")
print(as.data.frame(bmic_dist))

# BMI ranges by category (imputation 1)
bmi_check <- data_stacked %>%
  filter(.imp == 1) %>%
  group_by(bmic) %>%
  summarise(min_bmi = min(bmi), max_bmi = max(bmi),
            mean_bmi = mean(bmi), .groups = "drop")

message("\nBMI ranges by category (imputation 1):")
print(as.data.frame(bmi_check))

# =================================================================
# EXPORT
# =================================================================

write_xlsx(data_stacked, "data_stacked.xlsx")

message("\n========================================")
message("FILES GENERATED:")
message("  data_stacked.xlsx          (100 imputed datasets)")
message("  sleep_data_imp_10.xlsx     (single imputation)")
message("  Fig_S1_Missing_Pattern.tif")
message("  Fig_S2a_TracePlot_Physiological.tif")
message("  Fig_S2b_TracePlot_PSQI_Components.tif")
message("  Fig_S2c_DensityPlot_Physiological.tif")
message("  Fig_S2d_DensityPlot_PSQI_Components.tif")
message("========================================")

