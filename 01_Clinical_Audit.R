############################################################
# 01_Clinical_Audit.R
# Three-stage physiological data validation prior to
# multiple imputation
#
# Project: Biological and Anthropometric Determinants of
#          Sleep Quality Across the Adult Life Course
# Journal: Sleep Medicine
#
# References:
#   - Gallagher et al. (2000), Am J Clin Nutr
#   - Janssen et al. (2000), J Appl Physiol
#   - Watson et al. (1980), Am J Clin Nutr
#   - Pace & Rathbun (1945), J Biol Chem
#   - WHO (2011), Haemoglobin thresholds
#
# Input:  sleep_raw_original.csv
# Output: sleep_data_clean_for_imputation.csv
#         sleep_data_clean_for_imputation.xlsx
#         Audit_Report_Clinical_Flags.xlsx
############################################################

library(dplyr)
library(writexl)

# =================================================================
# WORKFLOW:
# Raw data → absolute ranges → sex-specific coherence →
# relational coherence → audit report → export
# =================================================================

# =================================================================
# STAGE 1 — ABSOLUTE RANGE VALIDATION
# =================================================================
# These thresholds are sex-independent. Any value outside these
# ranges is biologically implausible in living adults.

datos <- read.csv("sleep_raw_original.csv")

n_original <- nrow(datos)
message(sprintf("Observations loaded: %d", n_original))

# Compute BMI
datos <- datos %>%
  mutate(bmi = weight / ((height / 100)^2))

# Preserve original values for the audit report
datos <- datos %>%
  mutate(
    bf_original    = bf,
    mm_original    = mm,
    water_original = water,
    hb_original    = hb,
    bmi_original   = bmi
  )

# Apply absolute plausibility ranges
datos <- datos %>%
  mutate(
    hb    = ifelse(hb < 7  | hb > 22,  NA, hb),
    bf    = ifelse(bf < 3  | bf > 60,  NA, bf),
    mm    = ifelse(mm < 10 | mm > 70,  NA, mm),
    water = ifelse(water < 30 | water > 75, NA, water),
    bmi   = ifelse(bmi < 12 | bmi > 55, NA, bmi)
  )

# Stage 1 report
na_stage1 <- sapply(datos[, c("hb", "bf", "mm", "water", "bmi")],
                    function(x) sum(is.na(x)))
message("\n--- STAGE 1: Absolute range validation ---")
message("Missing values after absolute screening:")
print(na_stage1)

# =================================================================
# STAGE 2 — SEX-SPECIFIC PHYSIOLOGICAL COHERENCE
# =================================================================
# Reference ranges by sex derived from:
#   Body fat:    Gallagher et al. (2000), Am J Clin Nutr, 72(3):694-701
#   Muscle mass: Janssen et al. (2000), J Appl Physiol, 89(2):465-471
#   Body water:  Watson et al. (1980), Am J Clin Nutr, 33(1):27-39
#   Hemoglobin:  WHO (2011), Haemoglobin thresholds for anaemia

datos <- datos %>%
  mutate(
    # Body fat (%)
    # Female: 10% (essential fat) to 50% (extreme obesity)
    # Male:    5% (essential fat) to 45%
    flag_bf = case_when(
      sex == "female" & !is.na(bf) & (bf < 10 | bf > 50)  ~ 1,
      sex == "male"   & !is.na(bf) & (bf < 5  | bf > 45)  ~ 1,
      TRUE ~ 0
    ),

    # Skeletal muscle mass (kg)
    # Female: 15–50 kg
    # Male:   25–65 kg
    flag_mm = case_when(
      sex == "female" & !is.na(mm) & (mm < 15 | mm > 50)  ~ 1,
      sex == "male"   & !is.na(mm) & (mm < 25 | mm > 65)  ~ 1,
      TRUE ~ 0
    ),

    # Total body water (%)
    # Female: 35–60%
    # Male:   43–65%
    flag_water = case_when(
      sex == "female" & !is.na(water) & (water < 35 | water > 60) ~ 1,
      sex == "male"   & !is.na(water) & (water < 43 | water > 65) ~ 1,
      TRUE ~ 0
    ),

    # Hemoglobin (g/dL)
    # Female:  9.5–17.5 g/dL (includes mild anaemia)
    # Male:   10.5–19.0 g/dL
    flag_hb = case_when(
      sex == "female" & !is.na(hb) & (hb < 9.5  | hb > 17.5)  ~ 1,
      sex == "male"   & !is.na(hb) & (hb < 10.5 | hb > 19.0)  ~ 1,
      TRUE ~ 0
    )
  )

# Set flagged values to missing (variable-level, not case-level)
datos <- datos %>%
  mutate(
    bf    = ifelse(flag_bf == 1, NA, bf),
    mm    = ifelse(flag_mm == 1, NA, mm),
    water = ifelse(flag_water == 1, NA, water),
    hb    = ifelse(flag_hb == 1, NA, hb)
  )

# Stage 2 report
na_stage2 <- sapply(datos[, c("hb", "bf", "mm", "water")],
                    function(x) sum(is.na(x)))
message("\n--- STAGE 2: Sex-specific coherence ---")
message("Cumulative missing values (absolute + sex-specific):")
print(na_stage2)

flags_summary <- datos %>%
  summarise(
    flag_bf_total    = sum(flag_bf),
    flag_mm_total    = sum(flag_mm),
    flag_water_total = sum(flag_water),
    flag_hb_total    = sum(flag_hb)
  )
message("\nFlagged cases by variable:")
print(as.data.frame(flags_summary))

# =================================================================
# STAGE 3 — RELATIONAL COHERENCE (BETWEEN VARIABLES)
# =================================================================
# Detects inconsistencies in the relationship between variables,
# not in individual values.

datos <- datos %>%
  mutate(
    # 3.1 Expected total body water (%)
    # Based on: TBW ~ FFM x 0.73 (Pace & Rathbun, 1945)
    # FFM = (100 - BF%) / 100 x weight
    expected_water = case_when(
      sex == "female" & !is.na(bf) ~ 70 - (bf * 0.70),
      sex == "male"   & !is.na(bf) ~ 75 - (bf * 0.60),
      TRUE ~ NA_real_
    ),

    # 3.2 Relational flags
    # Flag: body fat + muscle mass > 90% of total body weight
    # (impossible — bone, organs, and other tissues account for ~15-20%)
    flag_rel_bf_mm = case_when(
      !is.na(bf) & !is.na(mm) & !is.na(weight) &
        ((bf / 100 * weight) + mm) / weight * 100 > 90 ~ 1,
      TRUE ~ 0
    ),

    # Flag: observed water deviates > 12 pp from expected
    flag_rel_water = case_when(
      !is.na(water) & !is.na(expected_water) &
        abs(water - expected_water) > 12 ~ 1,
      TRUE ~ 0
    ),

    # Flag: BMI–body fat discordance
    # Obese by BMI but athlete-level body fat (probable BIA error)
    # Underweight by BMI but obese-level body fat
    flag_rel_bmi_bf = case_when(
      !is.na(bmi) & !is.na(bf) & bmi > 30 &
        ((sex == "female" & bf < 18) | (sex == "male" & bf < 12)) ~ 1,
      !is.na(bmi) & !is.na(bf) & bmi < 20 & bf > 35 ~ 1,
      TRUE ~ 0
    ),

    # Combined relational flag
    flag_relational = ifelse(flag_rel_bf_mm == 1 |
                               flag_rel_water == 1 |
                               flag_rel_bmi_bf == 1, 1, 0)
  )

# Set to missing: only the most likely erroneous variable
datos <- datos %>%
  mutate(
    water = ifelse(flag_rel_water == 1, NA, water),
    bf    = ifelse(flag_rel_bmi_bf == 1, NA, bf),
    bf    = ifelse(flag_rel_bf_mm == 1, NA, bf),
    mm    = ifelse(flag_rel_bf_mm == 1, NA, mm)
  )

# Stage 3 report
na_stage3 <- sapply(datos[, c("hb", "bf", "mm", "water", "bmi")],
                    function(x) sum(is.na(x)))
message("\n--- STAGE 3: Relational coherence ---")
message("Cumulative missing values (absolute + sex + relational):")
print(na_stage3)

flags_rel_summary <- datos %>%
  summarise(
    flag_bf_mm  = sum(flag_rel_bf_mm),
    flag_water  = sum(flag_rel_water),
    flag_bmi_bf = sum(flag_rel_bmi_bf),
    flag_total  = sum(flag_relational)
  )
message("\nCases with relational inconsistencies:")
print(as.data.frame(flags_rel_summary))

# =================================================================
# AUDIT REPORT
# =================================================================

message("\n========================================")
message("COMPLETE AUDIT REPORT")
message("========================================")

# Global summary of missing values
na_final <- data.frame(
  Variable    = c("hb", "bf", "mm", "water", "bmi"),
  NAs_total   = sapply(datos[, c("hb", "bf", "mm", "water", "bmi")],
                       function(x) sum(is.na(x))),
  Pct_missing = round(sapply(datos[, c("hb", "bf", "mm", "water", "bmi")],
                             function(x) sum(is.na(x)) / length(x) * 100), 1)
)
message("\nFinal missing data summary:")
print(na_final, row.names = FALSE)

# Affected cases (detail)
affected_cases <- datos %>%
  filter(flag_bf == 1 | flag_mm == 1 | flag_water == 1 |
           flag_hb == 1 | flag_relational == 1) %>%
  dplyr::select(ide, sex, age,
                bf_original, bf, flag_bf,
                mm_original, mm, flag_mm,
                water_original, water, flag_water,
                hb_original, hb, flag_hb,
                flag_relational) %>%
  arrange(ide)

message(sprintf("\nTotal affected participants: %d of %d (%.1f%%)",
                nrow(affected_cases), n_original,
                nrow(affected_cases) / n_original * 100))

# Missing values by sex
comparison_sex <- datos %>%
  group_by(sex) %>%
  summarise(
    n        = n(),
    bf_NA    = sum(is.na(bf)),
    mm_NA    = sum(is.na(mm)),
    water_NA = sum(is.na(water)),
    hb_NA    = sum(is.na(hb)),
    .groups  = "drop"
  )
message("\nMissing values by sex:")
print(as.data.frame(comparison_sex))

# =================================================================
# QUANTIFY AUDIT IMPACT
# =================================================================
# Compare pre-existing missing values vs. audit-generated

datos_original <- read.csv("sleep_raw_original.csv")

na_original <- sapply(datos_original[, c("hb", "bf", "mm", "water")],
                      function(x) sum(is.na(x)))

datos_clean <- read.csv("sleep_data_clean_for_imputation.csv")
na_after    <- sapply(datos_clean[, c("hb", "bf", "mm", "water")],
                      function(x) sum(is.na(x)))

message("\n--- AUDIT IMPACT ---")
message(sprintf("Pre-existing missing values:    %d", sum(na_original)))
message(sprintf("Audit-generated missing values: %d", sum(na_after) - sum(na_original)))
message(sprintf("Total missing values:           %d", sum(na_after)))

# =================================================================
# EXPORT
# =================================================================

# Remove auxiliary columns before exporting
datos_export <- datos %>%
  dplyr::select(
    -bf_original, -mm_original, -water_original,
    -hb_original, -bmi_original,
    -flag_bf, -flag_mm, -flag_water, -flag_hb,
    -flag_rel_bf_mm, -flag_rel_water, -flag_rel_bmi_bf,
    -flag_relational, -expected_water
  )

write_xlsx(datos_export, "sleep_data_clean_for_imputation.xlsx")
write.csv(datos_export, "sleep_data_clean_for_imputation.csv", row.names = FALSE)

# Export audit report
write_xlsx(
  list(
    Summary        = na_final,
    Cases_Affected = as.data.frame(affected_cases),
    By_Sex         = as.data.frame(comparison_sex)
  ),
  "Audit_Report_Clinical_Flags.xlsx"
)

message("\n========================================")
message("FILES GENERATED:")
message("  sleep_data_clean_for_imputation.xlsx")
message("  sleep_data_clean_for_imputation.csv")
message("  Audit_Report_Clinical_Flags.xlsx")
message("========================================")

