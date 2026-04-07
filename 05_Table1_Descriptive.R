############################################################
# 05_Table1_Descriptive.R
# Demographic and anthropometric characteristics by sex
# Pre-imputation data (N = 289)
#
# Project: Biological and Anthropometric Determinants of
#          Sleep Quality Across the Adult Life Course
# Journal: Journal of Sleep Research
#
# Input:  sleep_data_clean_for_imputation2.csv
#         (output from 01_Clinical_Audit.R)
#
# Output: Table1_BySex.docx
#
# Note: Table 1 is stratified by sex (0% missing) rather
#       than PSQI (61% missing) to maximize sample size
#       and avoid selection bias.
############################################################

library(gtsummary)
library(dplyr)
library(flextable)

# =================================================================
# LOAD AND PREPARE DATA
# =================================================================

sleep_raw <- read.csv("sleep_data_clean_for_imputation.csv")

sleep_raw <- sleep_raw %>%
  mutate(
    # Ensure PSQI components are numeric
    across(Comp.1:Comp.7, as.numeric),

    # Compute derived variables
    psqi.total = Comp.1 + Comp.2 + Comp.3 + Comp.4 +
                 Comp.5 + Comp.6 + Comp.7,
    bmi = weight / (height / 100)^2,

    # BMI category (WHO; underweight merged with Normal)
    bmic = factor(
      case_when(
        bmi < 25 ~ "Normal",
        bmi < 30 ~ "Overweight",
        bmi >= 30 ~ "Obesity"
      ),
      levels = c("Normal", "Overweight", "Obesity")
    ),

    # Anaemia status (WHO thresholds by sex)
    hbc = factor(
      case_when(
        sex == "female" & hb < 12 ~ "Anemia",
        sex == "male"   & hb < 13 ~ "Anemia",
        !is.na(hb) ~ "Normal"
      ),
      levels = c("Normal", "Anemia")
    ),

    sex = factor(sex)
  )

# =================================================================
# BUILD TABLE 1
# =================================================================

tabla_final <- sleep_raw %>%
  select(sex, age, bmi, bmic, wc, bf, water,
         hb, hbc, psqi.total) %>%
  tbl_summary(
    by = sex,
    statistic = list(
      all_continuous()  ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    missing = "no",
    label = list(
      age        ~ "Age (years)",
      bmi        ~ "BMI (kg/m²)",
      bmic       ~ "BMI Category",
      wc         ~ "Waist Circumference (cm)",
      bf         ~ "Body Fat (%)",
      water      ~ "Total Body Water (%)",
      hb         ~ "Hemoglobin (g/dL)",
      hbc        ~ "Anemia Status",
      psqi.total ~ "PSQI Global Score"
    )
  ) %>%
  add_overall(last = FALSE) %>%
  add_p(
    test = list(
      all_continuous()  ~ "wilcox.test",
      all_categorical() ~ "chisq.test"
    )
  ) %>%
  add_n() %>%
  bold_labels() %>%
  bold_p(t = 0.05) %>%
  italicize_levels() %>%
  modify_header(
    label  = "**Variable**",
    stat_0 = "**Overall (N = {N})**"
  )

# Preview in console
tabla_final

# =================================================================
# EXPORT
# =================================================================

tabla_final %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = "Table1_BySex.docx")

message("\nFile saved: Table1_BySex.docx")

# Footnote for manuscript:
# "Data are median (Q1, Q3) or n (%). Missing data: body fat
#  (n = 138), muscle mass (n = 127), waist circumference (n = 20),
#  hemoglobin (n = 29). p-values from Wilcoxon rank-sum test
#  (continuous) and Pearson's chi-squared test (categorical)."
