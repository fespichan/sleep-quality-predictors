############################################################
# 09_Component_Screening.R
# Exploratory analysis: logistic regression for each of the
# 7 PSQI components (dichotomized 0 vs 1-3)
# Pooled without SMOTE using Rubin's rules
#
# Project: Biological and Anthropometric Determinants of
#          Sleep Quality Across the Adult Life Course
# Journal: Journal of Sleep Research
#
# Requires: data_stacked (from 02_Multiple_Imputation.R)
#
# Output: Table_S1_Component_Screening.xlsx
#         (3 sheets: All_Results, OR_Summary, P_Values)
############################################################

library(dplyr)
library(tidyr)
library(mice)
library(writexl)

# =================================================================
# SETUP
# =================================================================

predictors <- "sex + age + bmic + bf + hb"

components <- c("Comp.1", "Comp.2", "Comp.3", "Comp.4",
                "Comp.5", "Comp.6", "Comp.7")

comp_names <- c("Subjective Quality", "Sleep Latency",
                "Sleep Duration", "Sleep Efficiency",
                "Sleep Disturbances", "Medication Use",
                "Daytime Dysfunction")

n_imp <- max(data_stacked$.imp)
all_results <- data.frame()

# =================================================================
# COMPONENT-LEVEL LOGISTIC REGRESSION (RUBIN'S RULES)
# =================================================================

for (c in seq_along(components)) {

  comp <- components[c]
  name <- comp_names[c]

  message(sprintf("\n=== %s (%s) ===", comp, name))

  # Dichotomize: 0 = No problem vs 1-3 = Has problem
  data_stacked$outcome_temp <- factor(
    ifelse(data_stacked[[comp]] == 0, "No_problem", "Has_problem"),
    levels = c("No_problem", "Has_problem")
  )

  # Class distribution (first imputation)
  dist <- table(data_stacked$outcome_temp[data_stacked$.imp == 1])
  message(sprintf("  No_problem: %d | Has_problem: %d", dist[1], dist[2]))

  # Fit model on each imputed dataset
  model_list <- list()

  for (i in 1:n_imp) {
    capa_i <- data_stacked %>% filter(.imp == i) %>%
      dplyr::select(-.imp, -.id) %>%
      mutate(across(where(is.ordered),
                    ~factor(as.character(.x), ordered = FALSE))) %>%
      as.data.frame()

    formula_i <- as.formula(paste("outcome_temp ~", predictors))

    tryCatch({
      model_list[[i]] <- glm(formula_i, data = capa_i, family = "binomial")
    }, error = function(e) {
      model_list[[i]] <<- NULL
    })
  }

  # Remove failed models
  model_list <- Filter(Negate(is.null), model_list)

  if (length(model_list) >= 50) {
    pooled <- pool(as.mira(model_list))
    tabla  <- summary(pooled, exponentiate = TRUE, conf.int = TRUE)

    tabla$Component      <- comp
    tabla$Component_Name <- name
    all_results <- rbind(all_results, tabla)

    # Report significant predictors
    sig <- tabla %>% filter(p.value < 0.05 & term != "(Intercept)")
    if (nrow(sig) > 0) {
      message("  ** Significant predictors:")
      for (j in 1:nrow(sig)) {
        message(sprintf("    %s: OR = %.3f (%.3f-%.3f), p = %.4f",
                        sig$term[j], sig$estimate[j],
                        sig$`2.5 %`[j], sig$`97.5 %`[j], sig$p.value[j]))
      }
    } else {
      message("  No significant predictors")
    }
  } else {
    message("  ERROR: insufficient models converged")
  }
}

# Clean up
data_stacked$outcome_temp <- NULL

# =================================================================
# SUMMARY TABLES
# =================================================================

message("\n========================================")
message("SIGNIFICANCE SUMMARY (p-values)")
message("========================================")

summary_df <- all_results %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(Component_Name, term, estimate, p.value) %>%
  mutate(
    OR  = round(estimate, 3),
    p   = round(p.value, 4),
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  dplyr::select(Component_Name, term, OR, p, sig)

# Pivot: OR by component
table_or <- summary_df %>%
  dplyr::select(Component_Name, term, OR) %>%
  pivot_wider(names_from = Component_Name, values_from = OR)

# Pivot: p-values by component
table_p <- summary_df %>%
  mutate(p_sig = paste0(p, sig)) %>%
  dplyr::select(Component_Name, term, p_sig) %>%
  pivot_wider(names_from = Component_Name, values_from = p_sig)

message("\nOR by component:")
print(as.data.frame(table_or), row.names = FALSE)

message("\np-values by component:")
print(as.data.frame(table_p), row.names = FALSE)

# =================================================================
# EXPORT
# =================================================================

write_xlsx(
  list(
    All_Results = all_results,
    OR_Summary  = as.data.frame(table_or),
    P_Values    = as.data.frame(table_p)
  ),
  "Table_S1_Component_Screening.xlsx"
)

message("\nFile saved: Table_S1_Component_Screening.xlsx")
