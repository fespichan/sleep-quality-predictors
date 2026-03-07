############################################################
# 03_Reliability_Analysis.R
# Internal consistency of the PSQI across 100 multiply
# imputed datasets: Cronbach's alpha & McDonald's omega
#
# Project: Biological and Anthropometric Determinants of
#          Sleep Quality Across the Adult Life Course
# Journal: Sleep Medicine
#
# Requires: data_stacked (from 02_Multiple_Imputation.R)
#
# Output: Table_S2_PSQI_Reliability.xlsx
############################################################

library(psych)
library(dplyr)
library(writexl)

# =================================================================
# SETUP
# =================================================================

comp_vars <- c("Comp.1", "Comp.2", "Comp.3", "Comp.4",
               "Comp.5", "Comp.6", "Comp.7")

n_imp <- max(data_stacked$.imp)

# Storage vectors
alpha_vals  <- numeric(n_imp)
omega_total <- numeric(n_imp)
omega_h     <- numeric(n_imp)

message(sprintf("\nComputing reliability across %d imputed datasets...", n_imp))

# =================================================================
# COMPUTE RELIABILITY PER IMPUTATION
# =================================================================

for (i in 1:n_imp) {

  if (i %% 20 == 0) message(sprintf("  Dataset %d of %d", i, n_imp))

  capa_i <- data_stacked %>%
    filter(.imp == i) %>%
    dplyr::select(all_of(comp_vars)) %>%
    as.data.frame()

  # Cronbach's alpha
  alpha_i <- psych::alpha(capa_i, check.keys = FALSE)
  alpha_vals[i] <- alpha_i$total$raw_alpha

  # McDonald's omega (total and hierarchical)
  # Primary method: maximum likelihood (ml)
  # Fallback: principal axis (pa) if ml fails to converge
  tryCatch({
    omega_i <- psych::omega(capa_i, nfactors = 1, plot = FALSE,
                            fm = "ml", rotate = "none")
    omega_total[i] <- omega_i$omega.tot
    omega_h[i]     <- omega_i$omega_h
  }, error = function(e) {
    tryCatch({
      omega_i <- psych::omega(capa_i, nfactors = 1, plot = FALSE,
                              fm = "pa", rotate = "none")
      omega_total[i] <<- omega_i$omega.tot
      omega_h[i]     <<- omega_i$omega_h
    }, error = function(e2) {
      omega_total[i] <<- NA
      omega_h[i]     <<- NA
    })
  })
}

# =================================================================
# RESULTS
# =================================================================

message("\n========================================")
message("PSQI RELIABILITY (pooled across imputations)")
message("========================================")

reliability_summary <- data.frame(
  Metric = c("Cronbach's Alpha",
             "McDonald's Omega Total",
             "McDonald's Omega Hierarchical"),
  Mean = round(c(mean(alpha_vals),
                 mean(omega_total, na.rm = TRUE),
                 mean(omega_h, na.rm = TRUE)), 3),
  SD = round(c(sd(alpha_vals),
               sd(omega_total, na.rm = TRUE),
               sd(omega_h, na.rm = TRUE)), 3),
  Min = round(c(min(alpha_vals),
                min(omega_total, na.rm = TRUE),
                min(omega_h, na.rm = TRUE)), 3),
  Max = round(c(max(alpha_vals),
                max(omega_total, na.rm = TRUE),
                max(omega_h, na.rm = TRUE)), 3),
  N_valid = c(sum(!is.na(alpha_vals)),
              sum(!is.na(omega_total)),
              sum(!is.na(omega_h)))
)

print(reliability_summary, row.names = FALSE)

# =================================================================
# EXPORT
# =================================================================

write_xlsx(reliability_summary, "Table_S2_PSQI_Reliability.xlsx")

# Manuscript-ready text
message(sprintf("\nFor manuscript:"))
message(sprintf("  Cronbach's alpha = %.2f (SD = %.3f)",
                mean(alpha_vals), sd(alpha_vals)))
message(sprintf("  McDonald's omega total = %.2f (SD = %.3f)",
                mean(omega_total, na.rm = TRUE),
                sd(omega_total, na.rm = TRUE)))
message(sprintf("  Omega > alpha indicates heterogeneous factor loadings"))

message("\nFile saved: Table_S2_PSQI_Reliability.xlsx")

