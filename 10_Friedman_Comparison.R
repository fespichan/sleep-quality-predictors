# Friedman test
auc_matrix <- cbind(LR = auc_values_logit, RF = auc_values_rf, XGB = auc_values_xgb)
friedman.test(auc_matrix)

# Post-hoc pairwise
pairwise.wilcox.test(
  c(auc_values_logit, auc_values_rf, auc_values_xgb),
  rep(c("LR", "RF", "XGB"), each = 100),
  paired = TRUE, p.adjust.method = "bonferroni"
)

