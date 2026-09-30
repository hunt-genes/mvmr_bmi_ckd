# =============================================================================
# run_all.R  --  run every analysis in order. Prefer running scripts one at a
# time (see README); this is a convenience wrapper.
#   Rscript run_all.R
# =============================================================================
PROJECT_ROOT <- getwd()
scripts <- c(
  "R/01_tableS1_univariable_ivw.R",
  "R/02_tableS7_3_univariable_sensitivity.R",
  "R/03_tableS4_heterogeneity.R",
  "R/04_tableS5_gene_ontology_univariable.R",
  "R/05_tableS2_S6_mvmr.R",
  "R/06_tableS3_S7_4_mvmr_diagnostics.R",
  "R/07_leave_one_out.R"
)
for (s in scripts) {
  message("\n=== ", s, " ===")
  tryCatch(source(s, local = new.env()),
           error = function(e) message("  FAILED: ", conditionMessage(e)))
}
