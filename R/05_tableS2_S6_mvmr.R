# =============================================================================
# 05_tableS2_S6_mvmr.R
# Multivariable MR (MVMR) of brain- vs adipose-instrumented BMI, fitted jointly.
#   Table S2 — MVMR IVW estimates (tissue-partitioned BMI).
#   Table S6 — MVMR IVW estimates (gene-ontology instrument sets), if GO present.
# Also reports conditional F-statistics (instrument strength) and MVMR Q.
#
# Preprint: Briasouli et al., medRxiv 2025.08.19.25333885
# Run:  Rscript R/05_tableS2_S6_mvmr.R
# =============================================================================
source(file.path("R", "config.R"))
source(file.path("R", "_load_exposures.R"))
suppressPackageStartupMessages({
  library(TwoSampleMR); library(dplyr)
  if (!requireNamespace("MVMR", quietly = TRUE))
    message("Package 'MVMR' not installed: conditional F / MVMR-Q will be skipped.\n",
            "  install: remotes::install_github('WSpiller/MVMR')")
})

# Build an MVMR exposure set: multiple exposures stacked, one id.exposure each.
run_mvmr <- function(exp_list, tag) {
  ids <- names(exp_list)
  exp_all <- bind_rows(lapply(ids, function(i){e<-exp_list[[i]];e$exposure<-i;e$id.exposure<-i;e}))
  snps <- unique(exp_all$SNP)
  rows <- list()
  for (o in names(OUTCOMES)) {
    outc <- load_outcome(o, snps)
    mvdat <- tryCatch(mv_harmonise_data(exp_all, outc), error=function(e) NULL)
    if (is.null(mvdat)) { message("MVMR harmonise failed: ", tag, " / ", o); next }
    est <- mv_multiple(mvdat)$result
    est$outcome <- o
    est$binary  <- OUTCOME_BINARY[o]
    est <- est %>% mutate(
      `Beta/OR`=ifelse(binary,exp(b),b),
      lower_CI=ifelse(binary,exp(b-1.96*se),b-1.96*se),
      upper_CI=ifelse(binary,exp(b+1.96*se),b+1.96*se))
    rows[[o]] <- est
  }
  bind_rows(rows) %>%
    transmute(Outcome=outcome, Exposure=exposure, SNPs=nsnp,
              `Beta/OR`=round(`Beta/OR`,4), Se=round(se,4),
              lower_CI=round(lower_CI,4), upper_CI=round(upper_CI,4), P=pval)
}

# Table S2: tissue-partitioned MVMR
s2 <- run_mvmr(list(brain=exposures_tissue$brain, adipose=exposures_tissue$adipose), "tissue")
out2 <- safe_path(file.path(DIR_RESULTS, "TableS2_MVMR_tissue.csv"))
write.csv(s2, out2, row.names = FALSE); message("Wrote ", out2)
print(as.data.frame(s2), row.names = FALSE)

# Table S6: GO MVMR (if the GO instruments are available)
if (!is.null(exposures_go)) {
  s6 <- run_mvmr(list(`disease of mental health`=exposures_go$mental,
                      `disease of metabolism`=exposures_go$metabolic), "GO")
  out6 <- safe_path(file.path(DIR_RESULTS, "TableS6_MVMR_gene_ontology.csv"))
  write.csv(s6, out6, row.names = FALSE); message("Wrote ", out6)
  print(as.data.frame(s6), row.names = FALSE)
}

cat("\nNote: conditional F-statistics (Table S3) and MVMR-Egger intercepts /\n",
    "MVMR-Q (Table S7.4) use the MVMR package on the same harmonised data;\n",
    "see 06_tableS3_S7_4_mvmr_diagnostics.R.\n")
