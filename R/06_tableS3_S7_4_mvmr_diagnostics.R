# =============================================================================
# 06_tableS3_S7_4_mvmr_diagnostics.R
# Table S3   — conditional F-statistics (MVMR instrument strength).
# Table S7.4 — MVMR-Egger intercept (directional pleiotropy) and MVMR Q
#              (residual heterogeneity), per outcome, tissue-partitioned BMI.
#
# Requires the MVMR package (Sanderson/Spiller).
# Preprint: Briasouli et al., medRxiv 2025.08.19.25333885
# Run:  Rscript R/06_tableS3_S7_4_mvmr_diagnostics.R
# =============================================================================
source(file.path("R", "config.R"))
source(file.path("R", "_load_exposures.R"))
suppressPackageStartupMessages({library(TwoSampleMR); library(dplyr)})
if (!requireNamespace("MVMR", quietly = TRUE))
  stop("Install MVMR first: remotes::install_github('WSpiller/MVMR')")
library(MVMR)

exp_list <- list(brain=exposures_tissue$brain, adipose=exposures_tissue$adipose)
exp_all  <- bind_rows(lapply(names(exp_list), function(i){e<-exp_list[[i]];e$exposure<-i;e$id.exposure<-i;e}))
snps     <- unique(exp_all$SNP)

f_rows <- list(); d_rows <- list()
for (o in names(OUTCOMES)) {
  outc  <- load_outcome(o, snps)
  mvdat <- tryCatch(mv_harmonise_data(exp_all, outc), error=function(e) NULL)
  if (is.null(mvdat)) { message("skip ", o); next }
  # Assemble MVMR format: betas/ses for each exposure + outcome
  bx <- mvdat$exposure_beta; sx <- mvdat$exposure_se
  by <- mvdat$outcome_beta;  sy <- mvdat$outcome_se
  fmv <- format_mvmr(BXGs=bx, BYG=by, seBXGs=sx, seBYG=sy)
  sf  <- strength_mvmr(fmv, gencov = 0)          # conditional F
  pr  <- pleiotropy_mvmr(fmv, gencov = 0)        # MVMR-Egger intercept + Q
  f_rows[[o]] <- data.frame(Outcome=o, t(setNames(as.numeric(sf$exposure), paste0("condF_", colnames(sf)))))
  d_rows[[o]] <- data.frame(Outcome=o,
    Egger_intercept=unname(pr$pleiotropy_intercept[1]),
    Egger_intercept_SE=unname(pr$pleiotropy_se[1]),
    Egger_intercept_P=unname(pr$pleiotropy_pval[1]),
    Q=unname(pr$Qstat), Q_P=unname(pr$Qpval))
}
s3  <- bind_rows(f_rows); s74 <- bind_rows(d_rows)
o3  <- safe_path(file.path(DIR_RESULTS, "TableS3_conditional_F.csv"))
o74 <- safe_path(file.path(DIR_RESULTS, "TableS7_4_MVMR_egger_Q.csv"))
write.csv(s3, o3, row.names=FALSE); write.csv(s74, o74, row.names=FALSE)
message("Wrote ", o3, " and ", o74)
print(as.data.frame(s3), row.names=FALSE); print(as.data.frame(s74), row.names=FALSE)
