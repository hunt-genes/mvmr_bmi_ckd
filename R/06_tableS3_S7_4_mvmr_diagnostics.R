# =============================================================================
# 06_tableS3_S7_4_mvmr_diagnostics.R
# MVMR diagnostic statistics using the MVMR package:
#   Table S3  — conditional F-statistics (Sanderson-Windmeijer, with phenotypic
#               covariance estimated from the instrument beta matrix)
#   Table S7.4 — MVMR-Egger estimates, pleiotropy intercepts, MVMR-Q
#
# Uses the same tissue_stack_mvmr (BRAIN_ADIPOSE_WEIGHTED.xlsx betas) as
# script 05, to ensure consistency with Table S2.
#
# Preprint: Briasouli et al., medRxiv 2025.08.19.25333885
# Run:  Rscript R/06_tableS3_S7_4_mvmr_diagnostics.R
# =============================================================================
source(file.path("R", "config.R"))
source(file.path("R", "_load_exposures.R"))
suppressPackageStartupMessages({
  library(TwoSampleMR); library(dplyr)
  if (!requireNamespace("MVMR", quietly=TRUE))
    stop("Package 'MVMR' required. Install with: remotes::install_github('WSpiller/MVMR')")
  library(MVMR)
})

if (is.null(exposures_tissue_mvmr))
  stop("BRAIN_ADIPOSE_WEIGHTED.xlsx not found. See config.R.")

tissue_stack_mvmr <- bind_rows(
  exposures_tissue_mvmr$brain   %>% mutate(exposure="brain",   id.exposure="brain"),
  exposures_tissue_mvmr$adipose %>% mutate(exposure="adipose", id.exposure="adipose"))
all_snps <- unique(tissue_stack_mvmr$SNP)

s3_rows  <- list()
s74_rows <- list()

for (o in names(OUTCOMES)) {
  message("Diagnostics: ", o)
  outc  <- load_outcome(o, all_snps)
  mvdat <- tryCatch(mv_harmonise_data(tissue_stack_mvmr, outc),
                    error=function(e){ message("  harmonise failed"); NULL })
  if (is.null(mvdat) || nrow(mvdat$exposure_beta) < 3) next

  bx <- mvdat$exposure_beta[, c(1,2)]   # adipose, brain
  sx <- mvdat$exposure_se  [, c(1,2)]
  by <- mvdat$outcome_beta
  sy <- mvdat$outcome_se

  # format_mvmr for MVMR package
  fmv <- tryCatch(
    format_mvmr(BXGs=bx, BYG=by, seBXGs=sx, seBYG=sy, RSID=row.names(bx)),
    error=function(e){ message("  format_mvmr failed"); NULL })
  if (is.null(fmv)) next

  # ---- Table S3: conditional F (phenotypic covariance from instrument betas)
  pheno_mat <- cor(bx)
  pos <- tryCatch(phenocov_mvmr(pheno_mat, sx), error=function(e) NULL)
  if (!is.null(pos)) {
    fstats <- tryCatch(strength_mvmr(r_input=fmv, gencov=pos)[1, ],
                       error=function(e) NULL)
    if (!is.null(fstats)) {
      # columns: exposure1=adipose, exposure2=brain
      s3_rows[[o]] <- data.frame(
        Outcome       = o,
        BMI_brain     = round(as.numeric(fstats["exposure2"]), 3),
        BMI_adipose   = round(as.numeric(fstats["exposure1"]), 3),
        stringsAsFactors = FALSE)
      message("  condF brain=", round(as.numeric(fstats["exposure2"]),2),
              " adipose=", round(as.numeric(fstats["exposure1"]),2))
    }
  }

  # ---- Table S7.4: MVMR-Egger + MVMR-Q
  pr <- tryCatch(pleiotropy_mvmr(fmv, gencov=0), error=function(e) NULL)
  # Egger slope estimates
  # Build Egger model manually (add intercept column)
  bx_e <- cbind(bx, intercept=1)
  W    <- diag(1/sy^2)
  bx_et <- t(bx_e)
  beta_e <- tryCatch(solve(bx_et %*% W %*% bx_e) %*% bx_et %*% W %*% by,
                     error=function(e) NULL)
  if (!is.null(beta_e)) {
    fitted_e <- bx_e %*% beta_e
    Q_e      <- sum(((by - fitted_e)/sy)^2)
    df_e     <- nrow(bx) - ncol(bx_e)
    sigma2_e <- Q_e / df_e
    cov_e    <- sigma2_e * solve(bx_et %*% W %*% bx_e)
    se_e     <- sqrt(diag(cov_e))
    # MVMR-IVW estimates (for reporting alongside Egger)
    beta_ivw <- beta_e[1:2]   # adipose, brain
    se_ivw   <- se_e[1:2]
    int_b    <- beta_e[3]
    int_se   <- se_e[3]
    int_p    <- 2 * pnorm(-abs(int_b/int_se))
    bin      <- OUTCOME_BINARY[o]
    s74_rows[[o]] <- data.frame(
      Outcome          = o,
      N_variants       = nrow(bx),
      adipose_est      = ifelse(bin, round(exp(beta_ivw[1]),3), round(beta_ivw[1],4)),
      adipose_lo       = ifelse(bin, round(exp(beta_ivw[1]-1.96*se_ivw[1]),3),
                                     round(beta_ivw[1]-1.96*se_ivw[1],4)),
      adipose_hi       = ifelse(bin, round(exp(beta_ivw[1]+1.96*se_ivw[1]),3),
                                     round(beta_ivw[1]+1.96*se_ivw[1],4)),
      p_adipose        = round(2*pnorm(-abs(beta_ivw[1]/se_ivw[1])),4),
      brain_est        = ifelse(bin, round(exp(beta_ivw[2]),3), round(beta_ivw[2],4)),
      brain_lo         = ifelse(bin, round(exp(beta_ivw[2]-1.96*se_ivw[2]),3),
                                     round(beta_ivw[2]-1.96*se_ivw[2],4)),
      brain_hi         = ifelse(bin, round(exp(beta_ivw[2]+1.96*se_ivw[2]),3),
                                     round(beta_ivw[2]+1.96*se_ivw[2],4)),
      p_brain          = round(2*pnorm(-abs(beta_ivw[2]/se_ivw[2])),4),
      egger_intercept  = round(int_b, 6),
      egger_se         = round(int_se, 6),
      p_intercept      = round(int_p, 4),
      Q                = round(Q_e, 2),
      Q_df             = df_e,
      Q_pval           = round(pchisq(Q_e, df_e, lower.tail=FALSE), 5),
      stringsAsFactors = FALSE)
    message("  OK")
  }
}

# ---- Table S3 output -------------------------------------------------------
s3 <- bind_rows(s3_rows)
cat("\nTable S3: Conditional F-statistics\n")
print(as.data.frame(s3), row.names=FALSE)
out3 <- safe_path(file.path(DIR_RESULTS, "TableS3_conditional_F.csv"))
write.csv(s3, out3, row.names=FALSE)
message("Wrote ", out3)

# ---- Table S7.4 output -----------------------------------------------------
s74 <- bind_rows(s74_rows)
cat("\nTable S7.4: MVMR-Egger diagnostics\n")
print(as.data.frame(s74), row.names=FALSE)
out74 <- safe_path(file.path(DIR_RESULTS, "TableS7_4_MVMR_Egger.csv"))
write.csv(s74, out74, row.names=FALSE)
message("Wrote ", out74)
