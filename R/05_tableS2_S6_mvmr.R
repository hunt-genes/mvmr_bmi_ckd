# =============================================================================
# 05_tableS2_S6_mvmr.R
# Multivariable MR (MVMR) — tissue-partitioned and gene-ontology exposures
#   Table S2: MVMR IVW estimates (brain vs adipose tissue-specific betas)
#   Table S6: MVMR IVW estimates (mental vs metabolic GO instrument sets)
#
# Uses tissue-specific weighted betas from BRAIN_ADIPOSE_WEIGHTED.xlsx (S2)
# and wide-format zero-imputation for non-overlapping GO instruments (S6).
# Bypasses mv_multiple() (buggy in TwoSampleMR 0.7.x) with manual IVW-MVMR.
#
# Preprint: Briasouli et al., medRxiv 2025.08.19.25333885
# Run:  Rscript R/05_tableS2_S6_mvmr.R
# =============================================================================
source(file.path("R", "config.R"))
source(file.path("R", "_load_exposures.R"))
suppressPackageStartupMessages({library(TwoSampleMR); library(dplyr)})

# ---- manual IVW-MVMR (multiplicative random effects) -----------------------
# Bypasses mv_multiple() which has a dimnames bug in TwoSampleMR 0.7.x
ivw_mvmr_manual <- function(mvdat, outcome_label) {
  bx <- mvdat$exposure_beta   # SNPs x exposures
  by <- mvdat$outcome_beta
  sy <- mvdat$outcome_se
  W  <- diag(1/sy^2)
  bx_t <- t(bx)
  beta <- tryCatch(solve(bx_t %*% W %*% bx) %*% bx_t %*% W %*% by,
                   error=function(e){ message("  Singular matrix: ", outcome_label); return(NULL) })
  if (is.null(beta)) return(NULL)
  fitted <- bx %*% beta
  Q      <- sum(((by - fitted)/sy)^2)
  df     <- nrow(bx) - ncol(bx)
  sigma2 <- Q / df
  cov_b  <- sigma2 * solve(bx_t %*% W %*% bx)
  se_b   <- sqrt(diag(cov_b))
  pval   <- 2 * pnorm(-abs(as.numeric(beta)/se_b))
  data.frame(outcome=outcome_label, exposure=rownames(beta),
             nsnp=nrow(bx), b=as.numeric(beta), se=se_b, pval=pval,
             stringsAsFactors=FALSE)
}

fmt_mvmr_results <- function(rows) {
  bind_rows(rows) %>%
    mutate(binary = OUTCOME_BINARY[outcome],
           `Beta/OR` = ifelse(binary, exp(b), b),
           lower_CI  = ifelse(binary, exp(b - 1.96*se), b - 1.96*se),
           upper_CI  = ifelse(binary, exp(b + 1.96*se), b + 1.96*se)) %>%
    transmute(Outcome=outcome, Exposure=exposure, SNPs=nsnp,
              `Beta/OR`=round(`Beta/OR`,4), Se=round(se,5),
              lower_CI=round(lower_CI,4), upper_CI=round(upper_CI,4), P=pval)
}

# ============================================================
# Table S2: tissue-partitioned MVMR
# Uses tissue-specific weighted betas (beta_brain, beta_adipose)
# from BRAIN_ADIPOSE_WEIGHTED.xlsx — NOT the overall BMI beta.
# ============================================================
if (is.null(exposures_tissue_mvmr)) {
  stop("BRAIN_ADIPOSE_WEIGHTED.xlsx not found. ",
       "Place it in data/raw/ and ensure config.R points to it.")
}

tissue_stack_mvmr <- bind_rows(
  exposures_tissue_mvmr$brain   %>% mutate(exposure="brain",   id.exposure="brain"),
  exposures_tissue_mvmr$adipose %>% mutate(exposure="adipose", id.exposure="adipose"))

all_snps_mvmr <- unique(tissue_stack_mvmr$SNP)
s2_rows <- list()
for (o in names(OUTCOMES)) {
  message("MVMR tissue: ", o)
  outc  <- load_outcome(o, all_snps_mvmr)
  mvdat <- tryCatch(mv_harmonise_data(tissue_stack_mvmr, outc),
                    error=function(e){ message("  harmonise failed"); NULL })
  if (is.null(mvdat) || nrow(mvdat$exposure_beta) < 3) next
  res <- ivw_mvmr_manual(mvdat, o)
  if (!is.null(res)) s2_rows[[o]] <- res
}
s2 <- fmt_mvmr_results(s2_rows)
out2 <- safe_path(file.path(DIR_RESULTS, "TableS2_MVMR_tissue.csv"))
write.csv(s2, out2, row.names=FALSE)
message("Wrote ", out2)
print(as.data.frame(s2), row.names=FALSE)

# ============================================================
# Table S6: GO MVMR (wide-format zero-imputation)
# Mental and metabolic instruments are largely non-overlapping.
# Overlapping SNPs (in both sets) are removed from the metabolic
# set to avoid collinear columns. Non-overlapping SNPs are
# assigned zero effect for the pathway they do not belong to.
# ============================================================
if (!is.null(exposures_go)) {
  mental_exp    <- exposures_go$mental
  metabolic_exp <- exposures_go$metabolic

  # Remove overlapping SNPs from metabolic set
  overlap <- intersect(mental_exp$SNP, metabolic_exp$SNP)
  if (length(overlap) > 0) {
    message("Removing ", length(overlap), " overlapping SNPs from metabolic set: ",
            paste(overlap, collapse=", "))
    metabolic_exp <- metabolic_exp %>% filter(!SNP %in% overlap)
  }

  # Build wide df: one row per unique SNP, zero-impute missing pathway
  all_go_snps <- unique(c(mental_exp$SNP, metabolic_exp$SNP))
  mental_df   <- mental_exp   %>% select(SNP, beta_mental=beta.exposure,
                                         se_mental=se.exposure,
                                         ea=effect_allele.exposure,
                                         oa=other_allele.exposure,
                                         eaf=eaf.exposure)
  metabolic_df<- metabolic_exp %>% select(SNP, beta_metabolic=beta.exposure,
                                          se_metabolic=se.exposure)
  wide_df <- data.frame(SNP=all_go_snps) %>%
    left_join(mental_df,    by="SNP") %>%
    left_join(metabolic_df, by="SNP") %>%
    mutate(beta_mental    = ifelse(is.na(beta_mental),    0,   beta_mental),
           se_mental      = ifelse(is.na(se_mental),      1e6, se_mental),
           beta_metabolic = ifelse(is.na(beta_metabolic), 0,   beta_metabolic),
           se_metabolic   = ifelse(is.na(se_metabolic),   1e6, se_metabolic))

  message("GO wide df: ", nrow(wide_df), " SNPs (",
          sum(wide_df$beta_mental!=0), " mental, ",
          sum(wide_df$beta_metabolic!=0), " metabolic, ",
          sum(wide_df$beta_mental!=0 & wide_df$beta_metabolic!=0), " both)")

  all_snps_go <- wide_df$SNP
  s6_rows <- list()
  for (o in names(OUTCOMES)) {
    message("MVMR GO: ", o)
    outc <- load_outcome(o, all_snps_go)
    merged <- wide_df %>%
      inner_join(outc %>% select(SNP, by=beta.outcome, sy=se.outcome), by="SNP") %>%
      filter(!is.na(by), !is.na(sy), sy > 0)
    if (nrow(merged) < 5) { message("  too few SNPs"); next }
    bx <- as.matrix(merged[, c("beta_mental","beta_metabolic")])
    colnames(bx) <- c("disease of mental health","disease of metabolism")
    mvdat_wide <- list(exposure_beta=bx,
                       exposure_se=as.matrix(merged[, c("se_mental","se_metabolic")]),
                       outcome_beta=merged$by, outcome_se=merged$sy)
    res <- ivw_mvmr_manual(mvdat_wide, o)
    if (!is.null(res)) s6_rows[[o]] <- res
  }
  s6 <- fmt_mvmr_results(s6_rows)
  out6 <- safe_path(file.path(DIR_RESULTS, "TableS6_MVMR_gene_ontology.csv"))
  write.csv(s6, out6, row.names=FALSE)
  message("Wrote ", out6)
  print(as.data.frame(s6), row.names=FALSE)
} else {
  message("GO instruments not available; Table S6 skipped.")
}
