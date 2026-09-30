# =============================================================================
# 02_tableS7_3_univariable_sensitivity.R
# Table S7.3 — Univariable sensitivity analyses: IVW, weighted mode, MR-Egger,
#              weighted median, plus the MR-Egger intercept, per exposure/outcome.
#
# Preprint: Briasouli et al., medRxiv 2025.08.19.25333885
# Run:  Rscript R/02_tableS7_3_univariable_sensitivity.R
# =============================================================================
source(file.path("R", "config.R"))
source(file.path("R", "_load_exposures.R"))
suppressPackageStartupMessages({library(TwoSampleMR); library(dplyr)})

METHODS <- c("mr_ivw","mr_weighted_mode","mr_egger_regression","mr_weighted_median")
exposures <- list(adipose=exposures_tissue$adipose, brain=exposures_tissue$brain,
                  `total BMI`=exposure_total_bmi)
all_snps  <- unique(unlist(lapply(exposures, `[[`, "SNP")))
outcomes  <- setNames(lapply(names(OUTCOMES), load_outcome, snps = all_snps), names(OUTCOMES))

fmt_int <- function(b,se){d<-max(3,-floor(log10(se)));sprintf("%s (%s)",formatC(b,format="f",digits=d),formatC(se,format="f",digits=d))}
rows <- list()
for (e in names(exposures)) for (o in names(outcomes)) {
  dat <- harmonise_data(exposures[[e]], outcomes[[o]], action = HARMONISE_ACTION)
  dat <- dat[dat$mr_keep, ]; if (nrow(dat) < 3) next
  set.seed(SEED)
  r  <- mr(dat, method_list = METHODS)
  pl <- mr_pleiotropy_test(dat)
  r$exposure <- e; r$outcome <- o
  r$egger_intercept <- ifelse(r$method=="MR Egger", fmt_int(pl$egger_intercept, pl$se), NA)
  r$egger_intercept_P <- ifelse(r$method=="MR Egger", pl$pval, NA)
  rows[[paste(e,o)]] <- r
}
res <- bind_rows(rows) %>%
  mutate(binary=OUTCOME_BINARY[outcome],
         `Beta/OR`=ifelse(binary,exp(b),b),
         lower_CI=ifelse(binary,exp(b-1.96*se),b-1.96*se),
         upper_CI=ifelse(binary,exp(b+1.96*se),b+1.96*se)) %>%
  transmute(Exposure=exposure, Outcome=outcome, Method=method, SNPs=nsnp,
            `Beta/OR`=round(`Beta/OR`,4), Se=round(se,5),
            lower_CI=round(lower_CI,4), upper_CI=round(upper_CI,4), P=pval,
            `MR-Egger intercept (SE)`=egger_intercept, `Intercept P`=egger_intercept_P)

out <- safe_path(file.path(DIR_RESULTS, "TableS7_3_univariable_sensitivity.csv"))
write.csv(res, out, row.names = FALSE, na = "")
message("Wrote ", out); print(as.data.frame(res), row.names = FALSE)
