# =============================================================================
# 01_tableS1_univariable_ivw.R
# Table S1 — Univariable IVW MR of total BMI, brain- and adipose-instrumented
#            BMI on six kidney outcomes.
#
# Reproduces: manuscript Table S1 (also the IVW rows of Table S7.3).
# Preprint: Briasouli et al., medRxiv 2025.08.19.25333885
#           https://doi.org/10.1101/2025.08.19.25333885
#
# Run from the repository root:  Rscript R/01_tableS1_univariable_ivw.R
# =============================================================================
source(file.path("R", "config.R"))
source(file.path("R", "_load_exposures.R"))
suppressPackageStartupMessages({library(TwoSampleMR); library(dplyr)})

exposures <- list(`total BMI` = exposure_total_bmi,
                  brain       = exposures_tissue$brain,
                  adipose     = exposures_tissue$adipose)
all_snps  <- unique(unlist(lapply(exposures, `[[`, "SNP")))
outcomes  <- setNames(lapply(names(OUTCOMES), load_outcome, snps = all_snps), names(OUTCOMES))

rows <- list()
for (e in names(exposures)) for (o in names(outcomes)) {
  dat <- harmonise_data(exposures[[e]], outcomes[[o]], action = HARMONISE_ACTION)
  dat <- dat[dat$mr_keep, ]
  if (nrow(dat) < 2) next
  set.seed(SEED)
  r <- mr(dat, method_list = "mr_ivw")
  r$exposure <- e; r$outcome <- o
  rows[[paste(e,o)]] <- r
}
res <- bind_rows(rows) %>%
  mutate(binary = OUTCOME_BINARY[outcome],
         `Beta/OR` = ifelse(binary, exp(b), b),
         lower_CI  = ifelse(binary, exp(b-1.96*se), b-1.96*se),
         upper_CI  = ifelse(binary, exp(b+1.96*se), b+1.96*se)) %>%
  transmute(Outcome=outcome, Exposure=exposure, SNPs=nsnp,
            `Beta/OR`=round(`Beta/OR`,4), Se=round(se,4),
            lower_CI=round(lower_CI,4), upper_CI=round(upper_CI,4), P=pval)

out <- safe_path(file.path(DIR_RESULTS, "TableS1_univariable_IVW.csv"))
write.csv(res, out, row.names = FALSE)
message("Wrote ", out); print(as.data.frame(res), row.names = FALSE)
