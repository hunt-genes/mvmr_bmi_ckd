# =============================================================================
# 07_leave_one_out.R
# Leave-one-out IVW analysis (robustness of the BMI -> kidney estimates to any
# single variant). Writes one CSV per exposure/outcome and a summary.
# Preprint: Briasouli et al., medRxiv 2025.08.19.25333885
# Run:  Rscript R/07_leave_one_out.R
# =============================================================================
source(file.path("R", "config.R"))
source(file.path("R", "_load_exposures.R"))
suppressPackageStartupMessages({library(TwoSampleMR); library(dplyr)})

exposures <- list(`total BMI`=exposure_total_bmi, brain=exposures_tissue$brain,
                  adipose=exposures_tissue$adipose)
all_snps  <- unique(unlist(lapply(exposures, `[[`, "SNP")))
outcomes  <- setNames(lapply(names(OUTCOMES), load_outcome, snps = all_snps), names(OUTCOMES))

all <- list()
for (e in names(exposures)) for (o in names(outcomes)) {
  dat <- harmonise_data(exposures[[e]], outcomes[[o]], action = HARMONISE_ACTION)
  dat <- dat[dat$mr_keep, ]; if (nrow(dat) < 3) next
  loo <- mr_leaveoneout(dat); loo$exposure <- e; loo$outcome <- o
  all[[paste(e,o)]] <- loo
}
res <- bind_rows(all)
out <- safe_path(file.path(DIR_RESULTS, "leave_one_out_all.csv"))
write.csv(res, out, row.names = FALSE); message("Wrote ", out)
