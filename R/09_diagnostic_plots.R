# =============================================================================
# 09_diagnostic_plots.R
# Supplementary MR diagnostic plots for every exposure/outcome pair:
#   - leave-one-out (IVW dropping each SNP)
#   - single-SNP forest
#   - scatter (SNP-exposure vs SNP-outcome with method slopes)
#   - funnel
# One multi-page PDF per exposure. Uses TwoSampleMR's built-in plotters.
#
# Preprint: Briasouli et al., medRxiv 2025.08.19.25333885
# Run:  Rscript R/09_diagnostic_plots.R
# =============================================================================
source(file.path("R", "config.R"))
source(file.path("R", "_load_exposures.R"))
suppressPackageStartupMessages({library(TwoSampleMR); library(ggplot2)})

exposures <- list(`total BMI`=exposure_total_bmi, brain=exposures_tissue$brain,
                  adipose=exposures_tissue$adipose)
all_snps  <- unique(unlist(lapply(exposures, `[[`, "SNP")))
outcomes  <- setNames(lapply(names(OUTCOMES), load_outcome, snps = all_snps), names(OUTCOMES))

for (e in names(exposures)) {
  plots <- list()
  for (o in names(outcomes)) {
    dat <- harmonise_data(exposures[[e]], outcomes[[o]], action = HARMONISE_ACTION)
    dat <- dat[dat$mr_keep, ]; if (nrow(dat) < 3) next
    set.seed(SEED)
    res    <- mr(dat, method_list = c("mr_ivw","mr_egger_regression",
                                      "mr_weighted_median","mr_weighted_mode"))
    single <- mr_singlesnp(dat)
    loo    <- mr_leaveoneout(dat)
    plots <- c(plots,
               mr_scatter_plot(res, dat),
               mr_forest_plot(single),
               mr_leaveoneout_plot(loo),
               mr_funnel_plot(single))
  }
  if (!length(plots)) next
  out <- safe_path(file.path(DIR_FIGURES, paste0("diagnostics_", gsub("[^A-Za-z0-9]+","_",e), ".pdf")))
  pdf(out, width = 7, height = 6); for (p in plots) print(p); dev.off()
  message("Wrote ", out)
}

