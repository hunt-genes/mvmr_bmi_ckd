# =============================================================================
# 03_tableS4_heterogeneity.R
# Table S4 — Cochran's Q heterogeneity test (IVW) for each exposure/outcome.
# Preprint: Briasouli et al., medRxiv 2025.08.19.25333885
# Run:  Rscript R/03_tableS4_heterogeneity.R
# =============================================================================
source(file.path("R", "config.R"))
source(file.path("R", "_load_exposures.R"))
suppressPackageStartupMessages({library(TwoSampleMR); library(dplyr)})

exposures <- list(brain=exposures_tissue$brain, adipose=exposures_tissue$adipose,
                  `total BMI`=exposure_total_bmi)
all_snps  <- unique(unlist(lapply(exposures, `[[`, "SNP")))
outcomes  <- setNames(lapply(names(OUTCOMES), load_outcome, snps = all_snps), names(OUTCOMES))

fmt_p <- function(p){if(is.na(p))return(NA_character_);if(p==0)return("0");
  if(p>=0.001)return(as.character(signif(p,2)));e<-floor(log10(p));m<-signif(p/10^e,2)
  if(m>=10){m<-m/10;e<-e+1};sprintf("%sx10^%d",m,e)}
elab <- c(brain="BMI_brain", adipose="BMI_adipose", `total BMI`="BMI")

rows <- list()
for (e in names(exposures)) for (o in names(outcomes)) {
  dat <- harmonise_data(exposures[[e]], outcomes[[o]], action = HARMONISE_ACTION)
  dat <- dat[dat$mr_keep, ]; if (nrow(dat) < 3) next
  h <- mr_heterogeneity(dat, method_list = "mr_ivw")
  rows[[paste(e,o)]] <- data.frame(Outcome=o, Exposure=elab[[e]],
    Q=round(h$Q), Q_df=h$Q_df, Q_P=h$Q_pval,
    `Q (P-value)`=paste0(round(h$Q)," (",fmt_p(h$Q_pval),")"), check.names=FALSE)
}
res <- bind_rows(rows)
out <- safe_path(file.path(DIR_RESULTS, "TableS4_heterogeneity.csv"))
write.csv(res, out, row.names = FALSE, na = "")
message("Wrote ", out); print(as.data.frame(res), row.names = FALSE)
