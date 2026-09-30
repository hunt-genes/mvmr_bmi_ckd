# =============================================================================
# 04_tableS5_gene_ontology_univariable.R
# Table S5 — Univariable IVW of the "disease of mental health" and
#            "disease of metabolism" GO instrument sets on kidney outcomes.
# Preprint: Briasouli et al., medRxiv 2025.08.19.25333885
# Run:  Rscript R/04_tableS5_gene_ontology_univariable.R
# =============================================================================
source(file.path("R", "config.R"))
source(file.path("R", "_load_exposures.R"))
suppressPackageStartupMessages({library(TwoSampleMR); library(dplyr)})
if (is.null(exposures_go)) stop("GO instruments not available; check RAW$go_instruments in config.R")

exposures <- list(`disease of mental health`=exposures_go$mental,
                  `disease of metabolism`   =exposures_go$metabolic)
all_snps  <- unique(unlist(lapply(exposures, `[[`, "SNP")))
outcomes  <- setNames(lapply(names(OUTCOMES), load_outcome, snps = all_snps), names(OUTCOMES))

rows <- list()
for (e in names(exposures)) for (o in names(outcomes)) {
  dat <- harmonise_data(exposures[[e]], outcomes[[o]], action = HARMONISE_ACTION)
  dat <- dat[dat$mr_keep, ]; if (nrow(dat) < 2) next
  set.seed(SEED); r <- mr(dat, method_list = "mr_ivw")
  r$exposure <- e; r$outcome <- o; rows[[paste(e,o)]] <- r
}
res <- bind_rows(rows) %>%
  mutate(binary=OUTCOME_BINARY[outcome],
         `Beta/OR`=ifelse(binary,exp(b),b),
         lower_CI=ifelse(binary,exp(b-1.96*se),b-1.96*se),
         upper_CI=ifelse(binary,exp(b+1.96*se),b+1.96*se)) %>%
  transmute(Outcome=outcome, Exposure=exposure, SNPs=nsnp,
            `Beta/OR`=round(`Beta/OR`,4), Se=round(se,4),
            lower_CI=round(lower_CI,4), upper_CI=round(upper_CI,4), P=pval)
out <- safe_path(file.path(DIR_RESULTS, "TableS5_GO_univariable_IVW.csv"))
write.csv(res, out, row.names = FALSE)
message("Wrote ", out); print(as.data.frame(res), row.names = FALSE)
