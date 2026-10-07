# =============================================================================
# _load_exposures.R  --  build the BMI instrument sets used across scripts
# Not run on its own; sourced by the analysis scripts after config.R.
#
# All tissue instrument information comes from Table S7.1 (tissue_instruments),
# which contains both the overall BMI betas (for univariable MR) and the
# tissue-specific weighted betas beta_brain/beta_adipose (for MVMR).
#
# Provides:
#   exposures_tissue      — list(brain, adipose) using overall beta (univariable MR)
#   exposures_tissue_mvmr — list(brain, adipose) using tissue-specific betas (MVMR only)
#   exposure_total_bmi    — 915 total-BMI instruments
#   exposures_go          — list(mental, metabolic) gene-ontology instrument sets
# =============================================================================
suppressPackageStartupMessages({library(TwoSampleMR); library(readxl); library(dplyr)})

.read_tissue_file <- function() {
  f <- raw_path("tissue_instruments")
  if (!file.exists(f)) stop("Missing tissue instrument file: ", f)
  t1 <- if (grepl("\\.xlsx?$", f, ignore.case=TRUE))
    as.data.frame(read_excel(f, skip=2))   # skip title + blank row
  else read.csv(f, check.names=FALSE)
  # Rename if columns came in with skip
  if (!"SNP" %in% names(t1) && ncol(t1) >= 13) {
    names(t1)[1:13] <- c("SNP","CHR","effect_allele","other_allele","eaf",
                          "beta","se","pval","samplesize",
                          "Gene_adipose","PPA4_adipose","Gene_brain","PPA4_brain")
    if (ncol(t1) >= 17)
      names(t1)[14:17] <- c("beta_brain","se_brain","beta_adipose","se_adipose")
  }
  t1 <- t1[!is.na(t1$SNP), ]
  t1
}

## ---- univariable exposures (overall BMI beta) ------------------------------
.read_tissue <- function() {
  t1 <- .read_tissue_file()
  need <- c("SNP","effect_allele","other_allele","eaf","beta","se","pval","PPA4_adipose","PPA4_brain")
  if (!all(need %in% names(t1))) stop("tissue file missing columns: ",
        paste(setdiff(need, names(t1)), collapse=", "))
  mk <- function(sub, lab) {
    sub$Phenotype <- lab
    e <- format_data(sub, type="exposure", phenotype_col="Phenotype",
                     snp_col="SNP", beta_col="beta", se_col="se", eaf_col="eaf",
                     effect_allele_col="effect_allele", other_allele_col="other_allele",
                     pval_col="pval", samplesize_col=if("samplesize"%in%names(sub))"samplesize" else NULL)
    e$exposure <- lab; e$id.exposure <- lab; e
  }
  list(adipose = mk(t1[t1$PPA4_adipose > PPA_THRESHOLD, ], "adipose"),
       brain   = mk(t1[t1$PPA4_brain   > PPA_THRESHOLD, ], "brain"))
}

## ---- MVMR exposures (tissue-specific weighted betas) ----------------------
.read_tissue_mvmr <- function() {
  t1 <- .read_tissue_file()
  need <- c("SNP","effect_allele","other_allele","eaf","pval",
            "beta_brain","se_brain","beta_adipose","se_adipose")
  if (!all(need %in% names(t1))) {
    stop("Table S7.1 is missing tissue-specific beta columns (beta_brain, se_brain, ",
         "beta_adipose, se_adipose). These should be present in the published S7.1 file. ",
         "Missing: ", paste(setdiff(need, names(t1)), collapse=", "))
  }
  mk <- function(beta_col, se_col, lab) {
    sub <- t1[!is.na(t1[[beta_col]]), ]
    sub$Phenotype <- lab
    sub$beta_use  <- sub[[beta_col]]
    sub$se_use    <- sub[[se_col]]
    e <- format_data(sub, type="exposure", phenotype_col="Phenotype",
                     snp_col="SNP", beta_col="beta_use", se_col="se_use", eaf_col="eaf",
                     effect_allele_col="effect_allele", other_allele_col="other_allele",
                     pval_col="pval", samplesize_col=NULL)
    e$exposure <- lab; e$id.exposure <- lab; e
  }
  list(brain   = mk("beta_brain",   "se_brain",   "brain"),
       adipose = mk("beta_adipose", "se_adipose", "adipose"))
}

## ---- total BMI instruments -------------------------------------------------
.read_total_bmi <- function() {
  f <- raw_path("bmi_instruments")
  if (!file.exists(f)) stop("Missing BMI instrument file: ", f)
  d <- read.table(f, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  if (!all(c("SNP","beta.exposure") %in% names(d)))
    stop("bmi.effs.txt not in TwoSampleMR exposure format")
  d$exposure <- "total BMI"; d$id.exposure <- "total BMI"; d
}

## ---- gene-ontology instruments ---------------------------------------------
.read_go <- function() {
  f <- raw_path("go_instruments")
  if (!file.exists(f)) stop("Missing GO instrument file: ", f)
  raw <- read.csv(f, check.names=FALSE)
  is_true <- function(x) x %in% c(TRUE,"TRUE","True","true","T",1,"1")
  mk <- function(sub, lab) {
    sub$Phenotype <- lab; sub$.n <- NA_real_
    e <- format_data(sub, type="exposure", phenotype_col="Phenotype",
                     snp_col="SNP", beta_col="beta", se_col="se", eaf_col="AF",
                     effect_allele_col="ALT", other_allele_col="REF",
                     pval_col="pval", samplesize_col=".n")
    e$exposure <- lab; e$id.exposure <- lab; e
  }
  list(mental    = mk(raw[is_true(raw$mental_pathway), ],    "disease of mental health"),
       metabolic = mk(raw[is_true(raw$metabolic_pathway), ], "disease of metabolism"))
}

exposures_tissue      <- .read_tissue()
exposures_tissue_mvmr <- tryCatch(.read_tissue_mvmr(),
  error=function(e){ message("Tissue MVMR betas not loaded: ", conditionMessage(e)); NULL })
exposure_total_bmi    <- .read_total_bmi()
exposures_go          <- tryCatch(.read_go(),
  error=function(e){ message("GO instruments not loaded: ", conditionMessage(e)); NULL })

message(sprintf(
  "Instruments loaded — adipose: %d, brain: %d, total BMI: %d%s%s",
  nrow(exposures_tissue$adipose), nrow(exposures_tissue$brain),
  nrow(exposure_total_bmi),
  if (!is.null(exposures_go)) sprintf(", GO mental: %d, metabolic: %d",
    nrow(exposures_go$mental), nrow(exposures_go$metabolic)) else "",
  if (!is.null(exposures_tissue_mvmr)) sprintf(
    " | MVMR betas: brain=%d, adipose=%d",
    nrow(exposures_tissue_mvmr$brain), nrow(exposures_tissue_mvmr$adipose)) else
    " | WARNING: tissue MVMR betas not found in S7.1"))
