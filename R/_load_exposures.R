# =============================================================================
# _load_exposures.R  --  build the BMI instrument sets used across scripts
# Not run on its own; sourced by the analysis scripts after config.R.
# Provides: exposures_tissue (adipose, brain), exposure_total_bmi, exposures_go
# =============================================================================
suppressPackageStartupMessages({library(TwoSampleMR); library(readxl); library(dplyr)})

## tissue-partitioned instruments from Table S7.1 (PPA4 > threshold) ----------
.read_tissue <- function() {
  f <- raw_path("tissue_instruments")
  if (!file.exists(f)) stop("Missing tissue instrument file: ", f)
  t1 <- if (grepl("\\.xlsx?$", f, ignore.case=TRUE)) as.data.frame(read_excel(f)) else
        read.csv(f, check.names = FALSE)
  need <- c("SNP","effect_allele","other_allele","eaf","beta","se","pval","PPA4_adipose","PPA4_brain")
  if (!all(need %in% names(t1))) stop("tissue file missing columns: ",
        paste(setdiff(need, names(t1)), collapse=", "))
  mk <- function(sub, lab) {
    sub$Phenotype <- lab
    e <- format_data(sub, type="exposure", phenotype_col="Phenotype",
                     snp_col="SNP", beta_col="beta", se_col="se", eaf_col="eaf",
                     effect_allele_col="effect_allele", other_allele_col="other_allele",
                     pval_col="pval", samplesize_col=NULL)
    e$exposure <- lab; e$id.exposure <- lab; e
  }
  list(adipose = mk(t1[t1$PPA4_adipose > PPA_THRESHOLD, ], "adipose"),
       brain   = mk(t1[t1$PPA4_brain   > PPA_THRESHOLD, ], "brain"))
}

## total BMI instruments (S7.2) ----------------------------------------------
.read_total_bmi <- function() {
  f <- raw_path("bmi_instruments")
  if (!file.exists(f)) stop("Missing BMI instrument file: ", f)
  d <- read.table(f, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  if (!all(c("SNP","beta.exposure") %in% names(d)))
    stop("bmi.effs.txt not in TwoSampleMR exposure format")
  d$exposure <- "total BMI"; d$id.exposure <- "total BMI"; d
}

## gene-ontology instruments (mental / metabolic) ----------------------------
.read_go <- function() {
  f <- raw_path("go_instruments")
  if (!file.exists(f)) stop("Missing GO instrument file: ", f)
  raw <- if (grepl("\\.xlsx?$", f, ignore.case=TRUE)) as.data.frame(read_excel(f)) else
         read.csv(f, check.names = FALSE)
  # columns in the provided file: SNP, REF, ALT, AF, beta, se, pval, mental_pathway, metabolic_pathway
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

exposures_tissue   <- .read_tissue()
exposure_total_bmi <- .read_total_bmi()
exposures_go       <- tryCatch(.read_go(), error = function(e) { message("GO instruments not loaded: ", conditionMessage(e)); NULL })

message(sprintf("Instruments: adipose=%d, brain=%d, total BMI=%d%s",
        nrow(exposures_tissue$adipose), nrow(exposures_tissue$brain), nrow(exposure_total_bmi),
        if (!is.null(exposures_go)) sprintf(", GO mental=%d, GO metabolic=%d",
          nrow(exposures_go$mental), nrow(exposures_go$metabolic)) else ""))
