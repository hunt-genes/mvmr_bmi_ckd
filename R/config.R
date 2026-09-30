# =============================================================================
# config.R  --  shared paths, settings and outcome definitions
# =============================================================================
# Sourced by every analysis script. Edit ONLY this file to point the code at
# your local data. Nothing here runs an analysis.
#
# Repository layout
#   data/raw/        <- input GWAS summary statistics + instrument files (NOT in git)
#   data/processed/  <- harmonised objects written by the scripts
#   results/         <- output tables (CSV/XLSX), one per manuscript table
#   figures/         <- output plots
#   R/               <- analysis scripts (this folder)
# =============================================================================

## ---- 0. project root -------------------------------------------------------
# All paths are relative to the repository root. If you open the .Rproj or set
# the working directory to the repo root, these resolve automatically.
if (!exists("PROJECT_ROOT")) PROJECT_ROOT <- getwd()

DIR_RAW       <- file.path(PROJECT_ROOT, "data", "raw")
DIR_PROCESSED <- file.path(PROJECT_ROOT, "data", "processed")
DIR_RESULTS   <- file.path(PROJECT_ROOT, "results")
DIR_FIGURES   <- file.path(PROJECT_ROOT, "figures")
for (d in c(DIR_PROCESSED, DIR_RESULTS, DIR_FIGURES))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

## ---- 1. analysis settings --------------------------------------------------
HARMONISE_ACTION <- 1          # TwoSampleMR harmonise_data(action = 1): assume
                               # forward strand (all summary stats forward-strand)
SEED             <- 20250101   # for bootstrapped weighted median / mode SEs
PPA_THRESHOLD    <- 0.8        # tissue assignment: PPA4 > 0.8 (brain / adipose)

## ---- 2. input file names (in data/raw/) ------------------------------------
# >>> Rename the files on disk to these, OR change the names here to match. <<<
RAW <- list(
  bmi_instruments   = "bmi.effs.txt",                                  # 915 total-BMI instruments (S7.2)
  tissue_instruments= "Table_S7_1.csv",                                # brain/adipose SNPs w/ PPA4 (S7.1)
  go_instruments    = "filtered_metabolic_and_mental.csv",             # GO SNP sets (mental/metabolic)
  bun               = "metal_bun_meta_all1.TBL.map.annot.gc",          # BUN
  egfrcys           = "metal_eGFRcys_meta1.TBL.map.annot.gc",          # eGFR cystatin C
  egfrcrea          = "metal_eGFR_meta_ea1.TBL.map.annot.gc",          # eGFR creatinine
  egfrdecline       = "CKDGen_eGFR-decline_overall_adjDM.txt",         # annual eGFR decline
  ckd_wuttke        = "CKD_overall_EA_JW_20180223_nstud23.dbgap.txt",  # CKD (Wuttke et al.)
  ckd_pattaro       = "pattaro_gwas.txt"                               # CKD (Pattaro et al.)
)

## ---- 3. outcome column mappings --------------------------------------------
# Each outcome: file key + the columns needed by format_data(). `binary=TRUE`
# means the effect is a log-OR (reported as OR in the tables).
OUTCOMES <- list(
  BUN         = list(file="bun",        sep="\t", snp="RSID", beta="Effect", se="StdErr",
                     eaf="Freq1", ea="Allele1", oa="Allele2", p="P.value",   n="n",           binary=FALSE),
  eGFRcys     = list(file="egfrcys",    sep="\t", snp="RSID", beta="Effect", se="StdErr",
                     eaf="Freq1", ea="Allele1", oa="Allele2", p="P.value",   n="n",           binary=FALSE),
  eGFRcrea    = list(file="egfrcrea",   sep="\t", snp="RSID", beta="Effect", se="StdErr",
                     eaf="Freq1", ea="Allele1", oa="Allele2", p="P.value",   n="n",           binary=FALSE),
  eGFRdecline = list(file="egfrdecline",sep="\t", snp="RSID", beta="Effect", se="StdErr",
                     eaf="Freq1", ea="Allele1", oa="Allele2", p="P-value",   n="N_total_sum", binary=FALSE),
  `CKD Wuttke`= list(file="ckd_wuttke", sep=" ",  snp="RSID", beta="Effect", se="StdErr",
                     eaf="Freq1", ea="Allele1", oa="Allele2", p="P-value",   n="n_total_sum", binary=TRUE),
  `CKD Pattaro`=list(file="ckd_pattaro",sep=",",  snp="rsID", beta="beta",   se="se",
                     eaf="freqA1",ea="allele1",oa="allele2", p="pval",       n="N",           binary=TRUE)
)

## ---- 4. helpers ------------------------------------------------------------
raw_path <- function(key) file.path(DIR_RAW, RAW[[key]])

# never overwrite an output file: name.ext -> name_1.ext, name_2.ext, ...
safe_path <- function(path) {
  if (!file.exists(path)) return(path)
  d <- dirname(path); ext <- tools::file_ext(path)
  stem <- tools::file_path_sans_ext(basename(path)); i <- 1
  repeat {
    cand <- file.path(d, if (nzchar(ext)) sprintf("%s_%d.%s", stem, i, ext) else sprintf("%s_%d", stem, i))
    if (!file.exists(cand)) return(cand); i <- i + 1
  }
}

# read one outcome GWAS, keep only the instrument SNPs, return TwoSampleMR outcome df
load_outcome <- function(name, snps) {
  suppressPackageStartupMessages({library(data.table); library(TwoSampleMR)})
  s <- OUTCOMES[[name]]; f <- raw_path(s$file)
  if (!file.exists(f)) stop("Missing outcome file for ", name, ": ", f)
  mn <- make.names
  d <- as.data.frame(data.table::fread(f, sep = s$sep))
  names(d) <- mn(names(d))
  keep <- mn(c(s$snp, s$beta, s$se, s$ea, s$oa, s$p))
  if (!all(keep %in% names(d)))
    stop(name, ": missing columns ", paste(setdiff(keep, names(d)), collapse=", "))
  d <- d[d[[mn(s$snp)]] %in% snps, , drop = FALSE]
  d$Phenotype <- name
  o <- format_data(d, type="outcome", phenotype_col="Phenotype",
                   snp_col=mn(s$snp), beta_col=mn(s$beta), se_col=mn(s$se), eaf_col=mn(s$eaf),
                   effect_allele_col=mn(s$ea), other_allele_col=mn(s$oa),
                   pval_col=mn(s$p), samplesize_col=mn(s$n))
  o$outcome <- name; o$id.outcome <- name; o
}

OUTCOME_BINARY <- vapply(OUTCOMES, function(x) x$binary, logical(1))

message("config.R loaded. PROJECT_ROOT = ", PROJECT_ROOT)
