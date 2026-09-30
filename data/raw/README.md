# `data/raw/` — input data (not tracked by git)

The analysis scripts read every input from this folder. These files are **not**
included in the repository (they are large and access-controlled). Place them
here yourself, then either rename them to match `RAW` in `R/config.R`, or edit
`RAW` to match your filenames.

## Files expected

### BMI instruments
| Key in config | Expected file | What it is |
|---|---|---|
| `bmi_instruments` | `bmi.effs.txt` | 915 total-BMI instruments, TwoSampleMR exposure format (Table S7.2) |
| `tissue_instruments` | `Table_S7_1.csv` | brain/adipose SNPs with `PPA4_brain`, `PPA4_adipose` columns (Table S7.1) |
| `go_instruments` | `filtered_metabolic_and_mental.csv` | GO SNP sets; columns `SNP, REF, ALT, AF, beta, se, pval, mental_pathway, metabolic_pathway` |

### Outcome GWAS summary statistics
| Key | Expected file | Trait |
|---|---|---|
| `bun` | `metal_bun_meta_all1.TBL.map.annot.gc` | Blood urea nitrogen |
| `egfrcys` | `metal_eGFRcys_meta1.TBL.map.annot.gc` | eGFR (cystatin C) |
| `egfrcrea` | `metal_eGFR_meta_ea1.TBL.map.annot.gc` | eGFR (creatinine) |
| `egfrdecline` | `CKDGen_eGFR-decline_overall_adjDM.txt` | Annual eGFR decline |
| `ckd_wuttke` | `CKD_overall_EA_JW_20180223_nstud23.dbgap.txt` | CKD (Wuttke et al.) |
| `ckd_pattaro` | `pattaro_gwas.txt` | CKD (Pattaro et al.) |

## Sources
Outcome GWAS are from the CKDGen Consortium (Wuttke et al. 2019; Pattaro et al.
2016; Gorski et al. 2022) and Stanzick et al. 2021. BMI instruments derive from
the GIANT + UK Biobank meta-analysis (Yengo et al. 2018), tissue-partitioned as
in Leyden et al. 2025. See the manuscript Methods and Table 1 for full details
and accession numbers.
