# `data/raw/` — input data (not tracked by git)

Place the files listed below in this folder, then run the scripts in order
(`run_all.R` or one by one). Either rename your files to match the names in
the table, or edit the `RAW` list in `R/config.R` to match your filenames.

## Files expected

### Instrument files
| Key in config | Expected filename | What it is |
|---|---|---|
| `tissue_instruments` | `Table_S7_1.xlsx` | **Supplementary Table S7, sheet "Table 1"** — 183 SNPs with PPA4 scores, overall BMI betas, and tissue-specific weighted betas (`beta_brain`, `se_brain`, `beta_adipose`, `se_adipose`). This single file is used for both univariable MR (overall beta) and MVMR (tissue-specific betas). Download from the journal alongside the paper. |
| `bmi_instruments` | `bmi.effs.txt` | 915 total-BMI instruments in TwoSampleMR exposure format, from Yengo et al. 2018 (GIANT + UK Biobank). Available from the GIANT Consortium or MRC-IEU OpenGWAS. |
| `go_instruments` | `filtered_metabolic_and_mental.csv` | Gene-ontology instrument sets (mental health and metabolic pathway SNPs). Columns: `SNP, REF, ALT, AF, beta, se, pval, mental_pathway, metabolic_pathway`. Available from the authors on request. |

### Outcome GWAS summary statistics
All from the CKDGen Consortium: **https://ckdgen.imbi.uni-freiburg.de**

| Key | Expected filename | Trait | Reference |
|---|---|---|---|
| `bun` | `metal_bun_meta_all1.TBL.map.annot.gc` | Blood urea nitrogen | CKDGen |
| `egfrcys` | `metal_eGFRcys_meta1.TBL.map.annot.gc` | eGFR (cystatin C) | Wuttke et al. 2019 |
| `egfrcrea` | `metal_eGFR_meta_ea1.TBL.map.annot.gc` | eGFR (creatinine) | Wuttke et al. 2019 |
| `egfrdecline` | `CKDGen_eGFR-decline_overall_adjDM.txt` | Annual eGFR decline (BGZF compressed — `fread` handles automatically) | Gorski et al. 2022 |
| `ckd_wuttke` | `CKD_overall_EA_JW_20180223_nstud23.dbgap.txt` | CKD case/control — European ancestry, CKDGen only (**does not include UK Biobank**) | Wuttke et al. 2019 |
| `ckd_pattaro` | `pattaro_gwas.txt` | CKD case/control | Pattaro et al. 2016 |

## Notes
- **`Table_S7_1.xlsx` is the only non-public file you need to add manually.** Everything else is publicly downloadable from the links above.
- The `egfrdecline` file is BGZF-compressed despite the `.txt` extension. `data.table::fread()` decompresses it automatically — no manual gunzip needed.
- The CKD Wuttke file has space-separated columns (sep=" ") — this is handled in `config.R`.
- The Pattaro file is comma-separated (sep=",") with `rsID` as the SNP column — also handled in `config.R`.
