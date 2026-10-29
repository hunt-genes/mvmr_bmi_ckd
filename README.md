# Tissue-partitioned MR of BMI on kidney function

Code to reproduce the Mendelian randomization (MR) analyses in:

> Briasouli A, Leyden G, Andries A, Hallan SI, Richardson TG, Åsvold BO,
> Rasheed H, Brumpton BM. **Exploring the Causal Relationship Between Body Mass
> Index and Kidney Function: Tissue-Partitioned Mendelian Randomization.**
> https://www.mdpi.com/2073-4425/17/11/1339

The study partitions BMI genetic instruments into brain- and adipose-tissue
components (and gene-ontology pathway sets) and estimates their effects on six
kidney outcomes (BUN, eGFRcys, eGFRcrea, annual eGFR decline, and CKD from
Wuttke et al. and Pattaro et al.) using univariable and multivariable MR.

## Repository layout

```
.
├── R/                     # analysis scripts — one per table/analysis
│   ├── 00_install_dependencies.R
│   ├── config.R                 # <-- EDIT THIS: paths, settings, columns
│   ├── _load_exposures.R        # shared instrument loader (not run directly)
│   ├── 01_tableS1_univariable_ivw.R
│   ├── 02_tableS7_3_univariable_sensitivity.R
│   ├── 03_tableS4_heterogeneity.R
│   ├── 04_tableS5_gene_ontology_univariable.R
│   ├── 05_tableS2_S6_mvmr.R
│   ├── 06_tableS3_S7_4_mvmr_diagnostics.R
│   └── 07_leave_one_out.R
├── data/
│   ├── raw/               # <-- PUT INPUT DATA HERE (not tracked by git)
│   └── processed/         # harmonised objects written by the scripts
├── results/              # output tables (one CSV per manuscript table)
├── figures/              # output plots
├── run_all.R             # convenience: run every script in order
├── CITATION.cff
└── LICENSE
```

Code and input data are kept separate: **all inputs live in `data/raw/` and are
never committed** (see `.gitignore`). Scripts write only to `results/`,
`figures/` and `data/processed/`.

## Which script makes which table

| Script | Reproduces |
|---|---|
| `01_tableS1_univariable_ivw.R` | **Table S1** — univariable IVW (total / brain / adipose BMI) |
| `02_tableS7_3_univariable_sensitivity.R` | **Table S7.3** — IVW, weighted mode, MR-Egger, weighted median + Egger intercept |
| `03_tableS4_heterogeneity.R` | **Table S4** — Cochran's Q heterogeneity |
| `04_tableS5_gene_ontology_univariable.R` | **Table S5** — univariable IVW of the GO instrument sets |
| `05_tableS2_S6_mvmr.R` | **Table S2** (tissue MVMR) and **Table S6** (GO MVMR) |
| `06_tableS3_S7_4_mvmr_diagnostics.R` | **Table S3** (conditional F) and **Table S7.4** (MVMR-Egger + Q) |
| `07_leave_one_out.R` | leave-one-out robustness check |

## How to run

1. **Install dependencies** (once):
   ```r
   Rscript R/00_install_dependencies.R
   ```
   Requires R ≥ 4.1, plus `TwoSampleMR` (MRCIEU) and `MVMR` (WSpiller) from GitHub.

2. **Add your data** to `data/raw/` and check the filenames against `RAW` in
   `R/config.R`. See `data/raw/README.md` for the full list.

3. **Set the working directory to the repository root**, then run any script:
   ```r
   Rscript R/01_tableS1_univariable_ivw.R
   ```
   or run everything:
   ```r
   Rscript run_all.R
   ```

Each script prints its result table and writes a CSV to `results/`. Output files
are never overwritten — a numeric suffix is added if a file already exists.

## Analysis settings (in `R/config.R`)

- `HARMONISE_ACTION = 1` — `harmonise_data(action = 1)`, assuming forward-strand
  summary statistics.
- `SEED = 20250101` — fixes the bootstrap SEs of the weighted-median and
  weighted-mode estimators so results are reproducible.
- `PPA_THRESHOLD = 0.8` — tissue assignment (`PPA4 > 0.8`) for brain / adipose.

## Citing

Please cite the preprint (see `CITATION.cff`). This code is released under the
MIT License (see `LICENSE`).
