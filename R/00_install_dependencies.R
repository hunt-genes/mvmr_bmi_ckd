# =============================================================================
# 00_install_dependencies.R  --  install everything the scripts need.
# Run once:  Rscript R/00_install_dependencies.R
# =============================================================================
cran <- c("data.table","dplyr","readxl","remotes","openxlsx","ggplot2")
new  <- cran[!cran %in% rownames(installed.packages())]
if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")

# TwoSampleMR and MVMR are on GitHub
if (!requireNamespace("TwoSampleMR", quietly = TRUE))
  remotes::install_github("MRCIEU/TwoSampleMR")
if (!requireNamespace("MVMR", quietly = TRUE))
  remotes::install_github("WSpiller/MVMR")

message("Dependencies ready.")
