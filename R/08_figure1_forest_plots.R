# =============================================================================
# 08_figure1_forest_plots.R
# Figure 1 — Forest plots of MR estimates for total, brain- and adipose-
#            instrumented BMI across the six kidney outcomes.
#   Panel A: univariable MR, CKD outcomes (odds ratios)
#   Panel B: univariable MR, continuous outcomes (beta)
#   Panel C: MVMR, CKD outcomes (odds ratios)
#   Panel D: MVMR, continuous outcomes (beta)
#
# Reads the CSVs written by scripts 01 (Table S1) and 05 (Table S2), so the
# figure always matches the tables. Run those first.
#
# Preprint: Briasouli et al., medRxiv 2025.08.19.25333885
# Run:  Rscript R/08_figure1_forest_plots.R
# =============================================================================
source(file.path("R", "config.R"))
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(forcats)})
if (!requireNamespace("patchwork", quietly = TRUE))
  install.packages("patchwork", repos = "https://cloud.r-project.org")
library(patchwork)

CKD_OUTCOMES <- c("CKD Wuttke", "CKD Pattaro")
EXP_COLORS   <- c("total BMI"="#e31a1c", "brain"="#1f78b4", "adipose"="#33a02c")
EXP_LABELS   <- c("total BMI"="Total BMI", "brain"="Brain-BMI", "adipose"="Adipose-BMI")

# newest matching results file (scripts never overwrite, so pick the latest)
latest <- function(prefix) {
  f <- list.files(DIR_RESULTS, pattern = prefix, full.names = TRUE)
  if (!length(f)) stop("No results file matching '", prefix, "' in ", DIR_RESULTS,
                       ". Run the analysis script that makes it first.")
  f[order(file.mtime(f), decreasing = TRUE)][1]
}

prep <- function(path) {
  d <- read.csv(path, check.names = FALSE)
  names(d) <- sub("^Beta/OR$", "est", names(d))
  d %>% mutate(is_ckd = Outcome %in% CKD_OUTCOMES,
               Exposure = factor(Exposure, levels = names(EXP_COLORS)))
}

forest_panel <- function(d, ckd, title) {
  dd <- filter(d, is_ckd == ckd)
  ref <- if (ckd) 1 else 0
  xlab <- if (ckd) "Odds ratio (95% CI)" else "Beta (95% CI)"
  ggplot(dd, aes(x = fct_rev(factor(Outcome)), y = est,
                 ymin = lower_CI, ymax = upper_CI, color = Exposure)) +
    geom_hline(yintercept = ref, linetype = "dashed", color = "grey40") +
    geom_pointrange(position = position_dodge(width = 0.6), size = 0.6) +
    coord_flip() +
    scale_color_manual(values = EXP_COLORS, labels = EXP_LABELS) +
    labs(title = title, x = NULL, y = xlab, color = "Instrument") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}

uni <- prep(latest("TableS1_univariable_IVW"))
mv  <- prep(latest("TableS2_MVMR_tissue"))

pA <- forest_panel(uni, TRUE,  "A  Univariable MR — CKD")
pB <- forest_panel(uni, FALSE, "B  Univariable MR — continuous traits")
pC <- forest_panel(mv,  TRUE,  "C  MVMR — CKD")
pD <- forest_panel(mv,  FALSE, "D  MVMR — continuous traits")

fig1 <- (pA | pB) / (pC | pD) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

out_pdf <- safe_path(file.path(DIR_FIGURES, "Figure1_forest.pdf"))
out_png <- safe_path(file.path(DIR_FIGURES, "Figure1_forest.png"))
ggsave(out_pdf, fig1, width = 11, height = 8)
ggsave(out_png, fig1, width = 11, height = 8, dpi = 300)
message("Wrote ", out_pdf, " and ", out_png)
