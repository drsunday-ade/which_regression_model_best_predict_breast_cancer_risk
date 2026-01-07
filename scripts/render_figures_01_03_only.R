#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
options(warn = 1)

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(data_dir = getwd())
  if (length(args) == 0) return(out)
  for (a in args) if (grepl("^--data-dir=", a)) out$data_dir <- sub("^--data-dir=", "", a)
  out
}
ARGS <- parse_args()

DATA_DIR <- normalizePath(ARGS$data_dir, winslash = "/", mustWork = FALSE)
OUT_DIR  <- file.path(DATA_DIR, "outputs_wdbc_jds")
FIG_DIR  <- file.path(OUT_DIR, "figures")
TAB_DIR  <- file.path(OUT_DIR, "tables")
AUX_DIR  <- file.path(OUT_DIR, "aux")

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

req <- function(pkgs, repos="https://cloud.r-project.org") {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = repos, dependencies = TRUE)
    if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p, call. = FALSE)
  }
}
PKGS <- c("dplyr","tibble","readr","ggplot2","pROC","purrr","forcats","scales","viridisLite","jsonlite")
req(PKGS)

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(readr); library(ggplot2)
  library(pROC); library(purrr); library(forcats); library(scales); library(viridisLite)
  library(jsonlite)
})

clip_prob <- function(p, eps = 1e-6) pmin(pmax(p, eps), 1 - eps)

pr_curve_manual <- function(y, p) {
  y <- as.integer(y); p <- as.numeric(p)
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]; p <- p[ok]
  if (length(unique(y)) < 2) return(tibble(recall = numeric(0), precision = numeric(0)))
  o <- order(p, decreasing = TRUE)
  y <- y[o]
  tp <- cumsum(y == 1)
  fp <- cumsum(y == 0)
  pos <- sum(y == 1)
  recall <- tp / pos
  precision <- tp / (tp + fp)
  tibble(recall = c(0, recall), precision = c(1, precision))
}

calib_bins <- function(y, p, bins = 10) {
  p <- clip_prob(p)
  probs <- seq(0, 1, length.out = bins + 1)
  brks <- suppressWarnings(as.numeric(stats::quantile(p, probs = probs, na.rm = TRUE, type = 7)))
  brks <- sort(unique(brks))
  if (length(brks) < 4) brks <- seq(0, 1, length.out = bins + 1)
  
  tibble(y = y, p = p) %>%
    mutate(bin = cut(p, breaks = brks, include.lowest = TRUE, right = TRUE)) %>%
    group_by(bin) %>%
    summarise(n = n(), mean_p = mean(p), mean_y = mean(y), .groups = "drop")
}

# ---- load existing outputs ----
oof_path  <- file.path(AUX_DIR, "oof_predictions.csv")
perf_path <- file.path(TAB_DIR, "table_03_performance_summary.csv")
man_path  <- file.path(AUX_DIR, "run_manifest.json")

if (!file.exists(oof_path))  stop("Missing: ", oof_path, call. = FALSE)
if (!file.exists(perf_path)) stop("Missing: ", perf_path, call. = FALSE)
if (!file.exists(man_path))  stop("Missing: ", man_path, call. = FALSE)

oof  <- readr::read_csv(oof_path, show_col_types = FALSE)
perf <- readr::read_csv(perf_path, show_col_types = FALSE)
man  <- jsonlite::fromJSON(man_path)

BEST <- as.character(man$best_strategy)
top3 <- perf %>%
  arrange(desc(roc_auc_mean), log_loss_mean) %>%
  slice_head(n = 3) %>%
  pull(strategy) %>%
  as.character()

show_strats <- unique(c(BEST, top3))

oof_top <- oof %>%
  filter(strategy %in% show_strats) %>%
  mutate(strategy_label = as.character(strategy),
         y = as.integer(y),
         p = clip_prob(as.numeric(p)))

# Use base split + imap_dfr to avoid dplyr group_modify() grouping-variable restrictions
by_strat <- split(oof_top, oof_top$strategy_label)

DPI   <- 320
FIG_W <- 9.5
FIG_H <- 6.0

# ---- Figure 1: ROC ----
roc_df <- purrr::imap_dfr(by_strat, function(df, lab) {
  if (length(unique(df$y)) < 2) return(tibble())
  r  <- pROC::roc(response = df$y, predictor = df$p, quiet = TRUE, direction = "<")
  cc <- pROC::coords(r, x = "all", ret = c("specificity","sensitivity"), transpose = FALSE)
  tibble(
    specificity = as.numeric(cc[,"specificity"]),
    sensitivity = as.numeric(cc[,"sensitivity"]),
    strategy_label = lab
  )
})

p1 <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity, color = strategy_label)) +
  geom_path(linewidth = 1.15, alpha = 0.95) +
  geom_abline(linetype = "dashed") +
  coord_equal() +
  labs(
    title = "Figure 1. ROC curves (pooled out-of-fold predictions)",
    subtitle = paste0("Best strategy: ", BEST, " | Top comparators: ", paste(top3, collapse = ", ")),
    x = "False Positive Rate (1 − Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "Strategy"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 15),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(FIG_DIR, "figure_01_roc.png"), p1, width = FIG_W, height = FIG_H, dpi = DPI, bg = "white")

# ---- Figure 2: PR ----
pr_df <- purrr::imap_dfr(by_strat, function(df, lab) {
  pr_curve_manual(df$y, df$p) %>% mutate(strategy_label = lab)
})

p2 <- ggplot(pr_df, aes(x = recall, y = precision, color = strategy_label)) +
  geom_path(linewidth = 1.15, alpha = 0.95) +
  labs(
    title = "Figure 2. Precision–Recall curves (pooled out-of-fold predictions)",
    subtitle = "PR curves emphasize malignant-class performance (event = 1) and are prevalence-sensitive",
    x = "Recall (Sensitivity)",
    y = "Precision (PPV)",
    color = "Strategy"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 15),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(FIG_DIR, "figure_02_pr.png"), p2, width = FIG_W, height = FIG_H, dpi = DPI, bg = "white")

# ---- Figure 3: Calibration ----
calib_long <- purrr::imap_dfr(by_strat, function(df, lab) {
  calib_bins(df$y, df$p, bins = 10) %>% mutate(strategy_label = lab)
})

p3 <- ggplot(calib_long, aes(x = mean_p, y = mean_y, color = strategy_label)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(aes(size = n), alpha = 0.9) +
  geom_line(linewidth = 1.1, alpha = 0.95) +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  scale_size_continuous(range = c(2.5, 7)) +
  labs(
    title = "Figure 3. Calibration (reliability) curves (10 bins)",
    subtitle = "Bin-level mean predicted risk vs observed malignant frequency; marker size ∝ bin n",
    x = "Mean predicted probability",
    y = "Observed event rate",
    color = "Strategy",
    size = "Bin n"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 15),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(FIG_DIR, "figure_03_calibration.png"), p3, width = FIG_W, height = FIG_H, dpi = DPI, bg = "white")

cat("\nOK: Re-rendered Figure 1–3 into:\n", FIG_DIR, "\n\n", sep = "")
cat("  - figure_01_roc.png\n  - figure_02_pr.png\n  - figure_03_calibration.png\n")
