#!/usr/bin/env Rscript
# ==============================================================================
# WDBC (Breast Cancer Wisconsin Diagnostic) — End-to-End JDS Pipeline (R-ONLY)
# BULLETPROOF EDITION v4 (FINAL):
#   - Removes yardstick roc_curve/pr_curve API fragility by using:
#       * pROC for ROC curve points
#       * Manual PR curve construction (version-agnostic)
#   - Fixes group_modify() scoping: uses .y$strategy_label (group keys), not .x
#   - Keeps non-fatal plotting (logs warnings, continues)
#   - Retains your nested CV + calibration + stability + outputs structure
# ==============================================================================

options(stringsAsFactors = FALSE)
options(warn = 1)
Sys.setenv(CLI_NUM_COLORS = "1")
options(cli.num_colors = 1)

# ----------------------------- CLI PARSER --------------------------------------
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(data_dir = getwd(), quick = FALSE)
  if (length(args) == 0) return(out)
  
  for (a in args) {
    if (identical(a, "--quick")) out$quick <- TRUE
    if (grepl("^--data-dir=", a)) out$data_dir <- sub("^--data-dir=", "", a)
  }
  out
}
ARGS <- parse_args()

# ----------------------------- USER CONFIG -------------------------------------
DATA_DIR   <- normalizePath(ARGS$data_dir, winslash = "/", mustWork = FALSE)
DATA_FILE  <- file.path(DATA_DIR, "wdbc.data")
NAMES_FILE <- file.path(DATA_DIR, "wdbc.names")

OUT_DIR <- file.path(DATA_DIR, "outputs_wdbc_jds")
FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR <- file.path(OUT_DIR, "tables")
AUX_DIR <- file.path(OUT_DIR, "aux")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(AUX_DIR, showWarnings = FALSE, recursive = TRUE)

timestamp_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
LOG_FILE <- file.path(OUT_DIR, paste0("PIPELINE_LOG_", timestamp_tag, ".txt"))

log_line <- function(...) {
  msg <- paste0(..., collapse = "")
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}
safe_condition_message <- function(e) {
  if (inherits(e, "condition")) return(conditionMessage(e))
  as.character(e)
}
die <- function(msg) {
  log_line("FATAL: ", msg)
  stop(msg, call. = FALSE)
}

with_logged_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      log_line("WARNING: ", safe_condition_message(w))
      invokeRestart("muffleWarning")
    }
  )
}

# --------------------------- PACKAGE LOADER ------------------------------------
require_or_install <- function(pkgs, repos = "https://cloud.r-project.org") {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      cat("Installing missing R package: ", p, "\n", sep = "")
      tryCatch({
        install.packages(p, repos = repos, dependencies = TRUE)
      }, error = function(e) {
        die(sprintf("R package install failed for '%s': %s", p, safe_condition_message(e)))
      })
      if (!requireNamespace(p, quietly = TRUE)) {
        die(paste0("Package still not available after install: ", p))
      }
    }
  }
}

PKGS <- c(
  "dplyr","tibble","tidyr","purrr","readr","stringr","forcats","rlang",
  "ggplot2","scales","viridisLite",
  "glmnet","rsample","yardstick","pROC","jsonlite","lubridate"
)
require_or_install(PKGS)

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(purrr); library(readr)
  library(stringr); library(forcats); library(rlang)
  library(ggplot2); library(scales); library(viridisLite)
  library(glmnet); library(rsample); library(yardstick); library(pROC)
  library(jsonlite); library(lubridate)
})

log_line("WDBC JDS pipeline started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
log_line("DATA_DIR: ", DATA_DIR)
log_line("OUT_DIR:  ", OUT_DIR)
log_line("LOG_FILE: ", LOG_FILE)
if (ARGS$quick) log_line("MODE: QUICK (reduced repeats/bootstraps for deadline runs)")

# ------------------------ OPTIONAL AUTO-DOWNLOAD (UCI) --------------------------
download_if_missing <- function() {
  if (file.exists(DATA_FILE) && file.exists(NAMES_FILE)) return(invisible(TRUE))
  
  log_line("One or both input files are missing. Attempting UCI download...")
  
  uci_base  <- "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin"
  uci_data  <- paste0(uci_base, "/wdbc.data")
  uci_names <- paste0(uci_base, "/wdbc.names")
  
  tryCatch({
    if (!file.exists(DATA_FILE)) {
      download.file(uci_data, destfile = DATA_FILE, mode = "wb", quiet = TRUE)
      log_line("Downloaded: ", DATA_FILE)
    }
    if (!file.exists(NAMES_FILE)) {
      download.file(uci_names, destfile = NAMES_FILE, mode = "wb", quiet = TRUE)
      log_line("Downloaded: ", NAMES_FILE)
    }
  }, error = function(e) {
    log_line("Download attempt failed: ", safe_condition_message(e))
  })
  
  if (!file.exists(DATA_FILE)) die(paste0("Cannot find 'wdbc.data' at: ", DATA_FILE, " (and download failed)."))
  if (!file.exists(NAMES_FILE)) die(paste0("Cannot find 'wdbc.names' at: ", NAMES_FILE, " (and download failed)."))
  
  invisible(TRUE)
}

# ----------------------------- METRIC HELPERS ----------------------------------
clip_prob <- function(p, eps = 1e-6) pmin(pmax(p, eps), 1 - eps)

log_loss_vec <- function(y, p) {
  p <- clip_prob(p)
  -mean(y*log(p) + (1-y)*log(1-p))
}
brier_vec <- function(y, p) mean((y - p)^2)

roc_auc_vec <- function(y, p) {
  r <- pROC::roc(response = y, predictor = p, quiet = TRUE, direction = "<")
  as.numeric(pROC::auc(r))
}

pr_auc_vec <- function(y, p) {
  df <- tibble(truth = factor(y, levels = c(0,1)), estimate = p)
  as.numeric(yardstick::pr_auc(df, truth, estimate, event_level = "second")$.estimate)
}

calib_bins <- function(y, p, bins = 10) {
  p <- clip_prob(p)
  probs <- seq(0, 1, length.out = bins + 1)
  brks <- suppressWarnings(as.numeric(stats::quantile(p, probs = probs, na.rm = TRUE, type = 7)))
  brks <- unique(brks)
  if (length(brks) < 4) brks <- seq(0, 1, length.out = bins + 1)
  brks <- sort(unique(brks))
  if (length(brks) < 2) brks <- c(0, 1)
  
  tibble(y = y, p = p) %>%
    mutate(bin = cut(p, breaks = brks, include.lowest = TRUE, right = TRUE)) %>%
    group_by(bin) %>%
    summarise(n = n(), mean_p = mean(p), mean_y = mean(y), .groups = "drop") %>%
    mutate(abs_gap = abs(mean_y - mean_p),
           w_abs_gap = abs_gap * (n / sum(n)))
}
ece_score <- function(y, p, bins = 10) sum(calib_bins(y, p, bins = bins)$w_abs_gap)

calib_slope_intercept <- function(y, p) {
  p <- clip_prob(p)
  z <- qlogis(p)
  fit <- tryCatch(
    suppressWarnings(glm(y ~ z, family = binomial())),
    error = function(e) NULL
  )
  if (is.null(fit)) return(tibble(intercept = NA_real_, slope = NA_real_))
  co <- coef(fit)
  tibble(intercept = unname(co[1]), slope = unname(co[2]))
}

net_benefit <- function(y, p, thresh) {
  y <- as.integer(y)
  pred_pos <- p >= thresh
  tp <- sum(pred_pos & (y == 1))
  fp <- sum(pred_pos & (y == 0))
  n <- length(y)
  (tp / n) - (fp / n) * (thresh / (1 - thresh))
}

# --- Version-agnostic PR curve (manual) ---
pr_curve_manual <- function(y, p) {
  y <- as.integer(y)
  p <- as.numeric(p)
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
  
  # Add anchor point at recall=0, precision=1 for nicer plotting
  tibble(recall = c(0, recall), precision = c(1, precision))
}

safe_write_csv <- function(df, path) {
  tryCatch(readr::write_csv(df, path),
           error = function(e) die(paste0("Failed write_csv: ", path, " :: ", safe_condition_message(e))))
}

safe_write_html <- function(df, path, caption = NULL, digits = 4) {
  tryCatch({
    if (nrow(df) == 0 || ncol(df) == 0) {
      html <- paste0("<html><body><p>", ifelse(is.null(caption),"",caption),
                     "</p><p>(Empty table)</p></body></html>")
      writeLines(html, path)
      return(invisible(TRUE))
    }
    df2 <- df
    for (j in seq_len(ncol(df2))) {
      if (is.numeric(df2[[j]])) df2[[j]] <- formatC(df2[[j]], format = "f", digits = digits)
    }
    header <- paste0("<tr>", paste0("<th>", names(df2), "</th>", collapse=""), "</tr>")
    rows <- apply(df2, 1, function(r) paste0("<tr>", paste0("<td>", r, "</td>", collapse=""), "</tr>"))
    cap <- if (is.null(caption)) "" else paste0(
      "<caption style='caption-side:top;font-weight:bold;text-align:left;margin-bottom:8px;'>",
      caption, "</caption>"
    )
    html <- paste0(
      "<html><body>",
      "<table border='1' cellpadding='6' cellspacing='0' style='border-collapse:collapse;font-family:Arial;font-size:13px;'>",
      cap, "<thead>", header, "</thead><tbody>", paste0(rows, collapse=""), "</tbody></table>",
      "</body></html>"
    )
    writeLines(html, path)
    invisible(TRUE)
  }, error = function(e) die(paste0("Failed write_html: ", path, " :: ", safe_condition_message(e))))
}

mean_se <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(list(mean = NA_real_, se = NA_real_))
  m <- mean(x)
  se <- stats::sd(x) / sqrt(length(x))
  list(mean = m, se = se)
}

# ----------------------------- MODEL HELPERS -----------------------------------
make_folds <- function(dat, v, strata_col, seed) {
  set.seed(seed)
  rsample::vfold_cv(dat, v = v, strata = !!rlang::sym(strata_col))
}

safe_cv_glmnet <- function(xtr, ytr, alpha, nfolds, type.measure = "deviance") {
  nfolds <- as.integer(nfolds)
  nfolds <- max(2L, nfolds)
  tryCatch(
    with_logged_warnings(cv.glmnet(
      xtr, ytr,
      family = "binomial",
      alpha = alpha,
      nfolds = nfolds,
      type.measure = type.measure,
      standardize = TRUE,
      maxit = 100000
    )),
    error = function(e) die(paste0("cv.glmnet failed (alpha=", alpha, ", nfolds=", nfolds, "): ",
                                   safe_condition_message(e)))
  )
}

fit_best_glmnet <- function(xtr, ytr, model_key, inner_folds, alpha_grid, one_se) {
  pick_lambda <- function(cvfit) if (one_se) cvfit$lambda.1se else cvfit$lambda.min
  
  min_class <- min(sum(ytr == 0), sum(ytr == 1))
  inner_folds <- min(inner_folds, max(2L, min_class))
  if (inner_folds < 2) die("Training fold has <2 samples in a class; cannot do inner CV safely.")
  
  if (model_key == "ridge") {
    cvfit <- safe_cv_glmnet(xtr, ytr, alpha = 0, nfolds = inner_folds)
    lam <- pick_lambda(cvfit)
    dev <- cvfit$cvm[which.min(abs(cvfit$lambda - lam))]
    return(list(alpha = 0, lambda = lam, inner_deviance = as.numeric(dev)))
  }
  
  if (model_key == "lasso") {
    cvfit <- safe_cv_glmnet(xtr, ytr, alpha = 1, nfolds = inner_folds)
    lam <- pick_lambda(cvfit)
    dev <- cvfit$cvm[which.min(abs(cvfit$lambda - lam))]
    return(list(alpha = 1, lambda = lam, inner_deviance = as.numeric(dev)))
  }
  
  if (model_key == "enet") {
    best <- list(alpha = NA_real_, lambda = NA_real_, inner_deviance = Inf)
    for (a in alpha_grid) {
      cvfit <- safe_cv_glmnet(xtr, ytr, alpha = a, nfolds = inner_folds)
      lam <- pick_lambda(cvfit)
      dev <- cvfit$cvm[which.min(abs(cvfit$lambda - lam))]
      if (is.finite(dev) && dev < best$inner_deviance) best <- list(alpha = a, lambda = lam, inner_deviance = as.numeric(dev))
    }
    return(best)
  }
  
  die(paste0("Unknown model_key: ", model_key))
}

fit_glmnet_fixed <- function(x, y, alpha, lambda) {
  tryCatch(
    with_logged_warnings(glmnet(x, y, family = "binomial", alpha = alpha, lambda = lambda,
                                standardize = TRUE, maxit = 100000)),
    error = function(e) die(paste0("glmnet failed (alpha=", alpha, "): ", safe_condition_message(e)))
  )
}

predict_prob <- function(fit, xnew) as.numeric(predict(fit, newx = xnew, type = "response"))

train_oof_pred_fixed <- function(xtr, ytr, alpha, lambda, kfold, seed) {
  min_class <- min(sum(ytr == 0), sum(ytr == 1))
  kfold <- min(kfold, max(2L, min_class))
  if (kfold < 2) die("Calibration CV cannot run: too few samples per class in training fold.")
  
  dat <- tibble(idx = seq_along(ytr), y = factor(ytr, levels = c(0,1)))
  folds <- make_folds(dat, v = kfold, strata_col = "y", seed = seed)
  
  p_oof <- rep(NA_real_, length(ytr))
  for (i in seq_len(nrow(folds))) {
    asmt <- rsample::assessment(folds$splits[[i]])$idx
    anl  <- rsample::analysis(folds$splits[[i]])$idx
    fit <- fit_glmnet_fixed(xtr[anl,,drop=FALSE], ytr[anl], alpha = alpha, lambda = lambda)
    p_oof[asmt] <- predict_prob(fit, xtr[asmt,,drop=FALSE])
  }
  if (anyNA(p_oof)) die("Calibration OOF prediction produced NA values (unexpected).")
  clip_prob(p_oof)
}

fit_sigmoid_calibrator <- function(y, p_oof) {
  p_oof <- clip_prob(p_oof)
  z <- qlogis(p_oof)
  fit <- tryCatch(suppressWarnings(glm(y ~ z, family = binomial())), error = function(e) NULL)
  if (is.null(fit)) {
    log_line("WARNING: sigmoid calibrator glm failed; falling back to identity calibration.")
    return(NULL)
  }
  fit
}
apply_sigmoid_calibrator <- function(cal_fit, p) {
  if (is.null(cal_fit)) return(clip_prob(p))
  p <- clip_prob(p)
  z <- qlogis(p)
  as.numeric(predict(cal_fit, newdata = data.frame(z = z), type = "response"))
}

fit_isotonic_calibrator <- function(y, p_oof) {
  p_oof <- clip_prob(p_oof)
  ord <- order(p_oof)
  x <- p_oof[ord]; yy <- y[ord]
  iso <- tryCatch(isoreg(x, yy), error = function(e) NULL)
  if (is.null(iso)) {
    log_line("WARNING: isotonic isoreg failed; falling back to identity calibration.")
    return(NULL)
  }
  list(x = iso$x, yfit = pmin(pmax(iso$yf, 0), 1))
}
apply_isotonic_calibrator <- function(iso, p) {
  if (is.null(iso)) return(clip_prob(p))
  p <- clip_prob(p)
  as.numeric(approx(iso$x, iso$yfit, xout = p, rule = 2, ties = "ordered")$y)
}

coef_mask <- function(fit) {
  co <- as.matrix(coef(fit))
  as.numeric(abs(co[-1, 1]) > 1e-12)
}
jaccard <- function(a, b) {
  a <- as.logical(a); b <- as.logical(b)
  inter <- sum(a & b); uni <- sum(a | b)
  if (uni == 0) return(1)
  inter / uni
}

safe_plot <- function(expr, name) {
  tryCatch(expr, error = function(e) {
    log_line("WARNING: plot failed [", name, "]: ", safe_condition_message(e))
    invisible(NULL)
  })
}

# ------------------------------ MAIN PIPELINE ----------------------------------
main <- function() {
  with_logged_warnings({
    download_if_missing()
    
    # ------------------------ READ + CHECK wdbc.names ---------------------------
    names_lines <- readLines(NAMES_FILE, warn = FALSE, encoding = "UTF-8")
    if (length(names_lines) < 10) die("wdbc.names seems too short/unreadable. Replace with the official file.")
    
    attr_start <- which(str_detect(names_lines, regex("Attribute Information", ignore_case = TRUE)))
    if (length(attr_start) == 0) attr_start <- 1
    excerpt <- names_lines[attr_start:min(length(names_lines), attr_start + 120)]
    writeLines(excerpt, file.path(AUX_DIR, "wdbc_names_excerpt.txt"))
    
    has_diag <- any(str_detect(names_lines, regex("Diagnosis", ignore_case = TRUE)))
    has_30 <- any(str_detect(names_lines, regex("30", ignore_case = TRUE))) &&
      any(str_detect(names_lines, regex("features", ignore_case = TRUE)))
    log_line("wdbc.names checks: Diagnosis mentioned = ", has_diag, "; mentions 30 features = ", has_30)
    
    # ------------------------------ LOAD DATA -----------------------------------
    col_base <- c("radius","texture","perimeter","area","smoothness",
                  "compactness","concavity","concave_points","symmetry","fractal_dimension")
    feat_names <- c(paste0("mean_", col_base),
                    paste0("se_", col_base),
                    paste0("worst_", col_base))
    
    wdbc <- read.csv(DATA_FILE, header = FALSE, stringsAsFactors = FALSE)
    if (ncol(wdbc) != 32) die(sprintf("Expected 32 columns in wdbc.data, got %d. Wrong file or corrupted.", ncol(wdbc)))
    colnames(wdbc) <- c("id","diagnosis", feat_names)
    
    wdbc <- wdbc %>%
      mutate(
        diagnosis = factor(as.character(diagnosis), levels = c("B","M")),
        y = if_else(diagnosis == "M", 1L, 0L)
      )
    
    for (nm in feat_names) wdbc[[nm]] <- suppressWarnings(as.numeric(wdbc[[nm]]))
    if (any(!is.finite(as.matrix(wdbc[, feat_names])))) {
      die("Non-numeric or missing values detected in feature columns after coercion. Check wdbc.data formatting.")
    }
    
    # ----------------------------- ANALYSIS CONFIG ------------------------------
    OUTER_FOLDS <- 10
    N_REPEATS   <- 10
    INNER_FOLDS <- 5
    CALIB_CV    <- 3
    BOOT_B      <- 200
    STABILITY_PI_THR <- 0.75
    ONE_SE_RULE <- TRUE
    SEED        <- 20260106
    
    if (ARGS$quick) {
      N_REPEATS <- 5
      BOOT_B <- 100
      log_line("QUICK MODE applied: N_REPEATS=", N_REPEATS, " BOOT_B=", BOOT_B)
    }
    
    DPI   <- 320
    FIG_W <- 9.5
    FIG_H <- 6.0
    
    set.seed(SEED)
    
    log_line("CONFIG: OUTER_FOLDS=", OUTER_FOLDS,
             " N_REPEATS=", N_REPEATS,
             " INNER_FOLDS=", INNER_FOLDS,
             " CALIB_CV=", CALIB_CV,
             " BOOT_B=", BOOT_B,
             " PI_THR=", STABILITY_PI_THR,
             " ONE_SE_RULE=", ONE_SE_RULE,
             " SEED=", SEED)
    
    # ----------------------------- TABLE 1–2 -----------------------------------
    tab1 <- wdbc %>%
      summarise(
        n = n(),
        n_features = length(feat_names),
        n_benign = sum(y == 0),
        n_malignant = sum(y == 1),
        prevalence_malignant = mean(y == 1)
      )
    safe_write_csv(tab1, file.path(TAB_DIR, "table_01_data_summary.csv"))
    safe_write_html(tab1, file.path(TAB_DIR, "table_01_data_summary.html"),
                    caption = "Table 1. WDBC dataset summary.", digits = 6)
    
    tab2 <- wdbc %>%
      select(all_of(feat_names)) %>%
      summarise(across(everything(),
                       list(mean=~mean(.), sd=~sd(.), min=~min(.), max=~max(.)),
                       .names="{.col}__{.fn}")) %>%
      pivot_longer(everything(), names_to="feature_stat", values_to="value") %>%
      separate(feature_stat, into=c("feature","stat"), sep="__") %>%
      pivot_wider(names_from=stat, values_from=value) %>%
      arrange(feature)
    safe_write_csv(tab2, file.path(TAB_DIR, "table_02_feature_summary.csv"))
    safe_write_html(tab2, file.path(TAB_DIR, "table_02_feature_summary.html"),
                    caption = "Table 2. Feature summary (raw scale).", digits = 6)
    
    # ----------------------------- RESAMPLING CORE ------------------------------
    X <- as.matrix(wdbc[, feat_names])
    y <- as.integer(wdbc$y)
    n <- length(y)
    
    outer_dat <- tibble(row_id = 1:n, y = factor(y, levels=c(0,1)))
    alpha_grid <- c(0.10, 0.30, 0.50, 0.70, 0.90)
    
    metric_rows <- list()
    oof_rows <- list()
    stab_rows <- list()
    
    split_id <- 0
    total_splits <- OUTER_FOLDS * N_REPEATS
    
    log_line("Starting repeated outer CV: ", OUTER_FOLDS, "-fold x ", N_REPEATS,
             " repeats (", total_splits, " splits)")
    
    for (r in seq_len(N_REPEATS)) {
      outer_folds <- make_folds(outer_dat, v = OUTER_FOLDS, strata_col = "y", seed = SEED + 1000*r)
      
      for (k in seq_len(nrow(outer_folds))) {
        split_id <- split_id + 1
        
        tryCatch({
          test_idx  <- rsample::assessment(outer_folds$splits[[k]])$row_id
          train_idx <- rsample::analysis(outer_folds$splits[[k]])$row_id
          
          xtr <- X[train_idx,,drop=FALSE]; ytr <- y[train_idx]
          xte <- X[test_idx,,drop=FALSE];  yte <- y[test_idx]
          
          if (length(unique(yte)) < 2) {
            log_line("WARNING: degenerate test fold at split_id=", split_id, " (single class). Some metrics set to NA.")
          }
          
          for (model_key in c("ridge","lasso","enet")) {
            best <- fit_best_glmnet(
              xtr, ytr, model_key,
              inner_folds = INNER_FOLDS,
              alpha_grid = alpha_grid,
              one_se = ONE_SE_RULE
            )
            
            base_fit <- fit_glmnet_fixed(xtr, ytr, alpha = best$alpha, lambda = best$lambda)
            p_raw_te <- clip_prob(predict_prob(base_fit, xte))
            
            p_oof_tr <- train_oof_pred_fixed(
              xtr, ytr,
              alpha = best$alpha,
              lambda = best$lambda,
              kfold = CALIB_CV,
              seed = SEED + split_id
            )
            
            sig_fit <- fit_sigmoid_calibrator(ytr, p_oof_tr)
            iso_fit <- fit_isotonic_calibrator(ytr, p_oof_tr)
            
            p_sig_te <- clip_prob(apply_sigmoid_calibrator(sig_fit, p_raw_te))
            p_iso_te <- clip_prob(apply_isotonic_calibrator(iso_fit, p_raw_te))
            
            for (cal_key in c("raw","sigmoid","isotonic")) {
              strat <- paste0(model_key, "_", cal_key)
              p_te <- switch(cal_key, raw = p_raw_te, sigmoid = p_sig_te, isotonic = p_iso_te)
              
              roc_auc <- if (length(unique(yte)) < 2) NA_real_ else roc_auc_vec(yte, p_te)
              pr_auc  <- if (length(unique(yte)) < 2) NA_real_ else pr_auc_vec(yte, p_te)
              ll      <- log_loss_vec(yte, p_te)
              br      <- brier_vec(yte, p_te)
              
              metric_rows[[length(metric_rows)+1]] <- tibble(
                split_id = split_id,
                repeat_id = r,
                fold = k,
                strategy = strat,
                model = model_key,
                calibration = cal_key,
                best_inner_deviance = best$inner_deviance,
                alpha = best$alpha,
                lambda = best$lambda,
                roc_auc = roc_auc,
                pr_auc = pr_auc,
                log_loss = ll,
                brier = br,
                n_test = length(test_idx),
                n_pos_test = sum(yte == 1)
              )
              
              oof_rows[[length(oof_rows)+1]] <- tibble(
                split_id = split_id,
                strategy = strat,
                id = wdbc$id[test_idx],
                y = yte,
                p = p_te
              )
            }
            
            if (model_key %in% c("lasso","enet") && BOOT_B > 0) {
              sel_sets <- vector("list", BOOT_B)
              sel_counts <- rep(0L, ncol(xtr))
              
              for (b in seq_len(BOOT_B)) {
                boot_idx <- sample(seq_along(ytr), size = length(ytr), replace = TRUE)
                fit_b <- fit_glmnet_fixed(xtr[boot_idx,,drop=FALSE], ytr[boot_idx], alpha = best$alpha, lambda = best$lambda)
                m <- coef_mask(fit_b)
                sel_sets[[b]] <- as.logical(m)
                sel_counts <- sel_counts + as.integer(m)
              }
              
              sel_freq <- sel_counts / BOOT_B
              
              m_sub <- min(length(sel_sets), 50)
              sims <- c()
              if (m_sub >= 2) {
                for (i in 1:(m_sub-1)) for (j in (i+1):m_sub) sims <- c(sims, jaccard(sel_sets[[i]], sel_sets[[j]]))
              }
              mean_jacc <- if (length(sims) > 0) mean(sims) else 1
              mean_k <- mean(vapply(sel_sets, function(s) sum(s), numeric(1)))
              
              strat_label <- paste0(model_key, "_raw")
              
              stab_rows[[length(stab_rows)+1]] <- tibble(
                split_id = split_id,
                strategy = strat_label,
                model = model_key,
                feature = feat_names,
                selection_freq = as.numeric(sel_freq),
                stable_selected = as.integer(sel_freq >= STABILITY_PI_THR),
                mean_jaccard_subsample = as.numeric(mean_jacc),
                mean_k_selected = as.numeric(mean_k),
                alpha = best$alpha,
                lambda = best$lambda
              )
            }
          }
          
          if (split_id %% 5 == 0) log_line("Progress: completed split_id=", split_id, " / ", total_splits)
          
        }, error = function(e) {
          die(paste0("Split failed at split_id=", split_id, " (repeat_id=", r, ", fold=", k, "): ",
                     safe_condition_message(e)))
        })
      }
    }
    
    fold_metrics <- bind_rows(metric_rows)
    oof <- bind_rows(oof_rows)
    stab <- if (length(stab_rows) > 0) bind_rows(stab_rows) else tibble()
    
    if (nrow(fold_metrics) == 0 || nrow(oof) == 0) die("No CV results produced. Check earlier logs for failure context.")
    
    safe_write_csv(fold_metrics, file.path(AUX_DIR, "fold_metrics.csv"))
    safe_write_csv(oof, file.path(AUX_DIR, "oof_predictions.csv"))
    safe_write_csv(stab, file.path(AUX_DIR, "stability_selection.csv"))
    
    # ----------------------- CALIBRATION EXTRAS PER SPLIT -----------------------
    split_calib <- oof %>%
      group_by(split_id, strategy) %>%
      group_modify(~{
        cal <- calib_slope_intercept(.x$y, .x$p)
        tibble(
          ece_10 = ece_score(.x$y, .x$p, bins = 10),
          slope = cal$slope,
          intercept = cal$intercept
        )
      }) %>%
      ungroup()
    
    fold_metrics2 <- fold_metrics %>%
      left_join(split_calib, by = c("split_id","strategy"))
    
    # ----------------------- TABLE 3: PERFORMANCE SUMMARY -----------------------
    perf_summary <- fold_metrics2 %>%
      group_by(strategy, model, calibration) %>%
      summarise(
        n_splits = n(),
        roc = list(mean_se(roc_auc)),
        pr  = list(mean_se(pr_auc)),
        ll  = list(mean_se(log_loss)),
        br  = list(mean_se(brier)),
        sl  = list(mean_se(slope)),
        it  = list(mean_se(intercept)),
        ec  = list(mean_se(ece_10)),
        roc_auc_mean = roc[[1]]$mean, roc_auc_se = roc[[1]]$se,
        pr_auc_mean  = pr[[1]]$mean,  pr_auc_se  = pr[[1]]$se,
        log_loss_mean = ll[[1]]$mean, log_loss_se = ll[[1]]$se,
        brier_mean = br[[1]]$mean, brier_se = br[[1]]$se,
        slope_mean = sl[[1]]$mean, slope_se = sl[[1]]$se,
        intercept_mean = it[[1]]$mean, intercept_se = it[[1]]$se,
        ece10_mean = ec[[1]]$mean, ece10_se = ec[[1]]$se,
        .groups = "drop"
      ) %>%
      select(-roc,-pr,-ll,-br,-sl,-it,-ec) %>%
      arrange(desc(roc_auc_mean), log_loss_mean)
    
    safe_write_csv(perf_summary, file.path(TAB_DIR, "table_03_performance_summary.csv"))
    safe_write_html(perf_summary, file.path(TAB_DIR, "table_03_performance_summary.html"),
                    caption = "Table 3. Discrimination, calibration, and scoring-rule performance by strategy (mean ± SE across outer splits).",
                    digits = 6)
    
    # ----------------------- TABLE 4: STABILITY SUMMARY -------------------------
    stab_summary <- tibble()
    if (nrow(stab) > 0) {
      stab_summary <- stab %>%
        group_by(strategy, model, split_id) %>%
        summarise(
          mean_jaccard = first(mean_jaccard_subsample),
          mean_k_selected = first(mean_k_selected),
          k_stable = sum(selection_freq >= STABILITY_PI_THR, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        group_by(strategy, model) %>%
        summarise(
          mean_jaccard = mean(mean_jaccard, na.rm = TRUE),
          mean_k_selected = mean(mean_k_selected, na.rm = TRUE),
          mean_k_stable = mean(k_stable, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        arrange(desc(mean_jaccard), mean_k_selected)
      
      safe_write_csv(stab_summary, file.path(TAB_DIR, "table_04_stability_summary.csv"))
      safe_write_html(stab_summary, file.path(TAB_DIR, "table_04_stability_summary.html"),
                      caption = paste0("Table 4. Stability-selection summaries for sparse strategies (BOOT_B=", BOOT_B,
                                       ", PI_THR=", STABILITY_PI_THR, ")."),
                      digits = 6)
    }
    
    # ----------------------- SELECT BEST STRATEGY -------------------------------
    best_roc <- max(perf_summary$roc_auc_mean, na.rm = TRUE)
    best_row <- perf_summary %>% filter(roc_auc_mean == best_roc) %>% slice(1)
    best_roc_se <- best_row$roc_auc_se
    
    candidate <- perf_summary %>%
      mutate(within_1se = if (ONE_SE_RULE) roc_auc_mean >= (best_roc - best_roc_se) else TRUE) %>%
      filter(within_1se) %>%
      arrange(log_loss_mean, ece10_mean, brier_mean)
    
    candidate2 <- candidate %>%
      left_join(stab_summary %>% select(strategy, mean_jaccard, mean_k_selected), by = "strategy") %>%
      mutate(
        mean_jaccard = if_else(is.na(mean_jaccard), -Inf, mean_jaccard),
        mean_k_selected = if_else(is.na(mean_k_selected), Inf, mean_k_selected)
      ) %>%
      arrange(log_loss_mean, ece10_mean, desc(mean_jaccard), mean_k_selected)
    
    BEST_STRATEGY <- candidate2$strategy[1]
    log_line("Selected BEST_STRATEGY: ", BEST_STRATEGY)
    
    # ----------------------- TABLE 5: PAIRWISE VS BEST (OOF) --------------------
    oof_wide <- oof %>%
      mutate(strategy = as.character(strategy)) %>%
      group_by(split_id, id, y, strategy) %>%
      summarise(p = mean(p), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = strategy, values_from = p)
    
    if (!(BEST_STRATEGY %in% names(oof_wide))) {
      die(paste0("BEST_STRATEGY column not found in oof_wide: ", BEST_STRATEGY))
    }
    
    y_all <- oof_wide$y
    p_best <- oof_wide[[BEST_STRATEGY]]
    ok_best <- is.finite(y_all) & is.finite(p_best)
    if (length(unique(y_all[ok_best])) < 2) die("OOF pooled truth is degenerate; cannot compute pooled metrics.")
    
    roc_best <- roc_auc_vec(y_all[ok_best], p_best[ok_best])
    ll_best  <- log_loss_vec(y_all[ok_best], p_best[ok_best])
    
    other_strats <- setdiff(names(oof_wide), c("split_id","id","y", BEST_STRATEGY))
    
    pairwise <- tibble(strategy = other_strats) %>%
      mutate(
        roc_auc_best = roc_best,
        logloss_best = ll_best,
        roc_auc_this = map_dbl(strategy, ~{
          p_this <- oof_wide[[.x]]
          ok <- ok_best & is.finite(p_this)
          if (length(unique(y_all[ok])) < 2) return(NA_real_)
          roc_auc_vec(y_all[ok], p_this[ok])
        }),
        logloss_this = map_dbl(strategy, ~{
          p_this <- oof_wide[[.x]]
          ok <- ok_best & is.finite(p_this)
          log_loss_vec(y_all[ok], p_this[ok])
        }),
        roc_auc_diff = roc_auc_this - roc_auc_best,
        logloss_diff = logloss_this - logloss_best
      ) %>%
      arrange(logloss_diff)
    
    safe_write_csv(pairwise, file.path(TAB_DIR, "table_05_pairwise_vs_best.csv"))
    safe_write_html(pairwise, file.path(TAB_DIR, "table_05_pairwise_vs_best.html"),
                    caption = "Table 5. Pairwise performance deltas versus the selected best strategy (OOF pooled).",
                    digits = 6)
    
    # ----------------------- TABLE 6: FEATURE STABILITY PROFILE -----------------
    final_feat <- tibble(feature = feat_names, selection_freq_mean = NA_real_, stable = NA)
    if (nrow(stab) > 0 && (str_detect(BEST_STRATEGY, "^lasso") || str_detect(BEST_STRATEGY, "^enet"))) {
      best_raw <- paste0(str_extract(BEST_STRATEGY, "^(ridge|lasso|enet)"), "_raw")
      final_feat <- stab %>%
        filter(strategy == best_raw) %>%
        group_by(feature) %>%
        summarise(selection_freq_mean = mean(selection_freq, na.rm = TRUE), .groups = "drop") %>%
        mutate(stable = selection_freq_mean >= STABILITY_PI_THR) %>%
        arrange(desc(selection_freq_mean))
    }
    safe_write_csv(final_feat, file.path(TAB_DIR, "table_06_feature_stability_profile.csv"))
    safe_write_html(final_feat, file.path(TAB_DIR, "table_06_feature_stability_profile.html"),
                    caption = paste0("Table 6. Feature stability profile for best sparse family (PI_THR=", STABILITY_PI_THR, ")."),
                    digits = 6)
    
    # ------------------------------ FIGURES -------------------------------------
    top3 <- candidate2 %>% slice_head(n = 3) %>% pull(strategy)
    show_strats <- unique(c(BEST_STRATEGY, top3))
    
    oof_top <- oof %>%
      filter(strategy %in% show_strats) %>%
      mutate(
        strategy_label = as.character(strategy)
      )
    
    # Figure 1: ROC curves (pROC points; no yardstick API risk)
    safe_plot({
      roc_df <- oof_top %>%
        group_by(strategy_label) %>%
        group_modify(~{
          r <- pROC::roc(response = .x$y, predictor = .x$p, quiet = TRUE, direction = "<")
          cc <- pROC::coords(r, x = "all", ret = c("specificity","sensitivity"), transpose = FALSE)
          tibble(
            specificity = as.numeric(cc[,"specificity"]),
            sensitivity = as.numeric(cc[,"sensitivity"]),
            strategy_label = .y$strategy_label
          )
        }) %>%
        ungroup()
      
      p1 <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity, color = strategy_label)) +
        geom_path(linewidth = 1.15, alpha = 0.95) +
        geom_abline(linetype = "dashed") +
        coord_equal() +
        labs(
          title = "Figure 1. ROC curves (pooled out-of-fold predictions)",
          subtitle = paste0("Repeated nested CV: ", OUTER_FOLDS, "-fold × ", N_REPEATS, " repeats"),
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
    }, "figure_01_roc")
    
    # Figure 2: PR curves (manual PR curve; no yardstick API risk)
    safe_plot({
      pr_df <- oof_top %>%
        group_by(strategy_label) %>%
        group_modify(~{
          pr <- pr_curve_manual(.x$y, .x$p)
          pr %>% mutate(strategy_label = .y$strategy_label)
        }) %>%
        ungroup()
      
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
    }, "figure_02_pr")
    
    # Figure 3: Calibration curves (fix group_modify scoping: use .y keys)
    safe_plot({
      calib_long <- oof_top %>%
        group_by(strategy_label) %>%
        group_modify(~{
          bins <- calib_bins(.x$y, .x$p, bins = 10)
          bins %>% mutate(strategy_label = .y$strategy_label)
        }) %>%
        ungroup()
      
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
    }, "figure_03_calibration")
    
    # Figure 4: Decision curve analysis
    safe_plot({
      thr_grid <- seq(0.01, 0.80, by = 0.01)
      dca <- oof_top %>%
        group_by(strategy_label) %>%
        summarise(dca = list(tibble(
          threshold = thr_grid,
          net_benefit = map_dbl(thr_grid, ~net_benefit(y, p, .x))
        )), .groups = "drop") %>%
        unnest(dca)
      
      prev <- mean(wdbc$y == 1)
      dca_ref <- tibble(
        threshold = thr_grid,
        treat_none = 0,
        treat_all = prev - (1 - prev) * (threshold / (1 - threshold))
      ) %>%
        pivot_longer(cols = c(treat_none, treat_all), names_to = "ref", values_to = "net_benefit")
      
      p4 <- ggplot() +
        geom_line(data = dca, aes(x = threshold, y = net_benefit, color = strategy_label),
                  linewidth = 1.15, alpha = 0.95) +
        geom_line(data = dca_ref, aes(x = threshold, y = net_benefit, linetype = ref),
                  linewidth = 0.9) +
        scale_linetype_manual(values = c(treat_none = "dashed", treat_all = "dotdash")) +
        labs(
          title = "Figure 4. Decision curve analysis (net benefit across thresholds)",
          subtitle = "Clinical utility across decision thresholds; higher net benefit is better",
          x = "Decision threshold (predicted risk)",
          y = "Net benefit",
          color = "Strategy",
          linetype = "Reference"
        ) +
        theme_minimal(base_size = 13) +
        theme(plot.title = element_text(face = "bold", size = 15),
              axis.title = element_text(face = "bold"),
              legend.title = element_text(face = "bold"),
              panel.grid.minor = element_blank())
      
      ggsave(file.path(FIG_DIR, "figure_04_decision_curve.png"), p4, width = FIG_W, height = FIG_H, dpi = DPI, bg = "white")
    }, "figure_04_decision_curve")
    
    # Figure 5: Stability heatmap
    safe_plot({
      stab_heat <- tibble()
      if (nrow(stab) > 0) {
        stab_heat <- stab %>%
          filter(strategy %in% c("lasso_raw","enet_raw")) %>%
          group_by(strategy, feature) %>%
          summarise(selection_freq = mean(selection_freq, na.rm = TRUE), .groups = "drop") %>%
          group_by(strategy) %>%
          slice_max(order_by = selection_freq, n = 15, with_ties = FALSE) %>%
          ungroup() %>%
          mutate(strategy_label = strategy)
      }
      if (nrow(stab_heat) == 0) {
        stab_heat <- tibble(strategy_label = "(no stability available)", feature = feat_names[1:10], selection_freq = 0)
      }
      
      p5 <- ggplot(stab_heat, aes(x = strategy_label, y = fct_reorder(feature, selection_freq), fill = selection_freq)) +
        geom_tile(color = "white", linewidth = 0.35) +
        scale_fill_viridis_c(option = "C", end = 0.98, labels = percent_format(accuracy = 1)) +
        labs(
          title = "Figure 5. Feature selection stability (top features)",
          subtitle = paste0("Mean selection frequency across splits; BOOT_B=", BOOT_B, " per split (sparse families only)"),
          x = NULL, y = NULL, fill = "Selection\nfrequency"
        ) +
        theme_minimal(base_size = 13) +
        theme(plot.title = element_text(face = "bold", size = 15),
              axis.text.x = element_text(angle = 20, hjust = 1),
              legend.title = element_text(face = "bold"),
              panel.grid = element_blank())
      
      ggsave(file.path(FIG_DIR, "figure_05_stability_heatmap.png"), p5, width = FIG_W, height = FIG_H, dpi = DPI, bg = "white")
    }, "figure_05_stability_heatmap")
    
    # ------------------------------ MANIFEST + SUMMARY ---------------------------
    manifest <- list(
      timestamp_local = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      data_dir = DATA_DIR,
      out_dir = OUT_DIR,
      config = list(
        OUTER_FOLDS = OUTER_FOLDS, N_REPEATS = N_REPEATS, INNER_FOLDS = INNER_FOLDS,
        CALIB_CV = CALIB_CV, BOOT_B = BOOT_B, PI_THR = STABILITY_PI_THR,
        ONE_SE_RULE = ONE_SE_RULE, SEED = SEED, QUICK = ARGS$quick
      ),
      best_strategy = BEST_STRATEGY,
      r_version = R.version.string,
      platform = R.version$platform,
      session_info = capture.output(sessionInfo()),
      packages = lapply(PKGS, function(p) {
        v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
        list(name = p, version = v)
      })
    )
    
    writeLines(jsonlite::toJSON(manifest, pretty = TRUE, auto_unbox = TRUE),
               file.path(AUX_DIR, "run_manifest.json"))
    
    summary_txt <- paste0(
      "RUN SUMMARY — WDBC JDS Pipeline (R-ONLY)\n",
      "Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n",
      "DATA_DIR: ", DATA_DIR, "\n",
      "OUT_DIR:  ", OUT_DIR, "\n",
      "LOG_FILE: ", LOG_FILE, "\n\n",
      "BEST_STRATEGY: ", BEST_STRATEGY, "\n\n",
      "Key outputs:\n",
      "  Figures: ", FIG_DIR, "\n",
      "    figure_01_roc.png\n",
      "    figure_02_pr.png\n",
      "    figure_03_calibration.png\n",
      "    figure_04_decision_curve.png\n",
      "    figure_05_stability_heatmap.png\n\n",
      "  Tables:  ", TAB_DIR, "\n",
      "    table_01_data_summary.(csv/html)\n",
      "    table_02_feature_summary.(csv/html)\n",
      "    table_03_performance_summary.(csv/html)\n",
      "    table_04_stability_summary.(csv/html) [if sparse models]\n",
      "    table_05_pairwise_vs_best.(csv/html)\n",
      "    table_06_feature_stability_profile.(csv/html)\n\n",
      "  Aux:     ", AUX_DIR, "\n",
      "    oof_predictions.csv\n",
      "    fold_metrics.csv\n",
      "    stability_selection.csv\n",
      "    run_manifest.json\n",
      "    wdbc_names_excerpt.txt\n"
    )
    writeLines(summary_txt, file.path(OUT_DIR, "RUN_SUMMARY.txt"))
    
    log_line("Pipeline completed successfully.")
    log_line("See RUN_SUMMARY.txt and PIPELINE_LOG for details.")
  })
}

# ------------------------------ EXECUTE ----------------------------------------
tryCatch(
  main(),
  error = function(e) {
    log_line("FATAL (top-level): ", safe_condition_message(e))
    log_line("HINT: inspect LOG_FILE above. If packages updated mid-session, restart R and rerun.")
    if (!interactive()) quit(status = 1)
    stop(safe_condition_message(e), call. = FALSE)
  }
)
