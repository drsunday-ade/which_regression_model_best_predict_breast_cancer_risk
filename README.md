# which_regression_model_best_predict_breast_cancer_risk
Which regression model best predicts breast cancer—and can we trust its risk scores? End-to-end reproducible pipeline comparing ridge/lasso/elastic net with probability calibration, decision-curve analysis, and bootstrap stability selection using repeated nested cross-validation.

# Which Model Best Predicts Breast Cancer—and Can We Trust Its Risk Scores?
**A fully reproducible benchmark of penalized logistic regression (ridge, lasso, elastic net) with calibration and stability selection on the WDBC dataset**

**Author / Corresponding:** Sunday A. Adetunji  
Department of Biostatistics, Oregon State University  
ORCID: https://orcid.org/0000-0001-9321-9957  
Email: adetunjs@oregonstate.edu  

---

## 1. Project summary

This repository is the companion reproducibility package for a methods-focused study that compares **penalized logistic regression strategies** for binary risk prediction under **high collinearity** and **grouped feature structure**. Using the **Wisconsin Diagnostic Breast Cancer (WDBC)** dataset (569 breast masses; 30 continuous imaging-derived predictors; malignant prevalence 37.3%), we evaluate which modeling strategy provides the best **joint** performance across:

- **Discrimination** (ranking performance): ROC-AUC, PR-AUC  
- **Probability quality / calibration** (risk score reliability): log loss, Brier score, calibration curves, expected calibration error (ECE)  
- **Interpretability and robustness**: sparsity and **feature-selection stability** (bootstrap stability selection; mean Jaccard similarity)

The repository contains all **scripts, out-of-fold predictions, fold-level metrics, tables, and figures** needed to reproduce the Results, Discussion, and Reproducibility Appendix.

> **Important:** This is an educational/methodological benchmark using a public dataset. It is **not** a clinical decision support tool.

---

## 2. Key findings (high level)

Across repeated nested stratified cross-validation (10-fold outer CV × 10 repeats = 100 out-of-sample splits), the **elastic net + sigmoid (Platt) calibration** achieved the best overall trade-off between near-ceiling discrimination and improved probability reliability, while maintaining stable and interpretable selection under correlated predictors.

A concise snapshot (see `tables/table_03_performance_summary.*` and Figures 1–4):

- Best overall strategy: **enet_sigmoid**
- Strong discrimination and improved calibration vs. uncalibrated enet
- Stability: elastic net selections more reproducible than lasso (see Table 4 and Figure 5)

---

## 3. Repository structure

.
├── aux/
│ ├── fold_metrics.csv # Fold-level metrics for every strategy and outer split
│ ├── oof_predictions.csv # Pooled out-of-fold predictions (per sample × strategy)
│ ├── stability_selection.csv # Bootstrap stability-selection results (sparse families)
│ ├── run_manifest.json # Run metadata: seeds, versions, settings, file manifest
│ └── wdbc_names_excerpt.txt # Dataset field/feature documentation excerpt used by pipeline
│
├── data/
│ ├── wdbc.data # WDBC raw data file (as distributed)
│ └── wdbc.names # WDBC metadata / feature names (as distributed)
│
├── figures/
│ ├── figure_01_roc.png # Figure 1: ROC curves (pooled out-of-fold predictions)
│ ├── figure_02_pr.png # Figure 2: Precision–Recall curves (pooled OOF)
│ ├── figure_03_calibration.png # Figure 3: Calibration (reliability) curves (10 bins)
│ ├── figure_04_decision_curve.png # Figure 4: Decision curve analysis (net benefit)
│ └── figure_05_stability_heatmap.png # Figure 5: Feature selection stability heatmap
│
├── scripts/
│ ├── wdbc_jds_end_to_end_bulletproof.R # Main end-to-end pipeline: data → models → outputs
│ └── render_figures_01_03_only.R # Optional: re-render Figures 1–3 only
│
├── tables/
│ ├── table_01_data_summary.{html,csv} # Table 1: data summary
│ ├── table_02_feature_summary.{html,csv} # Table 2: feature summary
│ ├── table_03_performance_summary.{html,csv} # Table 3: performance summary
│ ├── table_04_stability_summary.{html,csv} # Table 4: stability summary
│ ├── table_05_pairwise_vs_best.{html,csv} # Table 5: paired comparisons vs best strategy
│ └── table_06_feature_stability_profile.{html,csv} # Table 6: feature-level stability profiles
│
├── methods_section.Rmd # Manuscript Methods section (ready to knit)
├── appendix_reproducibility.Rmd # Reproducibility appendix subsection (ready to knit)
├── PIPELINE_LOG_*.txt # Detailed run log (timestamped)
└── RUN_SUMMARY.txt # Human-readable run summary

swift
Copy code

---

## 4. Reproducibility: how to run

### 4.1. Requirements

- R ≥ 4.2 (recommended)
- R packages (typical set): `tidyverse`, `glmnet`, `rsample`, `yardstick`, `pROC`, `PRROC`,
  `ggplot2`, `dplyr`, `purrr`, `readr`, `tibble`, `jsonlite`, `knitr`, `rmarkdown`

Install packages (adjust to your environment):

```r
install.packages(c(
  "tidyverse","glmnet","rsample","yardstick","pROC","PRROC",
  "jsonlite","knitr","rmarkdown"
))
4.2. Run the full end-to-end pipeline (recommended)
From the repository root:

r
Copy code
source("scripts/wdbc_jds_end_to_end_bulletproof.R")
This script will:

Load and clean WDBC data (ID excluded; 30 standardized predictors)

Run repeated nested stratified CV with training-only preprocessing and calibration

Fit ridge/lasso/elastic net strategies with tuning and calibration variants

Run bootstrap stability selection for sparse families (as configured)

Generate:

Figures in figures/

Tables in tables/ (HTML + CSV)

Reproducibility artifacts in aux/ (OOF predictions, fold metrics, manifest)

4.3. Render/refresh select figures (optional)
r
Copy code
source("scripts/render_figures_01_03_only.R")
4.4. Knit manuscript components
Methods section:

r
Copy code
rmarkdown::render("methods_section.Rmd")
Reproducibility appendix:

r
Copy code
rmarkdown::render("appendix_reproducibility.Rmd")
5. What the aux/ folder is (and why it matters)
The aux/ directory contains analysis-grade, provenance-preserving outputs that support:

Auditability (exact split-level metrics and predictions)

Re-analysis (users can recompute any metric, add new metrics, or re-plot curves)

Reproducibility review (a reviewer can verify claims without re-running the full pipeline)

Recommended practice: keep aux/ in the public repository, but treat it as derived artifacts.
In the manuscript, reference it as a reproducibility supplement (see Section 8 below).

6. Outputs: tables and figures (what they represent)
Figures 1–2: Discrimination curves using pooled out-of-fold predictions across all outer splits.

Figure 3: Calibration (reliability) curves (10 bins), with bin size reflected by point size.

Figure 4: Decision curve analysis showing clinical net benefit across thresholds (vs treat-all / treat-none).

Figure 5: Feature selection stability heatmap comparing selection frequency across splits for sparse methods.

Tables provide structured, publication-ready summaries:

Table 1: Dataset composition and outcome prevalence

Table 2: Feature groups/structure and summary statistics

Table 3: Primary performance metrics by strategy (mean ± SE)

Table 4: Stability metrics by strategy (e.g., mean Jaccard similarity)

Table 5: Paired comparisons vs best strategy (OOF paired evaluation)

Table 6: Feature-level stability profile (selection frequencies)

7. How to cite this repository
If you use this code or outputs, please cite the repository and the associated manuscript.

Suggested citation (plain text):

Adetunji, S. A. (2026). Which Model Best Predicts Breast Cancer—and Can We Trust Its Risk Scores? A reproducible comparison of penalized logistic regression methods (code and reproducibility materials). GitHub repository.

BibTeX (edit as needed):

bibtex
Copy code
@misc{adetunji_wdbc_penalized_2026,
  author = {Adetunji, Sunday A.},
  title  = {Which Model Best Predicts Breast Cancer---and Can We Trust Its Risk Scores? A reproducducible comparison of penalized logistic regression methods},
  year   = {2026},
  howpublished = {GitHub repository},
  url    = {https://github.com/drsunday-ade/which_regression_model_best_predict_breast_cancer_risk}
}
8. Suggested manuscript reproducibility statement
A concise manuscript-ready statement (adapt to journal format):

“All code, derived outputs (out-of-fold predictions, fold-level metrics, and stability-selection artifacts), and manuscript-ready tables/figures are publicly available in the companion GitHub repository. The repository includes a run manifest capturing seeds, configuration, and software versions to enable exact reproduction of the reported results.”

9. License and data note
Code: add a LICENSE file appropriate to your preference (MIT/Apache-2.0 are common for academic code).

Data: WDBC is a public dataset distributed via the UCI Machine Learning Repository. This repository includes the original wdbc.data and wdbc.names files for reproducibility; please ensure downstream reuse complies with the dataset’s stated terms/attribution in the UCI source.

10. Contact
For questions, issues, or reproducibility reports:

Sunday A. Adetunji — adetunjs@oregonstate.edu

ORCID: https://orcid.org/0000-0001-9321-9957

perl
Copy code
::contentReference[oaicite:0]{index=0}
 ​:contentReference[oaicite:1]{index=1}​






