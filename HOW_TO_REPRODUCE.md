# How to Reproduce All Results

This document provides step-by-step instructions to reproduce every table, figure, and numerical result reported in the manuscript and Online Appendix.

---

## Prerequisites

### R Installation

R version 4.3.0 or later is required. Download from [https://cran.r-project.org/](https://cran.r-project.org/).

### Install Required Packages

Open R and run:

```r
install.packages(c(
  "cSEM", "rsample", "arrow", "dplyr", "purrr",
  "tidyr", "tibble", "readxl", "readr", "glue", "rlang",
  "future", "furrr", "parallelly", "progressr"
))
```

**Note:** Manuscript figures use base R graphics (no ggplot2 dependency). The `MASS` package (for `mvrnorm`) is included with base R.

### Working Directory

**Important:** All scripts use relative paths (e.g., `../outputs_main/`). You must `cd` into the appropriate `code_*/` directory before running each pipeline — **not** the repository root.

```bash
cd /path/to/TC-PVSS-PLSSEM/code_main       # for Pipeline A
cd /path/to/TC-PVSS-PLSSEM/code_figures     # for Pipeline C (figures)
cd /path/to/TC-PVSS-PLSSEM/code_simulation  # for Pipeline B (simulation)
```

Or in R:

```r
setwd("/path/to/TC-PVSS-PLSSEM/code_main")
```

---

## Part A: ACSI Empirical Analysis (Sections 4–5, Tables 4–6, Figures 1–5)

Scripts must be run in order. Each script reads outputs from the previous step.

### Step 0: Prepare Data

```bash
cd code_main
Rscript 00_prepare_data.R
```

**What it does:** Reads `inputs_main/ACSI Data 2015 ver2.xlsx`, recodes invalid values (>10) to NA, selects the 12 indicators used in the analysis, and saves the cleaned dataset.

**Output:** `outputs_main/acsi_clean.rds`

**Runtime:** ~2 seconds

### Step 1: Measurement Model Audit

```bash
Rscript 01_measurement_audit.R
```

**What it does:** Reads the cleaned data, applies listwise deletion for complete cases, and estimates the M0 baseline PLS-PM model using cSEM. Extracts indicator loadings, composite reliability (ρc), Cronbach's α, rho_A, AVE, HTMT, and Fornell-Larcker matrices. Reports threshold checks with pass/fail indicators.

**Outputs:**
- `outputs_main/measurement_audit_table_A0.csv` → Online Appendix Table A0
- `outputs_main/htmt_matrix.csv` → HTMT discriminant validity
- `outputs_main/fornell_larcker.csv` → Fornell-Larcker criterion
- Console output includes suggested text for manuscript Section 4.1

**Runtime:** ~10 seconds

**Note:** This script uses listwise deletion (cSEM requires complete data). The nested CV pipeline (Step 2) uses the full sample with per-fold mean imputation. Both sample sizes are reported for transparency.

### Step 2: Nested Cross-Validation with lm Benchmark

```bash
Rscript 02_nested_cv_search.R
```

**What it does:** Runs the full TC-PVSS nested cross-validation pipeline:
- Creates 50 outer splits (10-fold × 5 repeats, stratified by INDUSTRY)
- For each outer split: selects the best model via inner CV (5-fold × 2 repeats), refits on outer-train, predicts on outer-test
- Also fits M0 baseline and indicator-level lm benchmark per outer split
- Records per-case predictions, squared errors, and absolute errors for all three comparators

**Outputs:**
- `outputs_main/pred_errors_long.parquet` → Per-case prediction errors (all 41,195 predictions)
- `outputs_main/outer_results.parquet` → Per-split summary metrics

**Runtime:** ~2–4 hours (depending on hardware)

**Global seed:** `20260220` (set at line 33 of script)

### Step 3: Summarize, Test, and Visualize

```bash
Rscript 03_summarize_and_cvpattest.R
```

**What it does:**
- Computes pooled and per-split predictive performance (selected vs M0 vs lm)
- Computes model selection frequencies and edge inclusion probabilities
- Runs CVPAT-style paired tests at both split-level (primary, n=50) and per-case (supplementary)
- Generates by-industry performance breakdown
- Creates diagnostic figures (loss difference histograms, RMSE/MAE boxplots)

**Outputs:**
- `outputs_main/summary_tables.csv` → Table 5 source data
- `outputs_main/model_frequency.csv` → Figure 1 source data
- `outputs_main/performance_by_outer_split.csv` → Figures 3–4 source data
- `outputs_main/performance_by_industry.csv` → Table 6, Figure 5 source data
- `outputs_main/fig_loss_diff_mse_hist.png` → Diagnostic: MSE loss differences
- `outputs_main/fig_loss_diff_ae_hist.png` → Diagnostic: AE loss differences
- `outputs_main/fig_rmse_boxplot_outer.png` → Diagnostic: RMSE boxplot (selected vs M0 vs lm)
- `outputs_main/fig_mae_boxplot_outer.png` → Diagnostic: MAE boxplot (selected vs M0 vs lm)

**Runtime:** ~5 seconds

---

## Part B: Manuscript Figures (Figures 1–5, S1–S4)

### Main Text Figures

```bash
cd code_figures
Rscript 01_make_figures_main.R
```

**Reads:** CSV files from `../outputs_main/`

**Writes:** `../figures_main/Fig1_*.png`, `Fig1_*.tif`, ..., `Fig5_*.png`, `Fig5_*.tif` (600 DPI)

| Figure | Description | Source data |
|--------|-------------|-------------|
| Figure 1 | Model selection frequency across outer splits | `model_frequency.csv` (fallback from `model_frequency_by_outer_split.csv`) |
| Figure 2 | Optional-path inclusion probabilities | `summary_tables.csv` |
| Figure 3 | RMSE distribution (selected vs baseline) | `performance_by_outer_split.csv` |
| Figure 4 | MAE distribution (selected vs baseline) | `performance_by_outer_split.csv` |
| Figure 5 | RMSE improvement by industry | `performance_by_industry.csv` |

### Online Appendix Figures

```bash
Rscript 02_make_figures_simulation.R
```

**Reads:** CSV files from `../sim_outputs/`

**Writes:** `../figures_simulation/FigS1_*.png`, ..., `FigS4_*.tif` (600 DPI)

| Figure | Description | Source file |
|--------|-------------|-------------|
| Figure S1 | Absolute miscalibration by sample size | `Table_S2_optimism.csv` |
| Figure S2 | False-positive rate by sample size (modal evaluation) | `sim_main.csv` |
| Figure S3 | Stable-core correctness vs τ | `Table_S4_stable_core_tau.csv` |
| Figure S4 | Stable-core true RMSE vs τ | `Table_S4_stable_core_tau.csv` |

---

## Part C: Simulation Study (Section 6, Tables S1–S4, Figures S1–S4)

```bash
# WARNING: Full simulation requires ~8–24 hours depending on hardware
cd code_simulation
Rscript -e 'source("00_config.R"); source("01_core_functions.R"); source("02_run_simulation_parallel.R")'
Rscript 03_summarize_simulation.R
```

The simulation uses checkpointing via chunk parquet files in `sim_outputs/chunks/`. If interrupted, re-running will resume from the last completed chunk. All random seeds are fixed (global seed: `20260301`).

**Design:** 3 (N) × 2 (λ) × 2 (ρ) × 3 (Truth) × 200 replications = 7,200 total runs.

**Outputs (all in `../sim_outputs/`):**
- `sim_main.csv` / `.parquet` — Main results (7,200 rows)
- `sim_stable.csv` / `.parquet` — Stability-augmented results (τ sensitivity)
- `Table_S2_optimism.csv` — Optimism bias by condition (36 rows)
- `Table_S3_selection_accuracy.csv` — TPR/FPR/F1 by condition (36 rows)
- `Table_S4_stable_core_tau.csv` — τ sensitivity (720 rows)
- `run_log.txt` — Timestamps and run parameters

Pre-computed simulation outputs are included in the repository so that figures and tables can be reproduced without re-running the simulation.

---

## Mapping: Script Outputs → Manuscript Tables and Figures

| Manuscript element | Source script | Source file |
|--------------------|---------------|-------------|
| Table 4 (construct-indicator mapping) | — | Described in manuscript |
| Table 5 (pooled performance + CVPAT) | `03_summarize_and_cvpattest.R` | `outputs_main/summary_tables.csv` |
| Table 6 (industry-level performance) | `03_summarize_and_cvpattest.R` | `outputs_main/performance_by_industry.csv` |
| Table A0 (measurement audit) | `01_measurement_audit.R` | `outputs_main/measurement_audit_table_A0.csv` |
| Figures 1–5 | `01_make_figures_main.R` | `figures_main/` |
| Figures S1–S4 | `02_make_figures_simulation.R` | `figures_simulation/` |
| Table S2 (optimism bias) | `03_summarize_simulation.R` | `sim_outputs/Table_S2_optimism.csv` |
| Table S3 (selection accuracy) | `03_summarize_simulation.R` | `sim_outputs/Table_S3_selection_accuracy.csv` |
| Table S4 (stable-core τ) | `03_summarize_simulation.R` | `sim_outputs/Table_S4_stable_core_tau.csv` |

---

## Troubleshooting

**cSEM predict() output format varies across versions.** Script `02_nested_cv_search.R` includes a robust extraction function (`get_pred_vec()`) that handles multiple output formats. If prediction extraction fails, check your cSEM version (`packageVersion("cSEM")`); version ≥ 0.5.0 is required.

**Memory for parquet files.** The per-case prediction file (`pred_errors_long.parquet`) contains ~41,000 rows × 18 columns. The `arrow` package handles this efficiently, but systems with <4 GB RAM may need to increase available memory.

**Simulation checkpointing.** If the simulation is interrupted, simply re-run the same script. Completed cells are detected and skipped automatically.

---

## Contact

[Removed for double-blind review]
