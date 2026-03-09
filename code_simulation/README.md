# TC-PVSS Simulation B — FULL + Parallel + Checkpoints

This package is optimized for:
- **Parallel** execution on Windows using `future::multisession`
- **Checkpointing** (chunk parquet files written as the run proceeds)
- **Resume** support (skips chunks that already exist)
- Paper-ready summary tables for **optimism bias** and **selection accuracy**

## Run (exactly 4 sources)
```r
source("00_config.R")
source("01_core_functions.R")
source("02_run_simulation_parallel.R")
source("03_summarize_simulation.R")
```

## Output (written to ../sim_outputs/ when run from code_simulation/)
- sim_main.parquet / sim_main.csv
- sim_stable.parquet / sim_stable.csv
- Table_S2_optimism.csv
- Table_S3_selection_accuracy.csv
- Table_S4_stable_core_tau.csv
- run_log.txt (timestamps + parameters)
- chunks/ (checkpoint parquets, intermediate)

## FULL vs PILOT
Default is FULL (as requested). To run PILOT instead:
```r
Sys.setenv(TCPVSS_MODE="pilot")   # do this BEFORE source("00_config.R")
```

## Design notes

- Naïve CV uses the same fold count (10-fold) as the nested CV outer
  loop (NAIVE_V = OUTER_V), ensuring fair comparison.
- The loyalty construct (L) is measured by 3 reflective indicators
  (REPUR, HIGHPTOL, LOWPTOL) in the simulation DGP. This differs from
  the empirical illustration (single-item REPUR) to support a complete
  measurement model amenable to loading and reliability manipulation.
- Nested pipeline evaluation uses frequency-weighted true RMSE/MAE
  across all uniquely selected models, plus per-fold averaged TPR/FPR/F1.
  Modal-model metrics are retained as `_modal_` columns for comparison.

## Packages
```r
install.packages(c("dplyr","purrr","tidyr","rsample","cSEM","arrow",
                   "future","furrr","parallelly","progressr","glue","readr"))
```
`MASS` is used via `MASS::mvrnorm` (install.packages("MASS") if needed).
