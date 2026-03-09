# 02_run_simulation_parallel.R
# ============================================================
# FULL parallel runner + checkpoints + resume
# Writes chunk checkpoints to ../sim_outputs/chunks, then stitches.
#
# Nested pipeline evaluation uses:
#   - frequency-weighted true RMSE/MAE (pipeline-consistent)
#   - per-fold averaged TPR/FPR/F1
# Also retains modal-model metrics as _modal_ columns for comparison.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(purrr); library(tidyr); library(arrow)
  library(future); library(furrr); library(parallelly); library(progressr)
})

source("00_config.R")
source("01_core_functions.R")

dir.create("../sim_outputs", showWarnings = FALSE)
dir.create("../sim_outputs/chunks", recursive = TRUE, showWarnings = FALSE)

# Log file
log_path <- file.path("../sim_outputs", "run_log.txt")
cat("", file=log_path)  # reset (append later)

log_to_file <- function(...) {
  line <- paste0(timestamp(), " | ", paste0(..., collapse=""))
  cat(line, "\n")
  cat(line, "\n", file=log_path, append=TRUE)
}

# Build grid
grid <- tidyr::expand_grid(
  N = N_list,
  lambda = lambda_list,
  rho = rho_list,
  truth = truth_list,
  rep = seq_len(R_reps)
)

log_to_file("Grid size = ", nrow(grid), " cells; MODE=", MODE)
log_to_file("Budgets: OUTER_V=", OUTER_V, " OUTER_REPEATS=", OUTER_REPEATS,
            " INNER_V=", INNER_V, " INNER_REPEATS=", INNER_REPEATS)
log_to_file("NAIVE_V=", NAIVE_V, " NAIVE_REPEATS=", NAIVE_REPEATS)
log_to_file("N_test=", N_test, "; tau_list=", paste(tau_list, collapse=","))
log_to_file("CHUNK_SIZE=", CHUNK_SIZE, "; RESUME=", RESUME)

# Choose workers
workers <- min(parallelly::availableCores() - 2, MAX_WORKERS)
if (workers < 1) workers <- 1
log_to_file("Workers=", workers)

future::plan(future::multisession, workers=workers)
progressr::handlers("txtprogressbar")

# Helper: determine which chunk ids already exist
existing_chunks <- function() {
  files <- list.files("../sim_outputs/chunks", pattern="^main_\\d{6}_\\d{6}\\.parquet$", full.names=FALSE)
  if (length(files)==0) return(character(0))
  sub("^main_(\\d{6}_\\d{6})\\.parquet$", "\\1", files)
}

done_ids <- if (RESUME) existing_chunks() else character(0)

safe_one_cell <- function(N, lambda, rho, truth, rep) {
  train <- gen_acsi_data(N=N, lambda=lambda, rho=rho, truth_id=truth)
  test  <- gen_acsi_data(N=N_test, lambda=lambda, rho=rho, truth_id=truth)

  flags <- make_truth_flags(truth)
  truth_vec <- c(E_L=flags$E_L, Q_L=flags$Q_L, V_L=flags$V_L)

  tryCatch({
    # A) Naïve CV
    naive <- naive_cv_select(train, V=NAIVE_V, R=NAIVE_REPEATS, seed=1000+rep)
    cand_naive <- candidates %>% filter(model_id==naive$model_id) %>% slice(1)
    true_naive <- eval_model(train, test, cand_naive$model_syntax)

    # B) Nested CV
    nested <- nested_cv_select(train, OUTER_V=OUTER_V, OUTER_R=OUTER_REPEATS,
                               INNER_V=INNER_V, INNER_R=INNER_REPEATS, seed=2000+rep)

    # Frequency-weighted true error (pipeline-consistent)
    true_nested_fw <- nested_pipeline_true_error(train, test, nested$model_freq)

    # Modal-model true error (for robustness comparison)
    top_model <- nested$model_freq %>% arrange(desc(freq)) %>% slice(1) %>% pull(model_id)
    cand_nested_modal <- candidates %>% filter(model_id==top_model) %>% slice(1)
    true_nested_modal <- eval_model(train, test, cand_nested_modal$model_syntax)

    # C) Stable-core (tau sensitivity)
    stable_rows <- purrr::map_dfr(tau_list, function(tau) {
      core <- stable_core_select(nested$selected, tau=tau)
      true_core <- eval_model(train, test, core$model_syntax)
      tibble(
        tau=tau,
        core_model=core$model_id,
        true_rmse=true_core["rmse"],
        true_mae=true_core["mae"],
        pE=core$pE, pQ=core$pQ, pV=core$pV
      )
    })

    # Selection accuracy
    acc_naive <- acc_metrics(naive$E_L, naive$Q_L, naive$V_L, truth_vec)

    # Per-fold averaged accuracy (pipeline-consistent)
    acc_nested_pf <- per_fold_acc_metrics(nested$selected, truth_vec)

    # Modal-model accuracy (for robustness comparison)
    top_flags <- candidates %>% filter(model_id==top_model) %>% slice(1)
    acc_nested_modal <- acc_metrics(top_flags$E_L, top_flags$Q_L, top_flags$V_L, truth_vec)

    # Assemble output row
    main <- tibble(
      N=N, lambda=lambda, rho=rho, truth=truth, rep=rep,
      status="ok",

      # Naïve pipeline
      naive_model       = naive$model_id,
      naive_reported_rmse = naive$reported_rmse,
      naive_true_rmse   = true_naive["rmse"],
      naive_true_mae    = true_naive["mae"],
      naive_optimism_rmse = naive$reported_rmse - true_naive["rmse"],
      naive_TPR = acc_naive["TPR"],
      naive_FPR = acc_naive["FPR"],
      naive_F1  = acc_naive["F1"],

      # Nested pipeline — PRIMARY: frequency-weighted (pipeline-consistent)
      nested_top_model  = top_model,
      nested_reported_rmse = nested$reported_rmse,
      nested_true_rmse  = true_nested_fw["true_rmse"],
      nested_true_mae   = true_nested_fw["true_mae"],
      nested_optimism_rmse = nested$reported_rmse - true_nested_fw["true_rmse"],
      nested_TPR = acc_nested_pf["TPR"],
      nested_FPR = acc_nested_pf["FPR"],
      nested_F1  = acc_nested_pf["F1"],

      # Nested — SECONDARY: modal-model metrics (for robustness comparison)
      nested_modal_true_rmse = true_nested_modal["rmse"],
      nested_modal_true_mae  = true_nested_modal["mae"],
      nested_modal_optimism  = nested$reported_rmse - true_nested_modal["rmse"],
      nested_modal_TPR = acc_nested_modal["TPR"],
      nested_modal_FPR = acc_nested_modal["FPR"],
      nested_modal_F1  = acc_nested_modal["F1"]
    )

    stable <- stable_rows %>% mutate(N=N, lambda=lambda, rho=rho, truth=truth, rep=rep)

    list(main=main, stable=stable)
  }, error=function(e) {
    main <- tibble(
      N=N, lambda=lambda, rho=rho, truth=truth, rep=rep,
      status=paste0("error: ", conditionMessage(e))
    )
    stable <- tibble(N=N, lambda=lambda, rho=rho, truth=truth, rep=rep,
                     tau=NA_real_, core_model=NA_character_, true_rmse=NA_real_,
                     true_mae=NA_real_, pE=NA_real_, pQ=NA_real_, pV=NA_real_)
    list(main=main, stable=stable)
  })
}

# Chunked run with progress + resume
n <- nrow(grid)
with_progress({
  p <- progressr::progressor(along = 1:n)

  for (start in seq(1, n, by = CHUNK_SIZE)) {
    end <- min(start + CHUNK_SIZE - 1, n)
    part_id <- sprintf("%06d_%06d", start, end)

    if (RESUME && part_id %in% done_ids) {
      log_to_file("SKIP chunk ", part_id, " (already exists)")
      next
    }

    log_to_file("RUN chunk ", part_id, " (rows ", start, "..", end, " of ", n, ")")

    chunk <- grid[start:end, , drop=FALSE]

    res_chunk <- furrr::future_pmap(
      chunk,
      function(N, lambda, rho, truth, rep) { p(); safe_one_cell(N, lambda, rho, truth, rep) },
      .options = furrr::furrr_options(seed=TRUE, scheduling=2)
    )

    main_part   <- dplyr::bind_rows(purrr::map(res_chunk, "main"))
    stable_part <- dplyr::bind_rows(purrr::map(res_chunk, "stable"))

    arrow::write_parquet(main_part,   file.path("../sim_outputs/chunks", paste0("main_", part_id, ".parquet")))
    arrow::write_parquet(stable_part, file.path("../sim_outputs/chunks", paste0("stable_", part_id, ".parquet")))
  }
})

# Stitch chunks
log_to_file("Stitching chunks ...")
main_files   <- list.files("../sim_outputs/chunks", pattern="^main_.*\\.parquet$", full.names=TRUE)
stable_files <- list.files("../sim_outputs/chunks", pattern="^stable_.*\\.parquet$", full.names=TRUE)

main   <- dplyr::bind_rows(lapply(main_files, arrow::read_parquet))
stable <- dplyr::bind_rows(lapply(stable_files, arrow::read_parquet))

arrow::write_parquet(main,   "../sim_outputs/sim_main.parquet")
arrow::write_parquet(stable, "../sim_outputs/sim_stable.parquet")

write.csv(main,   "../sim_outputs/sim_main.csv", row.names=FALSE)
write.csv(stable, "../sim_outputs/sim_stable.csv", row.names=FALSE)

log_to_file("DONE. Wrote ../sim_outputs/sim_main.* and ../sim_outputs/sim_stable.*")

future::plan(future::sequential)
