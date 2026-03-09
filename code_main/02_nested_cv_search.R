# code_main/02_nested_cv_search.R
# ------------------------------------------------------------------
# TC-PVSS (Protocol v1) — Prediction-first model selection for PLS-SEM
# ACSI demo (2015 v2)
#
# Robust version (works across cSEM predict() output variants):
# - Does NOT depend on Prediction_metrics column names.
# - Extracts predictions from predict() output by searching in
#   Predictions_target / Prediction_target (and other variants), including list cases.
# - Suppresses benign warnings about disattenuation/benchmarks.
#
# Frozen v1:
# - Candidate models: 8 (ACSI backbone + optional E->L, Q->L, V->L)
# - Selection: min inner-CV RMSE(REPUR); tie-break MAE, then parsimony
# - Nested CV: Outer 10-fold × 5 repeats; Inner 5-fold × 2 repeats (stratified by INDUSTRY)
# - Mean imputation: training-only, applied to test (no leakage)
# - lm benchmark: indicator-level linear regression fitted per outer split (not a candidate)
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(cSEM)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(rsample)
  library(glue)
  library(arrow)
  library(rlang)
})

# -----------------------------
# Config
# -----------------------------
SEED_GLOBAL   <- 20260220
OUTER_V       <- 10
OUTER_REPEATS <- 5
INNER_V       <- 5
INNER_REPEATS <- 2
STRATA_VAR    <- "INDUSTRY"
TARGET_IND    <- "REPUR"

# -----------------------------
# Load cleaned data (from code_main/00_prepare_data.R)
# -----------------------------
acsi <- readRDS("../outputs_main/acsi_clean.rds") %>%
  mutate(.id = row_number())

INDICATORS <- c(
  "OVERALLX","CUSTOMX","WRONGX",
  "OVERALLQ","CUSTOMQ","WRONGQ",
  "PQ","QP",
  "SATIS","CONFIRM","IDEAL",
  "REPUR"
)

stopifnot(all(c(STRATA_VAR, ".id") %in% names(acsi)))
stopifnot(all(INDICATORS %in% names(acsi)))

# -----------------------------
# Mean imputation: train-only -> apply to test (no leakage)
# -----------------------------
impute_train_apply <- function(train_df, test_df, cols) {
  means <- vapply(train_df[cols], function(x) mean(as.numeric(x), na.rm = TRUE), numeric(1))
  train_imp <- train_df
  test_imp  <- test_df
  for (nm in cols) {
    train_imp[[nm]] <- as.numeric(train_imp[[nm]])
    test_imp[[nm]]  <- as.numeric(test_imp[[nm]])
    train_imp[[nm]][is.na(train_imp[[nm]])] <- means[[nm]]
    test_imp[[nm]][is.na(test_imp[[nm]])]   <- means[[nm]]
  }
  list(train=train_imp, test=test_imp)
}

# -----------------------------
# Model builder
# -----------------------------
make_model_syntax <- function(add_E_L = FALSE, add_Q_L = FALSE, add_V_L = FALSE) {
  rhs_L <- c("S",
             if (add_E_L) "E" else NULL,
             if (add_Q_L) "Q" else NULL,
             if (add_V_L) "V" else NULL)
  rhs_L <- paste(rhs_L, collapse = " + ")

  glue(
    "
    # Measurement (composites / Mode A style)
    E <~ OVERALLX + CUSTOMX + WRONGX
    Q <~ OVERALLQ + CUSTOMQ + WRONGQ
    V <~ PQ + QP
    S <~ SATIS + CONFIRM + IDEAL
    L <~ REPUR

    # Structural (ACSI backbone)
    Q ~ E
    V ~ E + Q
    S ~ E + Q + V
    L ~ {rhs_L}
    "
  )
}

candidates <- tibble(
  model_id = paste0("M", 0:7),
  E_L = c(FALSE, TRUE,  FALSE, FALSE, TRUE,  TRUE,  FALSE, TRUE),
  Q_L = c(FALSE, FALSE, TRUE,  FALSE, TRUE,  FALSE, TRUE,  TRUE),
  V_L = c(FALSE, FALSE, FALSE, TRUE,  FALSE, TRUE,  TRUE,  TRUE)
) %>%
  mutate(
    complexity   = as.integer(E_L) + as.integer(Q_L) + as.integer(V_L),
    model_syntax = pmap_chr(list(E_L, Q_L, V_L), ~as.character(make_model_syntax(..1, ..2, ..3)))
  )

# -----------------------------
# Robust extraction of predicted values for TARGET_IND from cSEM::predict() output
# Handles:
# - Predictions_target (matrix/data.frame or list)
# - Prediction_target (older)
# - other variants
# Selects the vector that matches expected length.
# -----------------------------
get_pred_vec <- function(pred_obj, target_name, expected_n = NULL) {

  want_len <- function(v) {
    if (is.null(expected_n)) return(TRUE)
    length(v) == expected_n
  }

  extract_from <- function(x) {
    if (is.null(x)) return(NULL)

    # numeric vector
    if (is.numeric(x) && length(x) > 1) {
      v <- as.numeric(x)
      if (want_len(v)) return(v)
    }

    # matrix / data.frame
    if (is.matrix(x) || is.data.frame(x)) {
      m <- as.matrix(x)
      if (!is.null(colnames(m)) && target_name %in% colnames(m)) {
        v <- as.numeric(m[, target_name])
        if (want_len(v)) return(v)
      }
      if (ncol(m) == 1) {
        v <- as.numeric(m[, 1])
        if (want_len(v)) return(v)
      }
    }

    # list (may contain matrices per construct/indicator)
    if (is.list(x)) {
      # direct named access
      if (!is.null(names(x)) && target_name %in% names(x)) {
        v <- extract_from(x[[target_name]])
        if (!is.null(v)) return(v)
      }
      # scan elements
      for (el in x) {
        v <- extract_from(el)
        if (!is.null(v)) return(v)
      }
    }

    NULL
  }

  preferred <- c(
    "Predictions_target", "Prediction_target",
    "Predictions", "Prediction",
    "predictions", "prediction"
  )

  for (nm in preferred) {
    if (nm %in% names(pred_obj)) {
      v <- extract_from(pred_obj[[nm]])
      if (!is.null(v)) return(v)
    }
  }

  # last resort: scan all entries
  for (nm in names(pred_obj)) {
    v <- extract_from(pred_obj[[nm]])
    if (!is.null(v)) return(v)
  }

  stop(
    "Could not extract predictions for ", target_name,
    " from cSEM::predict() output. names(pred_obj) = ",
    paste(names(pred_obj), collapse = ", ")
  )
}

# -----------------------------
# Inner CV scoring (manual)
# -----------------------------
inner_score_candidate <- function(train_df_raw, model_syntax, seed_inner) {
  set.seed(seed_inner)

  inner_folds <- vfold_cv(
    train_df_raw,
    v = INNER_V,
    repeats = INNER_REPEATS,
    strata = !!sym(STRATA_VAR)
  )

  fold_metrics <- map_dfr(inner_folds$splits, function(sp) {
    tr <- analysis(sp)
    te <- assessment(sp)

    imp <- impute_train_apply(tr, te, INDICATORS)
    tr_imp <- imp$train
    te_imp <- imp$test

    tryCatch({
      fit <- csem(.data = tr_imp[, INDICATORS], .model = model_syntax, .approach_weights = "PLS-PM")

      pr  <- suppressWarnings(
        predict(.object = fit, .approach_predict = "direct", .test_data = te_imp[, INDICATORS])
      )

      y    <- as.numeric(te_imp[[TARGET_IND]])
      yhat <- get_pred_vec(pr, TARGET_IND, expected_n = length(y))

      tibble(
        rmse = sqrt(mean((y - yhat)^2)),
        mae  = mean(abs(y - yhat))
      )
    }, error = function(e) {
      tibble(rmse = Inf, mae = Inf)
    })
  })

  list(
    rmse = mean(fold_metrics$rmse),
    mae  = mean(fold_metrics$mae)
  )
}

# -----------------------------
# Inner selection
# -----------------------------
select_model_inner <- function(train_df_raw, seed_base) {
  scores <- candidates %>%
    mutate(score = pmap(
      list(model_syntax, seq_len(n()) + seed_base),
      ~inner_score_candidate(train_df_raw, ..1, ..2)
    )) %>%
    mutate(rmse = map_dbl(score, "rmse"),
           mae  = map_dbl(score, "mae")) %>%
    arrange(rmse, mae, complexity)

  scores[1, ] %>%
    select(model_id, E_L, Q_L, V_L, complexity, rmse, mae, model_syntax)
}

# -----------------------------
# Outer evaluation per split
# -----------------------------
evaluate_on_outer_split <- function(split_obj, split_index) {
  train_raw <- analysis(split_obj)
  test_raw  <- assessment(split_obj)

  chosen <- select_model_inner(train_raw, seed_base = 10000 + split_index)

  imp <- impute_train_apply(train_raw, test_raw, INDICATORS)
  train_imp <- imp$train
  test_imp  <- imp$test

  fit_chosen <- csem(.data = train_imp[, INDICATORS], .model = chosen$model_syntax, .approach_weights = "PLS-PM")
  pr_chosen  <- suppressWarnings(
    predict(.object = fit_chosen, .approach_predict = "direct", .test_data = test_imp[, INDICATORS])
  )

  base <- candidates %>% filter(model_id == "M0") %>% slice(1)
  fit_base <- csem(.data = train_imp[, INDICATORS], .model = base$model_syntax, .approach_weights = "PLS-PM")
  pr_base  <- suppressWarnings(
    predict(.object = fit_base, .approach_predict = "direct", .test_data = test_imp[, INDICATORS])
  )

  # --- lm benchmark: regress TARGET on all other indicators ---
  # This is NOT a candidate model; it benchmarks whether the PLS structural
  # model adds value over a simple linear regression on manifest indicators.
  lm_predictors <- setdiff(INDICATORS, TARGET_IND)
  lm_formula    <- reformulate(lm_predictors, response = TARGET_IND)
  fit_lm        <- lm(lm_formula, data = train_imp)
  yhat_lm       <- stats::predict(fit_lm, newdata = test_imp)

  y <- as.numeric(test_imp[[TARGET_IND]])
  yhat_chosen <- get_pred_vec(pr_chosen, TARGET_IND, expected_n = length(y))
  yhat_base   <- get_pred_vec(pr_base, TARGET_IND, expected_n = length(y))

  err_chosen <- y - yhat_chosen
  err_base   <- y - yhat_base
  err_lm     <- y - yhat_lm

  tibble(
    outer_split = split_index,
    .row_id     = test_imp$.id,
    industry    = test_imp[[STRATA_VAR]],
    model_id    = chosen$model_id,
    E_L         = chosen$E_L,
    Q_L         = chosen$Q_L,
    V_L         = chosen$V_L,
    complexity  = chosen$complexity,
    actual      = y,
    pred_chosen = yhat_chosen,
    pred_base   = yhat_base,
    pred_lm     = as.numeric(yhat_lm),
    se_chosen   = err_chosen^2,
    ae_chosen   = abs(err_chosen),
    se_base     = err_base^2,
    ae_base     = abs(err_base),
    se_lm       = err_lm^2,
    ae_lm       = abs(err_lm)
  )
}

# -----------------------------
# Run outer CV
# -----------------------------
set.seed(SEED_GLOBAL)

outer_splits <- vfold_cv(
  acsi,
  v = OUTER_V,
  repeats = OUTER_REPEATS,
  strata = !!sym(STRATA_VAR)
)

res_long <- map2_dfr(
  outer_splits$splits,
  seq_along(outer_splits$splits),
  ~evaluate_on_outer_split(.x, split_index = .y)
)

dir.create("../outputs_main", showWarnings = FALSE)
write_parquet(res_long, "../outputs_main/pred_errors_long.parquet")

outer_summary <- res_long %>%
  group_by(outer_split, model_id, E_L, Q_L, V_L, complexity) %>%
  summarise(
    RMSE      = sqrt(mean(se_chosen, na.rm = TRUE)),
    MAE       = mean(ae_chosen, na.rm = TRUE),
    RMSE_base = sqrt(mean(se_base, na.rm = TRUE)),
    MAE_base  = mean(ae_base, na.rm = TRUE),
    RMSE_lm   = sqrt(mean(se_lm, na.rm = TRUE)),
    MAE_lm    = mean(ae_lm, na.rm = TRUE),
    n_test    = n(),
    .groups   = "drop"
  )

write_parquet(outer_summary, "../outputs_main/outer_results.parquet")

cat("Done.\nSaved:\n- ../outputs_main/pred_errors_long.parquet\n- ../outputs_main/outer_results.parquet\n")
