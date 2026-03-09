# 01_core_functions.R
# ============================================================
# Core DGP + model candidates + evaluation + three pipelines
# A) Naïve CV selection (optimistic)
# B) Nested CV selection (TC-PVSS)
# C) Stable-core TC-PVSS (tau)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(purrr); library(rsample); library(cSEM); library(glue)
})

# --------------------------
# Truth flags for optional edges into L
# --------------------------
make_truth_flags <- function(truth_id) {
  E_L <- truth_id %in% c("M1","M4","M5","M7")
  Q_L <- truth_id %in% c("M2","M4","M6","M7")
  V_L <- truth_id %in% c("M3","M5","M6","M7")
  list(E_L=E_L, Q_L=Q_L, V_L=V_L)
}

# --------------------------
# ACSI-like DGP parameters (tunable in sensitivity)
# --------------------------
beta_QE <- 0.60
beta_VE <- 0.30
beta_VQ <- 0.40
beta_SE <- 0.20
beta_SQ <- 0.30
beta_SV <- 0.40
beta_LS <- 0.50

beta_LE <- 0.15
beta_LQ <- 0.15
beta_LV <- 0.15

# --------------------------
# Generate dataset: latents + reflective indicators
# L is measured by 3 reflective indicators (REPUR, HIGHPTOL,
# LOWPTOL) in the simulation DGP. This differs from the empirical
# illustration (single-item REPUR) to support loading/reliability
# manipulation in the Monte Carlo design (see manuscript Section 6.1).
# --------------------------
gen_acsi_data <- function(N, lambda=0.70, rho=0.0, truth_id="M6") {
  flags <- make_truth_flags(truth_id)

  E <- rnorm(N, 0, 1)

  eps <- MASS::mvrnorm(N, mu=c(0,0), Sigma=matrix(c(1,rho,rho,1),2,2))
  eps_Q <- eps[,1]; eps_V <- eps[,2]

  Q <- beta_QE*E + eps_Q
  V <- beta_VE*E + beta_VQ*Q + eps_V
  S <- beta_SE*E + beta_SQ*Q + beta_SV*V + rnorm(N,0,1)

  L <- beta_LS*S +
    (if(flags$E_L) beta_LE*E else 0) +
    (if(flags$Q_L) beta_LQ*Q else 0) +
    (if(flags$V_L) beta_LV*V else 0) +
    rnorm(N,0,1)

  mk_ind <- function(lat) lambda*lat + sqrt(1-lambda^2)*rnorm(N,0,1)

  tibble::tibble(
    OVERALLX = mk_ind(E), CUSTOMX = mk_ind(E), WRONGX = mk_ind(E),
    OVERALLQ = mk_ind(Q), CUSTOMQ = mk_ind(Q), WRONGQ = mk_ind(Q),
    PQ = mk_ind(V), QP = mk_ind(V),
    SATIS = mk_ind(S), CONFIRM = mk_ind(S), IDEAL = mk_ind(S),
    REPUR    = mk_ind(L),
    HIGHPTOL = mk_ind(L),
    LOWPTOL  = mk_ind(L)
  )
}

INDICATORS <- c("OVERALLX","CUSTOMX","WRONGX",
                 "OVERALLQ","CUSTOMQ","WRONGQ",
                 "PQ","QP",
                 "SATIS","CONFIRM","IDEAL",
                 "REPUR","HIGHPTOL","LOWPTOL")

# --------------------------
# Candidate models M0..M7 (ACSI backbone + optional edges into L)
# L is measured by 3 reflective indicators
# --------------------------
make_model_syntax <- function(add_E_L=FALSE, add_Q_L=FALSE, add_V_L=FALSE) {
  rhs_L <- c("S",
             if(add_E_L) "E" else NULL,
             if(add_Q_L) "Q" else NULL,
             if(add_V_L) "V" else NULL)
  rhs_L <- paste(rhs_L, collapse=" + ")

  glue::glue("
    E <~ OVERALLX + CUSTOMX + WRONGX
    Q <~ OVERALLQ + CUSTOMQ + WRONGQ
    V <~ PQ + QP
    S <~ SATIS + CONFIRM + IDEAL
    L <~ REPUR + HIGHPTOL + LOWPTOL

    Q ~ E
    V ~ E + Q
    S ~ E + Q + V
    L ~ {rhs_L}
  ")
}

candidates <- tibble(
  model_id = paste0("M", 0:7),
  E_L = c(FALSE, TRUE,  FALSE, FALSE, TRUE,  TRUE,  FALSE, TRUE),
  Q_L = c(FALSE, FALSE, TRUE,  FALSE, TRUE,  FALSE, TRUE,  TRUE),
  V_L = c(FALSE, FALSE, FALSE, TRUE,  FALSE, TRUE,  TRUE,  TRUE)
) %>%
  mutate(
    model_syntax = purrr::pmap_chr(list(E_L,Q_L,V_L), ~as.character(make_model_syntax(..1,..2,..3))),
    complexity = as.integer(E_L)+as.integer(Q_L)+as.integer(V_L)
  )

# --------------------------
# Mean-impute train-only then apply to test (no leakage)
# --------------------------
impute_train_apply <- function(train_df, test_df) {
  means <- vapply(train_df[INDICATORS], function(x) mean(as.numeric(x), na.rm=TRUE), numeric(1))
  tr <- train_df; te <- test_df
  for (nm in INDICATORS) {
    tr[[nm]] <- as.numeric(tr[[nm]]); te[[nm]] <- as.numeric(te[[nm]])
    tr[[nm]][is.na(tr[[nm]])] <- means[[nm]]
    te[[nm]][is.na(te[[nm]])] <- means[[nm]]
  }
  list(train=tr, test=te)
}

# Robust extraction of predictions for REPUR across cSEM versions
get_pred <- function(pred_obj, target="REPUR") {
  x <- NULL
  if ("Predictions_target" %in% names(pred_obj)) x <- pred_obj$Predictions_target
  else if ("Prediction_target" %in% names(pred_obj)) x <- pred_obj$Prediction_target
  else if ("Predictions" %in% names(pred_obj)) x <- pred_obj$Predictions

  if (!is.null(x)) {
    if (is.matrix(x) || is.data.frame(x)) {
      m <- as.matrix(x)
      if (!is.null(colnames(m)) && target %in% colnames(m)) return(as.numeric(m[,target]))
      return(as.numeric(m[,1]))
    }
    if (is.list(x)) {
      for (el in x) {
        if (is.matrix(el) || is.data.frame(el)) {
          m <- as.matrix(el)
          if (!is.null(colnames(m)) && target %in% colnames(m)) return(as.numeric(m[,target]))
          if (ncol(m)==1) return(as.numeric(m[,1]))
        }
      }
    }
    if (is.numeric(x)) return(as.numeric(x))
  }

  for (nm in names(pred_obj)) {
    el <- pred_obj[[nm]]
    if (is.matrix(el) || is.data.frame(el)) {
      m <- as.matrix(el)
      if (!is.null(colnames(m)) && target %in% colnames(m)) return(as.numeric(m[,target]))
      if (ncol(m)==1) return(as.numeric(m[,1]))
    }
  }
  stop("Could not extract predictions. names(pred_obj)=", paste(names(pred_obj), collapse=", "))
}

# Fit+predict on heldout (single split)
eval_model <- function(train_df, test_df, model_syntax) {
  imp <- impute_train_apply(train_df, test_df)
  tr <- imp$train; te <- imp$test

  fit <- csem(.data=tr[,INDICATORS], .model=model_syntax, .approach_weights="PLS-PM")
  pr  <- suppressWarnings(predict(.object=fit, .approach_predict="direct", .test_data=te[,INDICATORS]))

  y <- as.numeric(te$REPUR)
  yhat <- get_pred(pr, "REPUR")

  c(
    rmse = sqrt(mean((y-yhat)^2)),
    mae  = mean(abs(y-yhat))
  )
}

# A) Naïve CV selection (single CV used for selection AND reported performance)
#    V and R are passed from config (NAIVE_V, NAIVE_REPEATS) to match
#    the resampling budget of the nested outer loop.
naive_cv_select <- function(data, V=10, R=1, seed=1) {
  set.seed(seed)
  folds <- vfold_cv(data, v=V, repeats=R)

  score_one <- function(model_syntax) {
    vals <- purrr::map_dbl(folds$splits, function(sp) {
      tr <- analysis(sp); te <- assessment(sp)
      eval_model(tr, te, model_syntax)["rmse"]
    })
    mean(vals)
  }

  scores <- candidates %>%
    mutate(cv_rmse = purrr::map_dbl(model_syntax, score_one)) %>%
    arrange(cv_rmse, complexity)

  best <- scores[1,]
  list(model_id=best$model_id, E_L=best$E_L, Q_L=best$Q_L, V_L=best$V_L, reported_rmse=best$cv_rmse)
}

# B) Nested CV selection (TC-PVSS): inner selects; outer reports unbiased error
#    Returns per-fold selection flags for pipeline-consistent evaluation.
nested_cv_select <- function(data, OUTER_V=10, OUTER_R=1, INNER_V=5, INNER_R=1, seed=1) {
  set.seed(seed)
  outer <- vfold_cv(data, v=OUTER_V, repeats=OUTER_R)

  sel <- vector("list", length(outer$splits))
  outer_rmse <- numeric(length(outer$splits))

  for (i in seq_along(outer$splits)) {
    tr_outer <- analysis(outer$splits[[i]])
    te_outer <- assessment(outer$splits[[i]])

    set.seed(seed + i)
    inner <- vfold_cv(tr_outer, v=INNER_V, repeats=INNER_R)

    score_one <- function(model_syntax) {
      vals <- purrr::map_dbl(inner$splits, function(sp) {
        tr <- analysis(sp); te <- assessment(sp)
        eval_model(tr, te, model_syntax)["rmse"]
      })
      mean(vals)
    }

    scores <- candidates %>%
      mutate(inner_rmse = purrr::map_dbl(model_syntax, score_one)) %>%
      arrange(inner_rmse, complexity)

    best <- scores[1,]

    outer_rmse[i] <- eval_model(tr_outer, te_outer, best$model_syntax)["rmse"]
    sel[[i]] <- best %>% select(model_id, E_L, Q_L, V_L, complexity)
  }

  sel_df <- bind_rows(sel)
  list(
    selected = sel_df,
    reported_rmse = mean(outer_rmse),
    model_freq = sel_df %>% count(model_id, sort=TRUE) %>% mutate(freq=n/sum(n))
  )
}

# Compute frequency-weighted true RMSE/MAE on independent test set.
# Each unique selected model is fit on the full training set and
# evaluated on the test set; results are weighted by selection
# frequency.  This matches the empirical paradigm where per-fold
# predictions from different models are pooled.
nested_pipeline_true_error <- function(train, test, model_freq_df) {
  metrics <- purrr::map_dfr(seq_len(nrow(model_freq_df)), function(j) {
    mid  <- model_freq_df$model_id[j]
    freq <- model_freq_df$freq[j]
    cand <- candidates %>% filter(model_id == mid) %>% slice(1)
    tr   <- eval_model(train, test, cand$model_syntax)
    tibble(model_id = mid, freq = freq, true_rmse = tr["rmse"], true_mae = tr["mae"])
  })
  c(
    true_rmse = sum(metrics$true_rmse * metrics$freq),
    true_mae  = sum(metrics$true_mae  * metrics$freq)
  )
}

# Per-fold averaged accuracy metrics for nested CV.
# Instead of computing TPR/FPR/F1 from the single modal model,
# average per-fold TPR/FPR/F1 across all outer folds.  This is
# consistent with the pipeline paradigm where each fold may
# select a different model.
per_fold_acc_metrics <- function(sel_df, truth_vec) {
  per_fold <- purrr::map_dfr(seq_len(nrow(sel_df)), function(j) {
    a <- acc_metrics(sel_df$E_L[j], sel_df$Q_L[j], sel_df$V_L[j], truth_vec)
    tibble(TPR = a["TPR"], FPR = a["FPR"], F1 = a["F1"])
  })
  c(
    TPR = mean(per_fold$TPR, na.rm = TRUE),
    FPR = mean(per_fold$FPR, na.rm = TRUE),
    F1  = mean(per_fold$F1,  na.rm = TRUE)
  )
}

# C) Stable-core derived from nested selections at tau
stable_core_select <- function(nested_sel_df, tau=0.7) {
  pE <- mean(nested_sel_df$E_L); pQ <- mean(nested_sel_df$Q_L); pV <- mean(nested_sel_df$V_L)
  E_L <- pE >= tau; Q_L <- pQ >= tau; V_L <- pV >= tau

  cand <- candidates %>% filter(E_L==!!E_L, Q_L==!!Q_L, V_L==!!V_L)
  if (nrow(cand)==0) {
    cand <- candidates %>%
      mutate(dist = abs(as.integer(E_L)-as.integer(!!E_L)) +
                   abs(as.integer(Q_L)-as.integer(!!Q_L)) +
                   abs(as.integer(V_L)-as.integer(!!V_L))) %>%
      arrange(dist, complexity) %>% slice(1)
  } else {
    cand <- cand %>% arrange(complexity) %>% slice(1)
  }

  list(
    model_id = cand$model_id[[1]],
    model_syntax = cand$model_syntax[[1]],
    E_L = E_L, Q_L = Q_L, V_L = V_L,
    pE = pE, pQ = pQ, pV = pV
  )
}

# Selection accuracy metrics (TPR/FPR/F1) for optional edges
acc_metrics <- function(sel_E, sel_Q, sel_V, truth_vec) {
  sel <- c(E_L=sel_E, Q_L=sel_Q, V_L=sel_V)
  TP <- sum(sel & truth_vec); FP <- sum(sel & !truth_vec); FN <- sum((!sel) & truth_vec)
  TPR <- ifelse(sum(truth_vec)>0, TP/sum(truth_vec), NA_real_)
  FPR <- ifelse(sum(!truth_vec)>0, FP/sum(!truth_vec), NA_real_)
  F1  <- ifelse((2*TP+FP+FN)>0, 2*TP/(2*TP+FP+FN), NA_real_)
  c(TPR=TPR, FPR=FPR, F1=F1)
}
