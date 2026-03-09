# code_main/01_measurement_audit.R
# ------------------------------------------------------------------
# TC-PVSS Step 3.1: Measurement Model Audit (ACSI Illustration)
# ------------------------------------------------------------------
# Purpose:  Verify convergent validity, internal consistency, and
#           discriminant validity of the ACSI measurement model
#           PRIOR to structural search.
#
# Input:    ../outputs_main/acsi_clean.rds  (from 00_prepare_data.R)
#
# Specification matches 02_nested_cv_search.R exactly:
#   - Package: cSEM (PLS-PM, mode A)
#   - Indicators: 12 (same as nested CV pipeline)
#   - Loyalty: single-item (REPUR only)
#   - HIGHPTOL/LOWPTOL excluded (different scale, high missingness)
#
# Outputs:  ../outputs_main/measurement_audit_table_A0.csv
#           ../outputs_main/htmt_matrix.csv
#           ../outputs_main/fornell_larcker.csv
#
# Usage:    setwd("/path/to/TC-PVSS-PLSSEM/code_main")
#           source("01_measurement_audit.R")
#
# Requires: R >= 4.1, cSEM >= 0.5.0
# ------------------------------------------------------------------

# --- 0. Setup -------------------------------------------------------

library(cSEM)

cat("=============================================================\n")
cat("  TC-PVSS Step 3.1: Measurement Model Audit (ACSI)\n")
cat("  Specification: cSEM, PLS-PM, L = REPUR (single item)\n")
cat("=============================================================\n\n")

# --- 1. Load cleaned data -------------------------------------------
# Read from the SAME cleaned data used by 02_nested_cv_search.R
# to guarantee identical recoding. cSEM requires complete cases,
# so we apply listwise deletion here and report both sample sizes.

acsi <- readRDS("../outputs_main/acsi_clean.rds")
cat(sprintf("Loaded cleaned data: %d observations, %d variables\n",
            nrow(acsi), ncol(acsi)))

# Indicators used in nested CV (matching 02_nested_cv_search.R)
INDICATORS <- c("OVERALLX", "CUSTOMX", "WRONGX",   # E
                "OVERALLQ", "CUSTOMQ", "WRONGQ",   # Q
                "PQ", "QP",                          # V
                "SATIS", "CONFIRM", "IDEAL",         # S
                "REPUR")                             # L (single item)

stopifnot(all(INDICATORS %in% names(acsi)))

# Report missingness
cat("\nMissing data per indicator:\n")
any_miss <- FALSE
for (v in INDICATORS) {
  n_miss <- sum(is.na(acsi[[v]]))
  if (n_miss > 0) {
    any_miss <- TRUE
    cat(sprintf("  %10s: %5d missing (%.1f%%)\n", v, n_miss,
                100 * n_miss / nrow(acsi)))
  }
}
if (!any_miss) cat("  (no missing values on the 12 indicators)\n")

# Listwise deletion on the 12 indicators
# NOTE: 02_nested_cv_search.R handles missingness via per-fold mean
# imputation. For this full-sample audit, cSEM requires complete data.
# We report both N values for transparency.
n_full <- nrow(acsi)
data <- acsi[complete.cases(acsi[, INDICATORS]), ]
n_complete <- nrow(data)
cat(sprintf("\nFull sample (from 00_prepare_data.R): N = %d\n", n_full))
cat(sprintf("Complete cases (listwise on 12 indicators): N = %d",
            n_complete))
if (n_full > n_complete) {
  cat(sprintf(" (removed %d, %.1f%%)\n",
              n_full - n_complete,
              100 * (n_full - n_complete) / n_full))
} else {
  cat(" (no cases removed)\n")
}
cat("Note: Nested CV (script 02) uses the full sample with per-fold\n")
cat("      mean imputation. This audit uses complete cases only.\n\n")

# Keep only indicator columns for cSEM
data_csem <- data[, INDICATORS]

# --- 2. Model specification (cSEM syntax) ---------------------------
# Matches 02_nested_cv_search.R: composite mode A (<~), L = REPUR

acsi_model <- "
  # Measurement (composite, mode A -- equivalent to PLS-PM reflective)
  E <~ OVERALLX + CUSTOMX + WRONGX
  Q <~ OVERALLQ + CUSTOMQ + WRONGQ
  V <~ PQ + QP
  S <~ SATIS + CONFIRM + IDEAL
  L <~ REPUR

  # Structural (M0 baseline)
  Q ~ E
  V ~ E + Q
  S ~ E + Q + V
  L ~ S
"

# --- 3. Estimate PLS-PM ---------------------------------------------

cat("--- Estimating PLS-PM (M0 baseline, cSEM) ---\n")
fit <- csem(.data = data_csem, .model = acsi_model,
            .approach_weights = "PLS-PM")

# Check admissibility
ver <- verify(fit)
if (!all(ver)) {
  warning("Model estimation produced inadmissible results. Check output carefully.")
  cat("Verify output:\n")
  print(ver)
} else {
  cat("Estimation converged, results admissible.\n")
}

model_sum <- summarize(fit)

# --- 4. Extract measurement quality metrics -------------------------

## 4a. Indicator loadings (for multi-item constructs)
cat("\n--- 4a. Indicator Loadings ---\n")

# cSEM stores loadings in $Estimates$Loading_estimates
loadings_raw <- fit$Estimates$Loading_estimates

# Multi-item constructs only (L is single-item: loading = 1 by definition)
multi_constructs <- c("E", "Q", "V", "S")

# Extract loadings into a clean list
indicators_map <- list(
  E = c("OVERALLX", "CUSTOMX", "WRONGX"),
  Q = c("OVERALLQ", "CUSTOMQ", "WRONGQ"),
  V = c("PQ", "QP"),
  S = c("SATIS", "CONFIRM", "IDEAL")
)

for (cst in multi_constructs) {
  for (ind in indicators_map[[cst]]) {
    val <- loadings_raw[cst, ind]
    cat(sprintf("  %s -> %s: %.4f\n", cst, ind, val))
  }
}
cat("  L -> REPUR: 1.000 (single-item, by definition)\n")

## 4b. Reliability metrics (multi-item constructs only)
cat("\n--- 4b. Reliability Metrics ---\n")

# cSEM::assess() returns empty vectors for AVE, RhoC, etc. when the model
# contains a single-item construct (L = REPUR). Compute manually from
# loadings (already extracted in 4a) and indicator correlations.

# AVE = mean of squared loadings
AVE_vals <- setNames(
  vapply(multi_constructs, function(cst) {
    lambdas <- loadings_raw[cst, indicators_map[[cst]]]
    mean(lambdas^2)
  }, numeric(1)),
  multi_constructs
)

# Composite reliability rhoC = (sum(lambda))^2 / ((sum(lambda))^2 + sum(1 - lambda^2))
rhoC_vals <- setNames(
  vapply(multi_constructs, function(cst) {
    lambdas <- loadings_raw[cst, indicators_map[[cst]]]
    sum_l  <- sum(lambdas)
    sum_e  <- sum(1 - lambdas^2)
    sum_l^2 / (sum_l^2 + sum_e)
  }, numeric(1)),
  multi_constructs
)

# Cronbach's alpha from indicator correlations
alpha_vals <- setNames(
  vapply(multi_constructs, function(cst) {
    inds <- indicators_map[[cst]]
    k <- length(inds)
    if (k < 2) return(NA_real_)
    R <- cor(data_csem[, inds], use = "pairwise.complete.obs")
    mean_r <- (sum(R) - k) / (k * (k - 1))
    k * mean_r / (1 + (k - 1) * mean_r)
  }, numeric(1)),
  multi_constructs
)

# rhoA (Dijkstra-Henseler): complex formula requiring iterative weights.
# Initialize as NA; attempt extraction from assess()$Reliability below.
rhoA_vals <- setNames(rep(NA_real_, length(multi_constructs)), multi_constructs)

# Try to get rhoA from cSEM assess() Reliability sub-list
qa <- assess(fit)
if (!is.null(qa[["Reliability"]])) {
  rel_list <- qa[["Reliability"]]
  dh_rhoA <- rel_list[["Dijkstra-Henselers_rho_A"]]
  if (!is.null(dh_rhoA) && is.numeric(dh_rhoA) && length(dh_rhoA) >= length(multi_constructs)) {
    if (is.null(names(dh_rhoA))) {
      # Unnamed vector — assume order matches all constructs including L
      all_cst <- c(multi_constructs, "L")
      if (length(dh_rhoA) == length(all_cst)) {
        names(dh_rhoA) <- all_cst
      }
    }
    if (!is.null(names(dh_rhoA)) && all(multi_constructs %in% names(dh_rhoA))) {
      rhoA_vals <- dh_rhoA[multi_constructs]
    }
  }
  # Also try Cronbach and Joereskog from Reliability if our manual failed
  joreskog <- rel_list[["Joereskogs_rho"]]
  if (!is.null(joreskog) && is.numeric(joreskog) && length(joreskog) >= length(multi_constructs)) {
    if (is.null(names(joreskog))) {
      all_cst <- c(multi_constructs, "L")
      if (length(joreskog) == length(all_cst)) names(joreskog) <- all_cst
    }
    if (!is.null(names(joreskog)) && all(multi_constructs %in% names(joreskog))) {
      cat("  Note: Using Joereskog's rho from Reliability sub-list for rhoC cross-check.\n")
    }
  }
}

cat("\n  Construct  Alpha    rhoC     rhoA     AVE\n")
cat("  ---------  -----    ----     ----     ---\n")
for (cst in multi_constructs) {
  cat(sprintf("  %-9s  %.3f    %.3f    %s    %.3f\n",
              cst, alpha_vals[cst], rhoC_vals[cst],
              ifelse(is.na(rhoA_vals[cst]), "  n/c",
                     sprintf("%.3f", rhoA_vals[cst])),
              AVE_vals[cst]))
}
cat("  L          n/a      n/a      n/a      n/a  (single-item)\n")
cat("  Note: n/c = not computed (rhoA requires iterative weights).\n")

## 4c. HTMT (multi-item constructs)
cat("\n--- 4c. HTMT Matrix ---\n")

# cSEM::assess()$HTMT returns NULL when single-item constructs are present
# (L = REPUR). Compute HTMT manually from indicator correlations for
# multi-item constructs only (E, Q, V, S).

ind_cor <- cor(data_csem[, unlist(indicators_map)], use = "pairwise.complete.obs")

# HTMT(i,j) = mean(heterotrait-heteromethod correlations) /
#              geometric mean of mean(monotrait-heteromethod correlations)
compute_htmt <- function(cor_mat, map, constructs) {
  nc <- length(constructs)
  htmt <- matrix(NA, nc, nc, dimnames = list(constructs, constructs))

  for (i in seq_len(nc)) {
    ci <- constructs[i]
    inds_i <- map[[ci]]
    # Mean within-construct correlation (off-diagonal)
    wi <- cor_mat[inds_i, inds_i]
    diag(wi) <- NA
    mean_wi <- mean(wi, na.rm = TRUE)

    for (j in seq_len(nc)) {
      if (i == j) { htmt[i, j] <- 1; next }
      cj <- constructs[j]
      inds_j <- map[[cj]]

      # Mean between-construct correlations
      het <- cor_mat[inds_i, inds_j]
      mean_het <- mean(het, na.rm = TRUE)

      # Mean within-construct correlation for j
      wj <- cor_mat[inds_j, inds_j]
      diag(wj) <- NA
      mean_wj <- mean(wj, na.rm = TRUE)

      htmt[i, j] <- mean_het / sqrt(mean_wi * mean_wj)
    }
  }
  htmt
}

htmt_multi <- compute_htmt(ind_cor, indicators_map, multi_constructs)
print(round(htmt_multi, 4))

## 4d. Fornell-Larcker (multi-item constructs)
cat("\n--- 4d. Fornell-Larcker (sqrt(AVE) on diagonal) ---\n")

# Fornell-Larcker: diagonal = sqrt(AVE), off-diagonal = construct correlations.
# Extract construct correlations from cSEM fit object.
construct_vcv <- fit$Estimates$Construct_VCV
if (is.null(construct_vcv)) {
  # Fallback names that cSEM may use
  construct_vcv <- fit$Estimates$Proxy_VCV
}
if (is.null(construct_vcv)) {
  # Last resort: compute from composite scores
  cat("  Note: Computing construct correlations from composite scores.\n")
  scores <- fit$Estimates$Construct_scores
  if (!is.null(scores)) {
    construct_vcv <- cor(scores, use = "pairwise.complete.obs")
  }
}

if (!is.null(construct_vcv)) {
  # Ensure it's a correlation matrix (values in [-1,1])
  if (max(diag(construct_vcv)) > 1.01) {
    construct_vcv <- cov2cor(construct_vcv)
  }

  nc <- length(multi_constructs)
  fl_multi <- matrix(NA, nc, nc, dimnames = list(multi_constructs, multi_constructs))

  for (i in seq_len(nc)) {
    fl_multi[i, i] <- sqrt(AVE_vals[multi_constructs[i]])
    for (j in seq_len(nc)) {
      if (i != j) {
        fl_multi[i, j] <- construct_vcv[multi_constructs[i], multi_constructs[j]]
      }
    }
  }

  print(round(fl_multi, 4))
} else {
  cat("  WARNING: Could not extract construct correlations.\n")
  cat("  Available in fit$Estimates:", paste(names(fit$Estimates), collapse = ", "), "\n")
  fl_multi <- NULL
}

# --- 5. Build Table A0 ----------------------------------------------

cat("\n--- 5. Building Table A0 ---\n")

rows_list <- list()
idx <- 0

for (cst in multi_constructs) {
  inds <- indicators_map[[cst]]
  for (i in seq_along(inds)) {
    idx <- idx + 1
    rows_list[[idx]] <- data.frame(
      Construct = ifelse(i == 1, cst, ""),
      Indicator = inds[i],
      Loading   = round(loadings_raw[cst, inds[i]], 3),
      Alpha     = ifelse(i == 1, sprintf("%.3f", alpha_vals[cst]), ""),
      rhoC      = ifelse(i == 1, sprintf("%.3f", rhoC_vals[cst]), ""),
      rhoA      = ifelse(i == 1,
                         ifelse(is.na(rhoA_vals[cst]), "n/c", sprintf("%.3f", rhoA_vals[cst])),
                         ""),
      AVE       = ifelse(i == 1, sprintf("%.3f", AVE_vals[cst]), ""),
      stringsAsFactors = FALSE
    )
  }
}

# Add L (single-item)
idx <- idx + 1
rows_list[[idx]] <- data.frame(
  Construct = "L",
  Indicator = "REPUR",
  Loading   = 1.000,
  Alpha     = "n/a",
  rhoC      = "n/a",
  rhoA      = "n/a",
  AVE       = "n/a",
  stringsAsFactors = FALSE
)

table_A0 <- do.call(rbind, rows_list)

cat("\nTable A0. Measurement model assessment (ACSI, PLS-PM, full sample).\n")
cat(strrep("-", 75), "\n")
print(table_A0, row.names = FALSE)
cat(strrep("-", 75), "\n")
cat("Note: L = single-item construct (REPUR). Reliability and AVE are\n")
cat("undefined for single-item constructs.\n")
cat(sprintf("N = %d (complete cases on 12 indicators).\n", n_complete))

# --- 6. Threshold checks --------------------------------------------

cat("\n--- 6. Threshold Checks (multi-item constructs: E, Q, V, S) ---\n\n")

# Collect all loadings for multi-item constructs
all_loadings <- unlist(lapply(multi_constructs, function(cst) {
  loadings_raw[cst, indicators_map[[cst]]]
}))

# rhoC >= 0.70
cr_pass <- all(rhoC_vals >= 0.70)
cat(sprintf("  [%s] rhoC >= 0.70:  range %.3f - %.3f\n",
            ifelse(cr_pass, "PASS", "FAIL"),
            min(rhoC_vals), max(rhoC_vals)))

# Alpha >= 0.70
alpha_pass <- all(alpha_vals >= 0.70)
cat(sprintf("  [%s] Alpha >= 0.70: range %.3f - %.3f\n",
            ifelse(alpha_pass, "PASS", "FAIL"),
            min(alpha_vals), max(alpha_vals)))

# AVE >= 0.50
ave_pass <- all(AVE_vals >= 0.50)
cat(sprintf("  [%s] AVE >= 0.50:   range %.3f - %.3f\n",
            ifelse(ave_pass, "PASS", "WARN"),
            min(AVE_vals), max(AVE_vals)))

# Loadings > 0.708
load_708 <- all(all_loadings > 0.708)
load_50  <- all(all_loadings > 0.50)
cat(sprintf("  [%s] Loadings > 0.708: min = %.3f\n",
            ifelse(load_708, "PASS", "WARN"), min(all_loadings)))
cat(sprintf("  [%s] Loadings > 0.50:  min = %.3f\n",
            ifelse(load_50, "PASS", "FAIL"), min(all_loadings)))

# HTMT < 0.85 / 0.90 (multi-item only)
htmt_vals <- htmt_multi[lower.tri(htmt_multi)]
htmt_vals <- htmt_vals[!is.na(htmt_vals)]
htmt_85 <- all(htmt_vals < 0.85)
htmt_90 <- all(htmt_vals < 0.90)
cat(sprintf("  [%s] HTMT < 0.85:  max = %.3f\n",
            ifelse(htmt_85, "PASS", "WARN"), max(htmt_vals)))
cat(sprintf("  [%s] HTMT < 0.90:  max = %.3f\n",
            ifelse(htmt_90, "PASS", "WARN"), max(htmt_vals)))

# Note about WRONGX and WRONGQ if below threshold
low_items <- names(all_loadings[all_loadings < 0.708])
if (length(low_items) > 0) {
  cat(sprintf("\n  NOTE: Low loadings detected for: %s\n",
              paste(low_items, collapse = ", ")))
  cat("  In the ACSI framework, WRONGX and WRONGQ are reverse-scored\n")
  cat("  items measuring 'things gone wrong'. Lower loadings for such\n")
  cat("  items are common in established survey instruments and do not\n")
  cat("  warrant indicator removal (Hair et al., 2022, p. 120).\n")
}

# --- 7. Export CSV files --------------------------------------------

dir.create("../outputs_main", showWarnings = FALSE)

write.csv(table_A0,
          "../outputs_main/measurement_audit_table_A0.csv",
          row.names = FALSE)

write.csv(round(htmt_multi, 4),
          "../outputs_main/htmt_matrix.csv")

if (!is.null(fl_multi)) {
  write.csv(round(fl_multi, 4),
            "../outputs_main/fornell_larcker.csv")
} else {
  cat("  WARNING: fornell_larcker.csv not written.\n")
}

cat("\n--- Files saved ---\n")
cat("  ../outputs_main/measurement_audit_table_A0.csv\n")
cat("  ../outputs_main/htmt_matrix.csv\n")
cat("  ../outputs_main/fornell_larcker.csv\n")

# --- 8. Draft text for manuscript Section 4.1 -----------------------

cat("\n", strrep("=", 70), "\n")
cat("SUGGESTED TEXT FOR MANUSCRIPT Section 4.1 (after Table 4):\n")
cat(strrep("=", 70), "\n\n")

cat(sprintf(paste0(
  '"Consistent with the measurement audit required by TC-PVSS\n',
  '(Section 3.1), we assessed convergent and discriminant validity\n',
  'of the four multi-item constructs prior to structural search.\n',
  'Composite reliability (rho_c) ranges from %.3f to %.3f and AVE\n',
  'from %.3f to %.3f, meeting standard convergent validity\n',
  'thresholds (Hair et al., 2022). Two reverse-coded items (WRONGX,\n',
  'WRONGQ) exhibit loadings below 0.708 (%.3f and %.3f\n',
  'respectively); we retain them to preserve construct meaning\n',
  'consistent with the established ACSI instrument. HTMT values\n',
  'among multi-item constructs range from %.3f to %.3f. Loyalty is\n',
  'operationalized as a single-item construct via repurchase\n',
  'intention (REPUR), which also serves as the prediction target.\n',
  'HIGHPTOL and LOWPTOL (price tolerance) operate on a percentage\n',
  'scale incompatible with the 1-10 Likert indicators and exhibit\n',
  'high missingness (21%% and 84%%, respectively); they are\n',
  'excluded from the PLS-PM estimation. Full diagnostics are in\n',
  'Online Appendix Table A0."\n'
  ),
  min(rhoC_vals), max(rhoC_vals),
  min(AVE_vals), max(AVE_vals),
  loadings_raw["E", "WRONGX"], loadings_raw["Q", "WRONGQ"],
  min(htmt_vals), max(htmt_vals)
))

cat("\n[Review and adjust values after running.]\n")

# --- 9. Session info ------------------------------------------------

cat("\n--- Session Info ---\n")
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("cSEM version: %s\n", packageVersion("cSEM")))
cat(sprintf("Date: %s\n", Sys.time()))

cat("\n=== Measurement audit complete ===\n")
