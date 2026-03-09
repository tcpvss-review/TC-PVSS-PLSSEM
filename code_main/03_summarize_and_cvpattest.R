# code_main/03_summarize_and_cvpattest.R
# ------------------------------------------------------------------
# TC-PVSS (Protocol v1) — Summarize nested CV results
# - Pooled and split-level predictive performance
# - Model selection frequencies and edge inclusion probabilities
# - CVPAT tests: split-level (primary) and per-case (supplementary)
# - lm benchmark comparison (selected vs indicator-level lm)
# - By-industry performance breakdown
# - Diagnostic figures (histograms, boxplots)
#
# Input:  ../outputs_main/pred_errors_long.parquet  (from 02_nested_cv_search.R)
# Output: ../outputs_main/summary_tables.csv
#         ../outputs_main/model_frequency.csv
#         ../outputs_main/performance_by_outer_split.csv
#         ../outputs_main/performance_by_industry.csv
#         ../outputs_main/fig_loss_diff_mse_hist.png
#         ../outputs_main/fig_loss_diff_ae_hist.png
#         ../outputs_main/fig_rmse_boxplot_outer.png
#         ../outputs_main/fig_mae_boxplot_outer.png
# ------------------------------------------------------------------

library(dplyr)
library(tibble)
library(arrow)
library(readr)

TARGET <- "REPUR"

res_long <- read_parquet("../outputs_main/pred_errors_long.parquet")

cat("=============================================================\n")
cat("  TC-PVSS: Summarize Nested CV Results\n")
cat("=============================================================\n\n")
cat(sprintf("Loaded %d prediction rows from %d outer splits.\n\n",
            nrow(res_long),
            n_distinct(res_long$outer_split)))

# =====================================================================
# 1) Overall pooled predictive performance (selected vs baseline)
# =====================================================================
perf_pooled <- res_long %>%
  summarise(
    RMSE_selected = sqrt(mean(se_chosen, na.rm = TRUE)),
    MAE_selected  = mean(ae_chosen, na.rm = TRUE),
    RMSE_baseline = sqrt(mean(se_base, na.rm = TRUE)),
    MAE_baseline  = mean(ae_base, na.rm = TRUE),
    RMSE_lm       = sqrt(mean(se_lm, na.rm = TRUE)),
    MAE_lm        = mean(ae_lm, na.rm = TRUE),
    n_predictions = n()
  )

cat("--- 1. Pooled performance ---\n")
print(as.data.frame(perf_pooled))

# =====================================================================
# 2) Per-split performance
# =====================================================================
perf_split <- res_long %>%
  group_by(outer_split) %>%
  summarise(
    model_id      = first(model_id),
    E_L           = first(E_L),
    Q_L           = first(Q_L),
    V_L           = first(V_L),
    complexity    = first(complexity),
    RMSE_selected = sqrt(mean(se_chosen, na.rm = TRUE)),
    MAE_selected  = mean(ae_chosen, na.rm = TRUE),
    RMSE_baseline = sqrt(mean(se_base, na.rm = TRUE)),
    MAE_baseline  = mean(ae_base, na.rm = TRUE),
    RMSE_lm       = sqrt(mean(se_lm, na.rm = TRUE)),
    MAE_lm        = mean(ae_lm, na.rm = TRUE),
    n_test        = n(),
    .groups       = "drop"
  )

write_csv(perf_split, "../outputs_main/performance_by_outer_split.csv")

# =====================================================================
# 3) Model selection frequency (unit = outer split, n = 50)
# =====================================================================
model_freq <- perf_split %>%
  count(model_id, sort = TRUE) %>%
  mutate(freq = n / sum(n))

cat("\n--- 3. Model selection frequency (per outer split) ---\n")
print(as.data.frame(model_freq))

write_csv(model_freq, "../outputs_main/model_frequency.csv")

# =====================================================================
# 4) Stability: edge inclusion probabilities (unit = outer split)
# =====================================================================
stability <- perf_split %>%
  summarise(
    n_splits = n(),
    p_E_L    = mean(E_L),
    p_Q_L    = mean(Q_L),
    p_V_L    = mean(V_L)
  )

cat("\n--- 4. Edge inclusion probabilities ---\n")
cat(sprintf("  E→L: %.2f (%d/%d splits)\n",
            stability$p_E_L,
            sum(perf_split$E_L),
            stability$n_splits))
cat(sprintf("  Q→L: %.2f (%d/%d splits)\n",
            stability$p_Q_L,
            sum(perf_split$Q_L),
            stability$n_splits))
cat(sprintf("  V→L: %.2f (%d/%d splits)\n",
            stability$p_V_L,
            sum(perf_split$V_L),
            stability$n_splits))

# =====================================================================
# 4b) lm benchmark: pooled performance
# =====================================================================
cat("\n--- 4b. lm benchmark (indicator-level linear regression) ---\n")
cat(sprintf("  lm RMSE: %.4f  |  lm MAE: %.4f\n",
            perf_pooled$RMSE_lm, perf_pooled$MAE_lm))
cat(sprintf("  Selected vs lm:  RMSE diff = %.4f  |  MAE diff = %.4f\n",
            perf_pooled$RMSE_selected - perf_pooled$RMSE_lm,
            perf_pooled$MAE_selected - perf_pooled$MAE_lm))
cat("  (Negative = selected model is better than lm)\n")

# =====================================================================
# 5) CVPAT tests
# =====================================================================

# Helper function
cvpat_test <- function(loss_selected, loss_baseline, label = "") {
  D    <- loss_selected - loss_baseline
  Dbar <- mean(D, na.rm = TRUE)
  S2   <- var(D, na.rm = TRUE)
  N    <- sum(!is.na(D))
  Tval <- Dbar / sqrt(S2 / N)
  pval <- 2 * pt(abs(Tval), df = N - 1, lower.tail = FALSE)
  tibble(test = label, N = N, mean_diff = Dbar, T_stat = Tval, p_value = pval)
}

# --- 5a. PRIMARY: Split-level paired test (n = 50, independent) ---
# Each outer split contributes one independent RMSE/MAE pair
cat("\n--- 5a. CVPAT: Split-level paired tests (PRIMARY, n = 50) ---\n")

# MSE-loss at split level
split_mse <- perf_split %>%
  mutate(mse_sel  = RMSE_selected^2,
         mse_base = RMSE_baseline^2,
         mse_lm   = RMSE_lm^2)
cvpat_split_mse <- cvpat_test(split_mse$mse_sel, split_mse$mse_base,
                              "split_MSE_vs_M0")

# MAE-loss at split level
cvpat_split_mae <- cvpat_test(perf_split$MAE_selected, perf_split$MAE_baseline,
                              "split_MAE_vs_M0")

# RMSE paired test (descriptive supplement)
cvpat_split_rmse <- cvpat_test(perf_split$RMSE_selected, perf_split$RMSE_baseline,
                               "split_RMSE_vs_M0")

# --- lm benchmark comparisons (split-level) ---
cvpat_split_mse_lm <- cvpat_test(split_mse$mse_sel, split_mse$mse_lm,
                                  "split_MSE_vs_lm")
cvpat_split_mae_lm <- cvpat_test(perf_split$MAE_selected, perf_split$MAE_lm,
                                  "split_MAE_vs_lm")

split_tests <- bind_rows(cvpat_split_mse, cvpat_split_mae, cvpat_split_rmse,
                          cvpat_split_mse_lm, cvpat_split_mae_lm)
print(as.data.frame(split_tests), digits = 4)

# --- 5b. SUPPLEMENTARY: Per-case pooled test (n = 41,195) ---
# NOTE: Cases repeat across 5 CV repeats, so observations are NOT
# independent. Standard errors are underestimated and p-values are
# anti-conservative. Reported for comparability with prior literature
# but split-level tests (5a) are the primary inferential evidence.
cat("\n--- 5b. CVPAT: Per-case pooled tests (SUPPLEMENTARY, n =",
    nrow(res_long), ") ---\n")
cat("NOTE: Per-case tests pool repeated predictions and are NOT independent.\n")
cat("      Use split-level tests (5a) as primary evidence.\n\n")

cvpat_case_mse <- cvpat_test(res_long$se_chosen, res_long$se_base,
                             "case_MSE_vs_M0")
cvpat_case_mae <- cvpat_test(res_long$ae_chosen, res_long$ae_base,
                             "case_MAE_vs_M0")

# --- lm benchmark comparisons (per-case) ---
cvpat_case_mse_lm <- cvpat_test(res_long$se_chosen, res_long$se_lm,
                                 "case_MSE_vs_lm")
cvpat_case_mae_lm <- cvpat_test(res_long$ae_chosen, res_long$ae_lm,
                                 "case_MAE_vs_lm")

case_tests <- bind_rows(cvpat_case_mse, cvpat_case_mae,
                         cvpat_case_mse_lm, cvpat_case_mae_lm)
print(as.data.frame(case_tests), digits = 4)

# =====================================================================
# 6) Export combined summary
# =====================================================================
all_tests <- bind_rows(split_tests, case_tests)

summary_out <- bind_rows(
  perf_pooled %>% mutate(section = "pooled_performance"),
  stability   %>% mutate(section = "stability"),
  all_tests   %>% mutate(section = "cvpat_tests")
)

write_csv(summary_out, "../outputs_main/summary_tables.csv")

# =====================================================================
# 7) By-industry performance (selected vs baseline vs lm)
# =====================================================================
cat("\n--- 7. Performance by industry ---\n")

by_industry <- res_long %>%
  group_by(industry) %>%
  summarise(
    RMSE_selected = sqrt(mean(se_chosen, na.rm = TRUE)),
    MAE_selected  = mean(ae_chosen, na.rm = TRUE),
    RMSE_baseline = sqrt(mean(se_base, na.rm = TRUE)),
    MAE_baseline  = mean(ae_base, na.rm = TRUE),
    RMSE_lm       = sqrt(mean(se_lm, na.rm = TRUE)),
    MAE_lm        = mean(ae_lm, na.rm = TRUE),
    n_predictions = n(),
    .groups       = "drop"
  )

write_csv(by_industry, "../outputs_main/performance_by_industry.csv")
print(as.data.frame(by_industry))

# =====================================================================
# 8) Diagnostic figures (base R -- no extra dependencies)
# =====================================================================
cat("\n--- 8. Generating diagnostic figures ---\n")

# 8a. Loss difference histograms (per-case, descriptive)
loss_diff <- res_long %>%
  mutate(d_mse = se_chosen - se_base,
         d_ae  = ae_chosen - ae_base)

png("../outputs_main/fig_loss_diff_mse_hist.png", width = 1200, height = 800, res = 150)
hist(loss_diff$d_mse, breaks = 60,
     main = "Loss difference (MSE): selected - baseline", xlab = "d_mse")
abline(v = 0, lty = 2, col = "red")
dev.off()

png("../outputs_main/fig_loss_diff_ae_hist.png", width = 1200, height = 800, res = 150)
hist(loss_diff$d_ae, breaks = 60,
     main = "Loss difference (AE): selected - baseline", xlab = "d_ae")
abline(v = 0, lty = 2, col = "red")
dev.off()

# 8b. Per-split RMSE/MAE boxplots (selected vs baseline vs lm)
png("../outputs_main/fig_rmse_boxplot_outer.png", width = 1200, height = 800, res = 150)
boxplot(perf_split$RMSE_selected, perf_split$RMSE_baseline, perf_split$RMSE_lm,
        names = c("Selected", "Baseline (M0)", "lm"),
        main = "RMSE across outer splits (n = 50)", ylab = "RMSE",
        col = c("#4DAF4A", "#377EB8", "#FF7F00"))
dev.off()

png("../outputs_main/fig_mae_boxplot_outer.png", width = 1200, height = 800, res = 150)
boxplot(perf_split$MAE_selected, perf_split$MAE_baseline, perf_split$MAE_lm,
        names = c("Selected", "Baseline (M0)", "lm"),
        main = "MAE across outer splits (n = 50)", ylab = "MAE",
        col = c("#4DAF4A", "#377EB8", "#FF7F00"))
dev.off()

cat("  Figures saved to ../outputs_main/fig_*.png\n")

# =====================================================================
# Done
# =====================================================================
cat("\n--- Files saved ---\n")
cat("  ../outputs_main/summary_tables.csv\n")
cat("  ../outputs_main/model_frequency.csv\n")
cat("  ../outputs_main/performance_by_outer_split.csv\n")
cat("  ../outputs_main/performance_by_industry.csv\n")
cat("  ../outputs_main/fig_loss_diff_mse_hist.png\n")
cat("  ../outputs_main/fig_loss_diff_ae_hist.png\n")
cat("  ../outputs_main/fig_rmse_boxplot_outer.png\n")
cat("  ../outputs_main/fig_mae_boxplot_outer.png\n")

cat("\n=== Done ===\n")
