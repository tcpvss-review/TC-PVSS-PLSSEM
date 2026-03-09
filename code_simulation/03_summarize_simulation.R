# 03_summarize_simulation.R
# ============================================================
# Paper-ready tables from sim_main/sim_stable
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(arrow); library(readr)
})

main <- read_parquet("../sim_outputs/sim_main.parquet")
stable <- read_parquet("../sim_outputs/sim_stable.parquet")

main_ok <- main %>% filter(status == "ok")

# Table S2: optimism bias (RMSE) + true RMSE/MAE
Table_S2 <- main_ok %>%
  group_by(N, lambda, rho, truth) %>%
  summarise(
    naive_optimism_rmse  = mean(naive_optimism_rmse,  na.rm=TRUE),
    nested_optimism_rmse = mean(nested_optimism_rmse, na.rm=TRUE),
    naive_true_rmse  = mean(naive_true_rmse,  na.rm=TRUE),
    nested_true_rmse = mean(nested_true_rmse, na.rm=TRUE),
    naive_true_mae   = mean(naive_true_mae,   na.rm=TRUE),
    nested_true_mae  = mean(nested_true_mae,  na.rm=TRUE),
    n = n(),
    .groups="drop"
  )
write_csv(Table_S2, "../sim_outputs/Table_S2_optimism.csv")

# Table S3: selection accuracy (optional edges)
Table_S3 <- main_ok %>%
  group_by(N, lambda, rho, truth) %>%
  summarise(
    naive_TPR = mean(naive_TPR, na.rm=TRUE),
    naive_FPR = mean(naive_FPR, na.rm=TRUE),
    naive_F1  = mean(naive_F1,  na.rm=TRUE),
    nested_TPR = mean(nested_TPR, na.rm=TRUE),
    nested_FPR = mean(nested_FPR, na.rm=TRUE),
    nested_F1  = mean(nested_F1,  na.rm=TRUE),
    n = n(),
    .groups="drop"
  )
write_csv(Table_S3, "../sim_outputs/Table_S3_selection_accuracy.csv")

# Table S4: stable-core tau sensitivity (true error on independent test set)
stable_ok <- stable %>% filter(!is.na(tau))
Table_S4 <- stable_ok %>%
  group_by(N, lambda, rho, truth, tau, core_model) %>%
  summarise(
    true_rmse = mean(true_rmse, na.rm=TRUE),
    true_mae  = mean(true_mae,  na.rm=TRUE),
    pE = mean(pE, na.rm=TRUE),
    pQ = mean(pQ, na.rm=TRUE),
    pV = mean(pV, na.rm=TRUE),
    n = n(),
    .groups="drop"
  )
write_csv(Table_S4, "../sim_outputs/Table_S4_stable_core_tau.csv")

message("Wrote:\n- ../sim_outputs/Table_S2_optimism.csv\n- ../sim_outputs/Table_S3_selection_accuracy.csv\n- ../sim_outputs/Table_S4_stable_core_tau.csv\n")
