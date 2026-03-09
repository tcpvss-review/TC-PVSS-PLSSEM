# 00_config.R
# ============================================================
# TC-PVSS Simulation B — FULL + Parallel + Checkpoints
# Default mode: FULL. Switch to pilot by:
#   Sys.setenv(TCPVSS_MODE="pilot")  # before sourcing this file
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(purrr); library(tidyr); library(rsample)
  library(cSEM); library(arrow); library(glue)
})

set.seed(20260301)

MODE <- tolower(Sys.getenv("TCPVSS_MODE", "full"))
if (!MODE %in% c("pilot","full")) MODE <- "full"
message("TCPVSS_MODE = ", MODE)

# --------------------------
# Factors + budgets
# --------------------------
if (MODE == "pilot") {
  # Fast sanity run
  N_list      <- c(200, 500)
  lambda_list <- c(0.70, 0.85)
  rho_list    <- c(0.0)
  truth_list  <- c("M0","M6","M7")

  R_reps      <- 30
  N_test      <- 2000

  OUTER_V       <- 10
  OUTER_REPEATS <- 1
  INNER_V       <- 5
  INNER_REPEATS <- 1
} else {
  # FULL (paper run)
  N_list      <- c(200, 500, 1000)
  lambda_list <- c(0.70, 0.85)
  rho_list    <- c(0.0, 0.4)
  truth_list  <- c("M0","M6","M7")

  R_reps      <- 200
  N_test      <- 5000

  OUTER_V       <- 10
  OUTER_REPEATS <- 1
  INNER_V       <- 5
  INNER_REPEATS <- 1
}

# Naïve CV uses the same fold count as the nested outer loop,
# ensuring any performance differences are attributable to the
# nesting structure, not to a resampling-budget asymmetry.
NAIVE_V       <- OUTER_V
NAIVE_REPEATS <- OUTER_REPEATS

# Stability thresholds for stable-core model
tau_list <- c(0.6, 0.7, 0.8)

TARGET_IND <- "REPUR"
optional_edges <- c("E_L","Q_L","V_L")

# --------------------------
# Parallel + checkpointing
# --------------------------
# Workers: uses up to 24, or availableCores-2, whichever smaller.
MAX_WORKERS <- 24
CHUNK_SIZE  <- 50     # 50–200 typical. 50 is safer for RAM/overhead balance.
RESUME      <- TRUE   # if TRUE: skip chunks already present in ../sim_outputs/chunks

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

# Log helper
log_line <- function(...) {
  cat(timestamp(), " | ", paste0(..., collapse=""), "\n", sep="")
}
