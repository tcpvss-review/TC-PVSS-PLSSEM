# ─────────────────────────────────────────────────────────────────────
# 02_make_figures_simulation.R
# Reproduce SIMULATION figures (Figures S1–S4) from full outputs.
#
# Data sources (all from ../sim_outputs/):
#   FigS1 — Table_S2_optimism.csv         (pipeline miscalibration)
#   FigS2 — sim_main.csv                  (modal-evaluation FPR)
#   FigS3 — Table_S4_stable_core_tau.csv  (stable-core correctness)
#   FigS4 — Table_S4_stable_core_tau.csv  (stable-core RMSE)
#
# Outputs: ../figures_simulation/ as 600-DPI PNG + TIFF
# ─────────────────────────────────────────────────────────────────────
rm(list = ls())

in_dir  <- file.path("..", "sim_outputs")
out_dir <- file.path("..", "figures_simulation")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Global palette ───────────────────────────────────────────────────
col_naive  <- "#2C6FAC"   # steel-blue   (naive CV)
col_nested <- "#D97A2B"   # burnt-orange (nested CV / TC-PVSS)
col_M0     <- "#2C6FAC"   # blue   (Truth = M0)
col_M6     <- "#D97A2B"   # orange (Truth = M6)
col_M7     <- "#3A9A5B"   # green  (Truth = M7)
col_grid   <- "#E8E8E8"   # subtle grid

pch_naive  <- 16  # filled circle
pch_nested <- 17  # filled triangle
pch_M0     <- 16
pch_M6     <- 17
pch_M7     <- 15  # filled square

# ── Helper: publication-grade par defaults ───────────────────────────
pub_par <- function(...) {
  par(family = "sans",
      mar    = c(4.2, 5.0, 2.5, 1.0),
      mgp    = c(3.0, 0.7, 0),
      tcl    = -0.3,
      las    = 1,
      cex.main = 1.05,
      cex.lab  = 0.95,
      cex.axis = 0.85,
      ...)
}

# ── Helper: save dual format ─────────────────────────────────────────
save_fig <- function(basename, width = 6, height = 4.5, expr) {
  for (ext in c("png", "tif")) {
    if (ext == "png") {
      png(file.path(out_dir, paste0(basename, ".png")),
          width = width, height = height, units = "in", res = 600,
          bg = "white")
    } else {
      tiff(file.path(out_dir, paste0(basename, ".tif")),
           width = width, height = height, units = "in", res = 600,
           compression = "lzw", bg = "white")
    }
    pub_par()
    eval(expr)
    dev.off()
  }
}

# ── Helper: subtle grid ──────────────────────────────────────────────
add_hgrid <- function(at) abline(h = at, col = col_grid, lwd = 0.5)

# ── Load data ────────────────────────────────────────────────────────
s2  <- read.csv(file.path(in_dir, "Table_S2_optimism.csv"))
s4  <- read.csv(file.path(in_dir, "Table_S4_stable_core_tau.csv"))
sim <- read.csv(file.path(in_dir, "sim_main.csv"))
sim_ok <- subset(sim, status == "ok")

# =====================================================================
#  Figure S2 — FPR by sample size (MODAL evaluation)
#  Manuscript Sec 6.4.2: "nested CV reduces FPR from 0.265 to 0.223"
#  Uses nested_modal_FPR from sim_main.csv (single most-frequent model)
# =====================================================================
fpr_by_N <- aggregate(
  cbind(naive_FPR, nested_modal_FPR) ~ N, data = sim_ok,
  FUN = function(x) mean(x, na.rm = TRUE)
)

save_fig("FigS2_FPR_by_N", expr = quote({
  yr <- range(c(fpr_by_N$naive_FPR, fpr_by_N$nested_modal_FPR))
  yr <- yr + c(-0.015, 0.015)
  plot(fpr_by_N$N, fpr_by_N$naive_FPR,
       type = "b", pch = pch_naive, col = col_naive, lwd = 2, cex = 1.3,
       xlab = "Sample size (N)", ylab = "False-positive rate (FPR)",
       main = "FPR by sample size (modal evaluation)",
       ylim = yr, xaxt = "n")
  axis(1, at = fpr_by_N$N)
  add_hgrid(seq(0.10, 0.40, 0.05))
  # Redraw over grid
  lines(fpr_by_N$N, fpr_by_N$naive_FPR,
        type = "b", pch = pch_naive, col = col_naive, lwd = 2, cex = 1.3)
  lines(fpr_by_N$N, fpr_by_N$nested_modal_FPR,
        type = "b", pch = pch_nested, col = col_nested, lwd = 2, cex = 1.3)
  # Value labels
  text(fpr_by_N$N, fpr_by_N$naive_FPR + 0.008,
       sprintf("%.3f", fpr_by_N$naive_FPR), cex = 0.65, col = col_naive)
  text(fpr_by_N$N, fpr_by_N$nested_modal_FPR - 0.008,
       sprintf("%.3f", fpr_by_N$nested_modal_FPR), cex = 0.65, col = col_nested)
  legend("topleft",
         legend = c(expression("Na\u00efve CV"),
                    "Nested CV (TC-PVSS, modal)"),
         col = c(col_naive, col_nested), lwd = 2,
         pch = c(pch_naive, pch_nested), pt.cex = 1.1,
         bty = "n", cex = 0.85)
  box(bty = "l")
}))

# =====================================================================
#  Figure S1 — Absolute miscalibration by sample size
#  |reported RMSE − true RMSE| averaged over conditions within each N
# =====================================================================
s2$abs_naive  <- abs(s2$naive_optimism_rmse)
s2$abs_nested <- abs(s2$nested_optimism_rmse)
bias_by_N <- aggregate(cbind(abs_naive, abs_nested) ~ N, data = s2, FUN = mean)

save_fig("FigS1_AbsMiscalibration_by_N", expr = quote({
  yr <- range(c(bias_by_N$abs_naive, bias_by_N$abs_nested))
  yr <- c(0, yr[2] * 1.15)
  plot(bias_by_N$N, bias_by_N$abs_naive,
       type = "b", pch = pch_naive, col = col_naive, lwd = 2, cex = 1.3,
       xlab = "Sample size (N)",
       ylab = expression("|reported " - " true| (RMSE)"),
       main = "Absolute miscalibration by sample size",
       ylim = yr, xaxt = "n")
  axis(1, at = bias_by_N$N)
  add_hgrid(seq(0, 0.025, 0.005))
  lines(bias_by_N$N, bias_by_N$abs_naive,
        type = "b", pch = pch_naive, col = col_naive, lwd = 2, cex = 1.3)
  lines(bias_by_N$N, bias_by_N$abs_nested,
        type = "b", pch = pch_nested, col = col_nested, lwd = 2, cex = 1.3)
  # Value labels
  text(bias_by_N$N, bias_by_N$abs_naive + yr[2] * 0.04,
       sprintf("%.4f", bias_by_N$abs_naive), cex = 0.65, col = col_naive)
  text(bias_by_N$N, bias_by_N$abs_nested - yr[2] * 0.04,
       sprintf("%.4f", bias_by_N$abs_nested), cex = 0.65, col = col_nested)
  legend("topright",
         legend = c(expression("Na\u00efve CV"),
                    "Nested CV (TC-PVSS)"),
         col = c(col_naive, col_nested), lwd = 2,
         pch = c(pch_naive, pch_nested), pt.cex = 1.1,
         bty = "n", cex = 0.85)
  box(bty = "l")
}))

# =====================================================================
#  Figure S3 — Stable-core correctness vs tau
#  P(selected core model = true DGP) weighted by cell n
# =====================================================================
truths <- sort(unique(s4$truth))
taus   <- sort(unique(s4$tau))
corr_df <- data.frame()
for (tr in truths) {
  for (ta in taus) {
    sub <- subset(s4, truth == tr & tau == ta)
    pc  <- sum((sub$core_model == sub$truth) * sub$n) / sum(sub$n)
    corr_df <- rbind(corr_df, data.frame(truth = tr, tau = ta, prop_correct = pc))
  }
}
col_truth <- c(col_M0, col_M6, col_M7)
pch_truth <- c(pch_M0, pch_M6, pch_M7)

save_fig("FigS3_StableCore_Correctness_vs_tau", expr = quote({
  plot(NULL, xlim = range(taus) + c(-0.02, 0.02),
       ylim = c(0, 1),
       xlab = expression("Stability threshold " * tau),
       ylab = "P(stable-core model = truth)",
       main = expression("Stable-core correctness vs " * tau),
       xaxt = "n")
  axis(1, at = taus)
  add_hgrid(seq(0.1, 0.9, 0.1))
  for (i in seq_along(truths)) {
    sub <- subset(corr_df, truth == truths[i])
    lines(sub$tau, sub$prop_correct,
          type = "b", pch = pch_truth[i], col = col_truth[i],
          lwd = 2, cex = 1.3)
  }
  legend("bottomleft",
         legend = paste0("Truth = ", truths),
         col = col_truth, lwd = 2,
         pch = pch_truth, pt.cex = 1.1,
         bty = "n", cex = 0.85)
  box(bty = "l")
}))

# =====================================================================
#  Figure S4 — Stable-core true RMSE vs tau
# =====================================================================
rmse_df <- aggregate(true_rmse ~ truth + tau, data = s4,
                     FUN = function(x) mean(x, na.rm = TRUE))

save_fig("FigS4_StableCore_RMSE_vs_tau", expr = quote({
  yr <- range(rmse_df$true_rmse)
  yr <- yr + c(-1, 1) * diff(yr) * 0.15
  plot(NULL, xlim = range(taus) + c(-0.02, 0.02),
       ylim = yr,
       xlab = expression("Stability threshold " * tau),
       ylab = "True RMSE (independent test set)",
       main = expression("Stable-core predictive performance vs " * tau),
       xaxt = "n")
  axis(1, at = taus)
  add_hgrid(pretty(yr, 8))
  for (i in seq_along(truths)) {
    sub <- subset(rmse_df, truth == truths[i])
    lines(sub$tau, sub$true_rmse,
          type = "b", pch = pch_truth[i], col = col_truth[i],
          lwd = 2, cex = 1.3)
  }
  legend("topright",
         legend = paste0("Truth = ", truths),
         col = col_truth, lwd = 2,
         pch = pch_truth, pt.cex = 1.1,
         bty = "n", cex = 0.85)
  box(bty = "l")
}))

cat("Done. Wrote simulation figures to:", out_dir, "\n")
