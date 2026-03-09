# code_figures/99_session_info.R
# ------------------------------------------------------------------
# Export R session information for reproducibility documentation.
# Output: ../logs/sessionInfo.txt
# ------------------------------------------------------------------

dir.create("../logs", showWarnings = FALSE)

sink("../logs/sessionInfo.txt")
cat("TC-PVSS Reproducibility Package — Session Information\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat(strrep("=", 60), "\n\n")
sessionInfo()
sink()

cat("Session info written to ../logs/sessionInfo.txt\n")
