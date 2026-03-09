# code_main/00_prepare_data.R
# ------------------------------------------------------------------
# TC-PVSS Step 0: Load and clean ACSI 2015 (v2) data
# ------------------------------------------------------------------
# - Read ACSI v2 XLSX from ../inputs_main/
# - Recode invalid codes to NA (per frozen rules)
# - Save cleaned data to ../outputs_main/acsi_clean.rds
#
# This script runs FIRST. All downstream scripts read acsi_clean.rds.
# ------------------------------------------------------------------

library(readxl)
library(dplyr)

in_path <- "../inputs_main/ACSI Data 2015 ver2.xlsx"

acsi_raw <- read_excel(in_path)

# --- Recode rules (frozen for Protocol v1) ---
recode_upper <- function(x, upper) {
  x <- as.numeric(x)
  x[x > upper] <- NA_real_
  x
}

scale_1_10 <- c("SATIS","CONFIRM","IDEAL",
                "OVERALLX","CUSTOMX","WRONGX",
                "OVERALLQ","CUSTOMQ","WRONGQ",
                "PQ","QP",
                "REPUR")

acsi <- acsi_raw %>%
  mutate(across(all_of(scale_1_10), ~recode_upper(.x, 10))) %>%
  mutate(COMP    = recode_upper(COMP, 1)) %>%
  mutate(HANDLE  = recode_upper(HANDLE, 10)) %>%
  mutate(HIGHPTOL= recode_upper(HIGHPTOL, 26)) %>%
  mutate(LOWPTOL = recode_upper(LOWPTOL, 26))

# Keep columns needed for main analysis (core) + strata
keep <- c("INDUSTRY","YEAR",
          "OVERALLX","CUSTOMX","WRONGX",
          "OVERALLQ","CUSTOMQ","WRONGQ",
          "PQ","QP",
          "SATIS","CONFIRM","IDEAL",
          "REPUR")

acsi_core <- acsi %>% select(all_of(keep))

dir.create("../outputs_main", showWarnings = FALSE)
saveRDS(acsi_core, "../outputs_main/acsi_clean.rds")

cat("Saved cleaned data:", nrow(acsi_core), "rows x", ncol(acsi_core), "cols -> ../outputs_main/acsi_clean.rds\n")
