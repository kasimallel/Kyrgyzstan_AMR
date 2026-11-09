# config.R -- paths and shared helpers. Every script sources this file; nothing else holds a path.
suppressPackageStartupMessages({library(readxl); library(dplyr); library(tidyr); library(stringr); library(purrr)})
PKG <- Sys.getenv("KGZ_PKG", unset = "")
if (PKG == "") PKG <- getwd()
PKG <- normalizePath(PKG, mustWork = TRUE)
DATA_DIR <- normalizePath(Sys.getenv("KGZ_DATA_DIR", unset = file.path(PKG, "data")), mustWork = FALSE)

# ---- input files (not distributed with the code; see README) ----
XLSX         <- file.path(DATA_DIR, "Database_with_confirmed_ID_AMR_survey_KYR_261026.xlsx")   # central-laboratory confirmed workbook
DHIS2_XLSX   <- file.path(DATA_DIR, "amr-survey_kgz_2026-08-25.xlsx")                          # DHIS2 extract (case report, cohort, discharge forms)
SAMPLING_CSV <- file.path(DATA_DIR, "NewSamplingFrameDec2022.csv")                             # PPS sampling frame (facility ApCL)
for (f in c(XLSX, DHIS2_XLSX, SAMPLING_CSV)) if (!file.exists(f)) stop("Input file not found: ", f, "\n  Set KGZ_DATA_DIR to the folder holding the three input files (see README.md).")

# ---- output layout (mirrors the original project) ----
ROOT <- PKG; V2 <- PKG
OUT <- file.path(PKG, "outputs")
for (d in c(OUT, file.path(OUT, "figures_TIFF_500dpi"), file.path(PKG, "legacy", "out"), file.path(PKG, "legacy", "figs"), file.path(PKG, "logs"))) dir.create(d, showWarnings = FALSE, recursive = TRUE)
rd <- function(sheet) read_excel(XLSX, sheet = sheet, guess_max = 1e6)   # guess_max avoids type-guess corruption of sparse columns
