#!/usr/bin/env Rscript
# run_all.R -- reproduces every analysis, table and figure of the Kyrgyzstan national BSI-AMR survey.
#
#   Rscript run_all.R                      # full run (~30-40 min; the models script dominates)
#   Rscript run_all.R --from=09            # restart from a step
#   Rscript run_all.R --only=12,13,14      # selected steps
#   KGZ_DATA_DIR=/path/to/inputs Rscript run_all.R
#
# Each step is run in its own clean R process; logs are written to logs/. The run stops at the first error.
args <- commandArgs(trailingOnly = TRUE)
this <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
PKG <- normalizePath(if (!is.na(this) && nzchar(this)) dirname(this) else getwd())
Sys.setenv(KGZ_PKG = PKG)
setwd(PKG)
steps <- c(                                   # step id, script, what it does
  "01" = "01_pathogens_ast.R",                # species naming (Species 1 -> 2 -> commensal field -> comment), fungal exclusion, AST long table
  "02" = "02_resistance_defs.R",              # harmonised S/I/R (I counted as resistant), MRSA inference, MDR/XDR/PDR, AWaRe, all-Gram-negative carbapenem/cephalosporin flags
  "03" = "03_cohort_weights.R",               # case-report cleaning, design weights, regional calibration (raking), isolate linkage
  "04" = "04_prevalence_incidence.R",         # CONSORT table, pathogen distribution, prevalence among isolates, incidence (1 in x admissions)
  "05" = "05_heatmap.R",                      # pathogen x antibiotic-class heatmap
  "06" = "06_extra_outputs.R",                # subgroup resistance, oblast incidence, Figures 2, 4-8
  "07" = "07_export_definitions.R",           # definitions tables (Magiorakos categories, AST codes, drugs)
  "08" = "08_make_legacy_lab.R",              # builds the confirmed-isolate table in the layout the models script expects
  "09" = "09_legacy_models.R",                # risk factors, mortality (IPW Cox MSM, Fine-Gray, exact pairs), LOS, stewardship, CRF value, variable reduction
  "11" = "11_composites_and_tiff.R",          # two-panel composites; PNG -> TIFF (500 dpi)
  "10" = "10_mortality_bars.R",               # Figure 9 and Supplementary crude-mortality bars
  "12" = "12_sexage_figures.R",               # sex/age-stratified mortality, LOS and prevalence figures
  "13" = "13_fig7_incidence_1_in_x.R",        # Figure 7 (calibrated incidence as 1 in x admissions)
  "14" = "14_collect_figures.R",              # collect final TIFFs into outputs/figures_TIFF_500dpi
  "15" = "15_headline_check.R")               # self-check: headline numbers must match the manuscript
sel <- names(steps)
from <- sub("^--from=", "", grep("^--from=", args, value = TRUE)); only <- sub("^--only=", "", grep("^--only=", args, value = TRUE))
if (length(only)) sel <- intersect(sel, strsplit(only, ",")[[1]]) else if (length(from)) sel <- sel[match(from, sel):length(sel)]
dir.create("logs", showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), "logs/sessionInfo.txt")
rscript <- file.path(R.home("bin"), "Rscript")
for (id in sel) {
  f <- steps[[id]]; t0 <- Sys.time(); cat(sprintf("[%s] %s ... ", id, f)); flush.console()
  status <- system2(rscript, c("--vanilla", shQuote(file.path(PKG, "R", f))), stdout = file.path(PKG, "logs", sub("\\.R$", ".log", f)), stderr = file.path(PKG, "logs", sub("\\.R$", ".log", f)))
  cat(sprintf("%s (%.1f min)\n", if (status == 0) "ok" else "FAILED", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (status != 0) stop("Step ", id, " failed; see logs/", sub("\\.R$", ".log", f))
}
cat("All requested steps finished.\n")
