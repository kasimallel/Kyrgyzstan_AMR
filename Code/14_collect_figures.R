# 14_collect_figures.R -- gather the final manuscript/supplement figures (all TIFF, 500 dpi) into outputs/figures_TIFF_500dpi
source(file.path(Sys.getenv("KGZ_PKG"), "R", "config.R"))
dest <- file.path(OUT, "figures_TIFF_500dpi")
from_outputs <- c("Fig2_pathogen_distribution", "Fig4_resistance_phenotypes", "Fig5_sex_age_phenotypes", "Fig6_prevalence_maps", "Fig8_incidence_maps", "05_heatmap_by_class")
from_legacy  <- c("Figure9_mortality_crude", "Figure10_mortality_KM_curves", "LOS_forest_all_mechanisms_bigfont", "LOS_forest_all_mechanisms_IPWonly", "LOS_forest_all_mechanisms",
                  "Mortality_propensity_overlap_all_mechanisms", "SuppFigure_mortality_crude_calibrated", "SuppFigure_mortality_KM_calibrated", "SuppFigure_mortality_adjusted_CoxMSM_IPW",
                  "SuppFigure_mortality_adjusted_FineGray_IPW", "Mortality_carbapenemR_gramneg_subgroup_KM", "Mortality_cephalosporinR_gramneg_subgroup_KM", "Mortality_3GCR_gramneg_subgroup_KM",
                  "Mortality_adjusted_curves_by_mechanism")
cp <- function(src, name) if (file.exists(src)) file.copy(src, file.path(dest, paste0(name, ".tiff")), overwrite = TRUE) else message("missing: ", src)
for (f in from_outputs) cp(file.path(OUT, paste0(f, ".tiff")), f)
for (f in from_legacy)  cp(file.path(PKG, "legacy", "figs", paste0(f, ".tiff")), f)
# the three sex/age figures and Figure 7 are written straight into the destination folder by scripts 12 and 13
cat("figures in", dest, ":", length(list.files(dest, pattern = "\\.tiff$")), "\n")
