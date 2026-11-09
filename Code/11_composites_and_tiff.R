# 11_composites_and_tiff.R -- (R port of the former Python step)
#  * two-panel composite figures (Figure 10; Supplementary calibrated KM)
#  * every PNG produced by the models script is also written as a 500-dpi LZW TIFF into legacy/figs
suppressPackageStartupMessages(library(magick))
source(file.path(Sys.getenv("KGZ_PKG"), "R", "config.R"))
OUTL <- file.path(PKG, "legacy", "out"); FIGL <- file.path(PKG, "legacy", "figs")
write_tiff <- function(img, path) image_write(img, path = path, format = "tiff", compression = "LZW", density = "500x500")
panel <- function(a, b, out) {
  A <- image_read(file.path(OUTL, a)); B <- image_read(file.path(OUTL, b))
  h <- min(image_info(A)$height, image_info(B)$height)
  A <- image_resize(A, paste0("x", h)); B <- image_resize(B, paste0("x", h))
  fs <- round(h * 0.06)
  A <- image_annotate(A, "A", gravity = "northwest", location = "+10+5", size = fs, color = "black")
  B <- image_annotate(B, "B", gravity = "northwest", location = "+10+5", size = fs, color = "black")
  write_tiff(image_background(image_append(c(A, B)), "white"), file.path(FIGL, paste0(out, ".tiff")))
}
panel("Mortality_KM_28day_by_group_bigfont.tiff", "Mortality_MDR_subgroup_KM_bigfont.tiff", "Figure10_mortality_KM_curves")
panel("Mortality_KM_28day_calibrated_weighted_primary.tiff", "Mortality_KM_28day_calibrated_weighted_MDR.tiff", "SuppFigure_mortality_KM_calibrated")
n <- 0
for (p in sort(list.files(OUTL, pattern = "\\.png$", full.names = TRUE))) {
  base <- sub("\\.png$", "", basename(p)); t <- file.path(OUTL, paste0(base, ".tiff"))
  if (file.exists(t)) { file.copy(t, file.path(FIGL, paste0(base, ".tiff")), overwrite = TRUE); next }
  write_tiff(image_background(image_read(p), "white"), file.path(FIGL, paste0(base, ".tiff"))); n <- n + 1
}
cat("converted", n, "png -> tiff; tiff files in legacy/figs:", length(list.files(FIGL, pattern = "\\.tiff$")), "\n")
