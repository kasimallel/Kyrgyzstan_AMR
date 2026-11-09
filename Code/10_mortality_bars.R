# 10_mortality_bars.R -- crude (Figure 9) and calibrated (Supplementary S16) 28-day mortality bar charts, same design as the previous version
source(file.path(Sys.getenv("KGZ_PKG"), "R", "config.R"))
suppressPackageStartupMessages(library(ggplot2))
LEG <- file.path(V2, "legacy/out"); FIG <- file.path(V2, "legacy/figs"); dir.create(FIG, showWarnings = FALSE)
cr <- readxl::read_excel(file.path(LEG, "Mortality_crude_by_group_all_mechanisms.xlsx"))
cw <- readxl::read_excel(file.path(LEG, "Mortality_crude_calibrated_weighted_all_mechanisms.xlsx"))
pairs <- list(c("Susceptible (any)", "Resistant (any)"), c("Non-MDR", "MDR"), c("Non-XDR", "XDR"), c("Carbapenem-S", "Carbapenem-R"), c("Cephalosporin-S", "Cephalosporin-R"), c("3GC-S", "3GCR"), c("MSSA", "MRSA"), c("VSE", "VRE"), c("Access-S", "Access-R"), c("Watch-S", "Watch-R"), c("Reserve-S", "Reserve-R"))
grp_type <- function(g) ifelse(g == "Uninfected", "Uninfected", ifelse(grepl("^(Susceptible|Non-|Carbapenem-S|Cephalosporin-S|3GC-S|MSSA|VSE|Access-S|Watch-S|Reserve-S)", g), "Susceptible", "Resistant"))
chi <- function(d, a, b, nn = "n_outcome_known", dd = "deaths_28day") { x <- d[d$group == a, ]; y <- d[d$group == b, ]
  if (!nrow(x) || !nrow(y)) return(c(NA, NA)); m <- matrix(c(x[[dd]], x[[nn]] - x[[dd]], y[[dd]], y[[nn]] - y[[dd]]), 2, byrow = TRUE)
  t <- suppressWarnings(chisq.test(m, correct = FALSE)); c(unname(t$statistic), t$p.value) }
theme_m <- theme_minimal(base_family = "Helvetica", base_size = 11) + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.x = element_line(colour = "#E6E6E6", linewidth = 0.3), legend.position = "top", axis.text.y = element_text(colour = "black"))
cols <- c(Uninfected = "#8C5E3C", Susceptible = "#4F9D69", Resistant = "#C9A227")
build <- function(d, pctcol, lab_fun, xlab, ci = FALSE, chi_fun) {
  d$grp <- factor(d$group, levels = rev(d$group)); d$type <- factor(grp_type(d$group), levels = names(cols)); d$y <- as.numeric(d$grp)
  d$xmax <- if (ci) d$calibrated_upper else d[[pctcol]]; d$lab <- lab_fun(d)
  p <- ggplot(d, aes(y = grp)) + geom_col(aes(x = .data[[pctcol]], fill = type), width = 0.72, colour = "black", linewidth = 0.3)
  if (ci) p <- p + geom_errorbar(aes(xmin = calibrated_lower, xmax = calibrated_upper), width = 0.25, orientation = "y", linewidth = 0.3)
  xr <- max(d$xmax, na.rm = TRUE)
  p <- p + geom_text(aes(x = xmax + 0.01 * xr, label = lab), hjust = 0, size = 3)
  for (pr in pairs) { a <- d[d$group == pr[1], ]; b <- d[d$group == pr[2], ]; if (!nrow(a) || !nrow(b)) next
    cs <- chi_fun(pr[1], pr[2]); if (any(is.na(cs))) next
    x0 <- max(a$xmax, b$xmax) + 0.30 * xr
    p <- p + annotate("segment", x = x0, xend = x0, y = a$y, yend = b$y, linewidth = 0.4) +
      annotate("text", x = x0 + 0.012 * xr, y = (a$y + b$y) / 2, label = sprintf("(χ²=%.2f, p=%s)", cs[1], formatC(cs[2], format = "f", digits = 3)), hjust = 0, size = 2.8, fontface = "italic", colour = "#333333") }
  p + scale_fill_manual(values = cols, name = NULL) + scale_x_continuous(limits = c(0, xr * 2.05), breaks = pretty(c(0, xr))) + labs(x = xlab, y = NULL) + theme_m
}
f9 <- build(cr, "mortality_pct_28day", function(d) sprintf("%d/%d  (%s%%)", d$deaths_28day, d$n_outcome_known, format(round(d$mortality_pct_28day, 1), nsmall = 0, trim = TRUE)), "Crude 28-day mortality (%)", FALSE, function(a, b) chi(cr, a, b))
ggsave(file.path(FIG, "Figure9_mortality_crude.tiff"), f9, width = 8.5, height = 8.5, dpi = 500, bg = "white", device = "tiff", compression = "lzw")
cw$n_outcome_known <- cw$n
f16 <- build(cw, "calibrated_pct", function(d) sprintf("%d/%d  (%s%%)", d$deaths_28day, d$n, format(round(d$calibrated_pct, 1), nsmall = 0, trim = TRUE)), "Calibrated 28-day mortality (%, 95% CI)", TRUE, function(a, b) chi(cw, a, b, "n", "deaths_28day"))
ggsave(file.path(FIG, "SuppFigure_mortality_crude_calibrated.tiff"), f16, width = 8.5, height = 8.5, dpi = 500, bg = "white", device = "tiff", compression = "lzw")
cat("done\n")
