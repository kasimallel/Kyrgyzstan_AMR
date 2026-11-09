# 11_sexage_figures.R -- readable re-draw of the sex/age-stratified crude mortality and LOS figures (figures only)
source(file.path(Sys.getenv("KGZ_PKG"), "R", "config.R"))
suppressPackageStartupMessages(library(ggplot2))
d <- readxl::read_excel(file.path(V2, "legacy/out/SexAge_stratified_mortality_LOS_all_mechanisms.xlsx"))
ord <- c("Uninfected","Susceptible (any)","Resistant (any)","Non-MDR","MDR","Non-XDR","XDR","Carbapenem-S","Carbapenem-R","Cephalosporin-S","Cephalosporin-R","3GC-S","3GCR","MSSA","MRSA","VSE","VRE","Access-S","Access-R","Watch-S","Watch-R","Reserve-S","Reserve-R")
d$group <- as.character(d$group); d$ypos <- match(d$group, rev(ord)); d$yc <- d$ypos + ifelse(d$sex == "Female", 0.21, -0.21)
cols <- c("Male, adult" = "#2C5F7C", "Male, paediatric" = "#A1BDD0", "Female, adult" = "#B0583C", "Female, paediatric" = "#E8BBA8")
yb <- setNames(seq_along(rev(ord)), rev(ord))
mk <- function(d, val, lo = NULL, hi = NULL, lab, xlab, xmax, file, ttl_note) {
  d <- d[!is.na(d[[val]]), ]; d$v <- d[[val]]; d$ad <- d$v * ifelse(is.na(d$adult_pct), 0, d$adult_pct) / 100
  seg <- rbind(data.frame(d[c("yc","sex")], x0 = 0, x1 = d$ad, f = paste0(d$sex, ", adult")), data.frame(d[c("yc","sex")], x0 = d$ad, x1 = d$v, f = paste0(d$sex, ", paediatric")))
  seg$f <- factor(seg$f, levels = names(cols))
  p <- ggplot() + geom_rect(data = seg, aes(xmin = x0, xmax = x1, ymin = yc - 0.19, ymax = yc + 0.19, fill = f), colour = "black", linewidth = 0.25)
  if (!is.null(lo)) { d$lo <- pmax(d[[lo]], 0); d$hi <- d[[hi]]; d$clip <- !is.na(d$hi) & d$hi > xmax; d$hi_p <- pmin(d$hi, xmax)
    p <- p + geom_segment(data = d[!is.na(d$hi), ], aes(x = lo, xend = hi_p, y = yc, yend = yc), linewidth = 0.35, colour = "grey15") +
      geom_segment(data = d[!is.na(d$hi) & !d$clip, ], aes(x = hi_p, xend = hi_p, y = yc - 0.1, yend = yc + 0.1), linewidth = 0.35, colour = "grey15") +
      geom_text(data = d[d$clip, ], aes(x = hi_p, y = yc, label = "›"), size = 4, colour = "grey15")
    d$xl <- ifelse(is.na(d$hi), d$v, pmin(d$hi, xmax)) } else d$xl <- d$v
  p <- p + geom_text(data = d, aes(x = xl + xmax * 0.012, y = yc, label = lab(d)), hjust = 0, size = 2.9, colour = "black") +
    scale_y_continuous(breaks = unname(yb), labels = names(yb), expand = expansion(add = 0.6)) +
    scale_x_continuous(limits = c(0, xmax * 1.28), breaks = pretty(c(0, xmax)), expand = expansion(mult = c(0, 0))) +
    scale_fill_manual(values = cols, name = NULL, drop = FALSE) + coord_cartesian(clip = "off") +
    labs(x = xlab, y = NULL, caption = ttl_note) +
    theme_minimal(base_family = "Helvetica", base_size = 11) +
    theme(legend.position = "top", axis.text.y = element_text(colour = "black", size = 10, margin = margin(r = 8)), axis.text.x = element_text(colour = "black", size = 9.5),
          panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = "#E4E4E4", linewidth = 0.3),
          axis.line.x = element_line(colour = "black", linewidth = 0.4), plot.caption = element_text(hjust = 0, size = 8.5, lineheight = 1.15, margin = margin(t = 10)), plot.caption.position = "plot",
          plot.margin = margin(8, 20, 8, 8), legend.text = element_text(size = 10))
  for (dir in c(file.path(V2, "outputs/figures_TIFF_500dpi"), file.path(V2, "legacy/figs")))
    ggsave(file.path(dir, paste0(file, ".tiff")), p, width = 10.5, height = 11.5, dpi = 500, bg = "white", device = "tiff", compression = "lzw")
}
note1 <- "Bar = crude (unadjusted) value for each sex; darker shade = adult (≥18 y) share of that bar, lighter shade = paediatric (neonates and children) share. Labels: deaths/patients with known 28-day outcome."
mk(d, "mortality_pct", lab = function(x) sprintf("%.1f%% (%d/%d)", x$mortality_pct, x$deaths_28day, x$n_outcome_known), xlab = "Crude 28-day mortality (%)", xmax = max(d$mortality_pct, na.rm = TRUE) * 1.02, file = "SexAge_mortality_by_mechanism",
   ttl_note = paste(strwrap(note1, 150), collapse = "\n"))
d$los_lo <- d$los_mean_days - 1.96 * d$los_se_days; d$los_hi <- d$los_mean_days + 1.96 * d$los_se_days
note2 <- "Bar = crude mean length of stay from enrolment (days) for each sex, with 95% CI; darker shade = adult (≥18 y) share of that bar, lighter shade = paediatric share. › = interval extends beyond the axis (small samples). Labels: mean (n)."
mk(d, "los_mean_days", "los_lo", "los_hi", lab = function(x) sprintf("%.1f (n=%d)", x$los_mean_days, x$n_los_known), xlab = "Crude mean length of stay (days)", xmax = 14, file = "SexAge_LOS_by_mechanism", ttl_note = paste(strwrap(note2, 150), collapse = "\n"))
dp <- readxl::read_excel(file.path(V2, "legacy/out/SexAge_stratified_prevalence_all_mechanisms.xlsx"))
ordp <- c("Any resistance","MDR","XDR","Carbapenem-R","Cephalosporin-R","3GCR","MRSA","Access-R","Watch-R","Reserve-R")
dp$group <- as.character(dp$group)
ord <- ordp; yb <- setNames(seq_along(rev(ord)), rev(ord)); dp$yc <- yb[dp$group] + ifelse(dp$sex == "Female", 0.21, -0.21)
note3 <- "Bar = crude prevalence of each resistance category among evaluable isolates of that sex (resistant/susceptible isolates with a result); darker shade = adult (≥18 y) share of that bar, lighter shade = paediatric (neonates and children) share. Labels: resistant/evaluable isolates."
mk2 <- function() {
  d <- dp[!is.na(dp$prevalence_pct), ]; d$v <- d$prevalence_pct; d$ad <- d$v * ifelse(is.na(d$adult_pct), 0, d$adult_pct) / 100; xmax <- 100
  seg <- rbind(data.frame(d[c("yc","sex")], x0 = 0, x1 = d$ad, f = paste0(d$sex, ", adult")), data.frame(d[c("yc","sex")], x0 = d$ad, x1 = d$v, f = paste0(d$sex, ", paediatric"))); seg$f <- factor(seg$f, levels = names(cols))
  p <- ggplot() + geom_rect(data = seg, aes(xmin = x0, xmax = x1, ymin = yc - 0.19, ymax = yc + 0.19, fill = f), colour = "black", linewidth = 0.25) +
    geom_text(data = d, aes(x = v + 1.2, y = yc, label = sprintf("%.1f%% (%d/%d)", prevalence_pct, n_resistant, n_tested)), hjust = 0, size = 3, colour = "black") +
    scale_y_continuous(breaks = unname(yb), labels = names(yb), expand = expansion(add = 0.6)) + scale_x_continuous(limits = c(0, 135), breaks = seq(0, 100, 25), expand = expansion(mult = c(0, 0))) +
    scale_fill_manual(values = cols, name = NULL, drop = FALSE) + coord_cartesian(clip = "off") + labs(x = "Crude prevalence (%)", y = NULL, caption = paste(strwrap(note3, 130), collapse = "\n")) +
    theme_minimal(base_family = "Helvetica", base_size = 11) +
    theme(legend.position = "top", axis.text.y = element_text(colour = "black", size = 10.5, margin = margin(r = 8)), axis.text.x = element_text(colour = "black", size = 9.5), panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = "#E4E4E4", linewidth = 0.3), axis.line.x = element_line(colour = "black", linewidth = 0.4), plot.caption = element_text(hjust = 0, size = 8.5, lineheight = 1.15, margin = margin(t = 10)), plot.caption.position = "plot",
          plot.margin = margin(8, 20, 8, 8), legend.text = element_text(size = 10))
  for (dir in c(file.path(V2, "outputs/figures_TIFF_500dpi"), file.path(V2, "legacy/figs"))) ggsave(file.path(dir, "SexAge_prevalence_by_mechanism.tiff"), p, width = 10.5, height = 7.5, dpi = 500, bg = "white", device = "tiff", compression = "lzw")
}
mk2()
cat("done\n")
