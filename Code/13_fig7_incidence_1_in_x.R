# 12_fig7_one_in_x.R -- Figure 7 redrawn with the x axis as "1 in x admissions"
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
source(file.path(Sys.getenv("KGZ_PKG"), "R", "config.R"))
r4 <- readRDS(file.path(OUT, "04_results.rds"))
lancet_blue <- "#2F5D8C"
theme_set(theme_minimal(base_family = "Helvetica", base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), plot.margin = margin(8, 14, 8, 8)))
m <- r4$merged %>% filter(grepl("^BSI due to|^All|Community|Hospital|^BSI with", indicator)) %>% filter(!is.na(one_in), n_cases > 0)
m$indicator <- factor(m$indicator, levels = rev(m$indicator))
m$lab <- sprintf("  %s (%s–%s)", formatC(m$one_in, format = "d", big.mark = ","), formatC(m$one_in_lo, format = "d", big.mark = ","), formatC(m$one_in_hi, format = "d", big.mark = ","))
f7 <- ggplot(m, aes(one_in, indicator)) +
  geom_pointrange(aes(xmin = one_in_lo, xmax = one_in_hi), colour = lancet_blue, size = 0.4) +
  geom_text(aes(x = one_in_hi, label = lab), hjust = 0, size = 3) +
  scale_x_continuous(trans = "log10", breaks = c(100, 300, 1000, 3000, 10000, 30000, 100000, 300000), labels = function(x) formatC(x, format = "d", big.mark = ","), expand = expansion(mult = c(0.05, 0.45))) +
  labs(x = "Calibrated incidence as of 1 in 'x' number of admissions (95% CI)", y = NULL) +
  theme(panel.grid.major.x = element_line(colour = "#E0E0E0", linewidth = 0.3), axis.line.x = element_line(colour = "black", linewidth = 0.4))
ggsave(file.path(OUT, "figures_TIFF_500dpi", "Fig7_incidence_by_category.tiff"), f7, width = 10, height = 7.2, dpi = 500, bg = "white", device = "tiff", compression = "lzw")
print(as.data.frame(m[, c("indicator", "one_in", "one_in_lo", "one_in_hi")]))
