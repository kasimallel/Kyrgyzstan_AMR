# 05_heatmap.R -- % resistant by pathogen x antibiotic class (template: Kyrgyzstan Heatmap by Class 271026_10h.docx)
source(file.path(Sys.getenv("KGZ_PKG"), "R", "config.R"))
suppressPackageStartupMessages({library(survey); library(ggplot2); library(writexl)})
options(survey.lonely.psu = "adjust")
res <- readRDS(file.path(OUT, "02_resistance.rds")); r4 <- readRDS(file.path(OUT, "04_results.rds"))
D <- r4$D; ast <- res$ast; tpl <- res$tpl; dagger <- res$dagger; cls <- res$cls
rows <- names(tpl)
prop_cal <- function(y, data) {
  d <- data; d$.y <- as.numeric(y); n <- nrow(d); k <- sum(d$.y == 1)
  des <- svydesign(ids = ~facility_id_clean, weights = ~calibrated_weight, data = d)
  r <- tryCatch(suppressWarnings(svyciprop(~I(.y == 1), des, method = "logit")), error = function(e) NULL)
  ci <- if (!is.null(r)) tryCatch(suppressWarnings(confint(r)), error = function(e) c(NA, NA)) else c(NA, NA)
  est <- if (!is.null(r)) 100 * as.numeric(r) else NA
  if (is.null(r) || any(is.na(ci)) || k %in% c(0, n)) { bt <- binom.test(k, n); ci <- bt$conf.int; if (is.na(est)) est <- 100 * k / n; exact <- TRUE } else exact <- FALSE
  c(est = est, lo = 100 * ci[1], hi = 100 * ci[2], exact = exact)
}
cells <- list()
for (rw in rows) {
  ids <- D$record_id[!is.na(D$hm_row) & D$hm_row == rw]; n_row <- length(ids)
  for (cl in cls) {
    ag <- tpl[[rw]][[cl]]; dag <- cl %in% dagger[[rw]]
    if (is.null(ag)) { cells[[length(cells) + 1]] <- tibble(row = rw, class = cl, agents = NA, status = "not_in_panel", n_row, n_tested = NA, n_R = NA, pct = NA, lo = NA, hi = NA, coverage = NA); next }
    a <- ast %>% filter(record_id %in% ids, agent %in% ag, !is.na(res_final))
    per <- a %>% group_by(record_id) %>% summarise(R = any(res_final == "R"), .groups = "drop")
    n_t <- nrow(per); n_r <- sum(per$R); cov <- ifelse(n_row > 0, n_t / n_row, 0)
    st <- if (dag) "dagger" else if (n_t == 0) "not_in_panel" else if (n_t < 10 || cov < 0.30) "suppressed" else "shown"
    est <- c(est = NA, lo = NA, hi = NA, exact = NA)
    if (n_t >= 2) est <- prop_cal(per$R, D[match(per$record_id, D$record_id), ])
    cells[[length(cells) + 1]] <- tibble(row = rw, class = cl, agents = paste(ag, collapse = ", "), status = st, n_row, n_tested = n_t, n_R = n_r, pct = unname(est["est"]), lo = unname(est["lo"]), hi = unname(est["hi"]), coverage = round(100 * cov), unweighted_pct = ifelse(n_t > 0, 100 * n_r / n_t, NA))
  }
}
H <- bind_rows(cells)
H$bin <- cut(H$pct, c(-Inf, 10, 25, 50, 75, Inf), right = FALSE, labels = c("0–<10%", "10–<25%", "25–<50%", "50–<75%", "≥75%"))
H$label <- case_when(H$status == "not_in_panel" ~ "–", H$status == "dagger" ~ "†", H$status == "suppressed" ~ sprintf("n<10 (%d/%d)", H$n_R, H$n_tested),
                     TRUE ~ sprintf("%.0f%% (%d/%d)", H$pct, H$n_R, H$n_tested))
H$label[H$status == "suppressed" & !is.na(H$n_tested) & H$n_tested >= 10] <- sprintf("cov<30%% (%d/%d)", H$n_R, H$n_tested)[H$status == "suppressed" & !is.na(H$n_tested) & H$n_tested >= 10]
write.csv(H, file.path(OUT, "05_heatmap_cells.csv"), row.names = FALSE)
write_xlsx(list(heatmap_cells = H), file.path(OUT, "05_heatmap_cells.xlsx"))
rown <- H %>% distinct(row, n_row)
cat("\nHeatmap row n:\n"); print(as.data.frame(rown))
cat("\nShown cells:\n"); print(as.data.frame(H %>% filter(status == "shown") %>% select(row, class, agents, n_tested, n_R, pct, lo, hi, coverage)), digits = 3)
cat("\nStatus counts:\n"); print(table(H$status))

# ---- figure ----------------------------------------------------------------------------------------------------------
pal <- c("0–<10%" = "#FEE5D9", "10–<25%" = "#FCAE91", "25–<50%" = "#FB6A4A", "50–<75%" = "#DE2D26", "≥75%" = "#A50F15")
H$fill <- ifelse(H$status == "shown", as.character(H$bin), ifelse(H$status == "suppressed", "Suppressed", ifelse(H$status == "dagger", "Not interpretable", "Not in panel")))
fillcols <- c(pal, Suppressed = "#F4F4F4", `Not interpretable` = "#CFCFCF", `Not in panel` = "#E7E6E6")
H$row_lab <- factor(paste0(H$row, " (n=", H$n_row, ")"), levels = rev(paste0(rows, " (n=", rown$n_row[match(rows, rown$row)], ")")))
H$class <- factor(H$class, levels = cls)
H$txt <- ifelse(H$status == "shown", sprintf("%.0f%%\n%d/%d", H$pct, H$n_R, H$n_tested), ifelse(H$status == "suppressed", sprintf("n<10\n%d/%d", H$n_R, H$n_tested), ifelse(H$status == "dagger", "†", "–")))
H$txtcol <- ifelse(H$fill %in% c("50–<75%", "≥75%", "25–<50%"), "white", "black")
g <- ggplot(H, aes(class, row_lab, fill = fill)) + geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = txt, colour = txtcol), size = 2.5, lineheight = 0.85, fontface = "plain", family = "Times") +
  scale_colour_identity() + scale_fill_manual(values = fillcols, breaks = names(fillcols), name = "Resistant (calibrated %)") +
  scale_x_discrete(position = "top", labels = function(x) stringr::str_wrap(x, 9)) + labs(x = NULL, y = NULL) + theme_minimal(base_family = "Times", base_size = 10) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 6, face = "bold", lineheight = 0.8), axis.text.y = element_text(hjust = 1, face = "italic", size = 9), legend.position = "bottom", legend.title = element_text(size = 8), legend.text = element_text(size = 8))
ggsave(file.path(OUT, "05_heatmap_by_class.tiff"), g, width = 15, height = 5.6, dpi = 500, bg = "white", compression = "lzw")
cat("saved figure\n")
