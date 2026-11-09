# 07_extra_outputs.R -- subgroup resistance among isolates, oblast table, figures (pathogen distribution, phenotypes, sex/age, incidence, maps)
source(file.path(Sys.getenv("KGZ_PKG"), "R", "config.R"))
suppressPackageStartupMessages({library(survey); library(ggplot2); library(writexl); library(sf); library(patchwork)})
options(survey.lonely.psu = "adjust")
r4 <- readRDS(file.path(OUT, "04_results.rds")); co <- readRDS(file.path(OUT, "03_cohort.rds"))
D <- r4$D; fac <- co$fac
src <- readLines(file.path(V2, "R/04_prevalence_incidence.R"))
# reuse helper functions from script 04 (prop_ci, rate_ci, fmt, pct_txt)
eval(parse(text = src[grep("^prop_ci <- function", src):(grep("^one_in <-", src))]))
ast_ok <- TRUE
fmtp <- function(o) sprintf("%s (%s–%s)", fmt(o$est), fmt(o$lo), fmt(o$hi))

# ---- subgroup resistance among isolates ------------------------------------------------------------------------
D$age_grp2 <- D$age_grp; D$ward2 <- D$ward_type
sg <- function(var, lev) {
  m <- !is.na(D[[var]]) & D[[var]] == lev
  rr <- function(ind, den) { o <- prop_ci(ind, m & den)
    if (var == "NameAU" && o$n_den >= 1) { bt <- binom.test(o$n_num, o$n_den); o$est <- 100 * o$n_num / o$n_den; o$lo <- 100 * bt$conf.int[1]; o$hi <- 100 * bt$conf.int[2] }   # oblast strata: 1-7 facilities, use exact unweighted CI
    o }
  a <- rr(D$R_any, !is.na(D$R_any)); b <- rr(D$mdr == "Yes", D$mdr %in% c("Yes","No")); c3 <- rr(D$ESBL_3GCR, !is.na(D$ESBL_3GCR)); cr <- rr(D$CR, !is.na(D$CR)); mr <- rr(D$MRSA, !is.na(D$MRSA))
  tibble(variable = var, level = lev, n_isolates = sum(m), share_of_isolates = fmtp(prop_ci(m)), any_R = fmtp(a), n_anyR = sprintf("%d/%d", a$n_num, a$n_den), MDR = fmtp(b), n_MDR = sprintf("%d/%d", b$n_num, b$n_den),
         `3GC-R` = fmtp(c3), `CR-GN` = fmtp(cr), MRSA = fmtp(mr))
}
subg <- bind_rows(lapply(list(c("infection_origin", "Community-onset"), c("infection_origin", "Hospital-onset"), c("age_grp", "Neonate"), c("age_grp", "Paediatric"), c("age_grp", "Adult"),
  c("sex", "Female"), c("sex", "Male"), c("ward_type", "ICU"), c("ward_type", "Non-ICU"), c("specialty", "Medical"), c("specialty", "Surgical"),
  c("NameAU", "Batken"), c("NameAU", "Bishkek"), c("NameAU", "Chui"), c("NameAU", "Jalal-Abad"), c("NameAU", "Naryn"), c("NameAU", "Osh"), c("NameAU", "Talas")), function(x) sg(x[1], x[2])))
subg$suppress_note <- ifelse(subg$n_isolates < 10, "n<10: imprecise", "")
print(as.data.frame(subg %>% select(variable, level, n_isolates, share_of_isolates, any_R, MDR, suppress_note)))

# ---- oblast incidence (pooled, exact Poisson) and resistance % -------------------------------------------------------------
cases_by <- D %>% group_by(facility_id_clean) %>% summarise(n_bsi = n(), n_anyR = sum(R_any %in% TRUE), n_mdr = sum(mdr %in% "Yes"), .groups = "drop")
ob <- fac %>% left_join(cases_by, by = c(fid = "facility_id_clean")) %>% mutate(across(c(n_bsi, n_anyR, n_mdr), ~ ifelse(is.na(.x), 0, .x))) %>%
  group_by(NameAU) %>% summarise(facilities = n(), admissions = sum(admissions), n_bsi = sum(n_bsi), n_anyR = sum(n_anyR), n_mdr = sum(n_mdr), .groups = "drop")
pc <- function(x, e) { lo <- ifelse(x == 0, 0, qchisq(0.025, 2 * x) / 2); hi <- qchisq(0.975, 2 * x + 2) / 2; cbind(est = 1000 * x / e, lo = 1000 * lo / e, hi = 1000 * hi / e) }
ob <- bind_cols(ob, as_tibble(pc(ob$n_bsi, ob$admissions)) %>% setNames(c("bsi_per1000", "bsi_lo", "bsi_hi")),
                as_tibble(pc(ob$n_anyR, ob$admissions)) %>% setNames(c("r_per1000", "r_lo", "r_hi")), as_tibble(pc(ob$n_mdr, ob$admissions)) %>% setNames(c("mdr_per1000", "mdr_lo", "mdr_hi")))
ob$one_in <- round(1000 / ob$bsi_per1000)
ob$isolates_pct_anyR <- ifelse(ob$n_bsi >= 10, round(100 * ob$n_anyR / ob$n_bsi, 1), NA)
print(as.data.frame(ob %>% select(NameAU, facilities, admissions, n_bsi, bsi_per1000, one_in, n_anyR, n_mdr)), digits = 3)

# ---- figures (TIFF, 500 dpi, LZW) -----------------------------------------------------------------------------------------------
sv <- function(name, plot, w, h) ggsave(file.path(OUT, paste0(name, ".tiff")), plot, width = w, height = h, dpi = 500, bg = "white", device = "tiff", compression = "lzw")
lancet <- c(blue = "#2F5D8C", red = "#B0413E", green = "#5E8C61", teal = "#5B8FA8", purple = "#7A5C8E", salmon = "#D9A38F", crimson = "#8C2F39", grey = "#A9B7C0", black = "#1B1919")   # muted Lancet-style palette
theme_set(theme_minimal(base_family = "Helvetica", base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), plot.margin = margin(8, 14, 8, 8)))

# Figure 2: pathogen distribution (counts and % of all isolates), coloured by organism group, no confidence intervals
dist <- r4$dist %>% filter(!is.na(pathogen)); N <- sum(dist$n)
grp_of <- function(x) case_when(grepl("^(Escherichia|Klebsiella|Enterobacter|Citrobacter|Proteus|Salmonella|Pantoea)", x) ~ "Enterobacterales",
  grepl("^(Acinetobacter|Pseudomonas|Stenotrophomonas|Moraxella)", x) ~ "Non-fermenting Gram-negatives",
  grepl("^Staphylococcus", x) ~ "Staphylococcus", grepl("^(Streptococcus|Viridans)", x) ~ "Streptococcus", grepl("^Enterococcus", x) ~ "Enterococcus",
  grepl("^(Sphingomonas|Aeromonas|Agrobacterium|Chromobacterium|Rodentibacter)", x) ~ "Other Gram-negatives", TRUE ~ "Other Gram-positives")
dist$grp <- grp_of(dist$pathogen)
big <- dist %>% filter(n >= 3); oth <- dist %>% filter(n < 3)
d2 <- bind_rows(big %>% select(pathogen, n, grp), tibble(pathogen = sprintf("Other species (%d species, <3 isolates each)", nrow(oth)), n = sum(oth$n), grp = "Other"))
d2$pct <- 100 * d2$n / N
d2$pathogen <- factor(d2$pathogen, levels = rev(c(d2$pathogen[d2$grp != "Other"][order(-d2$n[d2$grp != "Other"])], d2$pathogen[d2$grp == "Other"])))
gcol <- c(Enterobacterales = lancet[["blue"]], `Non-fermenting Gram-negatives` = lancet[["teal"]], `Other Gram-negatives` = lancet[["grey"]], Staphylococcus = lancet[["red"]], Streptococcus = lancet[["green"]], Enterococcus = lancet[["purple"]], `Other Gram-positives` = lancet[["salmon"]], Other = "#CFCFCF")
d2$grp <- factor(d2$grp, levels = names(gcol))
f2 <- ggplot(d2, aes(pct, pathogen, fill = grp)) + geom_col(width = 0.72) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", n, pct)), hjust = -0.12, size = 3.1, family = "Helvetica") +
  scale_fill_manual(values = gcol, breaks = setdiff(names(gcol), "Other"), name = NULL) + scale_x_continuous(expand = expansion(mult = c(0, 0.16)), breaks = seq(0, 30, 5)) +
  labs(x = sprintf("Percentage of culture-confirmed BSI isolates (n=%d)", N), y = NULL) +
  theme(axis.text.y = element_text(face = "italic", colour = "black"), axis.line.x = element_line(colour = "black", linewidth = 0.4), panel.grid.major.x = element_line(colour = "#E6E6E6", linewidth = 0.3),
        legend.position = "bottom", legend.text = element_text(size = 9), axis.title.x = element_text(margin = margin(t = 8)))
sv("Fig2_pathogen_distribution", f2, 9, 6.6)

# Figure 4: resistance phenotypes among isolates (calibrated)
l1 <- r4$level1 %>% filter(indicator %in% c(grep("Any resistance", r4$level1$indicator, value = TRUE), "Resistant to ≥1 Access agent", "Resistant to ≥1 Watch agent", "Resistant to ≥1 Reserve agent", grep("^MDR", r4$level1$indicator, value = TRUE), "XDR", "PDR",
        "3rd-gen cephalosporin-resistant Enterobacterales", grep("^Carbapenem", r4$level1$indicator, value = TRUE), "MRSA among S. aureus", "Vancomycin-resistant enterococci"))
l1$indicator <- c("Any resistance", "Access-resistant", "Watch-resistant", "Reserve-resistant", "MDR", "XDR", "PDR", "3GC-resistant Enterobacterales", "Carbapenem-resistant Gram-negatives", "MRSA (of S. aureus)", "VRE (of enterococci)")
l1$lab <- sprintf("%s (%d/%d)", l1$indicator, l1$n_resistant, l1$n_denominator); l1$lab <- factor(l1$lab, levels = rev(l1$lab))
f4 <- ggplot(l1, aes(calibrated_pct, lab)) + geom_col(fill = lancet[["blue"]], width = 0.65) + geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, orientation = "y", linewidth = 0.4) +
  geom_text(aes(x = pmin(hi, 100) + 1.5, label = sprintf("%.1f%%", calibrated_pct)), hjust = 0, size = 3.1) + scale_x_continuous(limits = c(0, 115), breaks = seq(0, 100, 25)) +
  labs(x = "Calibrated % of isolates (95% CI); n resistant/n evaluable in label", y = NULL) + theme(axis.line.x = element_line(colour = "black", linewidth = 0.4), panel.grid.major.x = element_line(colour = "#E6E6E6", linewidth = 0.3))
sv("Fig4_resistance_phenotypes", f4, 9.5, 5.2)

# Figure 5: crude prevalence by sex (bars) and age (adult dark / paediatric light segments), phenotype names on the y axis
phen <- list(`Any resistance` = list(D$R_any, !is.na(D$R_any)), MDR = list(D$mdr == "Yes", !is.na(D$mdr)), XDR = list(D$xdr == "Yes", !is.na(D$xdr)), `Carbapenem-R` = list(D$CR, !is.na(D$CR)),
  `Cephalosporin-R` = list(D$CephR, !is.na(D$CephR)), `3GCR` = list(D$ESBL_3GCR, !is.na(D$ESBL_3GCR)), MRSA = list(D$MRSA, !is.na(D$MRSA)), `Access-R` = list(D$R_Access, !is.na(D$R_Access)), `Watch-R` = list(D$R_Watch, !is.na(D$R_Watch)), `Reserve-R` = list(D$R_Reserve, !is.na(D$R_Reserve)))
D$age2 <- ifelse(D$age_grp == "Adult", "adult", ifelse(is.na(D$age_grp), NA, "paediatric"))
rows5 <- bind_rows(lapply(names(phen), function(nm) bind_rows(lapply(c("Male", "Female"), function(sx) {
  ind <- phen[[nm]][[1]]; den <- phen[[nm]][[2]] & D$sex %in% sx & !is.na(D$age2); N_s <- sum(den); K_s <- sum(ind[den] %in% TRUE)
  tibble(phenotype = nm, sex = sx, N = N_s, K = K_s, seg_adult = 100 * sum(ind[den & D$age2 == "adult"] %in% TRUE) / N_s, seg_paed = 100 * sum(ind[den & D$age2 == "paediatric"] %in% TRUE) / N_s) }))))
l5 <- rows5 %>% tidyr::pivot_longer(c(seg_adult, seg_paed), names_to = "seg", values_to = "pct") %>% mutate(fillgrp = paste0(sex, ", ", ifelse(seg == "seg_adult", "adult", "paediatric")))
l5$fillgrp <- factor(l5$fillgrp, levels = c("Male, adult", "Male, paediatric", "Female, adult", "Female, paediatric"))
l5$phenotype <- factor(l5$phenotype, levels = rev(names(phen))); l5$sex <- factor(l5$sex, levels = c("Female", "Male"))
tot5 <- rows5 %>% mutate(tot = 100 * K / N, phenotype = factor(phenotype, levels = rev(names(phen))), sex = factor(sex, levels = c("Female", "Male")))
f5 <- ggplot(l5, aes(pct, phenotype, fill = fillgrp, group = interaction(phenotype, sex))) + geom_col(aes(y = phenotype), position = position_dodge2(width = 0.9, reverse = TRUE), width = 0.9, colour = "black", linewidth = 0.25) +
  scale_fill_manual(values = c("Male, adult" = "#2C5F7C", "Male, paediatric" = "#A1BDD0", "Female, adult" = "#B0583C", "Female, paediatric" = "#E8BBA8"), name = NULL) + labs(x = "Crude prevalence (%)", y = NULL)
# dodge2 with stacked segments is awkward -> draw explicit rectangles
yb <- setNames(seq_along(rev(names(phen))), rev(names(phen)))
rect <- l5 %>% arrange(phenotype, sex, fillgrp) %>% group_by(phenotype, sex) %>% mutate(xmax = cumsum(pct), xmin = xmax - pct) %>% ungroup() %>%
  mutate(yc = yb[as.character(phenotype)] + ifelse(sex == "Female", 0.22, -0.22))
f5 <- ggplot(rect) + geom_rect(aes(xmin = xmin, xmax = xmax, ymin = yc - 0.2, ymax = yc + 0.2, fill = fillgrp), colour = "black", linewidth = 0.25) +
  geom_text(data = tot5 %>% mutate(yc = yb[as.character(phenotype)] + ifelse(sex == "Female", 0.22, -0.22)), aes(x = tot + 1.2, y = yc, label = sprintf("%.1f%% (%d/%d)", tot, K, N)), hjust = 0, size = 2.9) +
  scale_y_continuous(breaks = unname(yb), labels = names(yb), expand = expansion(add = 0.5)) + scale_x_continuous(limits = c(0, 135), breaks = seq(0, 100, 25)) +
  scale_fill_manual(values = c("Male, adult" = "#2C5F7C", "Male, paediatric" = "#A1BDD0", "Female, adult" = "#B0583C", "Female, paediatric" = "#E8BBA8"), name = NULL) +
  labs(x = "Crude prevalence (%)", y = NULL, caption = "Darker shade = adult proportion of that bar; lighter shade = paediatric (neonates and children).") +
  theme(legend.position = "top", axis.text.y = element_text(colour = "black", size = 10), panel.grid.major.x = element_line(colour = "#E0E0E0", linewidth = 0.3), axis.line.x = element_line(colour = "black", linewidth = 0.4), plot.caption = element_text(hjust = 1, size = 8.5))
sv("Fig5_sex_age_phenotypes", f5, 10, 8)

# Figure 7: incidence by category
m <- r4$merged %>% filter(grepl("^BSI due to|^All|Community|Hospital|^BSI with", indicator)) %>% filter(!is.na(per1000adm), n_cases > 0)
m$indicator <- factor(m$indicator, levels = rev(m$indicator))
f7 <- ggplot(m, aes(per1000adm, indicator)) + geom_pointrange(aes(xmin = adm_lo, xmax = adm_hi), colour = lancet[["blue"]], size = 0.4) +
  geom_text(aes(x = adm_hi, label = sprintf("  1 in %s", formatC(one_in, format = "d", big.mark = ","))), hjust = 0, size = 3) +
  scale_x_continuous(trans = "log10", expand = expansion(mult = c(0.05, 0.35))) + labs(x = "Calibrated incidence per 1,000 admissions (95% CI, log scale); label = 1 in x admissions", y = NULL)
sv("Fig7_incidence_by_category", f7, 10, 7.2)

# Figures 6 and 8: oblast maps, indicator printed under each oblast name
sh <- tryCatch(rnaturalearth::ne_states(country = "Kyrgyzstan", returnclass = "sf"), error = function(e) NULL)
if (!is.null(sh)) {
  sh$NameAU <- recode(sh$name_en, Chuy = "Chui"); mp <- left_join(sh, ob, by = "NameAU")
  mp$pctR <- ob$isolates_pct_anyR[match(mp$NameAU, ob$NameAU)]
  mp$pctM <- 100 * ob$n_mdr[match(mp$NameAU, ob$NameAU)] / pmax(ob$n_bsi[match(mp$NameAU, ob$NameAU)], 1); mp$pctM[ob$n_bsi[match(mp$NameAU, ob$NameAU)] < 10] <- NA
  suppressWarnings({ pts <- sf::st_point_on_surface(sf::st_zm(mp)) }); co <- sf::st_coordinates(pts); mp$lx <- co[, 1]; mp$ly <- co[, 2]
  mp$lx[mp$NameAU == "Bishkek"] <- 73.9; mp$ly[mp$NameAU == "Bishkek"] <- 43.3; mp$lx[mp$NameAU == "Chui"] <- 75.2; mp$ly[mp$NameAU == "Chui"] <- 42.55
  mk_map <- function(var, ttl, fmt_fun) { d <- mp; d$lab <- paste0(d$NameAU, "\n", fmt_fun(d[[var]]))
    ggplot(d) + geom_sf(aes(fill = .data[[var]]), colour = "white") + scale_fill_gradient(low = "#FEE5D9", high = "#A50F15", name = ttl, na.value = "grey85") +
      geom_text(aes(x = lx, y = ly, label = lab), size = 2.6, family = "Helvetica", lineheight = 0.9) + theme_void(base_family = "Helvetica") }
  na_txt <- function(f) function(v) ifelse(is.na(v), "n<10", f(v))
  A <- mk_map("bsi_per1000", "BSI per\n1,000 adm.", function(v) sprintf("%.2f", v)); B <- mk_map("r_per1000", "Resistant BSI\nper 1,000 adm.", function(v) sprintf("%.2f", v)); C <- mk_map("mdr_per1000", "MDR BSI per\n1,000 adm.", function(v) sprintf("%.2f", v))
  sv("Fig8_incidence_maps", A + B + C + plot_annotation(tag_levels = "A"), 14, 3.6)
  A <- mk_map("n_bsi", "Isolates (n)", function(v) sprintf("n=%d", v)); B <- mk_map("pctR", "% resistant\n(oblasts ≥10 isolates)", na_txt(function(v) sprintf("%.0f%%", v))); C <- mk_map("pctM", "% MDR\n(oblasts ≥10 isolates)", na_txt(function(v) sprintf("%.0f%%", v)))
  sv("Fig6_prevalence_maps", A + B + C + plot_annotation(tag_levels = "A"), 14, 3.6)
}
write_xlsx(list(subgroups = subg, oblast = ob), file.path(OUT, "07_subgroups_oblast.xlsx"))
saveRDS(list(subg = subg, ob = ob), file.path(OUT, "07_extra.rds"))
cat("done\n")
