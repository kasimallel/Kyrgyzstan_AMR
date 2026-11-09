# 04_prevalence_incidence.R -- calibrated prevalence among isolates (3 levels), calibrated incidence ("1 in x"), merged table
source(file.path(Sys.getenv("KGZ_PKG"), "R", "config.R"))
suppressPackageStartupMessages({library(survey); library(writexl)})
options(survey.lonely.psu = "adjust")
res <- readRDS(file.path(OUT, "02_resistance.rds")); co <- readRDS(file.path(OUT, "03_cohort.rds")); p1 <- readRDS(file.path(OUT, "01_pathogens_ast.rds"))
D <- co$link %>% left_join(res$iso %>% select(record_id, R_any, R_Access, R_Watch, R_Reserve, mdr, xdr, pdr, ESBL_3GCR, CR, CephR, MRSA, VR, mag_table, hm_row, n_agents), by = "record_id", suffix = c("", ".y"))
D <- D[, !grepl("\\.y$", names(D))]
ast <- res$ast; fac <- co$fac; sf <- co$sf
NAT_ADM <- sf$TotalInpatPop[1]
stopifnot(all(!is.na(D$calibrated_weight)))

# ---- helpers ---------------------------------------------------------------------------------------
prop_ci <- function(ind, mask = rep(TRUE, nrow(D)), data = D) {
  d <- data; d$.y <- as.numeric(ind); d$.m <- mask & !is.na(ind)
  n_den <- sum(d$.m); n_num <- sum(d$.y[d$.m] == 1)
  out <- list(n_num = n_num, n_den = n_den, unw = ifelse(n_den > 0, 100 * n_num / n_den, NA), est = NA_real_, lo = NA_real_, hi = NA_real_)
  if (n_den < 2) return(out)
  des <- subset(svydesign(ids = ~facility_id_clean, weights = ~calibrated_weight, data = d), .m)
  r <- tryCatch(suppressWarnings(svyciprop(~I(.y == 1), des, method = "logit")), error = function(e) NULL)
  ci <- if (!is.null(r)) tryCatch(suppressWarnings(confint(r)), error = function(e) c(NA, NA)) else c(NA, NA)
  if (is.null(r) || any(is.na(ci)) || out$n_num %in% c(0, n_den)) {
    r <- tryCatch(suppressWarnings(svyciprop(~I(.y == 1), des, method = "beta")), error = function(e) NULL)
    ci <- if (!is.null(r)) tryCatch(suppressWarnings(confint(r)), error = function(e) c(NA, NA)) else c(NA, NA) }
  if (!is.null(r)) { out$est <- 100 * as.numeric(r); out$lo <- 100 * ci[1]; out$hi <- 100 * ci[2] }
  if (!is.na(out$est) && (is.na(out$lo) || is.na(out$hi)) && n_den >= 2) {     # degenerate design variance at 0%/100%: exact (Clopper-Pearson) interval on counts
    bt <- binom.test(n_num, n_den); out$lo <- 100 * bt$conf.int[1]; out$hi <- 100 * bt$conf.int[2]; out$exact <- TRUE }
  if (is.na(out$est) && n_den >= 1) { out$est <- out$unw; bt <- binom.test(n_num, n_den); out$lo <- 100 * bt$conf.int[1]; out$hi <- 100 * bt$conf.int[2]; out$exact <- TRUE }
  out
}
fmt <- function(x, d = 1) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = d))
pct_txt <- function(o) if (o$n_den == 0) "—" else sprintf("%s (%s–%s)", fmt(o$est), fmt(o$lo), fmt(o$hi))

fac_base <- fac %>% transmute(fid, NameAU, admissions, patient_days, calibrated_weight)
rate_ci <- function(case_vec, exposure = "admissions", per = 1000) {   # case_vec: per-patient 0/1 aligned with D
  cnt <- tibble(fid = D$facility_id_clean, y = case_vec) %>% group_by(fid) %>% summarise(n_cases = sum(y, na.rm = TRUE), .groups = "drop")
  f <- fac_base %>% left_join(cnt, by = "fid") %>% mutate(n_cases = ifelse(is.na(n_cases), 0, n_cases), expo = .data[[exposure]], lexpo = log(expo))
  out <- tibble(n_cases = sum(f$n_cases), est = NA_real_, lo = NA_real_, hi = NA_real_)
  if (sum(f$n_cases) == 0) return(out)
  des <- svydesign(ids = ~fid, weights = ~calibrated_weight, data = f)
  m <- tryCatch(suppressWarnings(svyglm(n_cases ~ 1, design = des, family = quasipoisson(), offset = lexpo)), error = function(e) NULL)
  if (is.null(m)) return(out)
  ci <- tryCatch(exp(confint(m)[1, ]) * per, error = function(e) c(NA_real_, NA_real_))
  out$est <- unname(exp(coef(m)[1]) * per); out$lo <- ci[1]; out$hi <- ci[2]; out
}
one_in <- function(r) ifelse(is.na(r) | r <= 0, NA, round(1000 / r))

# ---- CONSORT ----------------------------------------------------------------------------------------------
cons <- p1$consort$by_label
pr <- p1$pathogen_rows
consort <- tibble(step = c("Central reference laboratory records", "  Contaminant", "  Inconclusive", "  Negative", "  Blank confirmation field",
                           "Confirmed true pathogen ('pathogen')", "  Removed: fungal (Candida, Cryptococcus, Malassezia)",
                           "Culture-confirmed bacterial BSI isolates analysed",
                           "  Species taken from Species 1", "  Species taken from Species 2 (fallback)", "  Species taken from the commensal-identification field", "  Species taken from the laboratory comment",
                           "  Polymicrobial (second, different pathogen in other sample; not analysed as separate isolate)"),
                  n = c(p1$consort$crl_rows, cons[["contaminant"]], cons[["inconclusive"]], cons[["negative"]], sum(cons[is.na(names(cons))]),
                        nrow(pr), sum(pr$excl_reason == "Fungal (removed)", na.rm = TRUE), sum(!is.na(pr$species)),
                        sum(pr$name_source == "Species 1", na.rm = TRUE), sum(pr$name_source == "Species 2 (fallback)", na.rm = TRUE),
                        sum(grepl("^Commensal", pr$name_source), na.rm = TRUE), sum(pr$name_source == "Laboratory comment", na.rm = TRUE), sum(pr$polymicrobial, na.rm = TRUE)))
print(as.data.frame(consort))

# ---- LEVEL 0: pathogen distribution -------------------------------------------------------------------------
grp_counts <- D %>% count(pathogen_group, sort = TRUE)
dist <- bind_rows(lapply(grp_counts$pathogen_group, function(g) { o <- prop_ci(D$pathogen_group == g); tibble(pathogen = g, n = o$n_num, pct_of_isolates = o$unw, cal_pct = o$est, lo = o$lo, hi = o$hi) }))

# ---- LEVEL 1: pathogen-agnostic resistance ---------------------------------------------------------------------
agn <- list(
  list("Any resistance (I or R to ≥1 interpretable agent)", D$R_any, !is.na(D$R_any)),
  list("Resistant to ≥1 Access agent", D$R_Access, !is.na(D$R_Access)),
  list("Resistant to ≥1 Watch agent", D$R_Watch, !is.na(D$R_Watch)),
  list("Resistant to ≥1 Reserve agent", D$R_Reserve, !is.na(D$R_Reserve)),
  list("MDR (Magiorakos) — organisms with a Table 1–5 definition", D$mdr == "Yes", D$mdr %in% c("Yes","No")),
  list("XDR", D$xdr == "Yes", D$xdr %in% c("Yes","No")),
  list("PDR", D$pdr == "Yes", D$pdr %in% c("Yes","No")),
  list("3rd-gen cephalosporin-resistant Enterobacterales", D$ESBL_3GCR, !is.na(D$ESBL_3GCR)),
  list("Carbapenem-resistant Gram-negatives", D$CR, !is.na(D$CR)),
  list("MRSA among S. aureus", D$MRSA, !is.na(D$MRSA)),
  list("Vancomycin-resistant enterococci", D$VR, !is.na(D$VR)))
level1 <- bind_rows(lapply(agn, function(a) { o <- prop_ci(a[[2]], a[[3]]); tibble(indicator = a[[1]], n_resistant = o$n_num, n_denominator = o$n_den, unweighted_pct = o$unw, calibrated_pct = o$est, lo = o$lo, hi = o$hi, ci_text = pct_txt(o)) }))

# ---- LEVEL 2: by pathogen -------------------------------------------------------------------------------------------
lvl2_groups <- grp_counts$pathogen_group[grp_counts$n >= 10]
level2 <- bind_rows(lapply(lvl2_groups, function(g) {
  m <- D$pathogen_group == g; mk <- function(ind, den = !is.na(ind)) prop_ci(ind, m & den)
  a <- mk(D$R_any); mdr <- mk(D$mdr == "Yes", D$mdr %in% c("Yes","No")); xdr <- mk(D$xdr == "Yes", D$xdr %in% c("Yes","No"))
  key <- switch(g, `Escherichia coli` = list("3GC-R", D$ESBL_3GCR), `Klebsiella pneumoniae` = list("3GC-R", D$ESBL_3GCR), `Acinetobacter baumannii` = list("Carbapenem-R", D$CR),
                `Staphylococcus aureus` = list("MRSA", D$MRSA), `Enterococcus faecium` = list("VRE", D$VR), `Enterococcus faecalis` = list("VRE", D$VR), NULL)
  k <- if (!is.null(key)) mk(key[[2]]) else NULL
  tibble(pathogen = g, n_isolates = sum(m), any_R = pct_txt(a), n_any_R = sprintf("%d/%d", a$n_num, a$n_den), MDR = pct_txt(mdr), n_MDR = sprintf("%d/%d", mdr$n_num, mdr$n_den),
         XDR = pct_txt(xdr), key_phenotype = if (!is.null(k)) key[[1]] else NA, key_pct = if (!is.null(k)) pct_txt(k) else NA, key_n = if (!is.null(k)) sprintf("%d/%d", k$n_num, k$n_den) else NA)
}))

# ---- LEVEL 3: drug-bug (organisms with >=10 isolates; agents with >=10 interpretable results) ------------------------
agents_lab <- setNames(res$drug$drug, res$drug$code)
db <- bind_rows(lapply(lvl2_groups, function(g) {
  ids <- D$record_id[D$pathogen_group == g]
  a <- ast %>% filter(record_id %in% ids, !is.na(res_final))
  bind_rows(lapply(unique(a$agent), function(ag) {
    r <- D %>% select(record_id) %>% left_join(a %>% filter(agent == ag) %>% select(record_id, res_final), by = "record_id")
    o <- prop_ci(r$res_final == "R", D$pathogen_group == g & !is.na(r$res_final))
    ni <- sum(a$agent == ag & a$res_final == "I")
    tibble(pathogen = g, agent = ag, antibiotic = agents_lab[[ag]], n_tested = o$n_den, n_R = o$n_num, n_I = ni, unweighted_pct = o$unw, calibrated_pct = o$est, lo = o$lo, hi = o$hi, ci_text = pct_txt(o), suppressed = o$n_den < 10)
  }))
})) %>% arrange(pathogen, desc(n_tested))

# ---- Incidence + merged table -----------------------------------------------------------------------------------------------
inc_ind <- function(y) { a <- rate_ci(y, "admissions"); d <- rate_ci(y, "patient_days"); tibble(n_cases = a$n_cases, per1000adm = a$est, adm_lo = a$lo, adm_hi = a$hi, one_in = one_in(a$est), one_in_lo = one_in(a$hi), one_in_hi = one_in(a$lo),
  national_cases = a$est / 1000 * NAT_ADM, per1000pd = d$est, pd_lo = d$lo, pd_hi = d$hi) }
rows <- list(
  list("All culture-confirmed bacterial BSI", rep(1, nrow(D)), NULL),
  list("Community-onset BSI", as.integer(D$infection_origin == "Community-onset"), NULL),
  list("Hospital-onset BSI", as.integer(D$infection_origin == "Hospital-onset"), NULL))
for (g in lvl2_groups) rows[[length(rows) + 1]] <- list(paste0("BSI due to ", g), as.integer(D$pathogen_group == g), as.integer(D$pathogen_group == g))
ph <- list(list("Any resistance", D$R_any), list("Access-resistant", D$R_Access), list("Watch-resistant", D$R_Watch), list("Reserve-resistant", D$R_Reserve),
           list("MDR", D$mdr == "Yes"), list("XDR", D$xdr == "Yes"), list("3GC-resistant Enterobacterales", D$ESBL_3GCR), list("Carbapenem-resistant Gram-negative", D$CR),
           list("MRSA", D$MRSA), list("VRE", D$VR))
for (x in ph) rows[[length(rows) + 1]] <- list(paste0("BSI with ", x[[1]]), as.integer(x[[2]] %in% TRUE), NULL)
merged <- bind_rows(lapply(rows, function(r) {
  inc <- inc_ind(r[[2]])
  prev <- if (r[[1]] == "All culture-confirmed bacterial BSI") prop_ci(rep(TRUE, nrow(D))) else {
    lab <- r[[1]]
    ind <- switch(lab, `Community-onset BSI` = D$infection_origin == "Community-onset", `Hospital-onset BSI` = D$infection_origin == "Hospital-onset",
      `BSI with Any resistance` = D$R_any, `BSI with Access-resistant` = D$R_Access, `BSI with Watch-resistant` = D$R_Watch, `BSI with Reserve-resistant` = D$R_Reserve,
      `BSI with MDR` = D$mdr == "Yes", `BSI with XDR` = D$xdr == "Yes", `BSI with 3GC-resistant Enterobacterales` = D$ESBL_3GCR, `BSI with Carbapenem-resistant Gram-negative` = D$CR,
      `BSI with MRSA` = D$MRSA, `BSI with VRE` = D$VR, D$pathogen_group == sub("^BSI due to ", "", lab))
    den <- switch(lab, `BSI with MDR` = D$mdr %in% c("Yes","No"), `BSI with XDR` = D$xdr %in% c("Yes","No"), `BSI with Any resistance` = !is.na(D$R_any), `BSI with Access-resistant` = !is.na(D$R_Access),
      `BSI with Watch-resistant` = !is.na(D$R_Watch), `BSI with Reserve-resistant` = !is.na(D$R_Reserve), `BSI with 3GC-resistant Enterobacterales` = !is.na(D$ESBL_3GCR), `BSI with Carbapenem-resistant Gram-negative` = !is.na(D$CR),
      `BSI with MRSA` = !is.na(D$MRSA), `BSI with VRE` = !is.na(D$VR), !is.na(ind))
    prop_ci(ind, den) }
  bind_cols(tibble(indicator = r[[1]], prevalence_among_isolates = pct_txt(prev), isolates_n_over_N = sprintf("%d/%d", prev$n_num, prev$n_den)), inc) }))
cat("\n== Merged prevalence (% of isolates) + incidence ==\n")
print(as.data.frame(merged %>% mutate(inc_txt = sprintf("%s (%s–%s)", fmt(per1000adm, 2), fmt(adm_lo, 2), fmt(adm_hi, 2))) %>% select(indicator, prevalence_among_isolates, isolates_n_over_N, n_cases, inc_txt, one_in, national_cases)), digits = 3)

write_xlsx(list(CONSORT = consort, Pathogen_distribution = dist, Level1_agnostic = level1, Level2_by_pathogen = level2, Level3_drug_bug = db, Prevalence_incidence_merged = merged),
           file.path(OUT, "04_prevalence_incidence_tables.xlsx"))
saveRDS(list(D = D, consort = consort, dist = dist, level1 = level1, level2 = level2, db = db, merged = merged), file.path(OUT, "04_results.rds"))
cat("\nLevel 1:\n"); print(as.data.frame(level1 %>% select(indicator, n_resistant, n_denominator, ci_text)))
cat("\nLevel 2:\n"); print(as.data.frame(level2))
