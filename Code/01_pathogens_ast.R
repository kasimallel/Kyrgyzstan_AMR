# 01_pathogens_ast.R -- true-pathogen sample, species naming chain, AST long table
source(file.path(Sys.getenv("KGZ_PKG"), "R", "config.R"))
crl <- rd("Central ref lab")
stopifnot("confrimed_ID" %in% names(crl))
crl$confirmed_ID <- crl$confrimed_ID

consort <- list(crl_rows = nrow(crl), by_label = table(crl$confirmed_ID, useNA = "ifany"))
p <- crl %>% filter(confirmed_ID == "pathogen")
cat("True-pathogen rows:", nrow(p), "\n")

# ---- name normalisation ------------------------------------------------------
FUNGI <- "candida|cryptococc|aspergill|malassezia|fung|yeast"
COns_pat <- "^staphylococcus (epidermidis|haemolyticus|hominis|lentus|capitis|warneri|saprophyticus)$"   # unspecified "Staphylococcus sp." is KEPT
VGS <- c("mitis","oralis","mutans","salivarius","sanguinis","parasanguinis","anginosus","constellatus","viridans","gordonii","intermedius")
norm_name <- function(x) {
  x <- str_squish(x)
  x <- ifelse(x %in% c("", "Other", "other"), NA, x)
  lx <- tolower(x)
  x <- case_when(
    str_detect(lx, "^klebsiella pneumoniae") ~ "Klebsiella pneumoniae",
    str_detect(lx, "^enterobacter cloacae") ~ "Enterobacter cloacae",
    str_detect(lx, "^salmonella enterica typhi") ~ "Salmonella Typhi",
    str_detect(lx, "^salmonella") ~ "Salmonella spp. (non-Typhi)",
    str_detect(lx, "^kocuria") ~ "Kocuria varians",
    str_detect(lx, "^escherichia sp") ~ "Escherichia coli",   # genus-level Escherichia = E. coli
    TRUE ~ x)
  x
}
class_name <- function(nm) {
  l <- tolower(nm)
  case_when(
    is.na(nm) ~ "none",
    str_detect(l, FUNGI) ~ "fungal",
    str_detect(l, COns_pat) ~ "CoNS",
    TRUE ~ "ok")
}
group_name <- function(nm) {
  l <- tolower(nm)
  case_when(
    is.na(nm) ~ NA_character_,
    l == "streptococcus pyogenes" ~ "Streptococcus pyogenes",
    l == "streptococcus agalactiae" ~ "Streptococcus agalactiae",
    l == "streptococcus pneumoniae" ~ "Streptococcus pneumoniae",
    str_detect(l, "^streptococcus (viridans|mitis|oralis|mutans|salivarius|sanguinis|parasanguinis|anginosus|constellatus)") ~ "Viridans group streptococci",
    TRUE ~ nm)
}
cand <- function(sp, oth) {   # species field; if "Other"/NA use free-text "Species, other"
  nm <- ifelse(!is.na(sp) & sp != "Other", sp, oth)
  norm_name(nm)
}
p$c1 <- cand(p[["Species 1: Specify the species"]], p[["Species 1: Species, other"]])
p$c2 <- cand(p[["Species 2: Specify the species"]], p[["Species 2: Species, other"]])
p$k1 <- class_name(p$c1); p$k2 <- class_name(p$c2)

# --- fallback 3: "Identification: If commensal (contamination)" fields; fallback 4: species named in the laboratory free-text comments
p$m1 <- norm_name(p[["Identification: If comensal (contamination) - Species 1:"]]); p$m2 <- norm_name(p[["Identification: If comensal (contamination) - Species 2:"]])
genera <- unique(c(sub(" .*", "", na.omit(c(p$c1, p$c2, p$m1, p$m2))), "Kocuria", "Malassezia", "Rodentibacter", "Streptococcus", "Enterobacter"))
stopw <- c("does","that","from","with","were","this","identified","negative","result","was","not","and","the")
cm_species <- function(txt) {
  if (is.na(txt)) return(NA_character_)
  for (g in genera) {
    m <- regmatches(txt, regexec(paste0("(?i)\\b", g, "\\b\\s+([a-z]{4,})"), txt, perl = TRUE))[[1]]
    if (length(m) == 2 && !tolower(m[2]) %in% stopw) return(paste(g, tolower(m[2])))
    if (grepl(paste0("(?i)\\b", g, "\\b"), txt, perl = TRUE)) return(paste(g, "sp."))
  }
  NA_character_
}
p$c5 <- vapply(p[["Comments: Comments"]], cm_species, character(1))
p$k3 <- class_name(p$m1); p$k4 <- class_name(p$m2); p$k5 <- class_name(p$c5)
# If a commensal field says "Other", the species is named in the comment
p$name_pool <- NA_character_
cands <- list(c1 = c("k1", "Species 1", 1L), c2 = c("k2", "Species 2 (fallback)", 2L), m1 = c("k3", "Commensal-ID field (Species 1)", 0L), m2 = c("k4", "Commensal-ID field (Species 2)", 0L), c5 = c("k5", "Laboratory comment", 0L))
p$species <- NA_character_; p$name_source <- NA_character_; p$sample_used <- NA_integer_
for (nm in names(cands)) {
  ok <- is.na(p$species) & p[[cands[[nm]][1]]] == "ok"
  p$species[ok] <- p[[nm]][ok]; p$name_source[ok] <- cands[[nm]][2]; p$sample_used[ok] <- as.integer(cands[[nm]][3])
}
p$excl_reason <- case_when(
  !is.na(p$species) ~ NA_character_,
  p$k1 == "fungal" | p$k2 == "fungal" | p$k3 == "fungal" | p$k4 == "fungal" | p$k5 == "fungal" ~ "Fungal (removed)",
  TRUE ~ "No species recorded in any field")
p$pathogen_group <- group_name(p$species)
# polymicrobial: the other sample also yields a valid but different pathogen
other_k <- ifelse(p$sample_used %in% 1L, p$k2, ifelse(p$sample_used %in% 2L, p$k1, "none")); other_c <- ifelse(p$sample_used %in% 1L, p$c2, p$c1)
p$polymicrobial <- !is.na(p$sample_used) & other_k == "ok" & group_name(other_c) != p$pathogen_group
p$other_species <- ifelse(p$polymicrobial, other_c, NA)

cat("\n== CONSORT-style counts ==\n")
print(consort$by_label)
cat("Pathogen rows:", nrow(p), "\n")
print(table(p$excl_reason, useNA = "ifany"))
cat("Analysis isolates:", sum(!is.na(p$sample_used)), "\n")
print(table(p$name_source, useNA = "ifany"))
cat("Polymicrobial (second valid, different pathogen):", sum(p$polymicrobial), "\n")
print(sort(table(p$pathogen_group), decreasing = TRUE))

# ---- AST long table for the sample used -----------------------------------------
ast_long <- map_dfr(1:2, function(s) {
  sub <- p %>% filter(sample_used == s)
  if (!nrow(sub)) return(NULL)
  getc <- function(df, nm) if (nm %in% names(df)) as.character(df[[nm]]) else rep(NA_character_, nrow(df))
  map_dfr(1:30, function(k) {
    pre <- sprintf("S%d - Antibiotic %d: ", s, k)
    tibble(record_id = sub[["Record ID"]], sample_used = s, slot = k,
           antibiotic_raw = as.character(sub[[paste0(pre, "Specify the antibiotic")]]),
           result_raw = as.character(sub[[paste0(pre, "AST results")]]),
           value_raw = as.character(sub[[paste0(pre, "Value (unit)")]]),
           new_abx_raw = getc(sub, paste0(pre, "Add new antibiotic")))
  })
}) %>% filter(!is.na(antibiotic_raw) | !is.na(result_raw))
cat("\nAST rows:", nrow(ast_long), "\n")
cat("\n== AST result vocabulary ==\n"); print(sort(table(ast_long$result_raw, useNA = "ifany"), decreasing = TRUE))
cat("\n== Antibiotic vocabulary ==\n"); print(sort(table(ast_long$antibiotic_raw, useNA = "ifany"), decreasing = TRUE))
cat("\n== 'Add new antibiotic' values ==\n"); print(table(ast_long$new_abx_raw, useNA = "ifany"))

keep <- c("Record ID","Case report record ID","Unique facility ID (from case report)","Unique ward ID (from case report)","Unique survey patient ID",
          "CRL ID","ID outcome: Conclusion","Comments: Comments","c1","c2","m1","m2","c5","k1","k2","sample_used","species","name_source","excl_reason","pathogen_group","polymicrobial","other_species")
saveRDS(list(pathogen_rows = p[, keep], ast_long = ast_long, consort = consort), file.path(OUT, "01_pathogens_ast.rds"))
write.csv(p[, keep], file.path(OUT, "01_pathogen_naming_audit.csv"), row.names = FALSE)
