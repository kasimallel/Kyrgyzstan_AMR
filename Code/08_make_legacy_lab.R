# 09_make_legacy_lab.R -- synthetic "confirmed lab" table in the legacy CentralRefLab layout (294 culture-confirmed isolates) + v2 flags
source(file.path(Sys.getenv("KGZ_PKG"), "R", "config.R"))
res <- readRDS(file.path(OUT, "02_resistance.rds")); co <- readRDS(file.path(OUT, "03_cohort.rds"))
iso <- res$iso; ast <- res$ast; drug <- res$drug
lk <- co$link %>% select(record_id, patient_id)
iso <- iso %>% left_join(lk, by = "record_id")
stopifnot(!anyNA(iso$patient_id), !anyDuplicated(iso$patient_id))
nm <- setNames(drug$drug, drug$code)
syn <- tibble(`Unique survey patient ID` = iso$patient_id, `ID outcome: Conclusion` = "MICROBIOLOGICALLY", `ID outcome: Date of test result` = NA_character_,
              `Species 1: Specify the species` = iso$species, `Species 2: Specify the species` = NA_character_)
a <- ast %>% filter(record_id %in% iso$record_id) %>% mutate(drugname = nm[agent]) %>% group_by(record_id) %>% mutate(slot = row_number()) %>% ungroup()
for (k in 1:28) {
  sl <- a %>% filter(slot == k) %>% select(record_id, drugname, res_final)
  m <- match(iso$record_id, sl$record_id)
  syn[[paste0("S1 - Antibiotic ", k, ": Specify the antibiotic")]] <- sl$drugname[m]
  syn[[paste0("S1 - Antibiotic ", k, ": AST results")]] <- sl$res_final[m]     # "R" (incl. I/I*) or "S" or NA
}
stopifnot(max(a$slot) <= 28)
flags <- iso %>% transmute(patient_id, mdr_v2 = as.integer(mdr == "Yes") + ifelse(is.na(mdr), NA, 0L), xdr_v2 = as.integer(xdr == "Yes") + ifelse(is.na(xdr), NA, 0L), pdr_v2 = as.integer(pdr == "Yes") + ifelse(is.na(pdr), NA, 0L),
  any_v2 = as.integer(R_any), access_v2 = as.integer(R_Access), watch_v2 = as.integer(R_Watch), reserve_v2 = as.integer(R_Reserve), n_cat_tested = n_cat_tested, n_cat_ns = n_cat_ns,
  species = species, pathogen_group = pathogen_group, polymicrobial = polymicrobial)
saveRDS(list(syn_lab = syn, flags = flags), file.path(V2, "legacy/syn_lab.rds"))
cat("synthetic lab rows:", nrow(syn), " with AST:", sum(!is.na(syn$`S1 - Antibiotic 1: AST results`)), "\n")
old <- readxl::read_excel(DHIS2_XLSX, sheet = "CaseReport", guess_max = 1e6)
cat("patient IDs also in old CaseReport extract:", sum(iso$patient_id %in% old[["Unique survey patient ID"]]), "of", nrow(iso), "\n")
