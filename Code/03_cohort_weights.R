# 03_cohort_weights.R -- Case report base, facility ID cleaning, origin, design weights, regional calibration, isolate linkage
source(file.path(Sys.getenv("KGZ_PKG"), "R", "config.R"))
suppressPackageStartupMessages(library(lubridate))
res <- readRDS(file.path(OUT, "02_resistance.rds")); iso <- res$iso

clean_facility_id <- function(x) { x <- toupper(trimws(x)); has <- grepl("^HF ?[0-9]{1,3}", x); d <- sub("^HF ?([0-9]{1,3}).*$", "\\1", x)
  out <- rep(NA_character_, length(x)); out[has] <- sprintf("HF%03d", as.integer(d[has])); out }

cr <- rd("Case report")
n_raw <- nrow(cr)
cr <- cr %>% rename(rec_id = `Record ID`, patient_id = `Unique survey patient ID`, facility_raw = `Unique facility ID`)
is_junk <- is.na(cr$patient_id) | grepl("DATABASE|TEST", cr$patient_id, ignore.case = TRUE)
cr <- cr[!is_junk, ]; n_after_junk <- nrow(cr)
dup <- duplicated(cr$patient_id); cr_dups <- sum(dup)
# keep first occurrence of each patient ID, but never drop a record that holds a confirmed isolate link
cr <- cr %>% mutate(has_iso = rec_id %in% iso$case_rid) %>% arrange(desc(has_iso)) %>% filter(!duplicated(patient_id)) %>% select(-has_iso)
cat(sprintf("Case report: %d raw; %d after TEST/blank removal; %d duplicate patient IDs removed -> %d patients\n", n_raw, n_after_junk, cr_dups, nrow(cr)))
cr$facility_id_clean <- clean_facility_id(cr$facility_raw)

sf <- read.csv(SAMPLING_CSV, stringsAsFactors = FALSE)
cr <- cr %>% left_join(sf %>% select(UniqueID, NameAU), by = c(facility_id_clean = "UniqueID"))

# --- demographics, ward, origin ---------------------------------------------------------
parse_dt <- function(x) suppressWarnings(parse_date_time(as.character(x), orders = c("ymd HMS","ymd HM","ymd","dmy HMS","dmy HM","dmy"), quiet = TRUE))
col_by <- function(pat) { i <- grep(pat, names(cr), fixed = FALSE); if (!length(i)) stop("column not found: ", pat); i[1] }
yes <- function(v) toupper(trimws(as.character(v))) %in% c("YES","TRUE","1")
no  <- function(v) toupper(trimws(as.character(v))) %in% c("NO","FALSE","0")
cr$sex <- case_when(toupper(cr$`Sex at birth`) == "FEMALE" ~ "Female", toupper(cr$`Sex at birth`) == "MALE" ~ "Male", TRUE ~ NA_character_)
cr$age_grp <- case_when(grepl("^Neonate", cr$`Age group`) ~ "Neonate", grepl("^Paediatric", cr$`Age group`) ~ "Paediatric", grepl("^Adult", cr$`Age group`) ~ "Adult", TRUE ~ NA_character_)
sp_code <- cr$`Speciality code whitin the ward that is appicable to the patient`
cr$ward_type <- case_when(grepl("^(AICU|NICU|PICU)|Intensive", sp_code) ~ "ICU", !is.na(sp_code) ~ "Non-ICU", TRUE ~ NA_character_)
cr$specialty <- case_when(grepl("^ASW|^PSW", sp_code) | grepl("Surg", sp_code) ~ "Surgical", grepl("^AMW|^PMW|^NMW", sp_code) ~ "Medical", TRUE ~ NA_character_)

orig_raw <- toupper(trimws(cr[["Likely infection origin: Please indicate likely origin of the infection"]]))
origin_clin <- case_when(orig_raw == "HOSPITAL" ~ "Hospital-onset", orig_raw == "COMMUNITY" ~ "Community-onset", TRUE ~ NA_character_)
adm <- parse_dt(cr[["Hospital admission: Hospital admission date and time"]]); smp <- parse_dt(cr[["Blood sampling 1: 1. Date (dd/mm/yy) and time (24h) of sampling"]])
hrs <- ifelse(!is.na(adm) & !is.na(smp) & adm <= smp, as.numeric(difftime(smp, adm, units = "hours")), NA)
origin_48 <- case_when(is.na(hrs) ~ NA_character_, hrs > 48 ~ "Hospital-onset", TRUE ~ "Community-onset")
prox_cols <- c(grep("transferred in from another hospital|Was the patient transferred from another hospital", names(cr), value = TRUE),
               grep("discharged form a hospital", names(cr), value = TRUE),
               grep("transferred from another ward", names(cr), value = TRUE),
               grep("admitted to ICU for more than 48", names(cr), value = TRUE),
               grep("did the patient have surgery on this admission", names(cr), value = TRUE, ignore.case = TRUE))
P <- sapply(cr[prox_cols], yes); Pn <- sapply(cr[prox_cols], no)
any_pos <- rowSums(P) > 0; all_neg <- rowSums(Pn) == length(prox_cols)
cr$infection_origin <- case_when(!is.na(origin_clin) ~ origin_clin, !is.na(origin_48) ~ origin_48, any_pos ~ "Hospital-onset", all_neg ~ "Community-onset", TRUE ~ NA_character_)
cr$origin_source <- case_when(!is.na(origin_clin) ~ "clinician", !is.na(origin_48) ~ "48h timing", any_pos ~ "proxy+", all_neg ~ "proxy all-neg", TRUE ~ "unclassified")
cat("\nInfection origin source:\n"); print(table(cr$origin_source)); cat("proxy columns used:", length(prox_cols), "\n")

# --- Charlson (adapted, as legacy) ---------------------------------------------------------
cm <- function(lbl) yes(cr[[grep(paste0("Comorbidities at the time of hospital admission: ", lbl), names(cr), fixed = FALSE)[1]]])
cr$charlson <- with(list(), 1*(cm("Congestive heart failure")|cm("Chronic pulmonary")|cm("Connective tissue")|cm("Dementia")|cm("Cerebrovascular")|cm("Peptic")) +
  0)
w1 <- c("Congestive heart failure","Chronic pulmonary disease","Connective tissue rheumatologic","Dementia","Cerebrovascular disease","Peptic ulcer disease")
cr$charlson <- rowSums(sapply(w1, cm)) + ifelse(cm("Diabetes with end organ damage"), 2L, ifelse(cm("Diabetes$"), 1L, 0L)) +
  ifelse(cm("Moderate or severe liver"), 3L, ifelse(cm("Mild liver"), 1L, 0L)) + ifelse(cm("Metastic solid"), 6L, ifelse(cm("Malignant neoplasm \\(active"), 2L, 0L)) +
  ifelse(cm("HIV"), 6L, 0L) + ifelse(cm("Hemi- or paraplegia"), 2L, 0L) + ifelse(cm("Renal disease"), 2L, 0L)
cr$charlson_grp <- ifelse(cr$charlson > 0, "Charlson >0", "Charlson 0")

# --- facility characteristics, weights, calibration ----------------------------------------
fac <- rd("Facility") %>% transmute(fid = gsub("\\s+", "", `Unique facility ID`), level_of_care = `Level of care`,
  admissions = as.numeric(`Number of inpatient admissions: Number of inpatient admissions`),
  patient_days = as.numeric(`Total of inpatient days of care: Total inpatient days of care`))
sf$p_select <- sf$ApCL / sum(sf$ApCL); sf$design_weight <- 1 / (85 * sf$p_select)
pop2024 <- tribble(~NameAU, ~pop_k, "Batken",583.4, "Bishkek",1294.8, "Chui",957.6, "Issyk-Kul",544.4, "Jalal-Abad",1335.8, "Naryn",312.1, "Osh",1393.6+463.1, "Talas",277.1)
cr <- cr %>% left_join(sf %>% select(UniqueID, design_weight), by = c(facility_id_clean = "UniqueID")) %>% left_join(fac %>% select(fid, level_of_care), by = c(facility_id_clean = "fid"))
reg <- cr %>% filter(!is.na(design_weight)) %>% group_by(NameAU) %>% summarise(wn = sum(design_weight), .groups = "drop") %>% left_join(pop2024, by = "NameAU") %>%
  mutate(rake = (pop_k / sum(pop2024$pop_k)) / (wn / sum(wn)))
cr <- cr %>% left_join(reg %>% select(NameAU, rake), by = "NameAU") %>% mutate(calibrated_weight = design_weight * rake)
fac_w <- sf %>% select(UniqueID, NameAU, design_weight) %>% left_join(reg %>% select(NameAU, rake), by = "NameAU") %>% mutate(calibrated_weight = design_weight * rake)
fac <- fac %>% left_join(fac_w, by = c(fid = "UniqueID")) %>% filter(!is.na(design_weight))
cat("\nFacilities with weights:", nrow(fac), "; patients with weight:", sum(!is.na(cr$calibrated_weight)), "of", nrow(cr), "\n")
print(reg)

# --- link isolates to patients ---------------------------------------------------------------------
link <- iso %>% left_join(cr %>% select(rec_id, patient_id, facility_id_clean, NameAU, sex, age_grp, ward_type, specialty, infection_origin, charlson_grp, level_of_care, design_weight, calibrated_weight), by = c("case_rid" = "rec_id"))
cat("\nIsolates linked to a Case report patient:", sum(!is.na(link$patient_id)), "of", nrow(link), "\n")
cat("Isolates with calibrated weight:", sum(!is.na(link$calibrated_weight)), "\n")
cat("Isolates per patient (patients with >1 isolate row):", sum(duplicated(link$patient_id[!is.na(link$patient_id)])), "\n")
print(table(link$infection_origin, useNA = "ifany")); print(table(link$NameAU, useNA = "ifany"))
saveRDS(list(cr = cr, fac = fac, link = link, reg = reg, pop2024 = pop2024, sf = sf,
             counts = list(n_raw = n_raw, n_after_junk = n_after_junk, dups = cr_dups, n_patients = nrow(cr))), file.path(OUT, "03_cohort.rds"))
