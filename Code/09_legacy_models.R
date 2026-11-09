#---------------------------------------------------------#
#title: "WHO Kyrgistan
#date: "August-2026"
#---------------------------------------------------------#

#Packages #########
options(error = function() { try(save.image("ws_error.RData")); q(status = 1) })
library(readr)
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(writexl)
library(purrr)
library(grid)
library(survey)
library(mice)
library(broom)
library(nnet)

#########

#Main paths: ######
source(file.path(Sys.getenv("KGZ_PKG"), "R", "config.R"))
setwd(file.path(PKG, "legacy", "out"))
######



# Import the Excel files: Datasets #######
# clean_patient_sheet(): applied to every patient-level raw sheet IMMEDIATELY
# on read, before any other processing touches it, so the cleanup propagates
# through the entire script (every downstream section, not just the ones
# that happened to deduplicate for their own purposes, e.g. SECTION 8's
# merge previously deduplicated CaseReport for the mortality/LOS/CRF-value
# analyses but left SECTIONS 2-7's prevalence/incidence/AMR/risk-factor/
# diagnostic-stewardship analyses running on the raw, undeduplicated
# extract -- exactly the inconsistency this fixes). Removes, in order:
#   1. rows with a missing patient ID (cannot be analysed or matched);
#   2. non-patient QA/test records -- any ID containing "DATABASE" or
#      "TEST" (case-insensitive), e.g. "DATABASE TEST", "test_180625",
#      "TEST C3", "HF067/W01/001 TEST" -- confirmed present across
#      CaseReport, Discharge, CohortEnrolment, and CentralRefLab;
#   3. duplicate IDs (keep first occurrence) -- confirmed present in
#      CaseReport (19), FollowUpD28 (8), and Discharge (10).
clean_patient_sheet <- function(df, id_col) {
  id <- df[[id_col]]
  is_junk <- is.na(id) | grepl("DATABASE|TEST", id, ignore.case = TRUE)
  df <- df[!is_junk, , drop = FALSE]
  df[!duplicated(df[[id_col]]), , drop = FALSE]
}

## --- guess_max = 1e6 on every read_excel() call, throughout this script ---
# CRITICAL: readxl's default column-type guessing samples only the first
# 1000 rows per column. Several sparse text columns (Pitt BSI Score
# components, admission diagnosis/primary-infection-site, antibiotic-history
# route/dates, Discharge's "Treatment history" fields) have their real
# values clustered beyond row 1000, so the default guess silently typed
# them as logical and NULLED OUT every real value below the sample window --
# not "empty in this extract", corrupted on read. Confirmed via a
# readxl-default-vs-guess_max=1e6 diff: 41 CaseReport columns and 48
# Discharge columns were affected (recovered values ranged 1-997 per
# column); CentralRefLab and CohortEnrolment were unaffected. guess_max=1e6
# forces a full-file scan for every read_excel() call below and at every
# other read_excel() call site in this script (grep 'read_excel(' to find
# them all if adding a new one).
Main_data_CaseReport <- read_excel(
  DHIS2_XLSX,
  sheet = "CaseReport", guess_max = 1e6
) %>% clean_patient_sheet("Unique survey patient ID")

Main_data_hospital <- read_excel(
  DHIS2_XLSX,
  sheet = "Facility", guess_max = 1e6)
#View(Main_data_hospital)

Main_data_hospitalWard <- read_excel(
  DHIS2_XLSX,
  sheet = "Facility__WardData", guess_max = 1e6)
#View(Main_data_hospitalWard)


Main_data_hospitalLab <- read_excel(
  DHIS2_XLSX,
  sheet = "CentralRefLab", guess_max = 1e6) %>% clean_patient_sheet("Unique survey patient ID")

Main_data_hospitalFollowup28 <- read_excel(
  DHIS2_XLSX,
  sheet = "FollowUpD28", guess_max = 1e6) %>% clean_patient_sheet("Unique survey patient ID #1")

Main_data_hospitalDischarge <- read_excel(
  DHIS2_XLSX,
  sheet = "Discharge", guess_max = 1e6) %>% clean_patient_sheet("Unique survey patient ID #1")

Main_data_hospitalCohortEnrolment <- read_excel(
  DHIS2_XLSX,
  sheet = "CohortEnrolment", guess_max = 1e6) %>% clean_patient_sheet("Unique survey patient ID (___/C3)")

syn_obj <- readRDS(file.path(PKG, "legacy", "syn_lab.rds")); syn_lab <- syn_obj$syn_lab; syn_flags <- syn_obj$flags
#########



# =============================================================================
# SECTION 1: Main_data_CaseReport
# =============================================================================
# Rebuilt from scratch for the amr-survey_kgz_2026-08-25.xlsx CaseReport
# sheet, whose form structure differs materially from the earlier
# Survey_Data_Pre_EQA_270726.xlsx (age-stratified severity signs split into
# separate "Adult patient"/"Paediatric patient"/"Neonate patient" domains;
# several domain labels renamed, e.g. "Diagnosis on admission..." split into
# "Primary reason for admission to hospital..." + "Clinically suspected
# primary site of infection...").

## --- date parsing helper (ISO 8601 strings, e.g. "2024-07-10T07:10:00.000Z") ---
parse_dt <- function(x) as.POSIXct(x, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC")

## --- Registry tables: (original column name -> new variable name) -------

# Signs that are asked identically across 2-3 age domains (Adult/Paediatric/
# Neonate) are collapsed into ONE integrated variable, `group` marks which
# ones get combined. A patient only ever has one age domain populated (rows
# with more than one non-NA age-domain value: 0, verified against source),
# so "true in ANY applicable column" is a safe combination rule.
integrated_map <- tribble(
  ~orig,                                                                                      ~new,                              ~group,
  "Adult patient (≥ 18 years): Meningitis",                                                    "sev_meningitis",                  "meningitis",
  "Paediatric patient (≥ 29 days and <18 years): Meningitis",                                   "sev_meningitis",                  "meningitis",
  "Neonate patient (≤ 28 days): Meningitis",                                                     "sev_meningitis",                  "meningitis",
  "Adult patient (≥ 18 years): Pneumonia",                                                       "sev_pneumonia",                   "pneumonia",
  "Paediatric patient (≥ 29 days and <18 years): Pneumonia",                                     "sev_pneumonia",                   "pneumonia",
  "Neonate patient (≤ 28 days): Pneumonia",                                                       "sev_pneumonia",                   "pneumonia",
  "Adult patient (≥ 18 years): Additional signs indicating possible BSI",                        "sev_additional_signs_bsi",        "additional_signs",
  "Paediatric patient (≥ 29 days and <18 years): Additional signs indicating possible BSI",      "sev_additional_signs_bsi",        "additional_signs",
  "Neonate patient (≤ 28 days): Additional signs indicating possible BSI",                        "sev_additional_signs_bsi",        "additional_signs",
  "Adult patient (≥ 18 years): Indication for parenteral antibiotics",                            "sev_parenteral_abx_indicated",    "parenteral_abx",
  "Paediatric patient (≥ 29 days and <18 years): Indication for parenteral antibiotics",          "sev_parenteral_abx_indicated",    "parenteral_abx",
  "Neonate patient (≤ 28 days): Indication for parenteral antibiotics",                            "sev_parenteral_abx_indicated",    "parenteral_abx",
  "Paediatric patient (≥ 29 days and <18 years): Convulsion",                                     "sev_convulsion",                  "convulsion",
  "Neonate patient (≤ 28 days): Convulsion",                                                        "sev_convulsion",                  "convulsion"
)

# Age-domain "Other"/"Please, specify" free-text columns are deliberately
# excluded (removed per request) -- not referenced anywhere below.

# Signs specific to one age domain only -- kept as separate variables,
# prefixed by age group. Source values are "true"/blank (checkbox style,
# no explicit false) -> 1 = true, 0 = anything else.
agespecific_map <- tribble(
  ~orig, ~new,
  "Adult patient (≥ 18 years): Case definition met with ≥ 2 of the following signs (qSOFA)", "adult_case_definition_qsofa",
  "Adult patient (≥ 18 years): Abnormal core temperature (<36°C OR >38°C)", "adult_abn_temp_36_38",
  "Adult patient (≥ 18 years): Hypotension: systolic blood pressure (SBP) ≤ 100 mmHg", "adult_hypotension_sbp100",
  "Adult patient (≥ 18 years): Altered mental status: Glasgow Coma Scale < 15", "adult_altered_mental_status_gcs15",
  "Adult patient (≥ 18 years): Tachypnoea: Respiratory rate (RR) ≥22/min", "adult_tachypnoea_ge22",
  "Paediatric patient (≥ 29 days and <18 years): Case definition met with ≥ 1 of the following signs (LqSOFA)", "paed_case_definition_lqsofa",
  "Paediatric patient (≥ 29 days and <18 years): Capillary refill time (CRT) ≥3 seconds", "paed_crt_ge3s",
  "Paediatric patient (≥ 29 days and <18 years): Heart rate (HR) >99th centile thresholds", "paed_hr_gt99th_centile",
  "Paediatric patient (≥ 29 days and <18 years): Respiratory Rate (RR) >99th centile thresholds", "paed_rr_gt99th_centile",
  "Paediatric patient (≥ 29 days and <18 years): Alert, Voice, Pain, Unresponsive (AVPU) scale shows signs of VPU", "paed_avpu_vpu",
  "Paediatric patient (≥ 29 days and <18 years): Abnormal core temperature (<36°C or >38.5°C tympanic OR <35.5°C or >38°C axillary)", "paed_abn_temp_tympanic_axillary",
  "Paediatric patient (≥ 29 days and <18 years): Unable to drink or breastfeed", "paed_unable_drink_breastfeed",
  "Paediatric patient (≥ 29 days and <18 years): Vomiting everything", "paed_vomiting_everything",
  "Paediatric patient (≥ 29 days and <18 years): Lethargic or unconscious", "paed_lethargic_unconscious",
  "Paediatric patient (≥ 29 days and <18 years): Severe malnutrition", "paed_severe_malnutrition",
  "Neonate patient (≤ 28 days): Case definition met with ≥ 2 of the following signs", "neo_case_definition_2signs",
  "Neonate patient (≤ 28 days): Abnormal core temperature (>37.5°C or <35.5°C)", "neo_abn_temp_37_5_35_5",
  "Neonate patient (≤ 28 days): Inappropriate tachycardia", "neo_inappropriate_tachycardia",
  "Neonate patient (≤ 28 days): Poor feeding", "neo_poor_feeding",
  "Neonate patient (≤ 28 days): Reduced movement; change in level of activity (sleepiness, irritability, lethargy, floppiness)", "neo_reduced_movement",
  "Neonate patient (≤ 28 days): Respiratory Rate (RR) >60/min", "neo_resp_rate_gt60",
  "Neonate patient (≤ 28 days): Respiratory distress (severe chest in-drawing or grunting or cyanosis)", "neo_resp_distress"
)

# Healthcare exposure + invasive-devices YES/NO fields -> 1/0/NA
yesno_map <- tribble(
  ~orig, ~new,
  "Healthcare exposure prior to this hospital admission: Was the patient transferred in from another hospital?", "transferred_other_hospital",
  "Healthcare exposure prior to this hospital admission: Was the patient discharged form a hospital (same or another one) in the 48h preceding admission to this hospital?", "discharged_48h_prior",
  "Healthcare exposure prior to this hospital admission: Was the patient transferred from another ward in the same hospital?", "transferred_other_ward",
  "Invasive devices: Did the patient have invasive devices (catheter, intubation, etc.) inserted during this hospital episode prior to onset of clinical signs of suspected BSI?", "invasive_device_used",
  "Invasive devices: Was the patient admitted to ICU for more than 48 hours (2 days) on this admission prior to onset of clinical signs of suspected BSI?", "icu_48h_prior",
  "Invasive devices: Did the patient have surgery on this admission prior to onset of clinical signs of suspected BSI?", "surgery_prior_to_bsi"
)

# Specific invasive devices -- true/blank -> 1/0, kept separate
device_map <- tribble(
  ~orig, ~new,
  "Invasive devices: If yes, select device - Central IV catheter", "device_central_iv",
  "Invasive devices: If yes, select device - Intubation / Mechanical ventilation", "device_intubation_mechvent",
  "Invasive devices: If yes, select device - Peripheral IV catheter", "device_peripheral_iv"
)

# Date fields -- rename, keep as datetime
date_map <- tribble(
  ~orig, ~new,
  "Invasive devices: Date and time of onset of clinical signs and/or symptoms of BSI", "bsi_onset_datetime",
  "Hospital admission: Hospital admission date and time", "hospital_admissionDate",
  "Hospital admission: Hospital admission date", "hospital_admission_date_only",
  "Blood sampling 1: 1. Date (dd/mm/yy) and time (24h) of sampling", "blood_sample1_date",
  "Blood sampling 2: 2. Date (dd/mm/yy) and time (24h) of sampling", "blood_sample2_date"
)

# Primary reason for admission -- true/blank -> 1/0
reason_map <- tribble(
  ~orig, ~new,
  "Primary reason for admission to hospital (diagnosis on admission):: Cardiovascular", "reason_cardiovascular",
  "Primary reason for admission to hospital (diagnosis on admission):: Connective tissue/rheumatologic", "reason_connective_tissue",
  "Primary reason for admission to hospital (diagnosis on admission):: Gastrointestinal", "reason_gastrointestinal",
  "Primary reason for admission to hospital (diagnosis on admission):: Orthopaedic (non-trauma)", "reason_orthopaedic",
  "Primary reason for admission to hospital (diagnosis on admission):: Pulmonary", "reason_pulmonary",
  "Primary reason for admission to hospital (diagnosis on admission):: Renal", "reason_renal",
  "Primary reason for admission to hospital (diagnosis on admission):: Trauma", "reason_trauma",
  "Primary reason for admission to hospital (diagnosis on admission):: Unknown", "reason_unknown",
  "Primary reason for admission to hospital (diagnosis on admission):: Endocrine/Metabolic", "reason_endocrine_metabolic",
  "Primary reason for admission to hospital (diagnosis on admission):: Gynaecological", "reason_gynaecological",
  "Primary reason for admission to hospital (diagnosis on admission):: Genitourinary", "reason_genitourinary",
  "Primary reason for admission to hospital (diagnosis on admission):: Haematological", "reason_haematological",
  "Primary reason for admission to hospital (diagnosis on admission):: Neurological", "reason_neurological",
  "Primary reason for admission to hospital (diagnosis on admission):: Infectious disease", "reason_infectious_disease",
  "Primary reason for admission to hospital (diagnosis on admission):: Oncologic", "reason_oncologic",
  "Primary reason for admission to hospital (diagnosis on admission):: Dermatological", "reason_dermatological"
)

# Clinically suspected primary site of infection -- true/blank -> 1/0
site_map <- tribble(
  ~orig, ~new,
  "Clinically suspected primary site of infection: Central nervous system", "site_cns",
  "Clinically suspected primary site of infection: ENT/Upper respiratory tract", "site_ent_upper_resp",
  "Clinically suspected primary site of infection: Lower respiratory tract", "site_lower_resp_tract",
  "Clinically suspected primary site of infection: Pneumonia", "site_pneumonia",
  "Clinically suspected primary site of infection: Urinary tract", "site_urinary_tract",
  "Clinically suspected primary site of infection: Skin and soft tissue", "site_skin_soft_tissue",
  "Clinically suspected primary site of infection: Gastrointestinal", "site_gastrointestinal",
  "Clinically suspected primary site of infection: Intra-abdominal", "site_intra_abdominal",
  "Clinically suspected primary site of infection: Necrotising enterocolitis", "site_necrotising_enterocolitis",
  "Clinically suspected primary site of infection: Bone and joint", "site_bone_joint",
  "Clinically suspected primary site of infection: Cardiovascular system", "site_cardiovascular",
  "Clinically suspected primary site of infection: Surgical site", "site_surgical_site",
  "Clinically suspected primary site of infection: Genital", "site_genital",
  "Clinically suspected primary site of infection: Eye", "site_eye",
  "Clinically suspected primary site of infection: Febrile neutropenia", "site_febrile_neutropenia"
)

# Comorbidities at time of hospital admission -- YES/NO -> 1/0/NA
comorb_map <- tribble(
  ~orig, ~new,
  "Comorbidities at the time of hospital admission: Congestive heart failure", "comorb_chf",
  "Comorbidities at the time of hospital admission: Chronic pulmonary disease", "comorb_copd",
  "Comorbidities at the time of hospital admission: Mild liver disease", "comorb_liver_mild",
  "Comorbidities at the time of hospital admission: Moderate or severe liver disease", "comorb_liver_severe",
  "Comorbidities at the time of hospital admission: Diabetes", "comorb_diabetes",
  "Comorbidities at the time of hospital admission: Diabetes with end organ damage", "comorb_diabetes_end_organ",
  "Comorbidities at the time of hospital admission: Connective tissue rheumatologic, autoimmune", "comorb_connective_tissue",
  "Comorbidities at the time of hospital admission: Dementia", "comorb_dementia",
  "Comorbidities at the time of hospital admission: Cerebrovascular disease", "comorb_cerebrovascular",
  "Comorbidities at the time of hospital admission: Hemi- or paraplegia", "comorb_hemi_paraplegia",
  "Comorbidities at the time of hospital admission: Peptic ulcer disease", "comorb_peptic_ulcer",
  "Comorbidities at the time of hospital admission: Malignant neoplasm (active, past 6 months)", "comorb_malignancy_active",
  "Comorbidities at the time of hospital admission: Mental health disorder", "comorb_mental_health",
  "Comorbidities at the time of hospital admission: Malignant neoplasm (remission, >6 months)", "comorb_malignancy_remission",
  "Comorbidities at the time of hospital admission: Metastic solid tumour", "comorb_metastatic_tumour",
  "Comorbidities at the time of hospital admission: HIV (not on ART)", "comorb_hiv_no_art",
  "Comorbidities at the time of hospital admission: HIV (on ART)", "comorb_hiv_on_art",
  "Comorbidities at the time of hospital admission: HIV <200 CD4 cell count", "comorb_hiv_cd4_low",
  "Comorbidities at the time of hospital admission: Malaria", "comorb_malaria",
  "Comorbidities at the time of hospital admission: Tuberculosis (active)", "comorb_tb_active",
  "Comorbidities at the time of hospital admission: Tuberculosis (previous)", "comorb_tb_previous",
  "Comorbidities at the time of hospital admission: Solid organ transplant", "comorb_solid_organ_transplant",
  "Comorbidities at the time of hospital admission: Congenital malformation", "comorb_congenital_malformation",
  "Comorbidities at the time of hospital admission: Malnutrition", "comorb_malnutrition",
  "Comorbidities at the time of hospital admission: Hematopoietic stem-cell transplantation (HSCT)", "comorb_hsct",
  "Comorbidities at the time of hospital admission: Renal disease", "comorb_renal_disease"
)

# Risk factors at admission -- YES/NO -> 1/0/NA
riskfx_map <- tribble(
  ~orig, ~new,
  "Risk factors at admission: Hospitalisation in the 3 months prior to this admission", "riskfx_hospitalisation_3m",
  "Risk factors at admission: Regular hospital contact (e.g., dialysis, cancer treatment) in the 3 months prior to this admission", "riskfx_hospital_contact_3m",
  "Risk factors at admission: Surgery in the 3 months before admission", "riskfx_surgery_3m"
)

# Occurrences in the 48h prior to blood sample collection -- YES/NO -> 1/0/NA
occur48h_map <- tribble(
  ~orig, ~new,
  "In the 48 hours prior to the date of blood sample collection, did any of the following occur?: Intravenous vasopressor agent required", "vasopressor_48h",
  "In the 48 hours prior to the date of blood sample collection, did any of the following occur?: Mechanical ventilation needed", "mech_vent_48h",
  "In the 48 hours prior to the date of blood sample collection, did any of the following occur?: Cardiac arrest", "cardiac_arrest_48h"
)

pitt_temp_col <- "Pitt BSI Score: patient observations on the date of blood sample collection (or up to 48 hours prior): Temperature (ºC) - most extreme value"

abx_hist_cols <- paste0(
  "History of antibiotic therapy NO: Indicate the drug (International Non-propietary Name [INN]) ", 1:5
)

abx_current_status_col <- "Antimicrobial use history: Is the patient currently under antibiotic therapy?"
abx_current_drug_cols <- c(
  "Antimicrobial use history: Indicate the drug (International Non-propietary Name [INN]) #1",
  "Antimicrobial use history: Indicate the drug (International Non-propietary Name [INN]) #2",
  "Antimicrobial use history: Indicate the drug (International Non-propietary Name [INN]) #3",
  "Antimicrobial use history: Indicate the drug (International Non-propietary Name [INN]) #4"
)

sample_id_map <- tribble(
  ~orig, ~new,
  "Blood sampling 1: 1. Sample ID", "blood_sample1_ID",
  "Blood sampling 2: 2. Sample ID", "blood_sample2_ID"
)

origin_col2 <- "Likely infection origin: Please indicate likely origin of the infection"

core_cols <- c(
  "Unique survey patient ID", "Unique facility ID", "Unique ward/department ID",
  "Enrolment date", "F1. Date of transfer to inpatient ward",
  "Date of birth", "Sex at birth", "Age group",
  "Speciality code whitin the ward that is appicable to the patient"
)

## --- Keep only the CaseReport variables needed downstream ----------------
Main_data_CaseReport <- Main_data_CaseReport %>%
  select(
    all_of(core_cols),
    all_of(unique(integrated_map$orig)),
    all_of(agespecific_map$orig),
    all_of(yesno_map$orig),
    all_of(device_map$orig),
    all_of(date_map$orig),
    all_of(reason_map$orig),
    all_of(site_map$orig),
    all_of(comorb_map$orig),
    all_of(riskfx_map$orig),
    all_of(occur48h_map$orig),
    all_of(pitt_temp_col),
    all_of(abx_hist_cols),
    all_of(abx_current_status_col),
    all_of(abx_current_drug_cols),
    all_of(sample_id_map$orig),
    all_of(origin_col2)
  )

## --- ID / core variables: short names, no spaces --------------------------
Main_data_CaseReport <- Main_data_CaseReport %>%
  rename(
    patient_id       = `Unique survey patient ID`,
    facility_id      = `Unique facility ID`,
    ward_id          = `Unique ward/department ID`,
    enrolment_date   = `Enrolment date`,
    transfer_date_F1 = `F1. Date of transfer to inpatient ward`,
    dob              = `Date of birth`,
    speciality_code  = `Speciality code whitin the ward that is appicable to the patient`
  )

## --- sex_at_birth: FEMALE = 1, MALE = 0 ------------------------------------
Main_data_CaseReport <- Main_data_CaseReport %>%
  mutate(
    sex_at_birth = case_when(
      toupper(trimws(`Sex at birth`)) == "FEMALE" ~ 1L,
      toupper(trimws(`Sex at birth`)) == "MALE"   ~ 0L,
      TRUE ~ NA_integer_
    )
  ) %>%
  select(-`Sex at birth`)

## --- age_group: 0 = NEONATE, 1 = PAEDIATRIC, 3 = ADULT ---------------------
# NOTE: coding follows exactly what was requested (0/1/3, skipping 2).
Main_data_CaseReport <- Main_data_CaseReport %>%
  mutate(
    age_group = case_when(
      toupper(trimws(`Age group`)) == "NEONATE"    ~ 0L,
      toupper(trimws(`Age group`)) == "PAEDIATRIC" ~ 1L,
      toupper(trimws(`Age group`)) == "ADULT"      ~ 3L,
      TRUE ~ NA_integer_
    )
  ) %>%
  select(-`Age group`)

## --- Dates: enrolment/transfer, days between, DOB + age_years -------------
Main_data_CaseReport <- Main_data_CaseReport %>%
  mutate(
    enrolment_date   = parse_dt(enrolment_date),
    transfer_date_F1 = parse_dt(transfer_date_F1),
    dob              = parse_dt(dob),
    # days_enrolment_to_transfer: F1 transfer date minus enrolment date, in
    # days with decimals (datetime-level precision, not just calendar days).
    days_enrolment_to_transfer = round(as.numeric(difftime(transfer_date_F1, enrolment_date, units = "days")), 2),
    # age_years: current chronological age as of when this script is RUN
    # (Sys.time() - dob), per request ("subtract now to that date") -- NOT
    # age at any fixed clinical event, so this value changes on every re-run.
    age_years = round(as.numeric(difftime(Sys.time(), dob, units = "days")) / 365.25, 2),
    # Sanity filter (finally applied -- see the DATA QUALITY NOTE this
    # replaces, previously flagged but never acted on): a small number of DOB
    # entries are unambiguous data-entry errors, not real patients -- e.g. a
    # birth year of "0004" (age ~2021 years) and a birth date in 2030 (a
    # negative age). NA'd out here rather than silently left in, since an
    # implausible age_years value was propagating into every downstream
    # age-based model, subgroup table, and this report's own new Section-1
    # descriptive table (as a nonsensical SD of ~48 years). Affects 5 of 9222
    # rows (4 aged >110, 1 aged <0); no other row is touched.
    age_years = if_else(age_years < -0.5 | age_years > 110, NA_real_, age_years)
  )

## --- Integrated age-domain signs: true in ANY applicable age-domain column
for (grp in unique(integrated_map$group)) {
  sub <- integrated_map %>% filter(group == grp)
  new_name <- unique(sub$new)
  flags <- Main_data_CaseReport %>%
    select(all_of(sub$orig)) %>%
    mutate(across(everything(), ~ tolower(trimws(.x)) == "true"))
  Main_data_CaseReport[[new_name]] <- as.integer(rowSums(flags, na.rm = TRUE) > 0)
}
Main_data_CaseReport <- Main_data_CaseReport %>% select(-all_of(integrated_map$orig))

## --- Age-specific-only signs: true/blank -> 1/0 ----------------------------
Main_data_CaseReport <- Main_data_CaseReport %>%
  mutate(across(all_of(agespecific_map$orig), ~ if_else(tolower(trimws(.x)) == "true", 1L, 0L, missing = 0L))) %>%
  rename(!!!setNames(agespecific_map$orig, agespecific_map$new))

## --- Healthcare exposure / invasive-devices YES/NO -> 1/0/NA --------------
Main_data_CaseReport <- Main_data_CaseReport %>%
  mutate(across(all_of(yesno_map$orig), ~ case_when(
    toupper(trimws(.x)) == "NO"  ~ 0L,
    toupper(trimws(.x)) == "YES" ~ 1L,
    TRUE ~ NA_integer_
  ))) %>%
  rename(!!!setNames(yesno_map$orig, yesno_map$new))

## --- Specific invasive devices: true/blank -> 1/0 --------------------------
Main_data_CaseReport <- Main_data_CaseReport %>%
  mutate(across(all_of(device_map$orig), ~ if_else(tolower(trimws(.x)) == "true", 1L, 0L, missing = 0L))) %>%
  rename(!!!setNames(device_map$orig, device_map$new))

## --- Date fields: rename + parse -------------------------------------------
Main_data_CaseReport <- Main_data_CaseReport %>%
  rename(!!!setNames(date_map$orig, date_map$new)) %>%
  mutate(
    bsi_onset_datetime           = parse_dt(bsi_onset_datetime),
    hospital_admissionDate       = parse_dt(hospital_admissionDate),
    hospital_admission_date_only = parse_dt(hospital_admission_date_only),
    blood_sample1_date           = parse_dt(blood_sample1_date),
    blood_sample2_date           = parse_dt(blood_sample2_date)
  )

## --- Primary reason for admission: true/blank -> 1/0 -----------------------
Main_data_CaseReport <- Main_data_CaseReport %>%
  mutate(across(all_of(reason_map$orig), ~ if_else(tolower(trimws(.x)) == "true", 1L, 0L, missing = 0L))) %>%
  rename(!!!setNames(reason_map$orig, reason_map$new))

## --- Clinically suspected primary site of infection: true/blank -> 1/0 ----
Main_data_CaseReport <- Main_data_CaseReport %>%
  mutate(across(all_of(site_map$orig), ~ if_else(tolower(trimws(.x)) == "true", 1L, 0L, missing = 0L))) %>%
  rename(!!!setNames(site_map$orig, site_map$new))

## --- Comorbidities: YES/NO -> 1/0/NA ---------------------------------------
Main_data_CaseReport <- Main_data_CaseReport %>%
  mutate(across(all_of(comorb_map$orig), ~ case_when(
    toupper(trimws(.x)) == "NO"  ~ 0L,
    toupper(trimws(.x)) == "YES" ~ 1L,
    TRUE ~ NA_integer_
  ))) %>%
  rename(!!!setNames(comorb_map$orig, comorb_map$new))

## --- Risk factors at admission: YES/NO -> 1/0/NA ---------------------------
Main_data_CaseReport <- Main_data_CaseReport %>%
  mutate(across(all_of(riskfx_map$orig), ~ case_when(
    toupper(trimws(.x)) == "NO"  ~ 0L,
    toupper(trimws(.x)) == "YES" ~ 1L,
    TRUE ~ NA_integer_
  ))) %>%
  rename(!!!setNames(riskfx_map$orig, riskfx_map$new))

## --- Occurrences in 48h prior to sample collection: YES/NO -> 1/0/NA ------
Main_data_CaseReport <- Main_data_CaseReport %>%
  mutate(across(all_of(occur48h_map$orig), ~ case_when(
    toupper(trimws(.x)) == "NO"  ~ 0L,
    toupper(trimws(.x)) == "YES" ~ 1L,
    TRUE ~ NA_integer_
  ))) %>%
  rename(!!!setNames(occur48h_map$orig, occur48h_map$new))

## --- Pitt bacteraemia score (PARTIAL) --------------------------------------
# Of the 10 raw Pitt BSI Score sub-components in this extract, 5 are 100%
# empty across all 9,247 rows (Respiratory rate, Heart rate, Systolic BP,
# Acute hypotensive event, Mental status) and are EXCLUDED per "remove not
# plausible variables". Temperature has real data for only ~452 patients
# (~4.9%). Mechanical ventilation and Cardiac arrest come from a DIFFERENT
# domain (48h-prior-to-sampling occurrences) but ARE part of the classic
# Pitt score and DO have real data (~485-488 patients each).
# pitt_score_partial therefore uses only 3 of the 5 canonical components:
#   Temperature <36C or >=40C -> 1 point (standard abnormal-temp threshold)
#   Mechanical ventilation needed -> 2 points (standard Pitt weight)
#   Cardiac arrest -> 4 points (standard Pitt weight)
# Hypotension and Mental status axes are NOT scored (no usable data).
# pitt_score_partial_n_components (0-3) records how many of the 3 available
# inputs were non-missing for that patient -- treat scores based on 0-1
# components with caution.
Main_data_CaseReport <- Main_data_CaseReport %>%
  mutate(
    pitt_temp_c = suppressWarnings(as.numeric(.data[[pitt_temp_col]])),
    pitt_temp_points = if_else(!is.na(pitt_temp_c) & (pitt_temp_c < 36 | pitt_temp_c >= 40), 1L, 0L, missing = 0L),
    pitt_score_partial = pitt_temp_points +
      if_else(mech_vent_48h == 1L, 2L, 0L, missing = 0L) +
      if_else(cardiac_arrest_48h == 1L, 4L, 0L, missing = 0L),
    pitt_score_partial_n_components = (!is.na(pitt_temp_c)) + (!is.na(mech_vent_48h)) + (!is.na(cardiac_arrest_48h))
  ) %>%
  select(-all_of(pitt_temp_col), -pitt_temp_c)

## --- Antibiotic name homologation (shared helper) --------------------------
# Source values are messy (mixed case, a Cyrillic С visually identical to
# Latin C, transliteration/typo variants e.g. "CEFOZALIN"/"cefosaline" for
# Cefazolin, "AMOKCILLIN" for Amoxicillin). Recognized variants map to one
# canonical English INN name; anything not confidently recognized is kept
# as "UNMAPPED: <original value>" rather than guessed.
homologate_antibiotic <- function(x) {
  y <- str_squish(toupper(trimws(x)))
  y <- str_replace_all(y, "С", "C")  # Cyrillic С (U+0421) looks identical to Latin C
  case_when(
    is.na(y) | y == "" ~ NA_character_,
    str_detect(y, "GENTAM")                            ~ "Gentamicin",
    str_detect(y, "CEFTRIAX")                          ~ "Ceftriaxone",
    str_detect(y, "CEFAZOL|CEFOZOL|CEFOZAL|CEFOSAL|CEPHOSAL") ~ "Cefazolin",
    str_detect(y, "CEFOTAX|CEFATAX")                   ~ "Cefotaxime",
    str_detect(y, "CEFTAZID")                          ~ "Ceftazidime",
    str_detect(y, "CEFEP|CEPIME|CEPEP|CEPEFIME")       ~ "Cefepime",
    str_detect(y, "CEFIXIME")                          ~ "Cefixime",
    str_detect(y, "CEFPROZIL")                         ~ "Cefprozil",
    str_detect(y, "CEFUROXIME")                        ~ "Cefuroxime",
    str_detect(y, "CIPROFLOX|CIPROX|CYPROFLOX")        ~ "Ciprofloxacin",
    str_detect(y, "LEVOFLOX")                          ~ "Levofloxacin",
    str_detect(y, "MEROPEN|MEROPIN")                   ~ "Meropenem",
    str_detect(y, "METR.*DAZ")                         ~ "Metronidazole",
    str_detect(y, "VANCOM")                            ~ "Vancomycin",
    str_detect(y, "AMIKAC|AMICAC")                     ~ "Amikacin",
    str_detect(y, "AMPIC")                             ~ "Ampicillin",
    str_detect(y, "^AMOX|^AMOK|^AMOC")                 ~ "Amoxicillin",
    str_detect(y, "AZITHRO|AZITRO")                    ~ "Azithromycin",
    str_detect(y, "CLARITHROMYCIN")                    ~ "Clarithromycin",
    str_detect(y, "LINEZOLID")                         ~ "Linezolid",
    str_detect(y, "^PENIC")                            ~ "Penicillin",
    str_detect(y, "CEFALOSPORIN|CEPHALOSPORIN")        ~ "Cephalosporin (class, unspecified)",
    str_detect(y, "AMINOGLYCOSIDE")                    ~ "Aminoglycoside (class, unspecified)",
    TRUE ~ paste0("UNMAPPED: ", str_to_title(x))
  )
}

## --- History of prior antibiotic therapy (5 slots): homologate, rename ----
for (i in 1:5) {
  Main_data_CaseReport[[paste0("abx_hist_", i)]] <- homologate_antibiotic(Main_data_CaseReport[[abx_hist_cols[i]]])
}
Main_data_CaseReport <- Main_data_CaseReport %>% select(-all_of(abx_hist_cols))

## --- Currently under antibiotic therapy: true = 1, false = 0 --------------
# Unrecorded status (field blank/unparseable) is coded 0 (not on antibiotics),
# not NA: a blank answer here is treated as an absence of the exposure rather
# than missing data, per explicit instruction, so this covariate now enters
# every downstream model complete (no MICE imputation needed for it).
Main_data_CaseReport <- Main_data_CaseReport %>%
  mutate(
    abx_current_therapy = case_when(
      tolower(trimws(.data[[abx_current_status_col]])) == "true"  ~ 1L,
      tolower(trimws(.data[[abx_current_status_col]])) == "false" ~ 0L,
      TRUE ~ 0L
    )
  ) %>%
  select(-all_of(abx_current_status_col))

## --- Current-admission antibiotic drug names (renamed) ---------------------
# CORRECTED DATA QUALITY NOTE: these columns were previously believed 100%
# empty -- that was wrong, an artifact of readxl's default guess_max=1000
# column-type guessing silently nulling out values past row 1000 (see the
# guess_max=1e6 fix at this script's read_excel() calls, SECTION 1). With
# the full-file read, "Indicate the drug ... #1" alone has 1295 non-missing
# values (of 9247). Renamed here for structural completeness, same as
# before; per-antibiotic DURATION and gap-to-blood-culture variables still
# aren't built for this domain (start/end-date coverage remains partial
# even after the fix), but the drug-name signal itself is real and is now
# used to build abx_history immediately below.
for (i in 1:4) {
  Main_data_CaseReport[[paste0("abx_current_", i)]] <- Main_data_CaseReport[[abx_current_drug_cols[i]]]
}

## --- abx_history: any historical antibiotic drug name recorded (#1 slot) --
# 1 = "Antimicrobial use history: Indicate the drug (INN) #1" has non-blank
# text; 0 = blank. Distinct from abx_current_therapy (on/off status) and
# abx_exposure_48h (SECTION 12, timing-based). Candidate covariate --
# compared against abx_current_therapy ("current_abx") for model inclusion;
# see the comparison note where mortality/risk-factor adjustment_covars is
# defined.
Main_data_CaseReport <- Main_data_CaseReport %>%
  mutate(
    abx_history = as.integer(!is.na(.data[[abx_current_drug_cols[1]]]) & trimws(.data[[abx_current_drug_cols[1]]]) != "")
  ) %>%
  select(-all_of(abx_current_drug_cols))

## --- Blood sample IDs -------------------------------------------------------
Main_data_CaseReport <- Main_data_CaseReport %>%
  rename(!!!setNames(sample_id_map$orig, sample_id_map$new))

## --- Likely infection origin: hospital = 1, community = 0 ------------------
Main_data_CaseReport <- Main_data_CaseReport %>%
  mutate(
    infection_origin_hosp = case_when(
      toupper(trimws(.data[[origin_col2]])) == "HOSPITAL"  ~ 1L,
      toupper(trimws(.data[[origin_col2]])) == "COMMUNITY" ~ 0L,
      TRUE ~ NA_integer_
    )
  ) %>%
  select(-all_of(origin_col2))

## --- ward_type_icu (ICU vs Non-ICU) and specialty_group (Medical vs Surgical) ---
# Reused unchanged from the prior file version -- same Global-PPS
# WardCategory-Specialty code scheme, same speciality_code column.
ward_code_lookup <- tribble(
  ~code,        ~ward_type_icu, ~specialty_group,
  "PMW-PMW",    "Non-ICU",      "Medical",
  "AMW-AMW",    "Non-ICU",      "Medical",
  "AMW-OBG",    "Non-ICU",      "Medical",
  "ASW-ASW",    "Non-ICU",      "Surgical",
  "NMW-NICU",   "ICU",          NA_character_,
  "AICU-AICU",  "ICU",          NA_character_,
  "AMW-ID",     "Non-ICU",      "Medical",
  "PMW-ID",     "Non-ICU",      "Medical",
  "AMW-NEU",    "Non-ICU",      "Medical",
  "AICU-SUR",   "ICU",          "Surgical",
  "AMW-P",      "Non-ICU",      "Medical",
  "ASW-ORT",    "Non-ICU",      "Surgical",
  "PMW-PICU",   "ICU",          NA_character_,
  "AMW-CAR",    "Non-ICU",      "Medical",
  "PMW-PSW",    "Non-ICU",      "Surgical",
  "ASW-URO",    "Non-ICU",      "Surgical",
  "ASW-ONCO",   "Non-ICU",      "Surgical",
  "NMW-NMW",    "Non-ICU",      "Medical",
  "ASW-ENT",    "Non-ICU",      "Surgical",
  "ASW-DIG",    "Non-ICU",      "Surgical",
  "AICU-MED",   "ICU",          "Medical",
  "ASW-NEU",    "Non-ICU",      "Surgical",
  "ASW-CV",     "Non-ICU",      "Surgical",
  "NMW-NSW",    "Non-ICU",      "Surgical",
  "AMW-DB",     "Non-ICU",      "Medical",
  "AICU-CAR",   "ICU",          "Medical",
  "AMW-HO",     "Non-ICU",      "Medical",
  "ASW-PLAS",   "Non-ICU",      "Surgical",
  "AMW-REN",    "Non-ICU",      "Medical",
  "PMW-HO",     "Non-ICU",      "Medical"
)

Main_data_CaseReport <- Main_data_CaseReport %>%
  left_join(ward_code_lookup, by = c("speciality_code" = "code"))

# ward_type_icu_num: homologated 1/0 dummy (1 = true ICU admission, 0 = non-
# ICU admission, NA = unmapped/missing speciality_code) built from the text
# ward_type_icu ("ICU"/"Non-ICU"). Needed anywhere a numeric indicator is
# required -- e.g. SECTION 13's microbiological-outcome-group profile table
# previously compared this CHARACTER variable directly against numeric 1
# (`ward_type_icu == 1`), which is never TRUE for a string, silently
# returning a uniform, meaningless 0% for every outcome group.
Main_data_CaseReport <- Main_data_CaseReport %>%
  mutate(ward_type_icu_num = case_when(
    ward_type_icu == "ICU" ~ 1L,
    ward_type_icu == "Non-ICU" ~ 0L,
    TRUE ~ NA_integer_
  ))

## --- hosp_comm_acquis48: hospital- vs community-acquired from timing ------
# Hours from hospital_admissionDate (has time-of-day) to blood_sample1_date.
# Only computed when admission is at/before sampling (biologically valid
# sequence); >48h -> 1 (hospital-acquired), <=48h -> 0 (community-acquired).
Main_data_CaseReport <- Main_data_CaseReport %>%
  mutate(
    hours_admission_to_sample1 = if_else(
      !is.na(blood_sample1_date) & !is.na(hospital_admissionDate) & hospital_admissionDate <= blood_sample1_date,
      as.numeric(difftime(blood_sample1_date, hospital_admissionDate, units = "hours")),
      NA_real_
    ),
    hosp_comm_acquis48 = case_when(
      is.na(hours_admission_to_sample1) ~ NA_integer_,
      hours_admission_to_sample1 > 48 ~ 1L,
      TRUE ~ 0L
    )
  )

## --- infection_origin_composite: hierarchical CAI/HAI, built to close the
## coverage gap in infection_origin_hosp (48.5%) and hosp_comm_acquis48
## (23.9%). Both are limited by the same root cause: hospital_admissionDate
## itself is only 24.9% complete (2,306/9,247) -- switching the timing
## variable to use bsi_onset_datetime instead of blood_sample1_date would
## make this WORSE, not better (bsi_onset_datetime is only 11.8% complete).
## Coverage is instead extended using five near-complete (96-99%) structural
## exposure variables that each independently imply substantial healthcare
## contact before BSI onset. Applied as a 3-tier hierarchy, most to least
## direct evidence; infection_origin_source records which tier produced each
## patient's value, so a more conservative analysis can restrict to tier 1-2
## only if preferred:
##   tier 1 "clinician"        infection_origin_hosp, where recorded (4,485; 48.5%)
##   tier 2 "timing_48h"       hosp_comm_acquis48, where tier 1 missing but this is
##                             computable (adds 120; +1.3pp)
##   tier 3 "proxy_positive"   any of {transferred_other_ward, icu_48h_prior,
##                             surgery_prior_to_bsi, transferred_other_hospital,
##                             discharged_48h_prior} == 1 -> hospital-acquired.
##                             Each of these independently implies the patient was
##                             already an inpatient (this facility or another) or in
##                             recent healthcare contact before BSI onset (adds 796)
##   tier 3 "proxy_all_negative"  all five proxies recorded and none positive ->
##                             community-acquired. WEAKEST tier: absence of these
##                             specific exposures is suggestive, not proof, of
##                             community origin (adds 3,649)
## 197 patients (2.1%) remain unclassifiable (missing tier 1 and 2, and at
## least one of the five proxies also missing). Total coverage: 97.9%.
infection_origin_proxy_vars <- c("transferred_other_ward", "icu_48h_prior", "surgery_prior_to_bsi",
                                  "transferred_other_hospital", "discharged_48h_prior")
infection_origin_proxy_df <- Main_data_CaseReport %>% select(all_of(infection_origin_proxy_vars))
infection_origin_proxy_any_positive <- rowSums(infection_origin_proxy_df == 1L, na.rm = TRUE) > 0
infection_origin_proxy_all_negative <- rowSums(is.na(infection_origin_proxy_df)) == 0 &
  rowSums(infection_origin_proxy_df == 1L, na.rm = TRUE) == 0

Main_data_CaseReport <- Main_data_CaseReport %>%
  mutate(
    infection_origin_composite = case_when(
      !is.na(infection_origin_hosp) ~ infection_origin_hosp,
      !is.na(hosp_comm_acquis48) ~ hosp_comm_acquis48,
      infection_origin_proxy_any_positive ~ 1L,
      infection_origin_proxy_all_negative ~ 0L,
      TRUE ~ NA_integer_
    ),
    infection_origin_source = case_when(
      !is.na(infection_origin_hosp) ~ "clinician",
      !is.na(hosp_comm_acquis48) ~ "timing_48h",
      infection_origin_proxy_any_positive ~ "proxy_positive",
      infection_origin_proxy_all_negative ~ "proxy_all_negative",
      TRUE ~ NA_character_
    )
  )
rm(infection_origin_proxy_df, infection_origin_proxy_any_positive, infection_origin_proxy_all_negative)

## --- infection_origin_composite_strict: composite WITHOUT tier 3 -----------
# Same variable, but stops after tier 2 (clinician + 48h timing) and leaves
# tier 3 (the structural-proxy inference) as NA rather than guessing. This is
# the conservative companion to infection_origin_composite: identical values
# wherever both are non-missing (tiers 1-2 are shared), narrower coverage
# (48.5% + 1.3pp = 49.8%, vs 97.9% for the full composite), but nothing in it
# rests on the "5 proxies all negative -> community-acquired" inference that
# is the full composite's weakest link. Derived directly from
# infection_origin_source so the two variables cannot drift apart.
Main_data_CaseReport <- Main_data_CaseReport %>%
  mutate(
    infection_origin_composite_strict = if_else(
      infection_origin_source %in% c("clinician", "timing_48h"),
      infection_origin_composite,
      NA_integer_
    )
  )



# =============================================================================
# SECTION 2: Main_data_hospitalLab
# =============================================================================
# Everything in this section builds lab-side data from Main_data_hospitalLab
# (CentralRefLab) only -- species, AST results for up to 28 antibiotics per
# isolate for both Species 1 and 2, and the AMR resistance summary variables
# derived from them -- assembled into one frame, lab_full. Main_data_CaseReport
# is only touched ONCE, at the very end of the section, via a single merge:
# every CaseReport row is kept as-is, and only patients whose CentralRefLab
# conclusion is confirmed get lab_full's columns attached (everyone else gets
# NA for those columns).
# NOTE: in this new file, the conclusion column is renamed to "ID outcome:
# Conclusion" and its confirmed-infection value is truncated to just
# "MICROBIOLOGICALLY" (was "MICROBIOLOGICALLY_CONFIRMED_INFECTION" in the
# previous file) -- updated below to match.

lab_confirmed_species <- syn_lab %>%
  select(
    patient_id = `Unique survey patient ID`,
    pathogen1 = `Species 1: Specify the species`,
    test_result_date = `ID outcome: Date of test result`,
    pathogen2 = `Species 2: Specify the species`
  ) %>%
  mutate(
    lab_match = 1L,
    test_result_date = parse_dt(test_result_date)
  )

## --- Species 1/Species 2 fallback correction --------------------------------
# The form's Species-2 slot is meant for a genuine second (polymicrobial)
# isolate, but 21 of 490 microbiologically-confirmed lab records have
# Species 1 left blank while Species 2 carries the organism name (17 of
# those 21 also have determinate S2 antibiotic-susceptibility results) --
# a data-entry pattern where the isolate was recorded in the wrong slot,
# not a genuine second organism. Left uncorrected, pathogen1_clean (built
# from pathogen1 below) shows NA for these 21 patients, which silently
# excludes them from pathogen distribution (Table 3), the AMR heatmap, and
# every organism-restricted resistance mechanism (carbapenem/cephalosporin/
# 3GCR restricted to Enterobacterales/P. aeruginosa/A. baumannii; MRSA
# restricted to S. aureus; VRE restricted to Enterococcus) via
# is_cre_pa_ab()/genus1 downstream -- even though 8/1/3 of these 21 belong
# to those organism groups respectively. The category-level resistance
# flags themselves (carb_BSI_positive etc., SECTION 2 below) are NOT
# affected by this bug: they already pool AST results across all 56 S1+S2
# slots regardless of which slot the isolate landed in.
# Fix: when Species 1 is blank AND Species 2 is present, treat Species 2 as
# the sole (primary) isolate, not a second one. Genuine polymicrobial cases
# (both slots populated) and the 5 cases where both are blank are untouched.
lab_confirmed_species <- lab_confirmed_species %>%
  mutate(
    .p1_missing = is.na(pathogen1) | trimws(pathogen1) == "",
    .p2_present = !is.na(pathogen2) & trimws(pathogen2) != "",
    # Carried through to lab_full (below) to also swap that isolate's AST
    # results (abx_s1_*/ast_s1_* <-> abx_s2_*/ast_s2_*): the organism's
    # actual antibiotic panel was tested under the S2 slots (that's where
    # the lab recorded it), so the organism-restricted resistance tables
    # (which read carb_resistance_s1/ceph_resistance_s1/etc specifically,
    # not the pooled S1+S2 category flags) would otherwise still show these
    # patients as untested despite now having a correct pathogen1_clean.
    pathogen_swapped = .p1_missing & .p2_present,
    pathogen1 = if_else(pathogen_swapped, pathogen2, pathogen1),
    pathogen2 = if_else(pathogen_swapped, NA_character_, pathogen2)
  ) %>%
  select(-.p1_missing, -.p2_present)

## --- Species 1 antibiotic susceptibility (AST) results, up to 28 slots ---
# Everything through the end of this SECTION builds on Main_data_hospitalLab
# only -- Main_data_CaseReport is not touched until the single merge at the
# very end of the section.
# The form allows up to 28 antibiotics tested per isolate (columns
# "S1 - Antibiotic N: Specify the antibiotic" / "S1 - Antibiotic N: AST
# results", N = 1..28); in this extract only 24 of the 28 slots exist as
# columns at all (20, 23, 26, 27 were never used) -- generated generically
# below so any slot populated in a future extract is picked up automatically.
# Antibiotic names here come from a closed dropdown (not free text like the
# antibiotic-history fields), so "homologate" just means stripping the
# parenthetical alternate name some options carry, e.g.
# "Trimethoprim-sulfamethoxazole (Co-trimoxazole)" -> "Trimethoprim-sulfamethoxazole".
# AST results: R and I -> 1, S -> 0, anything else (e.g. "U", blank) -> NA.
homologate_antibiotic_dropdown <- function(x) {
  str_remove(trimws(x), "\\s*\\([^)]*\\)\\s*$")
}

s1_ast_map <- tibble(n = 1:28) %>%
  mutate(
    spec_orig = paste0("S1 - Antibiotic ", n, ": Specify the antibiotic"),
    ast_orig  = paste0("S1 - Antibiotic ", n, ": AST results"),
    abx_new   = paste0("abx_s1_", n),
    ast_new   = paste0("ast_s1_", n)
  )

lab_ast_s1 <- syn_lab %>%
  select(
    patient_id = `Unique survey patient ID`,
    any_of(s1_ast_map$spec_orig),
    any_of(s1_ast_map$ast_orig)
  )

for (i in 1:28) {
  spec_o <- s1_ast_map$spec_orig[i]; ast_o <- s1_ast_map$ast_orig[i]
  abx_n  <- s1_ast_map$abx_new[i];   ast_n <- s1_ast_map$ast_new[i]
  lab_ast_s1[[abx_n]] <- if (spec_o %in% names(lab_ast_s1)) {
    homologate_antibiotic_dropdown(lab_ast_s1[[spec_o]])
  } else {
    NA_character_
  }
  lab_ast_s1[[ast_n]] <- if (ast_o %in% names(lab_ast_s1)) {
    case_when(
      toupper(trimws(lab_ast_s1[[ast_o]])) %in% c("R", "I") ~ 1L,
      toupper(trimws(lab_ast_s1[[ast_o]])) == "S"            ~ 0L,
      TRUE ~ NA_integer_
    )
  } else {
    NA_integer_
  }
}
lab_ast_s1 <- lab_ast_s1 %>% select(patient_id, starts_with("abx_s1_"), starts_with("ast_s1_"))

## --- Species 2 antibiotic susceptibility (AST) results, up to 28 slots ---
# Mirrors the Species 1 block above exactly, just against the "S2 -
# Antibiotic N: ..." columns (23 of 28 slots exist in this extract: 1-19,
# 21, 22, 25, 28 -- a different subset than Species 1's, handled the same
# generic way). Reuses homologate_antibiotic_dropdown() from the S1 block.
s2_ast_map <- tibble(n = 1:28) %>%
  mutate(
    spec_orig = paste0("S2 - Antibiotic ", n, ": Specify the antibiotic"),
    ast_orig  = paste0("S2 - Antibiotic ", n, ": AST results"),
    abx_new   = paste0("abx_s2_", n),
    ast_new   = paste0("ast_s2_", n)
  )

lab_ast_s2 <- syn_lab %>%
  select(
    patient_id = `Unique survey patient ID`,
    any_of(s2_ast_map$spec_orig),
    any_of(s2_ast_map$ast_orig)
  )

for (i in 1:28) {
  spec_o <- s2_ast_map$spec_orig[i]; ast_o <- s2_ast_map$ast_orig[i]
  abx_n  <- s2_ast_map$abx_new[i];   ast_n <- s2_ast_map$ast_new[i]
  lab_ast_s2[[abx_n]] <- if (spec_o %in% names(lab_ast_s2)) {
    homologate_antibiotic_dropdown(lab_ast_s2[[spec_o]])
  } else {
    NA_character_
  }
  lab_ast_s2[[ast_n]] <- if (ast_o %in% names(lab_ast_s2)) {
    case_when(
      toupper(trimws(lab_ast_s2[[ast_o]])) %in% c("R", "I") ~ 1L,
      toupper(trimws(lab_ast_s2[[ast_o]])) == "S"            ~ 0L,
      TRUE ~ NA_integer_
    )
  } else {
    NA_integer_
  }
}
lab_ast_s2 <- lab_ast_s2 %>% select(patient_id, starts_with("abx_s2_"), starts_with("ast_s2_"))

## --- Combine all lab-side pieces into one frame (still no CaseReport merge) ---
lab_full <- lab_confirmed_species %>%
  left_join(lab_ast_s1, by = "patient_id", na_matches = "never") %>%
  left_join(lab_ast_s2, by = "patient_id", na_matches = "never")

## --- Swap S1/S2 AST blocks for the pathogen_swapped rows ------------------
# Mirrors the pathogen1/pathogen2 swap above: for these patients the S2
# antibiotic/AST columns hold the ONLY isolate's real testing data, but
# every organism-restricted resistance computation downstream (SECTION 3,
# e.g. carb_resistance_s1) reads the S1-labelled columns specifically as
# "this isolate's panel" -- so the S1/S2 column blocks themselves must swap
# too, not just the species name, or these patients still show as untested
# for their own organism group despite pathogen1_clean now being correct.
s1_col_names <- c(paste0("abx_s1_", 1:28), paste0("ast_s1_", 1:28))
s2_col_names <- c(paste0("abx_s2_", 1:28), paste0("ast_s2_", 1:28))
swap_idx <- which(lab_full$pathogen_swapped)
if (length(swap_idx) > 0) {
  tmp <- lab_full[swap_idx, s1_col_names]
  lab_full[swap_idx, s1_col_names] <- lab_full[swap_idx, s2_col_names]
  lab_full[swap_idx, s2_col_names] <- tmp
}
lab_full <- lab_full %>% select(-pathogen_swapped)

## --- AMR resistance summary variables (from abx_s1_*/ast_s1_*, abx_s2_*/ast_s2_*) ---
# Built by scanning all 56 antibiotic/AST slot pairs (S1 1-28 + S2 1-28)
# together, patient by patient, directly on lab_full (still lab-side only).
# Three-state logic throughout: 1 = resistance found, 0 = tested but never
# resistant, NA = never tested at all (so "not resistant" is never
# conflated with "not tested" -- a patient with no AST data gets NA, not 0).
#
# CAVEAT on MDR/XDR: the antibiotic_category_map below is a single,
# organism-AGNOSTIC set of ~20 antimicrobial categories (aminoglycosides,
# carbapenems, cephalosporins, etc.), used as a practical simplification.
# The rigorous Magiorakos et al. (2012) MDR/XDR/PDR definitions are
# organism-specific (Enterobacterales, Pseudomonas, Acinetobacter,
# Staphylococcus, and Enterococcus each have a DIFFERENT category panel).
# This dataset doesn't currently carry a variable stratifying category
# counts by organism, so mdr_resistance/xdr_resistance here are a
# category-count proxy, not the organism-specific validated definitions --
# treat as indicative, not a substitute for organism-stratified MDR review.
abx_cols <- c(paste0("abx_s1_", 1:28), paste0("abx_s2_", 1:28))
ast_cols <- c(paste0("ast_s1_", 1:28), paste0("ast_s2_", 1:28))

abx_mat <- as.matrix(lab_full[, abx_cols])
ast_mat <- as.matrix(lab_full[, ast_cols])
storage.mode(ast_mat) <- "integer"

# resistance_any: 1 if ANY of the 56 ast_s1_*/ast_s2_* slots is resistant.
any_tested_amr    <- rowSums(!is.na(ast_mat)) > 0
any_resistant_amr <- rowSums(ast_mat == 1L, na.rm = TRUE) > 0
lab_full$resistance_any <- case_when(
  any_resistant_amr ~ 1L,
  any_tested_amr ~ 0L,
  TRUE ~ NA_integer_
)

# Antibiotic -> antimicrobial category, covering every distinct homologated
# name observed across all S1+S2 slots in this extract (verified: 0 unmapped).
antibiotic_category_map <- tribble(
  ~antibiotic, ~category,
  "Trimethoprim-sulfamethoxazole", "Folate pathway inhibitors",
  "Gentamicin", "Aminoglycosides",
  "Tigecycline", "Glycylcyclines",
  "Levofloxacin", "Fluoroquinolones",
  "Clindamycin", "Lincosamides",
  "Vancomycin", "Glycopeptides",
  "Linezolid", "Oxazolidinones",
  "Rifampicin", "Rifamycins",
  "Teicoplanin", "Glycopeptides",
  "Daptomycin", "Lipopeptides",
  "Tetracycline", "Tetracyclines",
  "Erythromycin", "Macrolides",
  "Oxacillin", "Penicillins",
  "Fusidic acid", "Fusidanes",
  "Amikacin", "Aminoglycosides",
  "Ciprofloxacin", "Fluoroquinolones",
  "Ampicillin", "Penicillins",
  "Ceftazidime", "Cephalosporins",
  "Meropenem", "Carbapenems",
  "Cefotaxime", "Cephalosporins",
  "Cefepime", "Cephalosporins",
  "Cefoxitin", "Cephalosporins",
  "Ertapenem", "Carbapenems",
  "Colistin", "Polymyxins",
  "Nitrofurantoin", "Nitrofurans",
  "Amoxicillin-clavulanic acid", "Beta-lactam/beta-lactamase-inhibitor combinations",
  "Netilmicin", "Aminoglycosides",
  "Fosfomycin", "Fosfomycins",
  "Benzylpenicillin", "Penicillins",
  "Ceftriaxone", "Cephalosporins",
  "Moxifloxacin", "Fluoroquinolones",
  "Ceftaroline", "Cephalosporins",
  "Chloramphenicol", "Phenicols",
  "Ampicillin-sulbactam", "Beta-lactam/beta-lactamase-inhibitor combinations",
  "Ceftazidime-avibactam", "Cephalosporins",
  "Fluconazole", "Antifungals (azoles)",
  "Imipenem", "Carbapenems",
  "Micafungin", "Antifungals (echinocandins)",
  "Voriconazole", "Antifungals (azoles)",
  "Amoxicillin", "Penicillins",
  "Amphotericin B", "Antifungals (polyenes)",
  "Anidulafungin", "Antifungals (echinocandins)",
  "Aztreonam", "Monobactams",
  "Cefuroxime", "Cephalosporins",
  "Itraconazole", "Antifungals (azoles)",
  "Meropenem-vaborbactam", "Carbapenems",
  "Nalidixic acid", "Quinolones (non-fluorinated)",
  "Piperacillin-tazobactam", "Beta-lactam/beta-lactamase-inhibitor combinations",
  "Posaconazole", "Antifungals (azoles)",
  "Tobramycin", "Aminoglycosides"
)

cat_lookup <- setNames(antibiotic_category_map$category, antibiotic_category_map$antibiotic)
cat_mat <- matrix(cat_lookup[abx_mat], nrow = nrow(abx_mat), ncol = ncol(abx_mat))

n_cat_tested <- integer(nrow(lab_full))
n_cat_resistant <- integer(nrow(lab_full))
for (i in seq_len(nrow(lab_full))) {
  cats_i <- cat_mat[i, ]
  ast_i  <- ast_mat[i, ]
  valid <- !is.na(cats_i) & !is.na(ast_i)
  if (any(valid)) {
    n_cat_tested[i] <- length(unique(cats_i[valid]))
    n_cat_resistant[i] <- length(unique(cats_i[valid & ast_i == 1L]))
  }
}
n_cat_susceptible_only <- n_cat_tested - n_cat_resistant

# mdr_resistance: MDR (Magiorakos-style, category-count proxy) = resistant
# to >=1 agent in >=3 antimicrobial categories.
lab_full$mdr_resistance <- case_when(
  n_cat_tested == 0 ~ NA_integer_,
  n_cat_resistant >= 3 ~ 1L,
  TRUE ~ 0L
)
# xdr_resistance: XDR = MDR AND susceptible in at most 2 of the categories
# actually tested (i.e. non-susceptible to all-but-<=2 tested categories).
lab_full$xdr_resistance <- case_when(
  n_cat_tested == 0 ~ NA_integer_,
  lab_full$mdr_resistance == 1L & n_cat_susceptible_only <= 2 ~ 1L,
  TRUE ~ 0L
)

# Shared helper for the specific-class flags below: 1 = resistant to
# >=1 antibiotic matching the target set, 0 = target class tested but never
# resistant, NA = target class never tested for this patient.
flag_resistance_to <- function(is_target_fn) {
  is_target_mat <- matrix(is_target_fn(abx_mat), nrow = nrow(abx_mat), ncol = ncol(abx_mat))
  tested    <- rowSums(is_target_mat & !is.na(ast_mat), na.rm = TRUE) > 0
  resistant <- rowSums(is_target_mat & ast_mat == 1L, na.rm = TRUE) > 0
  case_when(resistant ~ 1L, tested ~ 0L, TRUE ~ NA_integer_)
}

# carb_resistance: any carbapenem (matched by name containing "PENEM", so
# it also automatically covers any future carbapenem/BLI-combo added to
# the dropdown, e.g. doripenem, imipenem-relebactam).
lab_full$carb_resistance <- flag_resistance_to(
  function(x) !is.na(x) & str_detect(toupper(x), "PENEM")
)
# methicilin_resistance: oxacillin, methicillin, or cefoxitin (the standard
# phenotypic surrogate panel for mecA-mediated resistance).
lab_full$methicilin_resistance <- flag_resistance_to(
  function(x) !is.na(x) & toupper(x) %in% c("OXACILLIN", "METHICILLIN", "CEFOXITIN")
)
# cephalosporin_resis: any cephalosporin (matched by name starting with
# "CEF", covering all generations present: Cefoxitin, Cefuroxime,
# Cefotaxime, Ceftriaxone, Ceftazidime(-avibactam), Cefepime, Ceftaroline).
lab_full$cephalosporin_resis <- flag_resistance_to(
  function(x) !is.na(x) & str_detect(toupper(x), "^CEF")
)
# 3GC_resistance: THIRD-generation cephalosporins specifically (Ceftriaxone,
# Cefotaxime, Ceftazidime, and the ceftazidime-backbone combination
# ceftazidime-avibactam) -- narrower than cephalosporin_resis, which pools
# every generation (2nd: cefoxitin/cefuroxime; 4th: cefepime; 5th:
# ceftaroline are all EXCLUDED here). No cephalosporin name in this dropdown
# is 3rd-gen without also matching one of these four exact names, so this
# uses an exact-name list rather than a pattern.
# NOTE: variable name starts with a digit, which is not a valid bare R
# symbol -- assigned via [["3GC_resistance"]] rather than $3GC_resistance;
# reference it the same way, e.g. Main_data_CaseReport[["3GC_resistance"]].
lab_full[["3GC_resistance"]] <- flag_resistance_to(
  function(x) !is.na(x) & toupper(x) %in% c("CEFTRIAXONE", "CEFOTAXIME", "CEFTAZIDIME", "CEFTAZIDIME-AVIBACTAM")
)
# vancomycin_resistance: vancomycin specifically (not the wider glycopeptide
# class -- teicoplanin is excluded, per the literal request).
lab_full$vancomycin_resistance <- flag_resistance_to(
  function(x) !is.na(x) & toupper(x) == "VANCOMYCIN"
)

# ---------------------------------------------------------------------------
# AWaRe (Access/Watch/Reserve) classification, WHO 2023 AWaRe book /
# 2023 Model List of Essential Medicines antimicrobial annex. Every one of
# the 50 distinct antimicrobial names appearing anywhere in abx_s1_1:28 /
# abx_s2_1:28 is classified below; 7 are antifungals (Amphotericin B,
# Anidulafungin, Fluconazole, Itraconazole, Micafungin, Posaconazole,
# Voriconazole) and are explicitly OUT OF SCOPE for AWaRe, which classifies
# ANTIBACTERIALS only -- kept as their own category, not silently dropped,
# so a resistance-panel count that includes them is still traceable.
#
# CONFIDENCE FLAGS (agents/authors should verify against the primary WHO
# document before submission -- these are the classification's genuine
# edge cases, not a blanket disclaimer):
#   - Fosfomycin: route-dependent in the WHO list (oral = Access, for
#     uncomplicated cystitis; IV/systemic = Reserve). An AST panel run on a
#     BLOOD isolate is necessarily testing systemic exposure, so classified
#     Reserve here; flagged because the source data does not itself encode
#     route.
#   - Netilmicin, Tobramycin: aminoglycosides, but NOT individually named in
#     the core WHO AWaRe list (only Amikacin and Gentamicin are). Classified
#     Access here by direct pharmacological-class analogy to those two, not
#     from an explicit WHO listing.
#   - Nalidixic acid: a first-generation quinolone predating, and not
#     included in, the current WHO Essential Medicines List. Classified
#     Watch here by analogy to the fluoroquinolone class it precedes
#     (ciprofloxacin/levofloxacin/moxifloxacin), not from an explicit WHO
#     AWaRe listing.
#   - Ampicillin-sulbactam, Aztreonam: WHO AWaRe lists these; included here
#     with ordinary (not flagged) confidence, but noted since some regional
#     adaptations of AWaRe classify beta-lactam/inhibitor combinations
#     differently.
aware_classification <- tribble(
  ~antibiotic, ~aware_category,
  "Amikacin", "Access",
  "Amoxicillin", "Access",
  "Amoxicillin-clavulanic acid", "Access",
  "Ampicillin", "Access",
  "Ampicillin-sulbactam", "Access",
  "Benzylpenicillin", "Access",
  "Chloramphenicol", "Access",
  "Clindamycin", "Access",
  "Gentamicin", "Access",
  "Netilmicin", "Access",
  "Nitrofurantoin", "Access",
  "Oxacillin", "Access",
  "Streptomycin", "Watch",
  "Tetracycline", "Access",
  "Tobramycin", "Access",
  "Trimethoprim-sulfamethoxazole", "Access",
  "Cefepime", "Watch",
  "Cefotaxime", "Watch",
  "Cefoxitin", "Watch",
  "Ceftazidime", "Watch",
  "Ceftriaxone", "Watch",
  "Cefuroxime", "Watch",
  "Ciprofloxacin", "Watch",
  "Erythromycin", "Watch",
  "Ertapenem", "Watch",
  "Fusidic acid", "Watch",
  "Imipenem", "Watch",
  "Levofloxacin", "Watch",
  "Meropenem", "Watch",
  "Moxifloxacin", "Watch",
  "Nalidixic acid", "Watch",
  "Piperacillin-tazobactam", "Watch",
  "Rifampicin", "Watch",
  "Teicoplanin", "Watch",
  "Vancomycin", "Watch",
  "Aztreonam", "Reserve",
  "Ceftaroline", "Reserve",
  "Ceftazidime-avibactam", "Reserve",
  "Colistin", "Reserve",
  "Daptomycin", "Reserve",
  "Fosfomycin", "Reserve",
  "Linezolid", "Reserve",
  "Meropenem-vaborbactam", "Reserve",
  "Tigecycline", "Reserve",
  "Amphotericin B", "Antifungal (not AWaRe)",
  "Anidulafungin", "Antifungal (not AWaRe)",
  "Fluconazole", "Antifungal (not AWaRe)",
  "Itraconazole", "Antifungal (not AWaRe)",
  "Micafungin", "Antifungal (not AWaRe)",
  "Posaconazole", "Antifungal (not AWaRe)",
  "Voriconazole", "Antifungal (not AWaRe)",
)
# Every antimicrobial name actually appearing in the raw AST panel must be
# classified above -- if a future data refresh adds a new drug, this fails
# loudly rather than silently excluding it from the AWaRe breakdown.
abx_names_in_data <- unique(toupper(trimws(c(abx_mat))))
abx_names_in_data <- abx_names_in_data[!is.na(abx_names_in_data)]
stopifnot(all(abx_names_in_data %in% toupper(aware_classification$antibiotic)))
aware_lookup <- setNames(aware_classification$aware_category, toupper(aware_classification$antibiotic))

# access_resistance / watch_resistance / reserve_resistance: same
# flag_resistance_to() construction as carb_resistance/methicilin_resistance
# above -- 1 = resistant to >=1 tested antibiotic in that AWaRe category, 0 =
# that category was tested and never resistant, NA = never tested.
lab_full$access_resistance <- flag_resistance_to(
  function(x) !is.na(x) & aware_lookup[toupper(trimws(x))] == "Access"
)
lab_full$watch_resistance <- flag_resistance_to(
  function(x) !is.na(x) & aware_lookup[toupper(trimws(x))] == "Watch"
)
lab_full$reserve_resistance <- flag_resistance_to(
  function(x) !is.na(x) & aware_lookup[toupper(trimws(x))] == "Reserve"
)

# pdr_resistance: PDR (Magiorakos-style, same category-count proxy as
# mdr/xdr) = resistant in EVERY antimicrobial category actually tested (0
# categories still susceptible). Because this proxy has no organism-specific
# panel, an isolate tested against very few categories can spuriously
# qualify -- n_amr_categories_tested (kept below, not discarded like the
# other intermediates) makes that visible wherever pdr_resistance is
# reported, per Magiorakos et al. (2012)'s own "where feasible" caveat on PDR.
lab_full$pdr_resistance <- case_when(
  n_cat_tested == 0 ~ NA_integer_,
  lab_full$xdr_resistance == 1L & n_cat_susceptible_only == 0 ~ 1L,
  TRUE ~ 0L
)
lab_full$n_amr_categories_tested <- n_cat_tested
lab_full$n_amr_categories_resistant <- n_cat_resistant

## --- pathogen1_clean / pathogen2_clean: collapse nominotypical-subspecies duplication ---
# 4 of the 68 distinct pathogen1 values (and 1 of pathogen2's) repeat the
# species epithet as a trailing "subspecies" word, e.g. "Klebsiella
# pneumoniae pneumoniae" -- standard nomenclature for the nominotypical
# subspecies, taxonomically identical to plain "Klebsiella pneumoniae", not
# a distinct organism. Collapsed here. Genuine subspecies/serovars with a
# DIFFERENT trailing word (e.g. "Klebsiella pneumoniae ozaenae", "Salmonella
# enterica Typhi", "Leuconostoc mesenteroides cremoris") are real taxonomic
# distinctions and are left untouched. No spelling errors were found in the
# raw pathogen1/pathogen2 values themselves (all 68 + 45 distinct names
# checked against standard bacterial/fungal nomenclature).
# NOTE (flagged, not auto-corrected): "Micrococcus kristinae" (2 rows) and
# "Rothia kristinae" (2 rows) are very likely the SAME organism reported
# under two different, non-current genus names by different lab
# technicians (Micrococcus kristinae was reclassified to Kocuria kristinae,
# and separately to Rothia kristinae, over successive taxonomic revisions).
# Left as distinct values -- reconciling requires microbiology judgement
# this script shouldn't assume.
collapse_subspecies_repeat <- function(x) {
  out <- x
  words <- strsplit(x, "\\s+")
  has_dup <- !is.na(x) & vapply(words, function(w) length(w) >= 3 && w[length(w)] == w[length(w) - 1], logical(1))
  out[has_dup] <- vapply(words[has_dup], function(w) paste(w[-length(w)], collapse = " "), character(1))
  out
}
lab_full$pathogen1_clean <- collapse_subspecies_repeat(lab_full$pathogen1)
lab_full$pathogen2_clean <- collapse_subspecies_repeat(lab_full$pathogen2)

## --- pathogen_domain / major_bacterial_group: genus-level classification ---
# Genus = first word of pathogen1_clean. "Other" (the survey's own
# free-text catch-all, not a real taxon) maps to domain NA -- unlike a named
# organism, we genuinely don't know if it's bacterial or fungal, so it must
# NOT default to "Bacteria". This is what lets SECTION 3.6 restrict
# MDR/XDR/PDR (bacteria-only concepts per Magiorakos et al. 2012) away from
# the 3 fungal isolates (Candida albicans x2, Cryptococcus neoformans x1)
# without guessing at the "Other" rows.
genus_group_map <- tribble(
  ~genus,             ~major_bacterial_group,
  "Escherichia",       "Enterobacterales",
  "Klebsiella",        "Enterobacterales",
  "Enterobacter",      "Enterobacterales",
  "Citrobacter",       "Enterobacterales",
  "Proteus",           "Enterobacterales",
  "Serratia",          "Enterobacterales",
  "Salmonella",        "Enterobacterales",
  "Cronobacter",       "Enterobacterales",
  "Pantoea",           "Enterobacterales",
  "Pseudomonas",       "Non-fermenting Gram-negative bacilli",
  "Acinetobacter",     "Non-fermenting Gram-negative bacilli",
  "Stenotrophomonas",  "Non-fermenting Gram-negative bacilli",
  "Burkholderia",      "Non-fermenting Gram-negative bacilli",
  "Moraxella",         "Non-fermenting Gram-negative bacilli",
  "Aeromonas",         "Other Gram-negative",
  "Agrobacterium",     "Other Gram-negative",
  "Ochrobactrum",      "Other Gram-negative",
  "Sphingomonas",      "Other Gram-negative",
  "Pasteurella",       "Other Gram-negative",
  "Staphylococcus",    "Staphylococcus spp.",
  "Streptococcus",     "Streptococcus/Enterococcus spp.",
  "Enterococcus",      "Streptococcus/Enterococcus spp.",
  "Micrococcus",       "Other Gram-positive",
  "Rothia",            "Other Gram-positive",
  "Listeria",          "Other Gram-positive",
  "Aerococcus",        "Other Gram-positive",
  "Gemella",           "Other Gram-positive",
  "Granulicatella",    "Other Gram-positive",
  "Leuconostoc",       "Other Gram-positive",
  "Kocuria",           "Other Gram-positive",
  "Rodentibacter",     "Other Gram-negative",
  "Chromobacterium",   "Other Gram-negative",
  "Malassezia",        "Fungi",
  "Candida",           "Fungi",
  "Cryptococcus",      "Fungi",
  "Other",             "Unspecified/Other"
)
lab_full$genus1 <- sub("^(\\S+).*$", "\\1", lab_full$pathogen1_clean)
lab_full <- lab_full %>%
  left_join(genus_group_map, by = c(genus1 = "genus")) %>%
  mutate(
    pathogen_domain = case_when(
      major_bacterial_group == "Fungi" ~ "Fungi",
      major_bacterial_group == "Unspecified/Other" ~ NA_character_,
      !is.na(major_bacterial_group) ~ "Bacteria",
      TRUE ~ NA_character_
    )
  )


## --- revised-plan overrides: any resistance (I or R), AWaRe, MDR/XDR/PDR (Magiorakos categories; same criteria on tested classes for other species) ---
.m <- match(lab_full$patient_id, syn_flags$patient_id)
lab_full$resistance_any <- syn_flags$any_v2[.m]
lab_full$access_resistance <- syn_flags$access_v2[.m]; lab_full$watch_resistance <- syn_flags$watch_v2[.m]; lab_full$reserve_resistance <- syn_flags$reserve_v2[.m]
lab_full$mdr_resistance <- syn_flags$mdr_v2[.m]; lab_full$xdr_resistance <- syn_flags$xdr_v2[.m]; lab_full$pdr_resistance <- syn_flags$pdr_v2[.m]
lab_full$n_amr_categories_tested <- syn_flags$n_cat_tested[.m]; lab_full$n_amr_categories_resistant <- syn_flags$n_cat_ns[.m]
rm(abx_mat, ast_mat, cat_mat, cat_lookup, n_cat_tested, n_cat_resistant,
   n_cat_susceptible_only, any_tested_amr, any_resistant_amr)


## --- Single merge into Main_data_CaseReport (everything above is lab-side only) ---
Main_data_CaseReport <- Main_data_CaseReport %>%
  left_join(
    lab_full,
    by = "patient_id",
    na_matches = "never"
  ) %>%
  mutate(lab_match = if_else(is.na(lab_match), 0L, lab_match))



# =============================================================================
# SECTION 3: Hospital characteristics (Main_data_hospital)
# =============================================================================
# NOTE: "Hospital data: Total inpatient days of care" and "Hospital data:
# Number of inpatient admissions" from the previous file are renamed in this
# new file to "Total of inpatient days of care: Total inpatient days of
# care" and "Number of inpatient admissions: Number of inpatient
# admissions" -- updated below.
Hospital_characteristics <- Main_data_hospital %>%
  select(
    `Hospital ID`,
    `Level of care`,
    `Total of inpatient days of care: Total inpatient days of care`,
    `Number of inpatient beds: Number of inpatient beds - Total beds`,
    `Hospital workforce during the study period: Number of FTE medical doctors`,	
    `Hospital workforce during the study period: Number of FTE qualified nurses`,
    `Number of inpatient admissions: Number of inpatient admissions`
  ) %>%
  rename(
    Level_of_care = `Level of care`,
    Total_patientDays = `Total of inpatient days of care: Total inpatient days of care`,
    Total_beds= `Number of inpatient beds: Number of inpatient beds - Total beds`,
    NumberFTE_medics= `Hospital workforce during the study period: Number of FTE medical doctors`,	
    NumberFTE_nurses= `Hospital workforce during the study period: Number of FTE qualified nurses`,
    Total_inpatient_admissions = `Number of inpatient admissions: Number of inpatient admissions`
  ) %>%
  # `Hospital ID` == "HF 034" (a stray internal space) fails to match
  # sampling_frame's "HF034" everywhere downstream (SECTION 3.2's
  # design_weight join, SECTION 3.5's facility_rate_base) -- normalized here,
  # once, at source, rather than patched separately in each place that reads
  # Hospital_characteristics. All other 39 IDs are unaffected (no whitespace).
  mutate(`Hospital ID` = gsub("\\s+", "", `Hospital ID`))


# =============================================================================
# SECTION 3.1: ANALYSIS of PREVALENCE
# =============================================================================
# Prevalence of BSI among suspected cases: proportion of enrolled patients
# (indication for blood culture) with a confirmed BSI.

## --- BSI_positive: lab_full merged onto Main_data_CaseReport by patient_id ---
# 1 for every patient_id that matched a confirmed CentralRefLab record in
# lab_full, 0 for everyone else. na.omit() on lab_full$patient_id mirrors the
# na_matches = "never" used everywhere else in this script -- lab_full itself
# has 2 rows with a missing patient_id (from the source data), and without
# guarding against that, any CaseReport row with a missing patient_id would
# spuriously match them (NA %in% c(NA, ...) is TRUE in base R). This is
# numerically identical to lab_match (already merged in SECTION 2); kept as
# its own explicitly-named variable here since that's what was asked for this
# analysis.
Main_data_CaseReport$BSI_positive <- if_else(
  !is.na(Main_data_CaseReport$patient_id) & Main_data_CaseReport$patient_id %in% na.omit(lab_full$patient_id),
  1L, 0L
)

# confirmed_bsi_patient_ids: the base (no-suffix) patient IDs of every
# microbiologically confirmed BSI case. Matched-control records across the
# Followup28/Discharge/CohortEnrolment sheets carry the SAME base ID as the
# case they were matched to, with "/C3" appended (verified directly against
# the raw sheets: e.g. case "HF017/W02/005" pairs with control
# "HF017/W02/005/C3"). A record is therefore only a TRUE matched control if
# its ID both ends in "/C3" AND its base ID identifies a confirmed-BSI case
# here -- not merely "ends in /C3 or /S1" (the previous, looser rule), and
# not "every row in CohortEnrolment" (SECTION 7's previous unconditional
# matched_control=1L). This is used identically at all three matched-control
# construction sites below (SECTIONS 5, 6, 7) so a record's status cannot
# disagree across sheets.
confirmed_bsi_patient_ids <- Main_data_CaseReport$patient_id[Main_data_CaseReport$BSI_positive == 1]
is_true_matched_control <- function(patient_id) {
  if_else(
    str_detect(patient_id, "/C3$") & str_remove(patient_id, "/C3$") %in% confirmed_bsi_patient_ids,
    1L, 0L, missing = 0L
  )
}

# matched_control_broad / is_cohort3_enrolled(): a SECOND, deliberately
# looser matched-control criterion, used ONLY for the IPW/marginal-
# structural-model analyses (SECTION 9 mortality Cox models, all
# resistance-mechanism-specific arms, and the LOS IPW models) -- NOT for
# the exact 1:1 paired-difference/McNemar analysis, which still uses
# is_true_matched_control() above unchanged, since that analysis genuinely
# needs a verified case-control pair.
#
# is_true_matched_control() requires the record's OWN ID to end "/C3" AND
# its base ID to trace to a CONFIRMED-BSI CaseReport patient. Of 643
# CohortEnrolment records, only 464 (72.2%) satisfy that -- verified
# directly against the raw sheets: the other 179 split into 170 whose
# paired case DID enrol but was never confirmed BSI-positive (a genuine,
# expected feature of prospective enrolment -- pairing happens against a
# SUSPECTED case, before culture confirmation), 4 that use a Cyrillic С
# (U+0421) rather than Latin C in the "/C3" suffix (a data-entry/encoding
# slip that silently defeats the regex), 3 that carry a "/S1" suffix
# instead of "/C3", and 2 with no base-ID match in CaseReport at all
# (apparent ID typos, e.g. "HFO25" for "HF025").
#
# For an IPW/MSM model, unlike the 1:1 pair analysis, none of that matters:
# every one of these 643 patients was still genuinely enrolled and
# followed as an uninfected survey participant, regardless of whether
# their nominally-paired case panned out, and IPW is a population-level
# (not pair-level) method.
#
# is_cohort3_enrolled() was FIRST implemented as a literal-ID lookup
# against the CohortEnrolment sheet's own patient_id list -- this
# initially looked like the more robust choice, since it sidesteps the
# Cyrillic-С and "/S1" suffix anomalies entirely. It turned out NOT to be:
# Discharge/FollowUp28/CohortEnrolment are 3 SEPARATE data-entry instances
# for the same physical patient, and a small number of patients have a
# clean "/C3"-suffixed ID in Discharge/FollowUp28 but an anomalous one
# (Cyrillic С, or "/S1") specifically in the CohortEnrolment sheet itself
# -- e.g. patient "HF017/W09/017": Discharge/FollowUp28 carry
# "HF017/W09/017/C3" (correctly strict-matched, matched_control=1), but
# CohortEnrolment itself carries "HF017/W09/017/S1". A literal lookup
# against ONLY the CohortEnrolment sheet's own ID column misses this
# patient's Discharge/FollowUp28 row entirely (that exact string was
# never in CohortEnrolment's own list), silently DROPPING 7 genuinely
# strict-matched (matched_control==1) controls from the broad population
# -- verified directly: this dropped the exact 1:1 paired-analysis
# control pool from 467 to 460 even though is_true_matched_control()
# itself was never touched, because matched_control_broad==0 for these 7
# rows meant they fell out of mortality_group entirely.
#
# Fixed by mirroring is_true_matched_control()'s own design: a per-row
# suffix check on THAT row's own patient_id (no cross-sheet lookup),
# broadened to accept /C3 (Latin), the Cyrillic-С variant, and /S1 --
# exactly the union the codebase's own SECTION 5/6 comment already
# documented as "the case-control/secondary-sample suffixes" before the
# stricter is_true_matched_control() rule was introduced. This raises the
# "Uninfected control" population available to IPW models from 464
# towards 643 (exact count depends on Discharge/FollowUp28 outcome
# linkage, applied downstream, unchanged), while restoring the strict
# 1:1 pair analysis to its correct 467-control pool.
is_cohort3_enrolled <- function(patient_id) {
  as.integer(
    str_detect(patient_id, "/C3$") |
      str_detect(patient_id, "/С3$") |
      str_detect(patient_id, "/S1$")
  )
}

# ---------------------------------------------------------------------------
# CoNS (coagulase-negative staphylococci) reclassification sensitivity flag
# (Terms of Reference SECTION 2.6). Built from the RAW CentralRefLab species
# fields, joined across EVERY id_outcome_conclusion value (MICROBIOLOGICALLY,
# CONTAMINATION, UNCONCLUSIVE, ANAEROBIC) -- deliberately NOT from
# pathogen1_clean/lab_full (built below, in SECTION 2, from
# lab_confirmed_species, which is itself pre-filtered to
# `ID outcome: Conclusion` == "MICROBIOLOGICALLY" only). Building the CoNS
# flag from that pre-filtered object would silently miss every CoNS-species
# record the lab's own protocol-defined criteria had already classified as
# contamination -- exactly the omission to avoid here. Coagulase-negative
# Staphylococcus is operationalised as genus Staphylococcus, species other
# than S. aureus (the standard clinical-microbiology convention; a small
# number of coagulase-POSITIVE non-aureus species, e.g. S. pseudintermedius,
# are folded into this definition for simplicity -- 5 records survey-wide,
# noted as a limitation).
is_CoNS_species <- function(x) {
  !is.na(x) & str_starts(x, "Staphylococcus") & x != "Staphylococcus aureus"
}
cons_raw <- Main_data_hospitalLab %>%
  transmute(
    patient_id = `Unique survey patient ID`,
    id_outcome_conclusion = `ID outcome: Conclusion`,
    species1_clean = collapse_subspecies_repeat(`Species 1: Specify the species`),
    species2_clean = collapse_subspecies_repeat(`Species 2: Specify the species`),
    has_CoNS_species = is_CoNS_species(species1_clean) | is_CoNS_species(species2_clean),
    # contaminated_species1/2_clean: a SEPARATE pair of raw fields, only
    # populated when id_outcome_conclusion == "CONTAMINATION" -- the lab
    # records what actually grew even when it calls the culture a
    # contaminant, rather than leaving species1/2 blank. Used below to build
    # has_CoNS_species_if_contaminated, the basis of the NEW
    # "CoNS-included" sensitivity population (distinct from
    # species1/2_clean, which are NA for contaminated records).
    contaminated_species1_clean = collapse_subspecies_repeat(`Identification: If comensal (contamination) - Species 1:`),
    contaminated_species2_clean = collapse_subspecies_repeat(`Identification: If comensal (contamination) - Species 2:`),
    has_CoNS_species_if_contaminated = is_CoNS_species(contaminated_species1_clean) | is_CoNS_species(contaminated_species2_clean)
  ) %>%
  filter(!is.na(patient_id)) %>%
  distinct(patient_id, .keep_all = TRUE)

# Validation cross-tab: how the lab's own protocol-defined clinical/
# microbiological criteria already split CoNS-species records between
# confirmed and contamination, BEFORE any sensitivity reclassification --
# reported in the Results/Supplementary Material as context for how
# conservative the sensitivity flag below is relative to current practice.
cons_validation_crosstab <- cons_raw %>%
  filter(has_CoNS_species) %>%
  count(id_outcome_conclusion, name = "n") %>%
  mutate(pct = round(100 * n / sum(n), 1))

Main_data_CaseReport <- Main_data_CaseReport %>%
  left_join(
    cons_raw %>% select(patient_id, id_outcome_conclusion, has_CoNS_species, has_CoNS_species_if_contaminated),
    by = "patient_id"
  ) %>%
  mutate(
    has_CoNS_species = if_else(is.na(has_CoNS_species), FALSE, has_CoNS_species),
    has_CoNS_species_if_contaminated = if_else(is.na(has_CoNS_species_if_contaminated), FALSE, has_CoNS_species_if_contaminated)
  )

# BSI_positive_CoNS_sensitive ("CoNS-excluded" sensitivity population): the
# maximally conservative Section-2.6 bound -- every CURRENTLY confirmed BSI
# case with a CoNS species recorded is reclassified as contamination,
# regardless of the clinical/microbiological criteria the lab itself
# applied. Non-CoNS contamination (Streptococcus spp., Leuconostoc,
# Aeromonas, Enterococcus, Sphingomonas, "Other" -- see
# cons_validation_crosstab) is untouched: this flag only ever moves CoNS
# records, never redefines contamination more broadly.
Main_data_CaseReport <- Main_data_CaseReport %>%
  mutate(BSI_positive_CoNS_sensitive = if_else(BSI_positive == 1 & has_CoNS_species, 0L, BSI_positive))
n_CoNS_reclassified <- sum(Main_data_CaseReport$BSI_positive == 1 & Main_data_CaseReport$has_CoNS_species, na.rm = TRUE)

# BSI_positive_CoNS_included ("CoNS-included" sensitivity population,
# distinct from and in the OPPOSITE direction to BSI_positive_CoNS_sensitive
# above): the primary BSI_positive definition (unchanged; the lab's own
# id_outcome_conclusion == "MICROBIOLOGICALLY" criterion) is the STARTING
# point for both sensitivity populations, per instruction that the main
# analysis throughout this report always uses the unchanged primary
# population. This flag additionally counts as confirmed-BSI any record the
# lab classified as CONTAMINATION but for which a CoNS species was
# nonetheless recorded in the separate "if comensal (contamination)" species
# fields -- i.e. it treats CoNS specifically as plausibly a missed true
# pathogen rather than automatically a contaminant, the opposite
# adjustment to BSI_positive_CoNS_sensitive's "always treat CoNS as
# contamination" bound. The two sensitivity populations therefore bracket
# the primary estimate from below (CoNS_sensitive) and above
# (CoNS_included); see Section 2.6 for both, reported alongside the
# unchanged primary throughout every downstream section (prevalence,
# incidence, mortality, risk factors, LOS, etc.), never replacing it.
Main_data_CaseReport <- Main_data_CaseReport %>%
  mutate(
    BSI_positive_CoNS_included = if_else(
      BSI_positive == 0 & !is.na(id_outcome_conclusion) &
        id_outcome_conclusion == "CONTAMINATION" & has_CoNS_species_if_contaminated,
      1L, BSI_positive
    )
  )
n_CoNS_included_from_contamination <- sum(
  Main_data_CaseReport$BSI_positive == 0 &
    !is.na(Main_data_CaseReport$id_outcome_conclusion) &
    Main_data_CaseReport$id_outcome_conclusion == "CONTAMINATION" &
    Main_data_CaseReport$has_CoNS_species_if_contaminated,
  na.rm = TRUE
)

## --- facility_id_clean: DATA QUALITY FIX required before any facility-level analysis ---
# Main_data_CaseReport$facility_id (source: "Unique facility ID") is meant to
# be just a hospital code like "HF069", but ~175 of 9247 rows carry a
# malformed value: a ward suffix appended ("HF069/W06"), a trailing slash
# ("HF032/"), a case typo ("Hf014"), or a truncated/unpadded number ("HF17",
# "HF81"). facility_id_clean extracts the "HF" + digits prefix and re-pads it
# to 3 digits, recovering 24 of those 175 rows. The remaining ~151 are NOT a
# formatting bug: "HF034" (149 rows) and "HF033" (1 row) are facility codes
# that simply do not exist anywhere in Hospital_characteristics's 40
# facilities (a real gap in the source data, not something to paper over),
# and "W12" (1 row) has no recoverable facility code at all. All of these,
# plus the 6 rows with a genuinely blank facility_id, come through as NA and
# are excluded from prevalence_datasetBSI below -- they are NOT silently
# dropped from Main_data_CaseReport itself, only from the facility-level
# aggregation, since there is nowhere in Hospital_characteristics for them to
# attach to.
clean_facility_id <- function(x) {
  x <- toupper(trimws(x))
  has_match <- grepl("^HF[0-9]{1,3}", x)
  digits <- sub("^HF([0-9]{1,3}).*$", "\\1", x)
  out <- rep(NA_character_, length(x))
  out[has_match] <- sprintf("HF%03d", as.integer(digits[has_match]))
  out
}
Main_data_CaseReport$facility_id_clean <- clean_facility_id(Main_data_CaseReport$facility_id)

## --- sampling_frame / NameAU (oblast): read once, here, so region is
## available as a subgroup for EVERY prevalence/incidence analysis below,
## not just SECTION 3.3's calibration. SECTION 3.2 reuses this same object
## for design_weight; SECTION 3.3 reuses it for region_calibration -- read
## only once to avoid three divergent copies.
sampling_frame <- read.csv(
  SAMPLING_CSV,
  stringsAsFactors = FALSE
)
Main_data_CaseReport <- Main_data_CaseReport %>%
  left_join(
    sampling_frame %>% select(UniqueID, NameAU),
    by = c(facility_id_clean = "UniqueID")
  )

## --- prevalence_datasetBSI: one row per facility ---
# Base = Hospital_characteristics (all 40 facilities, even ones with 0
# matched patients); n_tested/n_BSI_positive accumulated from
# Main_data_CaseReport by facility_id_clean and left-joined on.
facility_bsi_summary <- Main_data_CaseReport %>%
  filter(!is.na(facility_id_clean)) %>%
  group_by(facility_id_clean) %>%
  summarise(n_tested = n(), n_BSI_positive = sum(BSI_positive, na.rm = TRUE), .groups = "drop")

prevalence_datasetBSI <- Hospital_characteristics %>%
  left_join(facility_bsi_summary, by = c(`Hospital ID` = "facility_id_clean")) %>%
  mutate(
    Total_inpatient_admissions = as.numeric(Total_inpatient_admissions),
    n_tested = if_else(is.na(n_tested), 0L, n_tested),
    n_BSI_positive = if_else(is.na(n_BSI_positive), 0L, n_BSI_positive)
  )

## --- Methods note (confidence intervals) -----------------------------------
# Proportions (BSI-positive patients / patients tested, overall and by
# subgroup) are reported with the Wilson score interval rather than the
# simple (Wald) normal-approximation interval: the Wald interval can extend
# below 0% or above 100% and performs poorly for proportions near either
# bound or with small denominators, both of which arise in the subgroup
# strata below; the Wilson interval avoids this and is the interval
# recommended in general epidemiological reporting guidance (Newcombe, Stat
# Med 1998) for exactly this reason. The facility-level rate (BSI-positive
# cases per 100 total inpatient admissions) is a count over an exposure
# (admissions), not a simple proportion of tested patients, so it is instead
# reported with the exact Poisson (Garwood) interval, obtained from the
# relationship between the Poisson and chi-squared distributions. All
# intervals are two-sided 95% intervals.
wilson_ci <- function(x, n, conf = 0.95) {
  z <- qnorm(1 - (1 - conf) / 2)
  p_hat <- x / n
  denom <- 1 + z^2 / n
  center <- p_hat + z^2 / (2 * n)
  adj <- z * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2))
  tibble(est = p_hat, lower = (center - adj) / denom, upper = (center + adj) / denom)
}

poisson_rate_ci <- function(x, exposure, conf = 0.95) {
  alpha <- 1 - conf
  lower <- if (x == 0) 0 else qchisq(alpha / 2, 2 * x) / 2
  upper <- qchisq(1 - alpha / 2, 2 * x + 2) / 2
  tibble(est = x / exposure, lower = lower / exposure, upper = upper / exposure)
}

## --- (1) Overall prevalence: two distinct denominators, reported separately ---
# proportion: BSI_positive / all tested patients (the ~9247-patient
# denominator, matching the definition given: "proportion of enrolled
# patients with confirmed BSI"), as a percentage.
n_tested_overall <- sum(!is.na(Main_data_CaseReport$BSI_positive))
n_bsi_overall <- sum(Main_data_CaseReport$BSI_positive, na.rm = TRUE)
prevalence_overall_proportion <- wilson_ci(n_bsi_overall, n_tested_overall) %>%
  mutate(across(everything(), ~ .x * 100)) %>%
  rename(prevalence_pct = est, ci_lower = lower, ci_upper = upper) %>%
  mutate(n_tested = n_tested_overall, n_BSI_positive = n_bsi_overall, .before = 1)

# rate per 100 admissions: uses Total_inpatient_admissions (a facility-level
# count of ALL hospital admissions, not just tested/suspected patients) as
# the denominator -- a genuinely different, much larger denominator than
# n_tested above, so this rate is on a different, much smaller scale. Only
# computable at the facility (or overall, summed across facilities) level,
# since Total_inpatient_admissions isn't itself broken down by patient
# characteristics -- see the note on subgroup prevalence below.
total_admissions_overall <- sum(prevalence_datasetBSI$Total_inpatient_admissions, na.rm = TRUE)
total_bsi_overall_facility <- sum(prevalence_datasetBSI$n_BSI_positive, na.rm = TRUE)
prevalence_overall_rate_per100admissions <- poisson_rate_ci(total_bsi_overall_facility, total_admissions_overall) %>%
  mutate(across(everything(), ~ .x * 100)) %>%
  rename(rate_per_100_admissions = est, ci_lower = lower, ci_upper = upper) %>%
  mutate(n_BSI_positive = total_bsi_overall_facility, total_admissions = total_admissions_overall, .before = 1)

## --- Charlson comorbidity index (adapted) -----------------------------------
# Built from the comorb_* flags available in Main_data_CaseReport, using
# standard Charlson (Charlson et al. 1987) category weights where a matching
# flag exists. Two of the ~17 standard categories have NO corresponding
# source variable and are therefore NOT represented at all: myocardial
# infarction and peripheral vascular disease. This is consequently a
# PARTIAL/adapted index, not a complete validated Charlson score -- treat
# absolute values with that caveat; the relative ranking across patients
# (used for the median split below) is still meaningful. Mutually-exclusive
# severity pairs (diabetes with/without end-organ damage; mild vs
# moderate/severe liver disease; non-metastatic vs metastatic malignancy) are
# resolved to the higher weight only, per standard Charlson coding; missing
# comorbidity flags contribute 0 (absent), consistent with how missing data
# is handled elsewhere in this script (e.g. the Pitt score).
charlson_diabetes_pts <- case_when(
  Main_data_CaseReport$comorb_diabetes_end_organ == 1L ~ 2L,
  Main_data_CaseReport$comorb_diabetes == 1L ~ 1L,
  TRUE ~ 0L
)
charlson_liver_pts <- case_when(
  Main_data_CaseReport$comorb_liver_severe == 1L ~ 3L,
  Main_data_CaseReport$comorb_liver_mild == 1L ~ 1L,
  TRUE ~ 0L
)
charlson_malignancy_pts <- case_when(
  Main_data_CaseReport$comorb_metastatic_tumour == 1L ~ 6L,
  Main_data_CaseReport$comorb_malignancy_active == 1L ~ 2L,
  TRUE ~ 0L
)
charlson_hiv_pts <- if_else(
  Main_data_CaseReport$comorb_hiv_on_art == 1L | Main_data_CaseReport$comorb_hiv_no_art == 1L |
    Main_data_CaseReport$comorb_hiv_cd4_low == 1L,
  6L, 0L, missing = 0L
)
charlson_weight1_vars <- c("comorb_chf", "comorb_cerebrovascular", "comorb_dementia", "comorb_copd",
                           "comorb_connective_tissue", "comorb_peptic_ulcer")
charlson_weight1_pts <- Main_data_CaseReport %>%
  select(all_of(charlson_weight1_vars)) %>%
  mutate(across(everything(), ~ if_else(.x == 1L, 1L, 0L, missing = 0L))) %>%
  rowSums()
charlson_hemi_pts  <- if_else(Main_data_CaseReport$comorb_hemi_paraplegia == 1L, 2L, 0L, missing = 0L)
charlson_renal_pts <- if_else(Main_data_CaseReport$comorb_renal_disease == 1L, 2L, 0L, missing = 0L)

Main_data_CaseReport$charlson_score <- charlson_weight1_pts + charlson_diabetes_pts + charlson_liver_pts +
  charlson_malignancy_pts + charlson_hiv_pts + charlson_hemi_pts + charlson_renal_pts

# charlson_group: a literal "<median / >=median" split is degenerate for
# this population -- the median charlson_score is 0 (73% of patients have no
# captured comorbidity at all), so ">=median" would be everyone and
# "<median" would be empty. Using "0 vs >0" instead, which is what a
# median-0 split collapses to in practice, and is the standard fallback in
# comorbidity literature when the index is this zero-inflated.
charlson_median <- median(Main_data_CaseReport$charlson_score, na.rm = TRUE)
Main_data_CaseReport$charlson_group <- if_else(
  Main_data_CaseReport$charlson_score > 0, "Charlson >0", "Charlson = 0"
)

rm(charlson_diabetes_pts, charlson_liver_pts, charlson_malignancy_pts, charlson_hiv_pts,
   charlson_weight1_pts, charlson_hemi_pts, charlson_renal_pts)

## --- (2) Subgroup prevalence (proportion-based, Wilson CI) ------------------
# Deliberately uses the proportion (BSI_positive / n tested), not the
# per-100-admissions rate, for every subgroup: Total_inpatient_admissions is
# only available as a whole-facility count, not broken down by infection
# origin/age group/sex/ward type/specialty/comorbidity burden, so a
# "per-100-admissions" rate cannot be correctly computed WITHIN these
# subgroups (there is no stratified admissions denominator to divide by).
subgroup_prevalence <- function(df, group_var) {
  df %>%
    filter(!is.na(.data[[group_var]])) %>%
    mutate(group = as.character(.data[[group_var]])) %>%
    group_by(group) %>%
    summarise(n_tested = n(), n_BSI_positive = sum(BSI_positive, na.rm = TRUE), .groups = "drop") %>%
    bind_cols(wilson_ci(.$n_BSI_positive, .$n_tested) %>% mutate(across(everything(), ~ .x * 100))) %>%
    rename(prevalence_pct = est, ci_lower = lower, ci_upper = upper) %>%
    mutate(variable = group_var, .before = 1)
}

# infection_origin_composite: 0 = COMMUNITY, 1 = HOSPITAL -- hierarchical 3-tier
#   composite (97.9% coverage), used here instead of the clinician-only
#   infection_origin_hosp (48.5% coverage); see SECTION 1 for construction
#   and infection_origin_composite_strict (tiers 1-2 only) as a sensitivity
#   check, kept in the data/codebook but not reported here.
# age_group: 0 = NEONATE, 1 = PAEDIATRIC, 3 = ADULT
# sex_at_birth: 0 = MALE, 1 = FEMALE
# ward_type_icu: "ICU" vs "Non-ICU"
# specialty_group: "Medical" vs "Surgical"
# charlson_group: "Charlson = 0" vs "Charlson >0" (see note above)
prevalence_by_subgroup <- bind_rows(
  subgroup_prevalence(Main_data_CaseReport, "infection_origin_composite"),
  subgroup_prevalence(Main_data_CaseReport, "age_group"),
  subgroup_prevalence(Main_data_CaseReport, "sex_at_birth"),
  subgroup_prevalence(Main_data_CaseReport, "ward_type_icu"),
  subgroup_prevalence(Main_data_CaseReport, "specialty_group"),
  subgroup_prevalence(Main_data_CaseReport, "charlson_group"),
  subgroup_prevalence(Main_data_CaseReport, "NameAU")
)



# =============================================================================
# SECTION 3.2: WEIGHTED PREVALENCE ANALYSIS (sampling design)
# =============================================================================
# A second, SEPARATE estimate of everything computed in SECTION 3.1, now
# weighted by the facility sampling design used to select which of Kyrgyzstan's
# 119 hospitals were surveyed at all. See:
#   <data>/KyrgyzstanSampling.R
#   <data>/NewSamplingFrameDec2022.csv
#
# WHAT THE SAMPLING SCRIPT DID: probability-proportional-to-size (PPS)
# sampling WITH replacement. From the 119-facility national frame, it drew
# n = 85 "slots" (`sample(sf$UniqueID, size = 85, replace = TRUE, prob =
# ApCL)`, seed 314), where each facility's chance of being drawn in any one
# slot is proportional to ApCL (its anticipated patient case-load in the
# frame -- a pre-survey size measure, not the Total_inpatient_admissions
# collected during the survey itself). Because it's WITH replacement, a
# large-caseload facility can be drawn more than once; each facility's
# CLr85 (its count of slots) is exactly what PPS sampling is designed to
# produce -- bigger facilities absorb proportionally more of the 85 slots.
# Reproducing that draw here (same seed, same n, same prob) yields exactly
# 40 distinct facilities, which match this survey's 40 facilities exactly
# (originally 39/40 identical at verification time; the 40th, "HF 034" vs
# "HF034", was a stray-space formatting difference in Hospital_characteristics's
# "Hospital ID", not a real discrepancy -- since fixed at source in SECTION 3).
# This is strong confirmation that this IS the design behind the facilities
# actually surveyed.
#
# WHAT A SAMPLING WEIGHT MEANS HERE: PPS-with-replacement is a size-based
# design specifically so that, downstream, an EQUAL-probability (self-
# weighting) sample of patients can be drawn within facilities -- a facility
# selected 7 times was meant to contribute proportionally more to the
# analysis than one selected once. The standard Horvitz-Thompson design
# weight for a PPS-WR unit is the inverse of its per-draw selection
# probability, scaled by the number of draws: w_i = 1 / (n * p_i), where
# p_i = ApCL_i / sum(ApCL) over the FULL 119-facility frame (not just the 40
# selected -- p_i is defined relative to the population that could have been
# drawn). Facilities with a small anticipated case-load (low p_i) get a
# LARGER weight per contributing patient, exactly compensating for their
# lower chance of being included at all; large facilities are already
# "weighted up" by having been drawn -- and hence contributed more patients
# -- more often, so their per-patient weight is correspondingly smaller.
# This produces the second, weighted prevalence estimate below.
#
# HOW THE ESTIMATE ITSELF IS COMPUTED: a naive weighted mean of facility
# proportions would understate uncertainty, because it ignores that BSI
# status is correlated within a facility (patients in the same hospital
# aren't independent draws) and that facilities were sampled with unequal
# probability. Both are handled properly by the survey package: a
# svydesign() with facility_id_clean as the primary sampling unit (cluster)
# and design_weight as the sampling weight, from which svyciprop()/svyby()
# give a design-based point estimate AND a design-based confidence interval
# via Taylor linearisation -- not just a reweighted Wilson interval. This
# implicitly assumes patients enrolled within a selected facility reasonably
# represent that facility (e.g. all/most eligible suspected-BSI cases during
# the study period were enrolled, rather than a further unknown sub-sample)
# -- reasonable for a facility-based surveillance design, but not verifiable
# from the sampling script alone, so treat this as the best available
# design-based estimate rather than a certainty.
#
# COVERAGE: design weights come directly from NewSamplingFrameDec2022.csv,
# not from Hospital_characteristics, so they are NOT affected by the
# "HF 034" formatting issue -- all 149 HF034 patients (excluded from
# prevalence_datasetBSI in SECTION 3.1) ARE included here. Only patients
# with no resolvable facility_id_clean (7 rows: the 6 with a blank
# facility_id plus the 1 unrecoverable "W12") are excluded, since there is
# no facility to assign them a design weight from.

## sampling_frame was already read in SECTION 3.1 (needed there for the
## NameAU/oblast subgroup); reused here, not re-read.
n_draws <- 85
sampling_frame$p_select <- sampling_frame$ApCL / sum(sampling_frame$ApCL)
sampling_frame$design_weight <- 1 / (n_draws * sampling_frame$p_select)

Main_data_CaseReport <- Main_data_CaseReport %>%
  left_join(
    sampling_frame %>% select(UniqueID, design_weight),
    by = c(facility_id_clean = "UniqueID")
  )

## --- survey design: facility = cluster (PSU), design_weight = sampling weight ---
# The PPS design-weight-only prevalence estimate is no longer separately
# reported anywhere in this script (per instruction: unweighted + calibrated
# only); bsi_survey_design/weighted_case_data are retained here purely as
# the base design object the CoNS-sensitive population below is defined
# against before calibrated_weight exists (SECTION 3.3).
weighted_case_data <- Main_data_CaseReport %>% filter(!is.na(design_weight))
bsi_survey_design <- svydesign(
  ids = ~facility_id_clean,
  weights = ~design_weight,
  data = weighted_case_data
)

# -----------------------------------------------------------------------
# SECTION 2.6: CoNS reclassification sensitivity analysis (Terms of
# Reference SECTION 2.6). Re-derives diagnostic yield, population
# prevalence, and resistance prevalence under BSI_positive_CoNS_sensitive
# (built above) instead of the primary BSI_positive, i.e. treating every
# CoNS-species confirmed case as contamination. weighted_case_data/
# bsi_survey_design already carry BSI_positive_CoNS_sensitive because it was
# added to Main_data_CaseReport before this design object was built.
# -----------------------------------------------------------------------
cons_sensitivity_unweighted <- Main_data_CaseReport %>%
  summarise(
    n_tested = n(),
    n_BSI_positive = sum(BSI_positive, na.rm = TRUE),
    n_BSI_positive_CoNS_sensitive = sum(BSI_positive_CoNS_sensitive, na.rm = TRUE),
    n_BSI_positive_CoNS_included = sum(BSI_positive_CoNS_included, na.rm = TRUE),
    n_reclassified_to_contamination = sum(BSI_positive == 1 & has_CoNS_species, na.rm = TRUE),
    n_reclassified_from_contamination = n_CoNS_included_from_contamination
  ) %>%
  mutate(
    pct_primary = round(100 * n_BSI_positive / n_tested, 2),
    pct_CoNS_sensitive = round(100 * n_BSI_positive_CoNS_sensitive / n_tested, 2),
    pct_CoNS_included = round(100 * n_BSI_positive_CoNS_included / n_tested, 2)
  )

# Two sensitivity populations, bracketing the unchanged primary estimate
# from below (CoNS_sensitive: all CoNS treated as contamination) and above
# (CoNS_included: contaminated-but-CoNS-species cultures counted as
# confirmed) -- see the is_CoNS_species()/BSI_positive_CoNS_included
# comment (above) for the full rationale. The primary analysis throughout
# this report always uses the unchanged BSI_positive; these two rows are
# reported as sensitivity checks, never as a replacement.
cons_prevalence_comparison <- tibble(
  definition = c(
    "Primary (protocol-defined CoNS classification)",
    "CoNS-excluded sensitivity (all confirmed CoNS reclassified as contamination)",
    "CoNS-included sensitivity (contaminated cultures with a CoNS species counted as confirmed BSI)"
  ),
  n_confirmed_BSI = c(
    cons_sensitivity_unweighted$n_BSI_positive,
    cons_sensitivity_unweighted$n_BSI_positive_CoNS_sensitive,
    cons_sensitivity_unweighted$n_BSI_positive_CoNS_included
  ),
  unweighted_pct = c(
    cons_sensitivity_unweighted$pct_primary,
    cons_sensitivity_unweighted$pct_CoNS_sensitive,
    cons_sensitivity_unweighted$pct_CoNS_included
  )
)
# NOTE: the resistance-among-sensitive-population rows (which need
# resistant_BSI_positive, not built until later in SECTION 3) and the
# calibrated-weight columns (calibrated_weight not built until SECTION 3.3,
# after this point) are added to this same workbook further down the
# script -- see the write_xlsx call immediately after resistant_BSI_positive
# is constructed.

## --- (2) Weighted subgroup prevalence, mirroring SECTION 3.1's subgroups ---
subgroup_prevalence_weighted <- function(design, group_var) {
  sub_design <- subset(design, !is.na(get(group_var)))
  res <- svyby(
    as.formula("~I(BSI_positive == 1)"),
    as.formula(paste0("~", group_var)),
    sub_design, svyciprop, method = "logit", vartype = "ci"
  )
  tibble(
    variable = group_var,
    group = as.character(res[[group_var]]),
    prevalence_pct = res[["I(BSI_positive == 1)"]] * 100,
    ci_lower = res[["ci_l"]] * 100,
    ci_upper = res[["ci_u"]] * 100
  )
}

rm(weighted_case_data)



# =============================================================================
# SECTION 3.3: CALIBRATED WEIGHTED PREVALENCE (regional population)
# =============================================================================
# A third estimate. SECTION 3.2's design_weight corrects for each facility's
# probability of selection under the PPS design -- but that is only correct
# ON AVERAGE across repeated random draws; this ONE realised draw of 40
# facilities still leaves real regional imbalance (checked directly below).
# This section CALIBRATES design_weight against actual 2024 regional
# population counts from an external, authoritative source: the Demographic
# Yearbook of the Kyrgyz Republic (National Statistical Committee, 2020-2024
# edition, released 2025-10-03), Table I.1.7.1 ("Resident population by
# oblast, at the start of the year"), https://stat.gov.kg. Values used
# (thousands, 2024): Batken 583.4, Bishkek (city) 1294.8, Chui 957.6,
# Issyk-Kul 544.4, Jalal-Abad 1335.8, Naryn 312.1, Osh oblast 1393.6 + Osh
# (city) 463.1 = 1856.7 (combined, since the sampling frame's region field
# does not distinguish Osh oblast from Osh city), Talas 277.1. These 8
# figures sum to exactly 7,161.9 thousand, the frame's published national
# total -- confirming the extraction is internally consistent.
#
# WHY THIS MATTERS HERE (verified, not assumed): comparing design_weight's
# regional shares against true population shares shows real, substantial
# imbalance in this particular draw -- Issyk-Kul (7.6% of the population)
# contributes only 0.7% of weighted representation (it drew just 1 of 40
# facilities); Naryn (4.4% of the population) contributes only 1.2%. Talas
# is correspondingly over-represented (10.5% weighted vs 3.9% true). This is
# exactly the kind of single-draw sampling variability that a valid design
# does not, by itself, guarantee away.
#
# METHOD: single-margin calibration (a direct generalisation of
# post-stratification -- with only one calibration variable, region, this
# has a closed-form solution and does not require iterative proportional
# fitting/raking). For each region r, a raking factor is computed as
# (true population share of r) / (current design-weighted share of r);
# every patient's design_weight is then multiplied by their region's raking
# factor. This exactly forces the weighted regional shares to match the true
# population shares (verified below), while leaving each patient's weight
# WITHIN their own facility/region untouched relative to their facility-mates.

pop2024_by_region <- tribble(
  ~NameAU, ~pop_2024_thousands,
  "Batken", 583.4,
  "Bishkek", 1294.8,
  "Chui", 957.6,
  "Issyk-Kul", 544.4,
  "Jalal-Abad", 1335.8,
  "Naryn", 312.1,
  "Osh", 1393.6 + 463.1,
  "Talas", 277.1
)

## NameAU was already joined in SECTION 3.1.
region_calibration <- Main_data_CaseReport %>%
  filter(!is.na(design_weight)) %>%
  group_by(NameAU) %>%
  summarise(weighted_n = sum(design_weight), .groups = "drop") %>%
  left_join(pop2024_by_region, by = "NameAU") %>%
  mutate(
    weighted_share = weighted_n / sum(weighted_n),
    true_pop_share = pop_2024_thousands / sum(pop2024_by_region$pop_2024_thousands),
    raking_factor = true_pop_share / weighted_share
  )

Main_data_CaseReport <- Main_data_CaseReport %>%
  left_join(region_calibration %>% select(NameAU, raking_factor), by = "NameAU") %>%
  mutate(calibrated_weight = design_weight * raking_factor)

## --- survey design using the calibrated weight ---
weighted_case_data_calib <- Main_data_CaseReport %>% filter(!is.na(calibrated_weight))
bsi_survey_design_calib <- svydesign(
  ids = ~facility_id_clean,
  weights = ~calibrated_weight,
  data = weighted_case_data_calib
)

## --- (1) Overall calibrated prevalence, for comparison against SECTION 3.1/3.2 ---
prevalence_overall_calibrated <- svyciprop(~I(BSI_positive == 1), bsi_survey_design_calib, method = "logit")
prevalence_overall_calibrated_ci <- confint(prevalence_overall_calibrated)

## --- (2) Calibrated subgroup prevalence, mirroring SECTION 3.1/3.2's subgroups ---
prevalence_by_subgroup_calibrated <- bind_rows(
  subgroup_prevalence_weighted(bsi_survey_design_calib, "infection_origin_composite"),
  subgroup_prevalence_weighted(bsi_survey_design_calib, "age_group"),
  subgroup_prevalence_weighted(bsi_survey_design_calib, "sex_at_birth"),
  subgroup_prevalence_weighted(bsi_survey_design_calib, "ward_type_icu"),
  subgroup_prevalence_weighted(bsi_survey_design_calib, "specialty_group"),
  subgroup_prevalence_weighted(bsi_survey_design_calib, "charlson_group"),
  subgroup_prevalence_weighted(bsi_survey_design_calib, "NameAU")
)

rm(weighted_case_data_calib)



# =============================================================================
# SECTION 3.4: COMBINED COMPARISON TABLE (unweighted / calibrated)
# =============================================================================
# Brings SECTION 3.1 (unweighted) and 3.3 (calibrated) into one table, overall
# and by subgroup, each estimate with its own 95% CI. Per explicit
# instruction, the intermediate PPS design-weighted-only estimate (SECTION
# 3.2, uncalibrated) is no longer reported in any output table throughout
# this script -- design_weight itself is still computed and still used
# (calibrated_weight = design_weight * raking_factor, SECTION 3.3), but its
# own uncalibrated svyciprop()/svyratio() estimate is not a separate reported
# column anywhere from this point on. All figures are proportions AMONG
# TESTED PATIENTS (BSI-positive / patients with a blood culture drawn) -- not
# a rate per population, and not the per-100-admissions rate from SECTION
# 3.1(1), which cannot be computed within subgroups (see that section's
# note). n_tested/n_BSI_positive are the raw (unweighted) counts in every
# row, shown once for reference -- calibrated Ns are not simple counts and
# are not reported here; the survey::svyby() design effectively re-weights
# the same tested patients, it does not test additional ones.
#
# Hospital- vs community-acquired uses infection_origin_composite (the
# hierarchical 3-tier variable, 97.9% coverage -- see SECTION 1), not the
# clinician-only infection_origin_hosp (48.5% coverage). infection_origin_hosp
# and the conservative infection_origin_composite_strict (tiers 1-2 only, no
# proxy inference) remain in the dataset/codebook for sensitivity checks;
# they are not repeated in this summary table.

combine_estimates <- function(group_var) {
  unw <- subgroup_prevalence(Main_data_CaseReport, group_var) %>%
    select(variable, group, n_tested, n_BSI_positive,
           unweighted_pct = prevalence_pct, unweighted_lower = ci_lower, unweighted_upper = ci_upper)
  cal <- subgroup_prevalence_weighted(bsi_survey_design_calib, group_var) %>%
    select(group, calibrated_pct = prevalence_pct, calibrated_lower = ci_lower, calibrated_upper = ci_upper)
  unw %>% left_join(cal, by = "group")
}

overall_row <- tibble(
  variable = "Overall", group = "All tested patients",
  n_tested = n_tested_overall, n_BSI_positive = n_bsi_overall,
  unweighted_pct = prevalence_overall_proportion$prevalence_pct,
  unweighted_lower = prevalence_overall_proportion$ci_lower,
  unweighted_upper = prevalence_overall_proportion$ci_upper,
  calibrated_pct = as.numeric(prevalence_overall_calibrated) * 100,
  calibrated_lower = prevalence_overall_calibrated_ci[1] * 100,
  calibrated_upper = prevalence_overall_calibrated_ci[2] * 100
)

prevalence_combined_table <- bind_rows(
  overall_row,
  combine_estimates("infection_origin_composite"),
  combine_estimates("age_group"),
  combine_estimates("sex_at_birth"),
  combine_estimates("ward_type_icu"),
  combine_estimates("specialty_group"),
  combine_estimates("charlson_group"),
  combine_estimates("NameAU")
) %>%
  mutate(across(where(is.numeric) & !c(n_tested, n_BSI_positive), ~ round(.x, 2)))

write_xlsx(prevalence_combined_table, "Prevalence_combined_table.xlsx")

rm(overall_row)



# =============================================================================
# SECTION 3.5: COMBINED COMPARISON TABLE -- RATE PER 100 ADMISSIONS
# =============================================================================
# Same three estimates (unweighted / design-weighted / calibrated) as SECTION
# 3.4, but on the per-100-admissions scale from SECTION 3.1(1) instead of the
# among-tested-patients proportion scale.
#
# WHY THIS IS DIFFERENT FROM SECTION 3.4, NOT JUST A UNIT CHANGE:
# Total_inpatient_admissions (the denominator) is only available as a whole-
# facility count -- it is NOT broken down by infection origin, age group,
# sex, ward type, specialty, or comorbidity burden (SECTION 3.1's note on
# this). There is therefore no subgroup-specific admissions denominator to
# divide by. What IS reported below for each subgroup LEVEL is that level's
# BSI cases divided by the SAME whole-system admissions total used for the
# overall row -- i.e. a cause-specific rate, exactly analogous to a cause-
# specific mortality rate per population: "of every 100 admissions to a
# surveyed facility, X resulted in a hospital-acquired BSI, Y in a
# community-acquired BSI". Levels of the SAME variable sum (approximately,
# up to missingness in that variable -- see SECTION 3.1) to the overall
# per-100-admissions rate. This is NOT a subgroup-specific risk (it does NOT
# mean "X surgical admissions out of 100 surgical admissions get a BSI") --
# for that, SECTION 3.4's among-tested-patients proportion is the correct
# quantity.
#
# ESTIMATION: unweighted uses the exact Poisson (Garwood) interval, as in
# SECTION 3.1. Weighted/calibrated use survey::svyratio() -- a ratio-of-
# weighted-totals estimator with a Taylor-linearised CI -- on a FACILITY-
# level design (one row per facility, id = Hospital ID, weight =
# design_weight/calibrated_weight from SECTIONS 3.2/3.3), which is the
# correct unit for this: Total_inpatient_admissions and n_BSI_positive are
# both facility-level counts, not patient-level.

facility_weights <- sampling_frame %>%
  select(UniqueID, NameAU, design_weight) %>%
  left_join(region_calibration %>% select(NameAU, raking_factor), by = "NameAU") %>%
  mutate(calibrated_weight = design_weight * raking_factor)

# All 40 surveyed facilities are included here (the "HF 034" vs "HF034"
# stray-space mismatch that previously dropped 1 facility is fixed at source
# in SECTION 3, so Hospital_characteristics's "Hospital ID" now matches
# sampling_frame's UniqueID for all 40). The filter below is a safety net,
# not an expected exclusion.
facility_rate_base <- prevalence_datasetBSI %>%
  select(`Hospital ID`, Total_inpatient_admissions, n_tested, n_BSI_positive) %>%
  left_join(facility_weights %>% select(UniqueID, design_weight, calibrated_weight),
            by = c(`Hospital ID` = "UniqueID")) %>%
  filter(!is.na(design_weight), !is.na(Total_inpatient_admissions))

facility_group_bsi <- function(group_var) {
  Main_data_CaseReport %>%
    filter(!is.na(facility_id_clean), !is.na(.data[[group_var]])) %>%
    mutate(group = as.character(.data[[group_var]])) %>%
    group_by(facility_id_clean, group) %>%
    summarise(n_bsi_group = sum(BSI_positive, na.rm = TRUE), .groups = "drop")
}

## --- (1) Overall rate per 100 admissions: unweighted (Poisson CI, from SECTION 3.1) + calibrated (svyratio) ---
## Per instruction, the intermediate PPS design-weighted-only estimate is no
## longer a reported column anywhere in this script -- see SECTION 3.4's note.
facility_design_c <- svydesign(ids = ~`Hospital ID`, weights = ~calibrated_weight, data = facility_rate_base)

ratio_overall_c <- svyratio(~n_BSI_positive, ~Total_inpatient_admissions, facility_design_c)
ci_overall_c <- confint(ratio_overall_c) * 100
unw_overall_rate <- poisson_rate_ci(sum(facility_rate_base$n_BSI_positive), sum(facility_rate_base$Total_inpatient_admissions)) %>%
  mutate(across(everything(), ~ .x * 100))

overall_rate_row <- tibble(
  variable = "Overall", group = "All admissions",
  n_BSI_positive = sum(facility_rate_base$n_BSI_positive), total_admissions = sum(facility_rate_base$Total_inpatient_admissions),
  unweighted_rate_per100 = unw_overall_rate$est,
  unweighted_lower = unw_overall_rate$lower,
  unweighted_upper = unw_overall_rate$upper,
  calibrated_rate_per100 = as.numeric(coef(ratio_overall_c)) * 100,
  calibrated_lower = ci_overall_c[1], calibrated_upper = ci_overall_c[2]
)

## --- (2) Subgroup rate per 100 admissions (cause-specific, shared denominator) ---
rate_per100_by_group <- function(group_var) {
  fg <- facility_group_bsi(group_var)
  levels_g <- sort(unique(fg$group))
  purrr::map_dfr(levels_g, function(lvl) {
    fac_df <- facility_rate_base %>%
      left_join(
        fg %>% filter(group == lvl) %>% select(`Hospital ID` = facility_id_clean, n_bsi_group),
        by = "Hospital ID"
      ) %>%
      mutate(n_bsi_group = if_else(is.na(n_bsi_group), 0L, n_bsi_group))

    unw <- poisson_rate_ci(sum(fac_df$n_bsi_group), sum(fac_df$Total_inpatient_admissions)) %>%
      mutate(across(everything(), ~ .x * 100))

    des_c <- svydesign(ids = ~`Hospital ID`, weights = ~calibrated_weight, data = fac_df)
    ratio_c <- svyratio(~n_bsi_group, ~Total_inpatient_admissions, des_c)
    ci_c <- confint(ratio_c) * 100

    tibble(
      variable = group_var, group = lvl,
      n_BSI_positive = sum(fac_df$n_bsi_group), total_admissions = sum(fac_df$Total_inpatient_admissions),
      unweighted_rate_per100 = unw$est, unweighted_lower = unw$lower, unweighted_upper = unw$upper,
      calibrated_rate_per100 = as.numeric(coef(ratio_c)) * 100, calibrated_lower = ci_c[1], calibrated_upper = ci_c[2]
    )
  })
}

rate_per100_combined_table <- bind_rows(
  overall_rate_row,
  rate_per100_by_group("infection_origin_composite"),
  rate_per100_by_group("age_group"),
  rate_per100_by_group("sex_at_birth"),
  rate_per100_by_group("ward_type_icu"),
  rate_per100_by_group("specialty_group"),
  rate_per100_by_group("charlson_group"),
  rate_per100_by_group("NameAU")
) %>%
  mutate(across(where(is.numeric) & !c(n_BSI_positive, total_admissions), ~ round(.x, 3)))

write_xlsx(rate_per100_combined_table, "Rate_per100admissions_combined_table.xlsx")

rm(overall_rate_row, facility_design_c, ratio_overall_c, ci_overall_c, unw_overall_rate)



# =============================================================================
# SECTION 3.6: ANTIMICROBIAL RESISTANCE (AMR) ANALYSIS
# =============================================================================
# Three parts: (3.6.1) population prevalence of drug-resistant BSI, with
# bacterial isolates further split into MDR/XDR/PDR; (3.6.2) pathogen
# distribution among confirmed BSI cases; (3.6.3) pathogen x antibiotic
# resistance heatmap for the top 20 most-tested antibiotics. All reuse
# resistance_any/mdr_resistance/xdr_resistance/pdr_resistance and
# pathogen1_clean/pathogen_domain/major_bacterial_group from SECTION 2, and
# BSI_positive/facility_id_clean/design_weight/calibrated_weight from
# SECTIONS 3.1-3.3.

# -----------------------------------------------------------------------
# 3.6.1: Population prevalence of drug-resistant BSI (any / MDR / XDR / PDR)
# -----------------------------------------------------------------------
# Same "population" framing as SECTION 3.4: denominator = all tested
# patients (BSI-negative patients contribute a structural 0, since a
# patient without BSI cannot have a drug-resistant BSI). MDR/XDR/PDR are
# Magiorakos et al. (2012) bacteria-specific concepts, so those three rows
# further restrict to patients whose confirmed isolate is bacterial
# (pathogen_domain == "Bacteria"; the 3 fungal isolates and the "Other"
# free-text catch-all -- domain genuinely unknown -- get NA, not 0, for
# these three rows only). NA also applies wherever AST data was too
# incomplete to evaluate resistance_any/mdr_resistance/etc. at all.
build_population_amr_flag <- function(resistance_col, restrict_bacterial = FALSE, bsi_vec = NULL) {
  bsi <- if (is.null(bsi_vec)) Main_data_CaseReport$BSI_positive else bsi_vec
  domain <- Main_data_CaseReport$pathogen_domain
  resistance_val <- Main_data_CaseReport[[resistance_col]]
  not_bacterial <- is.na(domain) | domain != "Bacteria"
  case_when(
    bsi == 0 ~ 0L,
    restrict_bacterial & bsi == 1 & not_bacterial ~ NA_integer_,
    bsi == 1 & resistance_val == 1 ~ 1L,
    bsi == 1 & resistance_val == 0 ~ 0L,
    bsi == 1 & is.na(resistance_val) ~ NA_integer_,
    TRUE ~ NA_integer_
  )
}
Main_data_CaseReport$resistant_BSI_positive <- build_population_amr_flag("resistance_any", restrict_bacterial = FALSE)

# resistant_BSI_positive_CoNS_included: the SAME resistance_any lookup, but
# keyed off BSI_positive_CoNS_included rather than the primary BSI_positive
# -- required because build_population_amr_flag()'s bsi==0~0L branch would
# otherwise force EVERY one of the 122 newly-included contaminated-with-CoNS
# records to "confirmed susceptible" (0), when their true resistance status
# is genuinely untested (their resistance_any field is NA, since no AST
# panel was run under the CONTAMINATION workflow) -- an artefact of keying
# on the wrong bsi flag, not a real finding. Used only for the CoNS-included
# sensitivity check below; the primary resistant_BSI_positive (and every
# other analysis in this script) is unaffected and still keys off the
# unchanged primary BSI_positive.
Main_data_CaseReport$resistant_BSI_positive_CoNS_included <- build_population_amr_flag(
  "resistance_any", restrict_bacterial = FALSE, bsi_vec = Main_data_CaseReport$BSI_positive_CoNS_included
)

# Same construction, every other mechanism -- built here (not just for "any
# resistance") so the full Risk-Factors mechanism-by-mechanism table
# (SECTION 10) can be repeated in full under the CoNS-included population,
# per instruction. Each mirrors its primary-population counterpart in
# SECTION 3.6.1 exactly, substituting BSI_positive_CoNS_included for
# BSI_positive as the bsi_vec.
Main_data_CaseReport$mdr_BSI_bacterial_CoNS_included <- build_population_amr_flag(
  "mdr_resistance", restrict_bacterial = TRUE, bsi_vec = Main_data_CaseReport$BSI_positive_CoNS_included)
Main_data_CaseReport$xdr_BSI_bacterial_CoNS_included <- build_population_amr_flag(
  "xdr_resistance", restrict_bacterial = TRUE, bsi_vec = Main_data_CaseReport$BSI_positive_CoNS_included)
Main_data_CaseReport$carb_BSI_positive_CoNS_included <- build_population_amr_flag(
  "carb_resistance", restrict_bacterial = FALSE, bsi_vec = Main_data_CaseReport$BSI_positive_CoNS_included)
Main_data_CaseReport$methicillin_BSI_positive_CoNS_included <- build_population_amr_flag(
  "methicilin_resistance", restrict_bacterial = FALSE, bsi_vec = Main_data_CaseReport$BSI_positive_CoNS_included)
Main_data_CaseReport$cephalosporin_BSI_positive_CoNS_included <- build_population_amr_flag(
  "cephalosporin_resis", restrict_bacterial = FALSE, bsi_vec = Main_data_CaseReport$BSI_positive_CoNS_included)
Main_data_CaseReport$c3GC_BSI_positive_CoNS_included <- build_population_amr_flag(
  "3GC_resistance", restrict_bacterial = FALSE, bsi_vec = Main_data_CaseReport$BSI_positive_CoNS_included)
Main_data_CaseReport$access_BSI_positive_CoNS_included <- build_population_amr_flag(
  "access_resistance", restrict_bacterial = FALSE, bsi_vec = Main_data_CaseReport$BSI_positive_CoNS_included)
Main_data_CaseReport$watch_BSI_positive_CoNS_included <- build_population_amr_flag(
  "watch_resistance", restrict_bacterial = FALSE, bsi_vec = Main_data_CaseReport$BSI_positive_CoNS_included)
Main_data_CaseReport$reserve_BSI_positive_CoNS_included <- build_population_amr_flag(
  "reserve_resistance", restrict_bacterial = FALSE, bsi_vec = Main_data_CaseReport$BSI_positive_CoNS_included)

# Completes the SECTION 2.6 CoNS sensitivity workbook (deferred from just
# after cons_prevalence_comparison was built, above, because this row needs
# resistant_BSI_positive, which does not exist until this point in the
# script). Resistance prevalence among confirmed BSI cases, recomputed on
# the CoNS-sensitive population (denominator restricted to cases retained
# as confirmed BSI after CoNS reclassification, with a determinate
# resistance result).
cons_resistance_sensitivity <- Main_data_CaseReport %>%
  filter(BSI_positive_CoNS_sensitive == 1, !is.na(resistant_BSI_positive)) %>%
  summarise(
    n = n(),
    n_resistant = sum(resistant_BSI_positive == 1, na.rm = TRUE),
    pct_resistant = round(100 * mean(resistant_BSI_positive == 1, na.rm = TRUE), 1)
  )
# Mirror of the row above, but for the CoNS-INCLUDED sensitivity population
# (contaminated-with-CoNS cultures counted as confirmed) -- uses
# resistant_BSI_positive_CoNS_included (above), NOT resistant_BSI_positive,
# specifically so the 122 newly-included records correctly fall out as NA
# (genuinely untested) rather than being forced to "confirmed susceptible"
# by the primary flag's own bsi==0 branch. This denominator is necessarily
# smaller than n_confirmed_BSI for this population -- reported as such, not
# imputed.
cons_included_resistance_sensitivity <- tibble(
    bsi_included = Main_data_CaseReport$BSI_positive_CoNS_included,
    res_included = Main_data_CaseReport$resistant_BSI_positive_CoNS_included
  ) %>%
  filter(bsi_included == 1, !is.na(res_included)) %>%
  summarise(
    n = n(),
    n_resistant = sum(res_included == 1, na.rm = TRUE),
    pct_resistant = round(100 * mean(res_included == 1, na.rm = TRUE), 1)
  )

# Calibrated-weight column for cons_prevalence_comparison (SECTION 3.3's
# calibrated_weight/bsi_survey_design_calib now exist, unlike at the point
# cons_prevalence_comparison was first built in SECTION 2.6).
cons_prevalence_calib <- svyciprop(~I(BSI_positive_CoNS_sensitive == 1), bsi_survey_design_calib, method = "logit")
cons_prevalence_calib_ci <- confint(cons_prevalence_calib)
cons_included_data_calib <- Main_data_CaseReport %>% filter(!is.na(calibrated_weight))
cons_included_design_calib <- svydesign(ids = ~facility_id_clean, weights = ~calibrated_weight, data = cons_included_data_calib)
cons_prevalence_included_calib <- svyciprop(~I(BSI_positive_CoNS_included == 1), cons_included_design_calib, method = "logit")
cons_prevalence_included_calib_ci <- confint(cons_prevalence_included_calib)
cons_prevalence_comparison <- cons_prevalence_comparison %>%
  mutate(
    calibrated_pct = c(
      round(prevalence_overall_calibrated * 100, 2),
      round(cons_prevalence_calib[1] * 100, 2),
      round(cons_prevalence_included_calib[1] * 100, 2)
    ),
    calibrated_ci_lower = c(
      round(prevalence_overall_calibrated_ci[1] * 100, 2),
      round(cons_prevalence_calib_ci[1] * 100, 2),
      round(cons_prevalence_included_calib_ci[1] * 100, 2)
    ),
    calibrated_ci_upper = c(
      round(prevalence_overall_calibrated_ci[2] * 100, 2),
      round(cons_prevalence_calib_ci[2] * 100, 2),
      round(cons_prevalence_included_calib_ci[2] * 100, 2)
    )
  )

write_xlsx(
  list(
    validation_crosstab_lab_protocol = cons_validation_crosstab,
    prevalence_comparison = cons_prevalence_comparison,
    resistance_among_CoNS_excluded_confirmed_BSI = cons_resistance_sensitivity,
    resistance_among_CoNS_included_confirmed_BSI = cons_included_resistance_sensitivity,
    n_records_reclassified = tibble(
      n_CoNS_confirmed_reclassified_to_contamination = cons_sensitivity_unweighted$n_reclassified_to_contamination,
      n_contaminated_CoNS_reclassified_to_confirmed = n_CoNS_included_from_contamination
    )
  ),
  "CoNS_reclassification_sensitivity_analysis.xlsx"
)

Main_data_CaseReport$mdr_BSI_bacterial <- build_population_amr_flag("mdr_resistance", restrict_bacterial = TRUE)
Main_data_CaseReport$xdr_BSI_bacterial <- build_population_amr_flag("xdr_resistance", restrict_bacterial = TRUE)
Main_data_CaseReport$pdr_BSI_bacterial <- build_population_amr_flag("pdr_resistance", restrict_bacterial = TRUE)

# Class-specific population flags, same construction as resistant_BSI_positive
# (NOT bacteria-restricted: unlike mdr/xdr/pdr's organism-agnostic category-
# count proxy, these five are already narrow, named-antibiotic-class flags
# from SECTION 2 -- restricting further would just duplicate what "was this
# class ever tested" already does).
Main_data_CaseReport$carb_BSI_positive <- build_population_amr_flag("carb_resistance", restrict_bacterial = FALSE)
Main_data_CaseReport$methicillin_BSI_positive <- build_population_amr_flag("methicilin_resistance", restrict_bacterial = FALSE)
Main_data_CaseReport$cephalosporin_BSI_positive <- build_population_amr_flag("cephalosporin_resis", restrict_bacterial = FALSE)
Main_data_CaseReport$c3GC_BSI_positive <- build_population_amr_flag("3GC_resistance", restrict_bacterial = FALSE)
Main_data_CaseReport$vancomycin_BSI_positive <- build_population_amr_flag("vancomycin_resistance", restrict_bacterial = FALSE)

# AWaRe category population flags, same construction and same "not
# bacteria-restricted" reasoning as the class-specific flags directly above
# (each AWaRe category is itself already a named-antibiotic-class flag from
# SECTION 2, so restricting further would just duplicate what "was this
# category ever tested" already does).
Main_data_CaseReport$access_BSI_positive <- build_population_amr_flag("access_resistance", restrict_bacterial = FALSE)
Main_data_CaseReport$watch_BSI_positive <- build_population_amr_flag("watch_resistance", restrict_bacterial = FALSE)
Main_data_CaseReport$reserve_BSI_positive <- build_population_amr_flag("reserve_resistance", restrict_bacterial = FALSE)

## --- Organism-restricted class-resistance: WHO GLASS priority pathogen-drug pairs ---
# carb_resistance/methicilin_resistance/3GC_resistance/vancomycin_resistance
# (and the *_BSI_positive flags above) pool BOTH isolates (S1+S2 AST slots)
# per patient, which is right for a patient-level "any resistant isolate"
# flag but wrong for an organism-specific rate: a patient's carbapenem
# result could come from their pathogen2 isolate, which may not even be one
# of the organisms this metric is meant to describe. The four rows below
# instead use ONLY the S1 (Species 1 / primary isolate) AST slots, paired
# with pathogen1_clean -- the same correctness rule as SECTION 3.6.3's
# heatmap -- and restrict to the WHO GLASS priority organism group each
# drug class is actually surveilled against: carbapenems and 3rd-generation
# cephalosporins in Enterobacterales/P. aeruginosa/A. baumannii; methicillin
# (oxacillin/methicillin/cefoxitin surrogate) in Staphylococcus; vancomycin
# in Enterococcus. Denominator = susceptible + resistant isolates WITHIN
# that organism group (i.e. that group's isolates where the class was
# actually tested) -- NOT all tested patients like the unrestricted rows
# above, which are kept alongside these for comparison.
abx_s1_cols_r <- paste0("abx_s1_", 1:28)
ast_s1_cols_r <- paste0("ast_s1_", 1:28)
flag_resistance_to_s1 <- function(is_target_fn) {
  abx_mat_s1 <- as.matrix(Main_data_CaseReport[, abx_s1_cols_r])
  ast_mat_s1 <- as.matrix(Main_data_CaseReport[, ast_s1_cols_r])
  storage.mode(ast_mat_s1) <- "integer"
  is_target_mat <- matrix(is_target_fn(abx_mat_s1), nrow = nrow(abx_mat_s1), ncol = ncol(abx_mat_s1))
  tested <- rowSums(is_target_mat & !is.na(ast_mat_s1), na.rm = TRUE) > 0
  resistant <- rowSums(is_target_mat & ast_mat_s1 == 1L, na.rm = TRUE) > 0
  case_when(resistant ~ 1L, tested ~ 0L, TRUE ~ NA_integer_)
}
carb_resistance_s1 <- flag_resistance_to_s1(function(x) !is.na(x) & str_detect(toupper(x), "PENEM"))
c3gc_resistance_s1 <- flag_resistance_to_s1(function(x) !is.na(x) & toupper(x) %in% c("CEFTRIAXONE", "CEFOTAXIME", "CEFTAZIDIME", "CEFTAZIDIME-AVIBACTAM"))
# cephalosporin_resistance_s1: any-generation cephalosporin (same "^CEF" name
# match as the BSI-wide cephalosporin_resis, SECTION 2), S1-isolate-only and
# organism-restricted -- the any-generation counterpart to c3gc_resistance_s1
# immediately above, added so "cephalosporin resistance (any generation)" has
# the same Enterobacterales/P. aeruginosa/A. baumannii-restricted denominator
# as 3rd-generation cephalosporin and carbapenem resistance, instead of only
# existing as a population-wide (all-bacteria) measure.
cephalosporin_resistance_s1 <- flag_resistance_to_s1(function(x) !is.na(x) & str_detect(toupper(x), "^CEF"))
methicillin_resistance_s1 <- flag_resistance_to_s1(function(x) !is.na(x) & toupper(x) %in% c("OXACILLIN", "METHICILLIN", "CEFOXITIN"))
vancomycin_resistance_s1 <- flag_resistance_to_s1(function(x) !is.na(x) & toupper(x) == "VANCOMYCIN")

# Also assigned as real Main_data_CaseReport columns (not just bare vectors)
# so they survive SECTION 8's full_join into Main_data_Combined/
# mortality_analytic -- needed for the CoNS-included Risk-Factors
# replication below, which reads organism-restricted resistance status
# directly off mortality_analytic rather than re-deriving it.
Main_data_CaseReport$carb_resistance_s1 <- carb_resistance_s1
Main_data_CaseReport$c3gc_resistance_s1 <- c3gc_resistance_s1
Main_data_CaseReport$cephalosporin_resistance_s1 <- cephalosporin_resistance_s1
Main_data_CaseReport$methicillin_resistance_s1 <- methicillin_resistance_s1
Main_data_CaseReport$vancomycin_resistance_s1 <- vancomycin_resistance_s1

# is_cre_pa_ab(): the WHO GLASS priority organism group carbapenem, any-
# generation cephalosporin, and 3rd-generation cephalosporin resistance are
# surveilled against -- Enterobacterales, Pseudomonas aeruginosa, and
# Acinetobacter baumannii SPECIFICALLY, not the broader "Non-fermenting
# Gram-negative bacilli" / "Other Gram-negative" major_bacterial_group
# buckets (SECTION 2's genus_group_map), which also sweep in organisms such
# as Stenotrophomonas maltophilia and Burkholderia spp. -- species with
# well-documented INTRINSIC (not acquired) carbapenem/cephalosporin
# resistance (e.g. Stenotrophomonas' chromosomal L1/L2 metallo-beta-
# lactamases). Pooling those species into a "Gram-negative resistance" rate
# would inflate it with resistance that reflects species identity, not an
# acquired resistance signal -- exactly what restricting to the named
# WHO GLASS priority pathogens avoids. Takes vectors (not a data frame) so
# it works identically on Main_data_CaseReport and any of its derived
# subsets (mortality_analytic, mortality_paired, matched-cohort risk-factor
# data, ...), which all carry major_bacterial_group and pathogen1_clean.
is_cre_pa_ab <- function(major_bacterial_group_vec, pathogen1_clean_vec) {
  # revised: ALL Gram-negative isolates (Enterobacterales, non-fermenting Gram-negative bacilli, other Gram-negative)
  (!is.na(major_bacterial_group_vec) & major_bacterial_group_vec %in% c("Enterobacterales", "Non-fermenting Gram-negative bacilli", "Other Gram-negative"))
}
cre_pa_ab_group <- is_cre_pa_ab(Main_data_CaseReport$major_bacterial_group, Main_data_CaseReport$pathogen1_clean)
# MRSA/methicillin resistance is restricted to Staphylococcus AUREUS
# specifically (pathogen1_clean), not the wider "Staphylococcus spp." major
# bacterial group -- that group includes coagulase-negative Staphylococcus
# species (S. epidermidis, S. hominis, S. lentus, etc.), whose methicillin
# resistance is a materially different clinical entity (MRSE/MRSH-type,
# not MRSA) from true MRSA. Corrected throughout this script (prevalence,
# incidence, risk factors, mortality, LOS) to this species-level
# restriction; previously conflated with the whole Staphylococcus genus.
staph_group <- !is.na(Main_data_CaseReport$pathogen1_clean) & Main_data_CaseReport$pathogen1_clean == "Staphylococcus aureus"
enterococcus_group <- !is.na(Main_data_CaseReport$genus1) & Main_data_CaseReport$genus1 == "Enterococcus"

## Per instruction, the intermediate PPS design-weighted-only estimate is no
## longer a reported column anywhere in this script (SECTION 3.4's note) --
## every *_row()/heatmap_*() function below reports unweighted + calibrated
## only. design_weight itself is retained (calibrated_weight is derived from
## it), just not reported on its own.
organism_restricted_resistance_row <- function(label, resistance_s1_vec, organism_filter) {
  eligible <- organism_filter & !is.na(resistance_s1_vec)
  n_valid <- sum(eligible)
  n_pos <- sum(eligible & resistance_s1_vec == 1L)
  unw <- wilson_ci(n_pos, n_valid) %>% mutate(across(everything(), ~ .x * 100))

  row_data <- Main_data_CaseReport %>%
    mutate(res_val = resistance_s1_vec, elig = eligible) %>%
    filter(!is.na(calibrated_weight), elig)
  cal <- svyciprop(~I(res_val == 1), svydesign(ids = ~facility_id_clean, weights = ~calibrated_weight, data = row_data), method = "logit")
  cal_ci <- confint(cal)

  tibble(
    measure = label, n_valid = n_valid, n_positive = n_pos,
    unweighted_pct = unw$est, unweighted_lower = unw$lower, unweighted_upper = unw$upper,
    calibrated_pct = as.numeric(cal) * 100, calibrated_lower = cal_ci[1] * 100, calibrated_upper = cal_ci[2] * 100
  )
}

amr_pop_data <- Main_data_CaseReport %>% filter(!is.na(calibrated_weight))
amr_survey_design_calib <- svydesign(ids = ~facility_id_clean, weights = ~calibrated_weight, data = amr_pop_data)

amr_prevalence_row <- function(label, var) {
  x <- Main_data_CaseReport[[var]]
  n_valid <- sum(!is.na(x))
  n_pos <- sum(x, na.rm = TRUE)
  unw <- wilson_ci(n_pos, n_valid) %>% mutate(across(everything(), ~ .x * 100))

  form <- as.formula(paste0("~I(`", var, "` == 1)"))
  cal <- svyciprop(form, subset(amr_survey_design_calib, !is.na(get(var))), method = "logit")
  cal_ci <- confint(cal)

  tibble(
    measure = label, n_valid = n_valid, n_positive = n_pos,
    unweighted_pct = unw$est, unweighted_lower = unw$lower, unweighted_upper = unw$upper,
    calibrated_pct = as.numeric(cal) * 100, calibrated_lower = cal_ci[1] * 100, calibrated_upper = cal_ci[2] * 100
  )
}

amr_prevalence_combined_table <- bind_rows(
  amr_prevalence_row("Any drug-resistant BSI", "resistant_BSI_positive"),
  amr_prevalence_row("MDR BSI (bacterial isolates)", "mdr_BSI_bacterial"),
  amr_prevalence_row("XDR BSI (bacterial isolates)", "xdr_BSI_bacterial"),
  amr_prevalence_row("PDR BSI (bacterial isolates)", "pdr_BSI_bacterial"),
  amr_prevalence_row("Carbapenem-resistant BSI", "carb_BSI_positive"),
  amr_prevalence_row("Methicillin-resistant BSI", "methicillin_BSI_positive"),
  amr_prevalence_row("Cephalosporin-resistant BSI (any generation)", "cephalosporin_BSI_positive"),
  amr_prevalence_row("3rd-generation cephalosporin-resistant BSI", "c3GC_BSI_positive"),
  amr_prevalence_row("Vancomycin-resistant BSI", "vancomycin_BSI_positive"),
  amr_prevalence_row("Access-antibiotic resistant BSI (AWaRe)", "access_BSI_positive"),
  amr_prevalence_row("Watch-antibiotic resistant BSI (AWaRe)", "watch_BSI_positive"),
  amr_prevalence_row("Reserve-antibiotic resistant BSI (AWaRe)", "reserve_BSI_positive"),
  organism_restricted_resistance_row("Carbapenem resistance (Enterobacterales/P. aeruginosa/A. baumannii isolates)", carb_resistance_s1, cre_pa_ab_group),
  organism_restricted_resistance_row("Methicillin resistance (S. aureus isolates)", methicillin_resistance_s1, staph_group),
  organism_restricted_resistance_row("Cephalosporin resistance, any generation (Enterobacterales/P. aeruginosa/A. baumannii isolates)", cephalosporin_resistance_s1, cre_pa_ab_group),
  organism_restricted_resistance_row("3rd-gen cephalosporin resistance (Enterobacterales/P. aeruginosa/A. baumannii isolates)", c3gc_resistance_s1, cre_pa_ab_group),
  organism_restricted_resistance_row("Vancomycin resistance (Enterococcus isolates)", vancomycin_resistance_s1, enterococcus_group)
) %>%
  mutate(across(where(is.numeric) & !c(n_valid, n_positive), ~ round(.x, 2)))

write_xlsx(amr_prevalence_combined_table, "AMR_prevalence_combined_table.xlsx")
write_xlsx(aware_classification, "AWaRe_antibiotic_classification.xlsx")

# -----------------------------------------------------------------------
# 3.6.1b: Resistance proportion AMONG BSI-POSITIVE cultures specifically
# (denominator = confirmed BSI cases with a determinate result for that
# measure, not all evaluable enrolled patients like 3.6.1 above, which
# treats BSI-negative patients as a structural zero). Built for the
# summary resistance-proportion-by-mechanism figure; MRSA, carbapenem,
# 3GC, and vancomycin resistance already have an organism-restricted "among
# BSI+" equivalent in amr_prevalence_combined_table above (the correct,
# clinically meaningful denominator for an organism-specific mechanism is
# that organism's own tested isolates, not all BSI+ isolates regardless of
# species) and are not recomputed here.
# -----------------------------------------------------------------------
amr_among_bsi_row <- function(label, var) {
  x <- Main_data_CaseReport[[var]]
  eligible <- Main_data_CaseReport$BSI_positive == 1 & !is.na(x)
  n_valid <- sum(eligible)
  n_pos <- sum(eligible & x == 1)
  unw <- wilson_ci(n_pos, n_valid) %>% mutate(across(everything(), ~ .x * 100))

  row_data <- Main_data_CaseReport %>%
    mutate(res_val = .data[[var]]) %>%
    filter(eligible, !is.na(calibrated_weight))
  cal <- svyciprop(~I(res_val == 1), svydesign(ids = ~facility_id_clean, weights = ~calibrated_weight, data = row_data), method = "logit")
  cal_ci <- confint(cal)

  tibble(
    measure = label, n_valid = n_valid, n_positive = n_pos,
    unweighted_pct = unw$est, unweighted_lower = unw$lower, unweighted_upper = unw$upper,
    calibrated_pct = as.numeric(cal) * 100, calibrated_lower = cal_ci[1] * 100, calibrated_upper = cal_ci[2] * 100
  )
}

amr_among_bsi_positive_table <- bind_rows(
  amr_among_bsi_row("Any drug-resistant BSI", "resistant_BSI_positive"),
  amr_among_bsi_row("MDR BSI (bacterial isolates)", "mdr_BSI_bacterial"),
  amr_among_bsi_row("XDR BSI (bacterial isolates)", "xdr_BSI_bacterial"),
  amr_among_bsi_row("Cephalosporin-resistant BSI (any generation)", "cephalosporin_BSI_positive"),
  amr_among_bsi_row("Access-antibiotic resistant BSI (AWaRe)", "access_BSI_positive"),
  amr_among_bsi_row("Watch-antibiotic resistant BSI (AWaRe)", "watch_BSI_positive"),
  amr_among_bsi_row("Reserve-antibiotic resistant BSI (AWaRe)", "reserve_BSI_positive")
) %>%
  mutate(across(where(is.numeric) & !c(n_valid, n_positive), ~ round(.x, 2)))

write_xlsx(amr_among_bsi_positive_table, "AMR_prevalence_among_BSIpositive_table.xlsx")

rm(amr_pop_data)

# -----------------------------------------------------------------------
# 3.6.2: Pathogen distribution among confirmed BSI cases
# -----------------------------------------------------------------------
# Uses pathogen1_clean as THE causative organism per BSI case (standard
# single-organism-per-episode reporting), not pathogen2 -- 131 of 464
# identified cases (28.2%) also had a second isolate (pathogen2_clean);
# see polymicrobial_bsi_rate below rather than double-counting cases across
# both pathogen columns. "Species not identified" (26 cases, 5.6%) is
# BSI-confirmed but has no recorded pathogen1 at all -- kept as its own row
# rather than silently dropped. Percentages are of confirmed BSI cases
# (BSI_positive == 1), NOT of all tested patients (a different, much larger
# denominator, per SECTION 3.4).
polymicrobial_bsi_rate <- Main_data_CaseReport %>%
  filter(BSI_positive == 1) %>%
  summarise(n = n(), n_poly = sum(!is.na(pathogen1_clean) & !is.na(pathogen2_clean))) %>%
  mutate(pct_polymicrobial = round(100 * n_poly / n, 1))

bsi_design_only_calib <- subset(bsi_survey_design_calib, BSI_positive == 1)
n_bsi_total <- sum(Main_data_CaseReport$BSI_positive == 1, na.rm = TRUE)

pathogen_dist_row <- function(level_value, level_type) {
  if (level_value == "Species not identified") {
    ind_full <- is.na(Main_data_CaseReport$pathogen1_clean)
    d_c <- update(bsi_design_only_calib, pth_ind = is.na(pathogen1_clean))
  } else if (level_type == "species") {
    ind_full <- !is.na(Main_data_CaseReport$pathogen1_clean) & Main_data_CaseReport$pathogen1_clean == level_value
    d_c <- update(bsi_design_only_calib, pth_ind = !is.na(pathogen1_clean) & pathogen1_clean == level_value)
  } else {
    ind_full <- !is.na(Main_data_CaseReport$major_bacterial_group) & Main_data_CaseReport$major_bacterial_group == level_value
    d_c <- update(bsi_design_only_calib, pth_ind = !is.na(major_bacterial_group) & major_bacterial_group == level_value)
  }
  n_this <- sum(Main_data_CaseReport$BSI_positive == 1 & ind_full, na.rm = TRUE)
  unw <- wilson_ci(n_this, n_bsi_total) %>% mutate(across(everything(), ~ .x * 100))
  cal <- svyciprop(~pth_ind, d_c, method = "logit"); cal_ci <- confint(cal)
  tibble(
    level = level_type, pathogen = level_value, n = n_this, n_bsi_total = n_bsi_total,
    unweighted_pct = unw$est, unweighted_lower = unw$lower, unweighted_upper = unw$upper,
    calibrated_pct = as.numeric(cal) * 100, calibrated_lower = cal_ci[1] * 100, calibrated_upper = cal_ci[2] * 100
  )
}

species_levels <- Main_data_CaseReport %>%
  filter(BSI_positive == 1, !is.na(pathogen1_clean)) %>%
  count(pathogen1_clean, sort = TRUE) %>%
  pull(pathogen1_clean)
group_levels <- Main_data_CaseReport %>%
  filter(BSI_positive == 1, !is.na(major_bacterial_group)) %>%
  count(major_bacterial_group, sort = TRUE) %>%
  pull(major_bacterial_group)

pathogen_distribution_species <- bind_rows(
  purrr::map_dfr(species_levels, pathogen_dist_row, level_type = "species"),
  pathogen_dist_row("Species not identified", "species")
) %>%
  arrange(desc(n)) %>%
  mutate(across(where(is.numeric) & !c(n, n_bsi_total), ~ round(.x, 2)))

pathogen_distribution_by_group <- purrr::map_dfr(group_levels, pathogen_dist_row, level_type = "group") %>%
  arrange(desc(n)) %>%
  mutate(across(where(is.numeric) & !c(n, n_bsi_total), ~ round(.x, 2)))

write_xlsx(
  list(by_species = pathogen_distribution_species, by_major_group = pathogen_distribution_by_group,
       polymicrobial_rate = polymicrobial_bsi_rate),
  "Pathogen_distribution_combined_table.xlsx"
)

# -----------------------------------------------------------------------
# 3.6.3: Pathogen x antibiotic resistance heatmap (top 20 antibiotics)
# -----------------------------------------------------------------------
# Uses ONLY the S1 (Species 1 / primary isolate) AST slots, paired with
# pathogen1_clean -- S2's AST results belong to a DIFFERENT isolate
# (pathogen2), so mixing them into a "pathogen x antibiotic" cell would
# pair some resistance results with the wrong organism. Fungal isolates are
# excluded (this AMR panel is a bacterial antibiogram). "Most relevant" is
# operationalised as the 20 antibiotics with the most S1 AST results
# actually recorded (a data-driven proxy for both lab-panel priority and
# clinical relevance) -- the 20th (Fusidic acid) has 118 tests, the 21st
# (Cefepime) has 117, so this cutoff is not sensitive to exactly where 20
# falls.
abx_s1_cols <- paste0("abx_s1_", 1:28)
ast_s1_cols <- paste0("ast_s1_", 1:28)
abx_s1_test_long <- purrr::map2_dfr(abx_s1_cols, ast_s1_cols, function(a_col, s_col) {
  tibble(antibiotic = Main_data_CaseReport[[a_col]], ast = Main_data_CaseReport[[s_col]])
}) %>% filter(!is.na(antibiotic), !is.na(ast))
top20_antibiotics <- abx_s1_test_long %>% count(antibiotic, sort = TRUE) %>% slice_head(n = 20) %>% pull(antibiotic)

get_primary_ast <- function(df, antibiotic_name) {
  result <- rep(NA_integer_, nrow(df))
  for (i in seq_along(abx_s1_cols)) {
    match_i <- !is.na(df[[abx_s1_cols[i]]]) & df[[abx_s1_cols[i]]] == antibiotic_name
    result[match_i] <- df[[ast_s1_cols[i]]][match_i]
  }
  result
}

amr_heatmap_data <- Main_data_CaseReport %>%
  filter(BSI_positive == 1, pathogen_domain == "Bacteria", !is.na(calibrated_weight))
heatmap_design_c <- svydesign(ids = ~facility_id_clean, weights = ~calibrated_weight, data = amr_heatmap_data)

species_n_ge10 <- amr_heatmap_data %>% count(pathogen1_clean, sort = TRUE) %>% filter(n >= 10) %>% pull(pathogen1_clean)
group_n_all <- amr_heatmap_data %>% distinct(major_bacterial_group) %>% pull(major_bacterial_group)
heatmap_rows <- tibble(
  row_label = c(species_n_ge10, group_n_all, "All bacterial isolates"),
  row_type = c(rep("Species (n>=10)", length(species_n_ge10)), rep("Major bacterial group", length(group_n_all)), "Overall")
)

# IMPRECISE flag: n_tested < 10, OR the calibrated 95% CI spans > 50
# percentage points, OR every tested isolate agreed (n_resistant == 0 or
# == n_tested) -- all three are transparent, documented thresholds, not a
# claim of adequate precision below them. The third catches a real gap in
# the first two: svyciprop()'s logit-based CI collapses to a degenerate
# zero-width interval at a 0%/100% point estimate (e.g. n_tested=10,
# n_resistant=10 -> calibrated 100% (100.0-100.0)) rather than reflecting
# the genuine uncertainty of a small sample -- true population resistance
# is essentially never exactly 0% or 100%, so a cell landing there is a
# small-sample artefact regardless of how many isolates were tested.
heatmap_cell <- function(row_label, row_type, antibiotic_name, ast_vec) {
  row_filter <- if (row_type == "Species (n>=10)") {
    amr_heatmap_data$pathogen1_clean == row_label
  } else if (row_type == "Major bacterial group") {
    amr_heatmap_data$major_bacterial_group == row_label
  } else {
    rep(TRUE, nrow(amr_heatmap_data))
  }
  valid <- row_filter & !is.na(ast_vec)
  n_tested <- sum(valid)
  n_resistant <- sum(valid & ast_vec == 1L)
  if (n_tested == 0) {
    return(tibble(
      pathogen = row_label, row_type = row_type, antibiotic = antibiotic_name,
      n_tested = 0L, n_resistant = 0L, unweighted_pct = NA_real_,
      calibrated_pct = NA_real_, calibrated_lower = NA_real_, calibrated_upper = NA_real_,
      imprecise = NA
    ))
  }
  unw_pct <- 100 * n_resistant / n_tested

  d_c <- update(heatmap_design_c, ast_val = ast_vec, row_ind = row_filter)
  cal_est <- tryCatch({
    cal <- svyciprop(~I(ast_val == 1), subset(d_c, row_ind & !is.na(ast_val)), method = "logit")
    c(as.numeric(cal), confint(cal)) * 100
  }, error = function(e) c(NA_real_, NA_real_, NA_real_))

  tibble(
    pathogen = row_label, row_type = row_type, antibiotic = antibiotic_name,
    n_tested = n_tested, n_resistant = n_resistant, unweighted_pct = unw_pct,
    calibrated_pct = cal_est[1], calibrated_lower = cal_est[2], calibrated_upper = cal_est[3],
    imprecise = n_tested < 10 | (cal_est[3] - cal_est[2]) > 50 | n_resistant == 0 | n_resistant == n_tested
  )
}

amr_heatmap_long <- purrr::map_dfr(top20_antibiotics, function(abx) {
  ast_vec <- get_primary_ast(amr_heatmap_data, abx)
  purrr::pmap_dfr(heatmap_rows, function(row_label, row_type) heatmap_cell(row_label, row_type, abx, ast_vec))
}) %>%
  mutate(
    aware_category = aware_lookup[toupper(trimws(antibiotic))],
    across(where(is.numeric) & !c(n_tested, n_resistant), ~ round(.x, 1))
  ) %>%
  relocate(aware_category, .after = antibiotic)

write_xlsx(amr_heatmap_long, "AMR_heatmap_pathogen_x_antibiotic.xlsx")

## --- heatmap figures: unweighted / calibrated %, both (imprecise cells hatched) ---
# PPS design-weighted-only reporting was removed per instruction; only the
# unweighted (raw) and calibrated versions are rendered, both built by the
# same generic build_amr_heatmap_plot() renderer -- same rows, columns,
# ordering, and imprecision hatching, so the two read side by side.
heatmap_plot_data <- amr_heatmap_long %>%
  mutate(
    pathogen = factor(pathogen, levels = rev(c(species_n_ge10, group_n_all, "All bacterial isolates"))),
    antibiotic = factor(antibiotic, levels = top20_antibiotics)
  )

build_amr_heatmap_plot <- function(data, value_col, legend_label, title_suffix) {
  # Cell label: raw ratio n_resistant/n_tested (e.g. "1/5") instead of a
  # rounded percentage -- fill colour still encodes value_col (unweighted
  # or calibrated %), so the ratio is a precision cue layered on top of the
  # same colour scale, not a second scale to read. A reader can see at a
  # glance how many isolates a cell's colour rests on without a separate
  # lookup in the underlying xlsx.
  ggplot(data, aes(x = antibiotic, y = pathogen, fill = .data[[value_col]])) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = if_else(is.na(.data[[value_col]]), "",
                                   sprintf("%d/%d", n_resistant, n_tested))),
              size = 2.3) +
    geom_tile(data = filter(data, imprecise), fill = NA, color = "grey30", linewidth = 0.5, linetype = "dashed") +
    scale_fill_gradient(low = "#f7fbff", high = "#b2182b", na.value = "grey90", limits = c(0, 100), name = legend_label) +
    labs(
      title = paste("Antimicrobial resistance by pathogen and antibiotic --", title_suffix),
      subtitle = "Each cell: n resistant / n tested. Fill = % resistant (colour scale). Dashed border = imprecise (n<10 tested, 95% CI width >50pp, or every tested isolate agreed). Bold rows = major bacterial group / overall.",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
}

amr_heatmap_plot_unweighted <- build_amr_heatmap_plot(heatmap_plot_data, "unweighted_pct", "% resistant\n(unweighted)", "unweighted (raw)")
ggsave("AMR_heatmap_pathogen_x_antibiotic_unweighted.png", amr_heatmap_plot_unweighted, width = 13, height = 8, dpi = 500)

amr_heatmap_plot_calibrated <- build_amr_heatmap_plot(heatmap_plot_data, "calibrated_pct", "% resistant\n(calibrated)", "calibrated to oblast population")
ggsave("AMR_heatmap_pathogen_x_antibiotic_calibrated.png", amr_heatmap_plot_calibrated, width = 13, height = 8, dpi = 500)

rm(abx_s1_test_long, bsi_design_only_calib)

# -----------------------------------------------------------------------
# 3.6.4-3.6.6: Same three analyses, RATE PER 100 ADMISSIONS (SECTION 3.5's
# approach) instead of proportion among tested patients
# -----------------------------------------------------------------------
# Same logic and same caveat as SECTION 3.5 vs 3.4: Total_inpatient_admissions
# is only a whole-facility count, not broken down by pathogen, resistance
# phenotype, or antibiotic. Every rate below is therefore a CAUSE-SPECIFIC
# rate sharing the same facility_rate_base admissions denominator (all 40
# facilities, per SECTION 3.5's Hospital ID fix) -- "of every 100
# admissions, X resulted in an MDR BSI" -- not a subgroup-specific risk.
# Levels of the same variable (e.g. every pathogen1_clean species) sum
# (up to missingness) back to that variable's own overall row.
#
# rate_per100_row() is the general form of SECTION 3.5's rate_per100_by_group():
# takes any 0/1/NA per-patient indicator (case-confirmed = 1), aggregates it
# to a facility count (NA treated as "not a confirmed case", consistent with
# how resistant_BSI_positive/mdr_BSI_bacterial/etc. themselves treat
# "AST too incomplete to tell" -- i.e. this is a conservative, confirmed-
# cases-only numerator, not an imputed one), 0-fills facilities with no
# matching patient, then reports unweighted (Poisson/Garwood CI) and
# weighted/calibrated (facility-level svyratio, as in SECTION 3.5) rates.
# All 20 antibiotics x ~24 pathogen/group rows in 3.6.6 means several
# hundred cells, most with a facility numerator of 0 admissions-wide --
# svyratio degenerates safely to a point CI of [0,0] in that case (handled
# explicitly below, not left to survey:: to error on a zero-variance ratio).
rate_per100_row <- function(indicator_vec, label) {
  df <- tibble(facility_id_clean = Main_data_CaseReport$facility_id_clean, ind = indicator_vec) %>%
    filter(!is.na(facility_id_clean)) %>%
    group_by(facility_id_clean) %>%
    summarise(n_num = sum(ind, na.rm = TRUE), .groups = "drop")

  fac_df <- facility_rate_base %>%
    left_join(df, by = c(`Hospital ID` = "facility_id_clean")) %>%
    mutate(n_num = if_else(is.na(n_num), 0L, as.integer(n_num)))

  unw <- poisson_rate_ci(sum(fac_df$n_num), sum(fac_df$Total_inpatient_admissions)) %>%
    mutate(across(everything(), ~ .x * 100))

  if (sum(fac_df$n_num) == 0) {
    cal_est <- c(0, 0, 0)
  } else {
    des_c <- svydesign(ids = ~`Hospital ID`, weights = ~calibrated_weight, data = fac_df)
    ratio_c <- svyratio(~n_num, ~Total_inpatient_admissions, des_c)
    ci_c <- confint(ratio_c) * 100
    cal_est <- c(as.numeric(coef(ratio_c)) * 100, ci_c[1], ci_c[2])
  }
  # svyratio()'s CI is a symmetric Taylor-linearised (Wald-type) interval,
  # unlike the Poisson/Garwood interval used for the unweighted rate -- for
  # small ratios (most cells here) it can extend below 0, which is not a
  # meaningful rate. Floored at 0, not re-derived; the upper bound is
  # untouched.
  cal_est[2] <- max(cal_est[2], 0)

  bind_cols(
    label,
    tibble(
      n_numerator = sum(fac_df$n_num), total_admissions = sum(fac_df$Total_inpatient_admissions),
      unweighted_rate_per100 = unw$est, unweighted_lower = unw$lower, unweighted_upper = unw$upper,
      calibrated_rate_per100 = cal_est[1], calibrated_lower = cal_est[2], calibrated_upper = cal_est[3]
    )
  )
}

## --- 3.6.4: AMR prevalence, rate per 100 admissions ---
amr_rate_per100_table <- bind_rows(
  rate_per100_row(Main_data_CaseReport$resistant_BSI_positive, tibble(measure = "Any drug-resistant BSI")),
  rate_per100_row(Main_data_CaseReport$mdr_BSI_bacterial, tibble(measure = "MDR BSI (bacterial isolates)")),
  rate_per100_row(Main_data_CaseReport$xdr_BSI_bacterial, tibble(measure = "XDR BSI (bacterial isolates)")),
  rate_per100_row(Main_data_CaseReport$pdr_BSI_bacterial, tibble(measure = "PDR BSI (bacterial isolates)")),
  rate_per100_row(Main_data_CaseReport$carb_BSI_positive, tibble(measure = "Carbapenem-resistant BSI")),
  rate_per100_row(Main_data_CaseReport$methicillin_BSI_positive, tibble(measure = "Methicillin-resistant BSI")),
  rate_per100_row(Main_data_CaseReport$cephalosporin_BSI_positive, tibble(measure = "Cephalosporin-resistant BSI (any generation)")),
  rate_per100_row(Main_data_CaseReport$c3GC_BSI_positive, tibble(measure = "3rd-generation cephalosporin-resistant BSI")),
  rate_per100_row(Main_data_CaseReport$vancomycin_BSI_positive, tibble(measure = "Vancomycin-resistant BSI")),
  rate_per100_row(Main_data_CaseReport$access_BSI_positive, tibble(measure = "Access-antibiotic resistant BSI (AWaRe)")),
  rate_per100_row(Main_data_CaseReport$watch_BSI_positive, tibble(measure = "Watch-antibiotic resistant BSI (AWaRe)")),
  rate_per100_row(Main_data_CaseReport$reserve_BSI_positive, tibble(measure = "Reserve-antibiotic resistant BSI (AWaRe)")),
  rate_per100_row(
    as.integer(cre_pa_ab_group & !is.na(carb_resistance_s1) & carb_resistance_s1 == 1L),
    tibble(measure = "Carbapenem resistance (Enterobacterales/P. aeruginosa/A. baumannii isolates)")
  ),
  rate_per100_row(
    as.integer(staph_group & !is.na(methicillin_resistance_s1) & methicillin_resistance_s1 == 1L),
    tibble(measure = "Methicillin resistance (S. aureus isolates)")
  ),
  rate_per100_row(
    as.integer(cre_pa_ab_group & !is.na(cephalosporin_resistance_s1) & cephalosporin_resistance_s1 == 1L),
    tibble(measure = "Cephalosporin resistance, any generation (Enterobacterales/P. aeruginosa/A. baumannii isolates)")
  ),
  rate_per100_row(
    as.integer(cre_pa_ab_group & !is.na(c3gc_resistance_s1) & c3gc_resistance_s1 == 1L),
    tibble(measure = "3rd-gen cephalosporin resistance (Enterobacterales/P. aeruginosa/A. baumannii isolates)")
  ),
  rate_per100_row(
    as.integer(enterococcus_group & !is.na(vancomycin_resistance_s1) & vancomycin_resistance_s1 == 1L),
    tibble(measure = "Vancomycin resistance (Enterococcus isolates)")
  )
) %>%
  mutate(across(where(is.numeric) & !c(n_numerator, total_admissions), ~ round(.x, 4)))
write_xlsx(amr_rate_per100_table, "AMR_rate_per100admissions_combined_table.xlsx")

## --- 3.6.5: Pathogen distribution, rate per 100 admissions ---
pathogen_rate_row <- function(level_value, level_type) {
  if (level_value == "Species not identified") {
    ind <- as.integer(Main_data_CaseReport$BSI_positive == 1 & is.na(Main_data_CaseReport$pathogen1_clean))
  } else if (level_type == "species") {
    ind <- as.integer(Main_data_CaseReport$BSI_positive == 1 &
      !is.na(Main_data_CaseReport$pathogen1_clean) & Main_data_CaseReport$pathogen1_clean == level_value)
  } else {
    ind <- as.integer(Main_data_CaseReport$BSI_positive == 1 &
      !is.na(Main_data_CaseReport$major_bacterial_group) & Main_data_CaseReport$major_bacterial_group == level_value)
  }
  rate_per100_row(ind, tibble(level = level_type, pathogen = level_value))
}

pathogen_rate_species <- bind_rows(
  purrr::map_dfr(species_levels, pathogen_rate_row, level_type = "species"),
  pathogen_rate_row("Species not identified", "species")
) %>%
  arrange(desc(n_numerator)) %>%
  mutate(across(where(is.numeric) & !c(n_numerator, total_admissions), ~ round(.x, 4)))

pathogen_rate_by_group <- purrr::map_dfr(group_levels, pathogen_rate_row, level_type = "group") %>%
  arrange(desc(n_numerator)) %>%
  mutate(across(where(is.numeric) & !c(n_numerator, total_admissions), ~ round(.x, 4)))

write_xlsx(
  list(by_species = pathogen_rate_species, by_major_group = pathogen_rate_by_group),
  "Pathogen_rate_per100admissions_combined_table.xlsx"
)

## --- 3.6.6: AMR heatmap, rate per 100 admissions ---
heatmap_rate_cell <- function(row_label, row_type, antibiotic_name, ast_vec_full) {
  row_filter <- if (row_type == "Species (n>=10)") {
    !is.na(Main_data_CaseReport$pathogen1_clean) & Main_data_CaseReport$pathogen1_clean == row_label
  } else if (row_type == "Major bacterial group") {
    !is.na(Main_data_CaseReport$major_bacterial_group) & Main_data_CaseReport$major_bacterial_group == row_label
  } else {
    !is.na(Main_data_CaseReport$pathogen_domain) & Main_data_CaseReport$pathogen_domain == "Bacteria"
  }
  ind <- as.integer(row_filter & !is.na(ast_vec_full) & ast_vec_full == 1L)
  rate_per100_row(ind, tibble(pathogen = row_label, row_type = row_type, antibiotic = antibiotic_name))
}

amr_heatmap_rate_long <- purrr::map_dfr(top20_antibiotics, function(abx) {
  ast_vec_full <- get_primary_ast(Main_data_CaseReport, abx)
  purrr::pmap_dfr(heatmap_rows, function(row_label, row_type) heatmap_rate_cell(row_label, row_type, abx, ast_vec_full))
}) %>%
  mutate(
    aware_category = aware_lookup[toupper(trimws(antibiotic))],
    across(where(is.numeric) & !c(n_numerator, total_admissions), ~ round(.x, 4))
  ) %>%
  relocate(aware_category, .after = antibiotic)

write_xlsx(amr_heatmap_rate_long, "AMR_heatmap_rate_per100admissions_pathogen_x_antibiotic.xlsx")

## --- heatmap figure (weighted rate per 100 admissions) ---
# Same rows/columns/ordering as SECTION 3.6.3's percentage heatmap, so the
# two figures read side by side; no "imprecise" hatching here since
# rate_per100_row() doesn't compute the CI-width flag that svyciprop-based
# heatmap_cell() does (that threshold was defined in percentage points and
# doesn't translate to this scale) -- read cell counts (n_numerator, in the
# underlying xlsx) before trusting a rate built on very few admissions-wide
# cases.
heatmap_rate_plot_data <- amr_heatmap_rate_long %>%
  mutate(
    pathogen = factor(pathogen, levels = rev(c(species_n_ge10, group_n_all, "All bacterial isolates"))),
    antibiotic = factor(antibiotic, levels = top20_antibiotics)
  )
build_amr_heatmap_rate_plot <- function(data, value_col, legend_label, title_suffix) {
  ggplot(data, aes(x = antibiotic, y = pathogen, fill = .data[[value_col]])) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = if_else(n_numerator == 0, "", sprintf("%.3f", .data[[value_col]]))), size = 2.4) +
    scale_fill_gradient(low = "#f7fbff", high = "#b2182b", na.value = "grey90", name = legend_label) +
    labs(
      title = paste("AMR rate per 100 admissions, by pathogen and antibiotic --", title_suffix),
      subtitle = "Cause-specific rate vs. the whole-facility admissions denominator (SECTION 3.5/3.6.6) -- not a subgroup-specific risk. Bold rows = major bacterial group / overall.",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
}

amr_heatmap_rate_plot_unweighted <- build_amr_heatmap_rate_plot(heatmap_rate_plot_data, "unweighted_rate_per100", "Resistant cases\nper 100 admissions\n(unweighted)", "unweighted (raw)")
ggsave("AMR_heatmap_rate_per100admissions_pathogen_x_antibiotic_unweighted.png", amr_heatmap_rate_plot_unweighted, width = 13, height = 8, dpi = 500)

amr_heatmap_rate_plot_calibrated <- build_amr_heatmap_rate_plot(heatmap_rate_plot_data, "calibrated_rate_per100", "Resistant cases\nper 100 admissions\n(calibrated)", "calibrated to oblast population")
ggsave("AMR_heatmap_rate_per100admissions_pathogen_x_antibiotic_calibrated.png", amr_heatmap_rate_plot_calibrated, width = 13, height = 8, dpi = 500)



# =============================================================================
# SECTION 3.7: MASTER SUMMARY TABLE
# =============================================================================
# Consolidates every SECTION 3.4-3.6 table into one workbook, one sheet per
# domain, each sheet carrying BOTH scales (prevalence among tested AND rate
# per 100 admissions) and both estimates (unweighted / calibrated) side by
# side per row, joined on the row's own identifying
# key. Nothing is recomputed here -- this section only joins tables already
# built above, so any number here traces back to its SECTION 3.4-3.6 source.
# 4 sheets:
#   BSI_prevalence_and_rate      -- overall + subgroups (SECTIONS 3.4 + 3.5)
#   AMR_prevalence_and_rate      -- resistance measures, overall only, incl.
#                                    the organism-restricted rows (SECTION 3.6.1 + 3.6.4)
#   Pathogen_distribution_species -- causative-organism breakdown, species level (3.6.2 + 3.6.5)
#   Pathogen_distribution_group   -- same, major bacterial group level (3.6.2 + 3.6.5)
# A "_rate_" suffix always means per 100 admissions; a bare "_pct_"/"_lower"/
# "_upper" without that suffix is the prevalence-among-tested scale -- the
# two are NOT the same denominator (see SECTION 3.5's note) and must not be
# compared to each other directly.

# NOTE: the prevalence and rate tables both have bare "unweighted_lower" /
# "unweighted_upper" / "calibrated_lower" / "calibrated_upper" columns (one
# CI on the % scale, one on the rate scale) that are NOT part of the join
# key -- left_join() would silently suffix these .x/.y without an explicit
# rename first, which is exactly the kind of ambiguity this table exists to
# avoid.
# The "Overall" row's own group label differs between the two source tables
# ("All tested patients" in prevalence_combined_table vs "All admissions" in
# rate_per100_combined_table) -- joining on the raw group column would drop
# that row's rate_ columns entirely. join_key normalizes both to "Overall"
# for that one row; every subgroup row still joins on its real group label.
master_bsi_table <- prevalence_combined_table %>%
  mutate(join_key = if_else(variable == "Overall", "Overall", group)) %>%
  left_join(
    rate_per100_combined_table %>%
      select(-n_BSI_positive) %>%
      mutate(join_key = if_else(variable == "Overall", "Overall", group)) %>%
      select(-group) %>%
      rename_with(~ paste0("rate_", .x), -c(variable, join_key)),
    by = c("variable", "join_key")
  ) %>%
  select(-join_key)

master_amr_table <- amr_prevalence_combined_table %>%
  left_join(
    amr_rate_per100_table %>%
      select(-total_admissions) %>%
      rename_with(~ paste0("rate_", .x), -measure),
    by = "measure"
  )

master_pathogen_species_table <- pathogen_distribution_species %>%
  rename(n_bsi_cases = n) %>%
  left_join(
    pathogen_rate_species %>% select(level, pathogen, n_numerator, unweighted_rate_per100, unweighted_lower, unweighted_upper,
                                      calibrated_rate_per100, calibrated_lower, calibrated_upper) %>%
      rename_with(~ paste0("rate_", .x), -c(level, pathogen)),
    by = c("level", "pathogen")
  )

master_pathogen_group_table <- pathogen_distribution_by_group %>%
  rename(n_bsi_cases = n) %>%
  left_join(
    pathogen_rate_by_group %>% select(level, pathogen, n_numerator, unweighted_rate_per100, unweighted_lower, unweighted_upper,
                                       calibrated_rate_per100, calibrated_lower, calibrated_upper) %>%
      rename_with(~ paste0("rate_", .x), -c(level, pathogen)),
    by = c("level", "pathogen")
  )

write_xlsx(
  list(
    BSI_prevalence_and_rate = master_bsi_table,
    AMR_prevalence_and_rate = master_amr_table,
    Pathogen_distribution_species = master_pathogen_species_table,
    Pathogen_distribution_group = master_pathogen_group_table
  ),
  "Master_summary_table.xlsx"
)





# =============================================================================
# SECTION 4: ANALYSIS of INCIDENCE
# =============================================================================
# Incidence describes the RATE at which BSI occurs within the hospital
# population over time, using observed denominators (admissions,
# inpatient-days) rather than the tested-patient population used for
# prevalence (SECTIONS 3.1-3.7) -- a complementary, not competing, measure.
#
#   Any confirmed BSI (overall)   -- per 1,000 inpatient-days (primary) AND
#                                     per 1,000 admissions (secondary)
#   Community-acquired BSI        -- per 1,000 admissions
#   Hospital-acquired BSI         -- per 1,000 inpatient-days
#
# METHOD: survey-weighted Poisson regression with log(exposure) as an
# offset and design-based standard errors -- svyglm(n_cases ~ 1, design,
# family = quasipoisson(), offset = log_exposure), exponentiated. quasipoisson
# (not poisson) is used because svyglm's own documentation recommends it for
# design-based count models: it does not force the model-based Poisson
# mean-variance assumption, deferring instead to the design's sandwich/
# Taylor-linearised variance. The regression unit is the FACILITY (n_cases,
# admissions and inpatient-days are all facility-level counts, exactly as in
# SECTION 3.5), not the patient.
#
# TWO estimates throughout, not three: "unweighted" (facility cluster,
# weight = 1) and "weighted" (facility cluster, weight = calibrated_weight,
# SECTION 3.3's final weight) -- design_weight-only is skipped here since
# SECTIONS 3.2-3.3 already established that calibrated_weight is the
# preferred final weight.

## --- facility_incidence_base: one row per facility, both denominators ---
facility_incidence_base <- Hospital_characteristics %>%
  select(`Hospital ID`, Total_inpatient_admissions, Total_patientDays) %>%
  mutate(
    Total_inpatient_admissions = as.numeric(Total_inpatient_admissions),
    Total_patientDays = as.numeric(Total_patientDays)
  ) %>%
  left_join(facility_weights %>% select(UniqueID, NameAU, design_weight, calibrated_weight),
            by = c(`Hospital ID` = "UniqueID")) %>%
  filter(!is.na(design_weight), !is.na(Total_inpatient_admissions), !is.na(Total_patientDays))

## --- incidence_row(): facility-level svyglm(quasipoisson, offset), 2 estimates ---
# restrict_oblast, when given, restricts BOTH numerator (via the caller's
# indicator_vec, already oblast-filtered) AND denominator to that oblast's
# own facilities -- a genuine region-specific incidence rate. Without it,
# the denominator is the whole 40-facility base regardless of the
# numerator's subgroup (the same cause-specific-rate framing as SECTION
# 3.5/3.6.4-6: no facility-level breakdown exists for age/sex/ward/
# specialty/Charlson/infection-origin, so those subgroups' "incidence" is a
# share of the SAME national exposure, not a subgroup-specific exposure).
incidence_row <- function(indicator_vec, exposure_var, label, per = 1000, restrict_oblast = NULL) {
  base <- facility_incidence_base
  if (!is.null(restrict_oblast)) base <- base %>% filter(NameAU == restrict_oblast)

  case_counts <- tibble(facility_id_clean = Main_data_CaseReport$facility_id_clean, ind = indicator_vec) %>%
    filter(!is.na(facility_id_clean)) %>%
    group_by(facility_id_clean) %>%
    summarise(n_cases = sum(ind, na.rm = TRUE), .groups = "drop")

  fac_df <- base %>%
    left_join(case_counts, by = c(`Hospital ID` = "facility_id_clean")) %>%
    mutate(
      n_cases = if_else(is.na(n_cases), 0L, as.integer(n_cases)),
      exposure = .data[[exposure_var]],
      log_exposure = log(exposure),
      w_unweighted = 1
    )

  fit_rate <- function(weight_var) {
    if (nrow(fac_df) < 2) return(c(NA_real_, NA_real_, NA_real_))
    des <- svydesign(ids = ~`Hospital ID`, weights = reformulate(weight_var), data = fac_df)
    m <- tryCatch(
      suppressWarnings(svyglm(n_cases ~ 1, design = des, family = quasipoisson(), offset = log_exposure)),
      error = function(e) NULL
    )
    if (is.null(m)) return(c(NA_real_, NA_real_, NA_real_))
    est <- unname(exp(coef(m)[1]) * per)
    ci <- tryCatch(exp(confint(m)[1, ]) * per, error = function(e) c(NA_real_, NA_real_))
    c(est, unname(ci[1]), unname(ci[2]))
  }

  unw <- fit_rate("w_unweighted")
  cal <- fit_rate("calibrated_weight")

  tibble(
    measure = label, exposure_type = exposure_var,
    oblast = if (is.null(restrict_oblast)) "All" else restrict_oblast,
    n_facilities = nrow(fac_df), n_cases = sum(fac_df$n_cases), total_exposure = sum(fac_df$exposure),
    unweighted_rate_per1000 = unw[1], unweighted_lower = unw[2], unweighted_upper = unw[3],
    calibrated_rate_per1000 = cal[1], calibrated_lower = cal[2], calibrated_upper = cal[3]
  )
}

# -----------------------------------------------------------------------
# 4.1: Primary incidence estimates
# -----------------------------------------------------------------------
incidence_primary_table <- bind_rows(
  incidence_row(Main_data_CaseReport$BSI_positive, "Total_patientDays", "Any confirmed BSI (overall)"),
  incidence_row(Main_data_CaseReport$BSI_positive, "Total_inpatient_admissions", "Any confirmed BSI (overall)"),
  incidence_row(as.integer(Main_data_CaseReport$BSI_positive == 1 & Main_data_CaseReport$infection_origin_composite == 0),
                "Total_inpatient_admissions", "Community-acquired BSI"),
  incidence_row(as.integer(Main_data_CaseReport$BSI_positive == 1 & Main_data_CaseReport$infection_origin_composite == 1),
                "Total_patientDays", "Hospital-acquired BSI"),
  incidence_row(Main_data_CaseReport$resistant_BSI_positive, "Total_patientDays", "Any drug-resistant BSI"),
  incidence_row(Main_data_CaseReport$resistant_BSI_positive, "Total_inpatient_admissions", "Any drug-resistant BSI"),
  incidence_row(Main_data_CaseReport$mdr_BSI_bacterial, "Total_patientDays", "MDR BSI (bacterial isolates)"),
  incidence_row(Main_data_CaseReport$mdr_BSI_bacterial, "Total_inpatient_admissions", "MDR BSI (bacterial isolates)"),
  incidence_row(Main_data_CaseReport$xdr_BSI_bacterial, "Total_patientDays", "XDR BSI (bacterial isolates)"),
  incidence_row(Main_data_CaseReport$xdr_BSI_bacterial, "Total_inpatient_admissions", "XDR BSI (bacterial isolates)"),
  incidence_row(Main_data_CaseReport$pdr_BSI_bacterial, "Total_patientDays", "PDR BSI (bacterial isolates)"),
  incidence_row(Main_data_CaseReport$pdr_BSI_bacterial, "Total_inpatient_admissions", "PDR BSI (bacterial isolates)")
) %>%
  mutate(across(where(is.numeric) & !c(n_facilities, n_cases, total_exposure), ~ round(.x, 4)))

# -----------------------------------------------------------------------
# 4.2: Incidence by pathogen -- bacterial family / fungal, and by species
# -----------------------------------------------------------------------
incidence_by_pathogen_group <- purrr::map_dfr(group_levels, function(g) {
  ind <- as.integer(Main_data_CaseReport$BSI_positive == 1 &
    !is.na(Main_data_CaseReport$major_bacterial_group) & Main_data_CaseReport$major_bacterial_group == g)
  bind_rows(
    incidence_row(ind, "Total_patientDays", g),
    incidence_row(ind, "Total_inpatient_admissions", g)
  )
}) %>%
  mutate(level = "major_bacterial_group", .before = 1) %>%
  mutate(across(where(is.numeric) & !c(n_facilities, n_cases, total_exposure), ~ round(.x, 4)))

incidence_by_species <- purrr::map_dfr(species_n_ge10, function(sp) {
  ind <- as.integer(Main_data_CaseReport$BSI_positive == 1 &
    !is.na(Main_data_CaseReport$pathogen1_clean) & Main_data_CaseReport$pathogen1_clean == sp)
  bind_rows(
    incidence_row(ind, "Total_patientDays", sp),
    incidence_row(ind, "Total_inpatient_admissions", sp)
  )
}) %>%
  mutate(level = "species (n>=10 bacterial isolates)", .before = 1) %>%
  mutate(across(where(is.numeric) & !c(n_facilities, n_cases, total_exposure), ~ round(.x, 4)))

# -----------------------------------------------------------------------
# 4.3: Incidence by antibiotic class (never by individual antibiotic)
# -----------------------------------------------------------------------
incidence_by_amr_class <- bind_rows(
  incidence_row(Main_data_CaseReport$carb_BSI_positive, "Total_patientDays", "Carbapenem-resistant BSI"),
  incidence_row(Main_data_CaseReport$carb_BSI_positive, "Total_inpatient_admissions", "Carbapenem-resistant BSI"),
  incidence_row(Main_data_CaseReport$methicillin_BSI_positive, "Total_patientDays", "Methicillin-resistant BSI"),
  incidence_row(Main_data_CaseReport$methicillin_BSI_positive, "Total_inpatient_admissions", "Methicillin-resistant BSI"),
  incidence_row(Main_data_CaseReport$cephalosporin_BSI_positive, "Total_patientDays", "Cephalosporin-resistant BSI (any generation)"),
  incidence_row(Main_data_CaseReport$cephalosporin_BSI_positive, "Total_inpatient_admissions", "Cephalosporin-resistant BSI (any generation)"),
  incidence_row(Main_data_CaseReport$c3GC_BSI_positive, "Total_patientDays", "3rd-generation cephalosporin-resistant BSI"),
  incidence_row(Main_data_CaseReport$c3GC_BSI_positive, "Total_inpatient_admissions", "3rd-generation cephalosporin-resistant BSI"),
  incidence_row(Main_data_CaseReport$vancomycin_BSI_positive, "Total_patientDays", "Vancomycin-resistant BSI"),
  incidence_row(Main_data_CaseReport$vancomycin_BSI_positive, "Total_inpatient_admissions", "Vancomycin-resistant BSI"),
  incidence_row(Main_data_CaseReport$access_BSI_positive, "Total_patientDays", "Access-antibiotic resistant BSI (AWaRe)"),
  incidence_row(Main_data_CaseReport$access_BSI_positive, "Total_inpatient_admissions", "Access-antibiotic resistant BSI (AWaRe)"),
  incidence_row(Main_data_CaseReport$watch_BSI_positive, "Total_patientDays", "Watch-antibiotic resistant BSI (AWaRe)"),
  incidence_row(Main_data_CaseReport$watch_BSI_positive, "Total_inpatient_admissions", "Watch-antibiotic resistant BSI (AWaRe)"),
  incidence_row(Main_data_CaseReport$reserve_BSI_positive, "Total_patientDays", "Reserve-antibiotic resistant BSI (AWaRe)"),
  incidence_row(Main_data_CaseReport$reserve_BSI_positive, "Total_inpatient_admissions", "Reserve-antibiotic resistant BSI (AWaRe)")
) %>%
  mutate(across(where(is.numeric) & !c(n_facilities, n_cases, total_exposure), ~ round(.x, 4)))

# -----------------------------------------------------------------------
# 4.4: Incidence by subgroup -- mirrors SECTIONS 3.1-3.4's subgroups,
# oblast (NameAU) included per request, with a genuine region-specific
# denominator (restrict_oblast) rather than the shared national one used
# for the other subgroups (see incidence_row()'s comment above).
# -----------------------------------------------------------------------
subgroup_levels_list <- list(
  age_group = as.character(sort(unique(na.omit(Main_data_CaseReport$age_group)))),
  sex_at_birth = as.character(sort(unique(na.omit(Main_data_CaseReport$sex_at_birth)))),
  ward_type_icu = sort(unique(na.omit(Main_data_CaseReport$ward_type_icu))),
  specialty_group = sort(unique(na.omit(Main_data_CaseReport$specialty_group))),
  charlson_group = sort(unique(na.omit(Main_data_CaseReport$charlson_group))),
  infection_origin_composite = as.character(sort(unique(na.omit(Main_data_CaseReport$infection_origin_composite))))
)

incidence_by_subgroup_national <- purrr::imap_dfr(subgroup_levels_list, function(levels_vec, group_var) {
  purrr::map_dfr(levels_vec, function(lv) {
    ind <- as.integer(Main_data_CaseReport$BSI_positive == 1 &
      !is.na(Main_data_CaseReport[[group_var]]) & as.character(Main_data_CaseReport[[group_var]]) == lv)
    bind_rows(
      incidence_row(ind, "Total_patientDays", lv),
      incidence_row(ind, "Total_inpatient_admissions", lv)
    ) %>% mutate(variable = group_var, .before = 1)
  })
}) %>% rename(group = measure)

oblast_levels <- sort(unique(na.omit(Main_data_CaseReport$NameAU)))
incidence_by_oblast <- purrr::map_dfr(oblast_levels, function(ob) {
  ind <- as.integer(Main_data_CaseReport$BSI_positive == 1 &
    !is.na(Main_data_CaseReport$NameAU) & Main_data_CaseReport$NameAU == ob)
  bind_rows(
    incidence_row(ind, "Total_patientDays", ob, restrict_oblast = ob),
    incidence_row(ind, "Total_inpatient_admissions", ob, restrict_oblast = ob)
  )
}) %>%
  mutate(variable = "NameAU", .before = 1) %>%
  rename(group = measure)

incidence_by_subgroup <- bind_rows(incidence_by_subgroup_national, incidence_by_oblast) %>%
  mutate(across(where(is.numeric) & !c(n_facilities, n_cases, total_exposure), ~ round(.x, 4)))

# =========================================================================
# MAP DATA: oblast-level calibrated prevalence (per 100 admissions) and
# incidence (per 1000 admissions) for BSI+, any drug resistance, and MDR --
# generalises rate_per100_by_group()/facility_group_bsi() (SECTION 3.5,
# above) and incidence_row()'s existing per-oblast loop (this section,
# above) from BSI+-only to an arbitrary indicator column, for choropleth
# map figures by oblast.
cat("=== OBLAST x MECHANISM PREVALENCE + INCIDENCE (fast, for maps) ===\n")

facility_group_indicator <- function(indicator_var) {
  Main_data_CaseReport %>%
    filter(!is.na(facility_id_clean), !is.na(NameAU)) %>%
    group_by(facility_id_clean) %>%
    summarise(n_ind = sum(.data[[indicator_var]], na.rm = TRUE), .groups = "drop")
}

oblast_rate_per100 <- function(indicator_var, label) {
  fg <- facility_group_indicator(indicator_var)
  fac_df <- facility_rate_base %>%
    left_join(fg, by = c(`Hospital ID` = "facility_id_clean")) %>%
    mutate(n_ind = if_else(is.na(n_ind), 0L, n_ind)) %>%
    left_join(Main_data_CaseReport %>% distinct(facility_id_clean, NameAU), by = c(`Hospital ID` = "facility_id_clean"))
  purrr::map_dfr(sort(unique(na.omit(fac_df$NameAU))), function(ob) {
    sub <- fac_df %>% filter(NameAU == ob)
    # Some oblasts (e.g. Issyk-Kul, Naryn) have only 1-2 sampled facilities --
    # too few PSUs for svydesign()'s variance machinery. Fall back to a
    # simple calibrated-weighted ratio point estimate for the map (no CI
    # needed for a choropleth fill value).
    rate <- tryCatch({
      des_c <- svydesign(ids = ~`Hospital ID`, weights = ~calibrated_weight, data = sub)
      as.numeric(coef(svyratio(~n_ind, ~Total_inpatient_admissions, des_c))) * 100
    }, error = function(e) {
      100 * sum(sub$n_ind * sub$calibrated_weight) / sum(sub$Total_inpatient_admissions * sub$calibrated_weight)
    })
    tibble(oblast = ob, indicator = label, n_ind = sum(sub$n_ind), total_admissions = sum(sub$Total_inpatient_admissions),
           calibrated_rate_per100 = rate)
  })
}

oblast_prevalence_map_data <- bind_rows(
  oblast_rate_per100("BSI_positive", "Confirmed BSI"),
  oblast_rate_per100("resistant_BSI_positive", "Any drug resistance"),
  oblast_rate_per100("mdr_BSI_bacterial", "MDR")
)
write_xlsx(oblast_prevalence_map_data, "Oblast_map_prevalence_per100_calibrated.xlsx")

oblast_incidence_map_data <- purrr::map_dfr(oblast_levels, function(ob) {
  bind_rows(
    incidence_row(as.integer(Main_data_CaseReport$BSI_positive == 1 & Main_data_CaseReport$NameAU == ob),
                  "Total_inpatient_admissions", "Confirmed BSI", restrict_oblast = ob),
    incidence_row(as.integer(Main_data_CaseReport$resistant_BSI_positive == 1 & Main_data_CaseReport$NameAU == ob),
                  "Total_inpatient_admissions", "Any drug resistance", restrict_oblast = ob),
    incidence_row(as.integer(Main_data_CaseReport$mdr_BSI_bacterial == 1 & Main_data_CaseReport$NameAU == ob),
                  "Total_inpatient_admissions", "MDR", restrict_oblast = ob)
  )
}) %>% rename(indicator = measure)

# 1-facility oblasts (Issyk-Kul, Naryn) return NA from incidence_row()'s
# svyglm-based variance machinery (insufficient PSUs); backfill a simple
# calibrated-weighted point rate for the map's fill value only (no CI).
oblast_incidence_map_data <- oblast_incidence_map_data %>%
  mutate(calibrated_rate_per1000 = if_else(
    is.na(calibrated_rate_per1000) & n_facilities <= 1,
    1000 * n_cases / total_exposure,
    calibrated_rate_per1000
  ))
write_xlsx(oblast_incidence_map_data, "Oblast_map_incidence_per1000admissions_calibrated.xlsx")

cat("OBLAST MAP DATA SAVED\n")

# -----------------------------------------------------------------------
# 4.4b: CoNS sensitivity -- overall incidence recomputed under the two
# alternate CoNS-treatment populations (SECTION 2.6), reported alongside
# the unchanged-primary rows above, never replacing them. Overall level
# only (not repeated for every subgroup) -- consistent with how SECTION
# 2.6 itself reports the CoNS sensitivity check at the overall level, not
# resplit by every prevalence subgroup.
incidence_cons_sensitivity <- bind_rows(
  incidence_row(Main_data_CaseReport$BSI_positive, "Total_patientDays", "Primary (unchanged)"),
  incidence_row(Main_data_CaseReport$BSI_positive, "Total_inpatient_admissions", "Primary (unchanged)"),
  incidence_row(Main_data_CaseReport$BSI_positive_CoNS_sensitive, "Total_patientDays", "CoNS-excluded sensitivity"),
  incidence_row(Main_data_CaseReport$BSI_positive_CoNS_sensitive, "Total_inpatient_admissions", "CoNS-excluded sensitivity"),
  incidence_row(Main_data_CaseReport$BSI_positive_CoNS_included, "Total_patientDays", "CoNS-included sensitivity"),
  incidence_row(Main_data_CaseReport$BSI_positive_CoNS_included, "Total_inpatient_admissions", "CoNS-included sensitivity")
) %>%
  mutate(across(where(is.numeric) & !c(n_facilities, n_cases, total_exposure), ~ round(.x, 4)))

# -----------------------------------------------------------------------
# 4.5: write combined workbook
# -----------------------------------------------------------------------
write_xlsx(
  list(
    Primary = incidence_primary_table,
    By_pathogen_group = incidence_by_pathogen_group,
    By_species = incidence_by_species,
    By_AMR_class = incidence_by_amr_class,
    By_subgroup = incidence_by_subgroup,
    CoNS_sensitivity = incidence_cons_sensitivity
  ),
  "Incidence_combined_tables.xlsx"
)



# =============================================================================
# SECTION 5: Main_data_hospitalFollowup28
# =============================================================================
# mortality30day: "Status at 28-days review: Status at review" is a 5-level
# field -- ALIVE1/ALIVE2 (two flavours of "alive", pooled), DEAD, MORIBUND,
# and UNABLE (status could not be determined). 1 = DEAD or MORIBUND, 0 =
# ALIVE1 or ALIVE2, NA for UNABLE and for the 3 blank rows -- UNABLE is a
# genuine "don't know", not a 0, and is not conflated with it.
#
# Date_of_admission5enro / Date_of_discharge / LOS_days: Date_of_admission5enro
# sources from "Follow-up details: Date of enrollment in the survey" (as
# specified), EXCEPT where that produces a negative LOS_days -- there,
# Date_of_admission5enro is overwritten with "Follow-up details: Date of
# admission (dd/mm/yyyy)" (the sheet's own, separate admission-date field)
# and LOS_days is recomputed from that instead, per explicit follow-up
# instruction. Before this fallback, LOS_days had 31 negative rows and a
# max of 737,830 days (a garbled enrollment date); the fallback field is
# not itself guaranteed clean (10 of its own rows produce a negative LOS),
# so a small number of negative LOS_days values can still remain -- these
# are now a genuine source-data problem in the admission-date field itself,
# not an artifact of which field was chosen.
#
# in_hospital_mort: mortality30day wherever the patient is NOT confirmed to
# be at home at day 28 (the "Home" checkbox column is "true"/NA, not
# TRUE/FALSE -- NA is treated as "different from true", i.e. not confirmed
# home, per the literal instruction) -- 0 whenever they ARE confirmed home,
# regardless of mortality30day (a patient confirmed home cannot be an
# in-hospital death). NA propagates from mortality30day itself (UNABLE/blank).
#
# matched_control: patient_id ends in "/C3" or "/S1" (case-control /
# secondary-sample suffixes on the ward/sequence ID) -> 1, else 0 (including
# the 1 row with a blank patient_id, per "0 otherwise"). NOTE: 2 rows carry a
# suffix that looks like the same intent but doesn't match this exact
# pattern -- "045C3" (missing the "/") and "С3" (Cyrillic С, not Latin C) --
# left as 0 rather than guessed at; flagged here for awareness.
followup_status_col <- "Status at 28-days review: Status at review"
followup_home_col <- "Location of the patient at the time of the 28 days follow-up: Home"

Main_data_hospitalFollowup28_clean <- Main_data_hospitalFollowup28 %>%
  mutate(
    patient_id = `Unique survey patient ID #1`,
    mortality30day = case_when(
      .data[[followup_status_col]] %in% c("DEAD", "MORIBUND") ~ 1L,
      .data[[followup_status_col]] %in% c("ALIVE1", "ALIVE2") ~ 0L,
      TRUE ~ NA_integer_
    ),
    Date_of_admission5enro = parse_dt(`Follow-up details: Date of enrollment in the survey`),
    Date_of_admission_form = parse_dt(`Follow-up details: Date of admission (dd/mm/yyyy)`),
    Date_of_discharge = parse_dt(`Follow-up details: Date of discharge (dd/mm/yyyy)`),
    LOS_days_pre_fallback = as.numeric(difftime(Date_of_discharge, Date_of_admission5enro, units = "days")),
    Date_of_admission5enro = if_else(
      !is.na(LOS_days_pre_fallback) & LOS_days_pre_fallback < 0,
      Date_of_admission_form,
      Date_of_admission5enro
    ),
    LOS_days = as.numeric(difftime(Date_of_discharge, Date_of_admission5enro, units = "days")),
    # 7 rows are still implausible even after the fallback above: 5 stay
    # negative (the fallback admission-date field has its own bad entries
    # for these specific rows), and 2 have a garbled YEAR in one source
    # date (e.g. "0225-04-09", "0004-11-15"), producing a LOS in the
    # hundreds of thousands of days. Set to NA rather than left as a
    # negative or six-figure value -- per explicit follow-up instruction.
    LOS_days = if_else(LOS_days < 0 | LOS_days > 1000, NA_real_, LOS_days),
    in_hospital_mort = if_else(
      is.na(.data[[followup_home_col]]) | .data[[followup_home_col]] != "true",
      mortality30day,
      0L
    ),
    matched_control = is_true_matched_control(patient_id),
    matched_control_broad = is_cohort3_enrolled(patient_id)
  ) %>%
  # Winsorised at the 99th percentile: any LOS_days above the value that
  # bounds the bottom 99% of (remaining, non-NA) values is capped down to
  # that value -- caps the long upper tail without discarding those rows
  # (they keep every other variable; only LOS_days itself is capped).
  mutate(LOS_days = pmin(LOS_days, quantile(LOS_days, probs = 0.99, na.rm = TRUE))) %>%
  select(patient_id, mortality30day, Date_of_admission5enro, Date_of_discharge, LOS_days, in_hospital_mort, matched_control, matched_control_broad)




# =============================================================================
# SECTION 6: Main_data_hospitalDischarge
# =============================================================================
# patient_id / matched_control: identical construction to SECTION 5 (0 NA
# patient_ids in this sheet, unlike FollowUpD28's 1).
#
# ICU_entry: "Has the patient stayed in the ICU?" is "true"/"false"/NA
# (lowercase strings, 260 false / 959 true / 56 NA) -> 1/0/NA as specified;
# NA (status not captured) is left NA, not folded into either side.
#
# date_ICU_admission / date_ICU_discharge: parsed with parse_dt() (same
# ISO-8601 "...T18:00:00.000Z" source format as every other date field in
# this script).
#
# ICU_LOS_days: unlike SECTION 5's LOS_days, the raw source field ("1.
# Duration (in days)") has ZERO negative values in this extract (checked
# directly: min 0, max 71, mean 9.6) -- so no fallback-field swap is needed
# or built; the negative-value guard below is a safety net for future data,
# not something that currently fires. Capped at the 99th percentile
# (winsorised), same treatment as LOS_days.
#
# Complications index (any_complication / complications_count): the 20
# "Non infectious complications: ..." checkbox fields are YES/NO/UNK/NA
# (not true/NA), plus 1 free-text field ("Please, specify other non
# infectious complications", only 22/1275 non-blank) treated as positive
# whenever non-blank, per the same list. UNK (status genuinely unknown,
# ~2-3% of rows per field) is NA per field, not folded into 0 with NO --
# consistent with every other "don't know" field in this script.
# complications_count sums only the CONFIRMED (YES/NO) fields per patient
# (na.rm), so one UNK field doesn't blank out an otherwise-known count;
# any_complication/complications_count are NA only when ALL 20 checkbox
# fields are UNK/blank for that patient (the free-text field alone can't
# resolve that, since blank there just means nothing extra was specified).
discharge_icu_stayed_col <- "ICU Ward stays during hospital episode (from hospital admission to discharge/outcome): Has the patient stayed in the ICU?"
discharge_icu_los_col <- "ICU Ward stays during hospital episode (from hospital admission to discharge/outcome): 1. Duration (in days)"
discharge_icu_adm_col <- "ICU Ward stays during hospital episode (from hospital admission to discharge/outcome): Date of admission (dd/mm/yyyy)"
discharge_icu_disc_col <- "ICU Ward stays during hospital episode (from hospital admission to discharge/outcome): Date of discharge/outcome (dd/mm/yyyy)"

complication_cols <- c(
  "Non infectious complications: Acute renal injury",
  "Non infectious complications: Renal replacement therapy",
  "Non infectious complications: Acute respiratory diseases",
  "Non infectious complications: Bleeding",
  "Non infectious complications: Cardiac arrest",
  "Non infectious complications: Cardiac arrhytmia",
  "Non infectious complications: Cardiomyopathy",
  "Non infectious complications: Deep vein thrombosis",
  "Non infectious complications: Liver dysfunction",
  "Non infectious complications: Mental health disorder",
  "Non infectious complications: Pancreatitis",
  "Non infectious complications: Pulmonary embolism",
  "Non infectious complications: Seizure",
  "Non infectious complications: Shock",
  "Non infectious complications: Stroke: intracerebral hemorrhage",
  "Non infectious complications: Stroke: ischemic stroke",
  "Non infectious complications: Transfusion",
  "Non infectious complications: Other: physical impairment",
  "Non infectious complications: Other: cognitive impairment",
  "Non infectious complications: Other non infectious complications"
)
complication_free_text_col <- "Non infectious complications: Please, specify other non infectious complications"

complication_flags <- Main_data_hospitalDischarge %>%
  transmute(across(all_of(complication_cols), ~ case_when(
    .x == "YES" ~ 1L,
    .x == "NO" ~ 0L,
    TRUE ~ NA_integer_ # UNK or blank/NA -- genuinely unknown, not assumed absent
  )))
complication_free_text_flag <- if_else(
  !is.na(Main_data_hospitalDischarge[[complication_free_text_col]]) &
    trimws(Main_data_hospitalDischarge[[complication_free_text_col]]) != "",
  1L, 0L
)
# A patient is only "fully unknown" if every one of the 20 checkbox fields
# is NA -- a blank free-text field contributes no positive evidence but
# isn't itself missing data (it just means nothing extra was specified).
complication_all_unknown <- rowSums(!is.na(complication_flags)) == 0

Main_data_hospitalDischarge_clean <- Main_data_hospitalDischarge %>%
  mutate(
    patient_id = `Unique survey patient ID #1`,
    matched_control = is_true_matched_control(patient_id),
    matched_control_broad = is_cohort3_enrolled(patient_id),
    ICU_entry = case_when(
      .data[[discharge_icu_stayed_col]] == "true" ~ 1L,
      .data[[discharge_icu_stayed_col]] == "false" ~ 0L,
      TRUE ~ NA_integer_
    ),
    date_ICU_admission = parse_dt(.data[[discharge_icu_adm_col]]),
    date_ICU_discharge = parse_dt(.data[[discharge_icu_disc_col]]),
    ICU_LOS_days = suppressWarnings(as.numeric(.data[[discharge_icu_los_col]])),
    ICU_LOS_days = if_else(ICU_LOS_days < 0, NA_real_, ICU_LOS_days),
    complications_count = if_else(
      complication_all_unknown,
      NA_real_,
      rowSums(complication_flags, na.rm = TRUE) + complication_free_text_flag
    ),
    any_complication = if_else(complication_all_unknown, NA_integer_, if_else(complications_count > 0, 1L, 0L))
  ) %>%
  mutate(ICU_LOS_days = pmin(ICU_LOS_days, quantile(ICU_LOS_days, probs = 0.99, na.rm = TRUE))) %>%
  select(patient_id, matched_control, matched_control_broad, ICU_entry, date_ICU_admission, date_ICU_discharge, ICU_LOS_days,
         any_complication, complications_count)

rm(complication_flags, complication_free_text_flag, complication_all_unknown)





# =============================================================================
# SECTION 7: Main_data_hospitalCohortEnrolment
# =============================================================================
# CohortEnrolment (644 rows, 92 columns) is a MUCH shorter form than
# CaseReport (9247 rows, 243 columns): it captures baseline demographics,
# comorbidities, primary-admission-reason and risk factors, but NONE of the
# BSI-specific clinical workup (no severity signs, no blood culture
# sample/date fields, no invasive-devices/healthcare-exposure checklist, no
# suspected infection site, no Pitt score components, no antibiotic AST or
# free-text history). This is expected -- these are the survey's matched
# controls, who were never worked up for suspected BSI.
#
# Every variable below is built to be DIRECTLY COMPARABLE to its
# Main_data_CaseReport namesake: same variable name, same 1/0/NA coding,
# same derivation logic (age_group excepted -- see its own note). Where
# CohortEnrolment's source column wording differs from CaseReport's (most
# of them do, even for the same underlying question), the difference is
# resolved via an explicit alias below, verified 100% matched against the
# actual column names before use.
#
# NOT reproduced here because CohortEnrolment has no corresponding source
# field at all: yesno_map (prior healthcare exposure), device_map, site_map,
# integrated_map/agespecific_map (severity signs), date_map (BSI onset/
# admission/blood-sample dates), sample_id_map, occur48h_map, abx_hist_cols,
# Pitt score, infection_origin_hosp/hosp_comm_acquis48/infection_origin_composite,
# transfer_date_F1/days_enrolment_to_transfer.
#
# FLAGGED (not fixed here -- out of scope for this section, changes
# SECTION 1): two fields -- "Immunosuppression in the 48 hours prior to
# blood sample collection: Cytotoxic chemotherapy <=6 months" and
# "...Prednisone therapy >=10 mg/day" -- exist with IDENTICAL names in BOTH
# CaseReport and CohortEnrolment, but SECTION 1 does not currently code
# either one into Main_data_CaseReport. They're captured below as
# cohort-only bonus variables (immunosup48h_chemo/immunosup48h_prednisone);
# if wanted as genuinely comparable variables, SECTION 1 would need the
# same two columns added to its own core_cols/mapping.

## --- alias tables: cohort's source wording -> the SAME `new` names SECTION 1 uses ---
comorb_map_cohort <- comorb_map %>%
  mutate(orig_cohort = sub(
    "^Comorbidities at the time of hospital admission: ", "Comorbidities: ", orig
  ))
comorb_map_cohort$orig_cohort[comorb_map_cohort$new == "comorb_connective_tissue"] <-
  "Comorbidities: Connective tissue, rheumatologic, autoimmune"
comorb_map_cohort$orig_cohort[comorb_map_cohort$new == "comorb_metastatic_tumour"] <-
  "Comorbidities: Metastatic solid tumour"
comorb_map_cohort$orig_cohort[comorb_map_cohort$new == "comorb_hsct"] <-
  "Comorbidities: Hematopoietic stem-cell transplantation"

reason_map_cohort <- reason_map %>% mutate(orig_cohort = orig)
reason_map_cohort$orig_cohort[reason_map_cohort$new == "reason_connective_tissue"] <-
  "Primary reason for admission to hospital (diagnosis on admission):: Connective tissue / rheumatologic"
reason_map_cohort$orig_cohort[reason_map_cohort$new == "reason_unknown"] <-
  "Primary reason for admission to hospital (diagnosis on admission):: Unknown primary reason"

riskfx_map_cohort <- riskfx_map %>% mutate(orig_cohort = orig) # exact match, no aliasing needed

cohort_abx_current_status_col <- "Antimicrobial use history: Is the patient currently on antibiotics?"

stopifnot(
  all(comorb_map_cohort$orig_cohort %in% names(Main_data_hospitalCohortEnrolment)),
  all(reason_map_cohort$orig_cohort %in% names(Main_data_hospitalCohortEnrolment)),
  all(riskfx_map_cohort$orig_cohort %in% names(Main_data_hospitalCohortEnrolment)),
  cohort_abx_current_status_col %in% names(Main_data_hospitalCohortEnrolment)
)

## --- speciality_code: cohort uses a DIFFERENT code convention than
## ward_code_lookup (SECTION 1) -- specialty-first with "_" ("CAR_AICU"), or
## a bare ward-type/specialty token ("AICU", "NICU", "PSW") -- instead of
## CaseReport's "WARDTYPE-SPECIALTY" with "-" ("AICU-CAR"). Reverse-and-join
## for the "_" form; for a bare token, self-pair it if it's already a full
## ward-type code (AICU/AMW/ASW/PMW/NMW -> e.g. "ASW-ASW"), else prefix it
## with the ward-type its own initial implies (N.. -> "NMW-..", P.. ->
## "PMW-.."; both conventions observed in the Global-PPS scheme this form
## reuses). Verified: 641/641 non-blank cohort codes resolve to an existing
## ward_code_lookup row under this transform (3 genuinely blank codes -> NA).
transform_cohort_ward_code <- function(x) {
  ward_types <- c("AICU", "AMW", "ASW", "PMW", "NMW")
  out <- rep(NA_character_, length(x))
  has_us <- !is.na(x) & grepl("_", x)
  parts <- strsplit(x[has_us], "_")
  out[has_us] <- vapply(parts, function(p) paste0(p[2], "-", p[1]), character(1))
  no_us <- !is.na(x) & !grepl("_", x)
  bare <- x[no_us]
  bare_out <- rep(NA_character_, length(bare))
  is_wardtype <- bare %in% ward_types
  bare_out[is_wardtype] <- paste0(bare[is_wardtype], "-", bare[is_wardtype])
  starts_n <- !is_wardtype & grepl("^N", bare)
  bare_out[starts_n] <- paste0("NMW-", bare[starts_n])
  starts_p <- !is_wardtype & grepl("^P", bare)
  bare_out[starts_p] <- paste0("PMW-", bare[starts_p])
  out[no_us] <- bare_out
  out
}

## --- bonus cohort-only fields (no SECTION 1 equivalent -- see note above) ---
s4_map <- tribble(
  ~orig, ~new,
  "S4 Risk factors and comorbidities: Cardiac arrest", "s4_cardiac_arrest",
  "S4 Risk factors and comorbidities: Active malignancy", "s4_active_malignancy",
  "S4 Risk factors and comorbidities: Intravenous vasopressor agent required", "s4_iv_vasopressor",
  "S4 Risk factors and comorbidities: Mechanical ventilation needed", "s4_mech_vent",
  "S4 Risk factors and comorbidities: Other immunosuppression", "s4_other_immunosuppression",
  "S4 Risk factors and comorbidities: Solid organ transplantation", "s4_solid_organ_transplant",
  "S4 Risk factors and comorbidities: Child C cirrhosis", "s4_child_c_cirrhosis",
  "S4 Risk factors and comorbidities: Hematopoietic stem-cell transplantation", "s4_hsct",
  "S4 Risk factors and comorbidities: Neutropenia ≤500/mm3", "s4_neutropenia",
  "S4 Risk factors and comorbidities: End-stage renal disease requiring dialysis", "s4_esrd_dialysis",
  "S4 Risk factors and comorbidities: Insulin-dependent diabetes mellitus", "s4_insulin_dependent_dm",
  "S4 Risk factors and comorbidities: HIV (CD4 ≤ 200/mm3)", "s4_hiv_cd4_low"
)
stopifnot(all(s4_map$orig %in% names(Main_data_hospitalCohortEnrolment)))

s4_flags <- Main_data_hospitalCohortEnrolment %>%
  transmute(across(all_of(s4_map$orig), ~ case_when(
    toupper(trimws(.x)) == "NO" ~ 0L,
    toupper(trimws(.x)) == "YES" ~ 1L,
    TRUE ~ NA_integer_
  ))) %>%
  rename(!!!setNames(s4_map$orig, s4_map$new))

comorb_flags_cohort <- Main_data_hospitalCohortEnrolment %>%
  transmute(across(all_of(comorb_map_cohort$orig_cohort), ~ case_when(
    toupper(trimws(.x)) == "NO" ~ 0L,
    toupper(trimws(.x)) == "YES" ~ 1L,
    TRUE ~ NA_integer_
  ))) %>%
  rename(!!!setNames(comorb_map_cohort$orig_cohort, comorb_map_cohort$new))

riskfx_flags_cohort <- Main_data_hospitalCohortEnrolment %>%
  transmute(across(all_of(riskfx_map_cohort$orig_cohort), ~ case_when(
    toupper(trimws(.x)) == "NO" ~ 0L,
    toupper(trimws(.x)) == "YES" ~ 1L,
    TRUE ~ NA_integer_
  ))) %>%
  rename(!!!setNames(riskfx_map_cohort$orig_cohort, riskfx_map_cohort$new))

reason_flags_cohort <- Main_data_hospitalCohortEnrolment %>%
  transmute(across(all_of(reason_map_cohort$orig_cohort), ~ if_else(tolower(trimws(.x)) == "true", 1L, 0L, missing = 0L))) %>%
  rename(!!!setNames(reason_map_cohort$orig_cohort, reason_map_cohort$new))

Main_data_hospitalCohortEnrolment_clean <- Main_data_hospitalCohortEnrolment %>%
  transmute(
    patient_id = `Unique survey patient ID (___/C3)`,
    # Previously matched_control = 1L unconditionally for every row in this
    # sheet (i.e. "enrolled in CohortEnrolment" was treated as sufficient on
    # its own). Corrected to the same is_true_matched_control() criterion
    # used in SECTIONS 5/6: only TRUE if the ID ends "/C3" AND its base ID
    # identifies a confirmed-BSI case. Of 644 CohortEnrolment records (before
    # this fix), only 466 satisfy this -- the remaining ~178 have no base-ID
    # match to a confirmed-BSI CaseReport patient (most plausibly enrolled
    # against a suspected-BSI case that was not ultimately microbiologically
    # confirmed) and are consequently NOT retained as matched controls.
    matched_control = is_true_matched_control(patient_id),
    matched_control_broad = is_cohort3_enrolled(patient_id),
    facility_id = `Unique facility ID`,
    ward_id = `Unique ward/department ID identifying BSI and enrolling the patient`,
    speciality_code = transform_cohort_ward_code(`Specialty code of admitting ward that is applicable to the patient`),
    enrolment_date = parse_dt(`Enrolment date`),
    dob = parse_dt(`Date of Birth (dd/mm/yyyy)`),
    age_years = round(as.numeric(difftime(Sys.time(), dob, units = "days")) / 365.25, 2),
    # Same implausible-DOB sanity filter as SECTION 1's CaseReport age_years
    # (mistyped birth years producing multi-hundred-year or negative ages).
    age_years = if_else(age_years < -0.5 | age_years > 110, NA_real_, age_years),
    # age_group: NOT sourced from a direct "Age group" field (CohortEnrolment
    # has none) -- derived using the SAME NEONATE(<=28 days)/PAEDIATRIC(29
    # days-<18y)/ADULT(>=18y) cutoffs the form's own age-domain labels use
    # elsewhere, unlike CaseReport's age_group which comes straight from a
    # pre-computed source field. Same 0/1/3 coding. Deliberately uses age AT
    # ENROLMENT (enrolment_date - dob), NOT age_years (current age as of
    # whenever this script is run) -- using age_years here would silently
    # reclassify every neonate/young infant enrolled in 2024-2025 as older
    # once enough real time has passed (verified: doing it that way produces
    # ZERO neonates). age_years itself is kept as "current age" only,
    # matching SECTION 1's own established (if similarly quirky) convention
    # for CaseReport's age_years -- it is not used for any classification.
    age_group = case_when(
      is.na(dob) | is.na(enrolment_date) ~ NA_integer_,
      as.numeric(difftime(enrolment_date, dob, units = "days")) <= 28 ~ 0L,
      as.numeric(difftime(enrolment_date, dob, units = "days")) < 18 * 365.25 ~ 1L,
      TRUE ~ 3L
    ),
    sex_at_birth = case_when(
      toupper(trimws(`Sex at birth`)) == "FEMALE" ~ 1L,
      toupper(trimws(`Sex at birth`)) == "MALE" ~ 0L,
      TRUE ~ NA_integer_
    ),
    # Unrecorded status coded 0 (not on antibiotics), not NA -- see the
    # matching comment at abx_current_therapy's first definition (SECTION 1).
    abx_current_therapy = case_when(
      tolower(trimws(.data[[cohort_abx_current_status_col]])) == "true" ~ 1L,
      tolower(trimws(.data[[cohort_abx_current_status_col]])) == "false" ~ 0L,
      TRUE ~ 0L
    ),
    # bonus, cohort-only (see note above)
    immunosup48h_chemo = case_when(
      toupper(trimws(`Immunosuppression in the 48 hours prior to blood sample collection: Cytotoxic chemotherapy ≤ 6 months`)) == "NO" ~ 0L,
      toupper(trimws(`Immunosuppression in the 48 hours prior to blood sample collection: Cytotoxic chemotherapy ≤ 6 months`)) == "YES" ~ 1L,
      TRUE ~ NA_integer_
    ),
    immunosup48h_prednisone = case_when(
      toupper(trimws(`Immunosuppression in the 48 hours prior to blood sample collection: Prednisone therapy ≥ 10 mg/day`)) == "NO" ~ 0L,
      toupper(trimws(`Immunosuppression in the 48 hours prior to blood sample collection: Prednisone therapy ≥ 10 mg/day`)) == "YES" ~ 1L,
      TRUE ~ NA_integer_
    ),
    s4_n_hospital_episodes = suppressWarnings(as.integer(`S4 Risk factors and comorbidities: Number of hospital episodes`))
  ) %>%
  left_join(ward_code_lookup, by = c("speciality_code" = "code")) %>%
  bind_cols(comorb_flags_cohort, riskfx_flags_cohort, reason_flags_cohort, s4_flags)

rm(comorb_map_cohort, reason_map_cohort, riskfx_map_cohort, comorb_flags_cohort,
   riskfx_flags_cohort, reason_flags_cohort, s4_flags, s4_map, cohort_abx_current_status_col)




# =============================================================================
# SECTION 8: Merging datasets
# =============================================================================
# Combines all four patient-level datasets by patient_id into one wide
# dataset, Main_data_Combined -- one row per distinct patient across the
# UNION of everyone who appears in any of the four (verified: 9,912 distinct
# patients from 9,247 CaseReport + 1,275 Discharge + 1,261 FollowUp28 + 644
# CohortEnrolment). Structurally, Discharge/FollowUp28 are outcome-tracking
# forms covering roughly HALF case-patients and half matched controls each
# (~630/1270 found in CaseReport, ~604/1270 found in CohortEnrolment,
# confirming the "CaseReport only covers about half" observation);
# CohortEnrolment (matched controls) is NEVER found in CaseReport at all
# (verified 0/644) -- controls were never worked up for suspected BSI.
#
# COLUMN-NAME COLLISIONS: SECTION 7 deliberately gave
# Main_data_hospitalCohortEnrolment_clean the same 57 column names as
# Main_data_CaseReport (facility_id, dob, comorb_*, reason_*, riskfx_*,
# sex_at_birth, age_group, etc.) specifically so the two cohorts could be
# compared. Because a patient can only ever come from ONE side of that pair
# (verified above), those 57 are COALESCED after the join into one clean
# column each, rather than left as duplicated .x/.y pairs. matched_control
# itself is shared by Discharge/FollowUp28/CohortEnrolment (all three
# outcome-tracking sources) and is reconciled the same way, defaulting to 0
# for the CaseReport-only patients (who are, by construction, never
# controls).
#
# RESISTANCE METRICS: per explicit instruction, the named resistance/
# pathogen-found metrics get 0 (not NA) specifically for patients who could
# NOT be matched to CaseReport/lab data at all -- distinct from a patient
# who IS in CaseReport but tested BSI-negative (those already carry their
# own correct value: NA for the raw isolate-level flags -- never tested --
# and 0 for the population-level *_BSI_positive flags, per SECTION 3.6).
# Every OTHER CaseReport-sourced column (severity signs, Pitt score,
# antibiotic history, etc.) is left as a genuine NA when unmatched, per
# instruction.
resistance_metric_cols <- c(
  "resistance_any", "mdr_resistance", "xdr_resistance", "pdr_resistance",
  "carb_resistance", "methicilin_resistance", "cephalosporin_resis", "3GC_resistance", "vancomycin_resistance",
  "resistant_BSI_positive", "mdr_BSI_bacterial", "xdr_BSI_bacterial", "pdr_BSI_bacterial",
  "carb_BSI_positive", "methicillin_BSI_positive", "cephalosporin_BSI_positive", "c3GC_BSI_positive", "vancomycin_BSI_positive",
  "BSI_positive", "lab_match"
)
stopifnot(all(resistance_metric_cols %in% names(Main_data_CaseReport)))

overlap_cols <- intersect(names(Main_data_CaseReport), names(Main_data_hospitalCohortEnrolment_clean))
overlap_cols <- setdiff(overlap_cols, "patient_id")

## --- de-duplicate each source by patient_id (keep first row) ---
# A small number of patient_ids repeat within a single source (20 in
# CaseReport -- including one obvious junk row, "DATABASE TEST", not a real
# patient -- 10 in Discharge, 8 in FollowUp28, 0 in CohortEnrolment).
# Left un-deduplicated, a full_join on a non-unique key fans out into extra
# rows rather than merging cleanly; kept first occurrence per patient_id.
dedupe_by_patient <- function(df) df %>% distinct(patient_id, .keep_all = TRUE)

discharge_dedup <- dedupe_by_patient(Main_data_hospitalDischarge_clean)
followup_dedup <- dedupe_by_patient(Main_data_hospitalFollowup28_clean)
cohort_dedup <- dedupe_by_patient(Main_data_hospitalCohortEnrolment_clean)
casereport_dedup <- dedupe_by_patient(Main_data_CaseReport) %>% mutate(.in_casereport = TRUE)

## --- reconcile matched_control (and its broad IPW-population counterpart,
## matched_control_broad) across the 3 sources that carry them ---
matched_control_roster <- bind_rows(
  discharge_dedup %>% select(patient_id, matched_control, matched_control_broad),
  followup_dedup %>% select(patient_id, matched_control, matched_control_broad),
  cohort_dedup %>% select(patient_id, matched_control, matched_control_broad)
) %>%
  group_by(patient_id) %>%
  summarise(
    matched_control = max(matched_control, na.rm = TRUE),
    matched_control_broad = max(matched_control_broad, na.rm = TRUE),
    .groups = "drop"
  )

discharge_dedup <- discharge_dedup %>% select(-matched_control, -matched_control_broad)
followup_dedup <- followup_dedup %>% select(-matched_control, -matched_control_broad)
cohort_dedup <- cohort_dedup %>% select(-matched_control, -matched_control_broad)

## --- combine the 3 non-CaseReport sources (no overlapping columns among
## them besides patient_id/matched_control, already pulled out above) ---
combined_outcomes <- discharge_dedup %>%
  full_join(followup_dedup, by = "patient_id") %>%
  full_join(cohort_dedup, by = "patient_id")

## --- full join CaseReport (all 287 columns), suffixing the 57 overlapping ones ---
Main_data_Combined <- combined_outcomes %>%
  full_join(casereport_dedup, by = "patient_id", suffix = c("", "_case"))

for (cn in overlap_cols) {
  case_col <- paste0(cn, "_case")
  if (case_col %in% names(Main_data_Combined)) {
    Main_data_Combined[[cn]] <- coalesce(Main_data_Combined[[cn]], Main_data_Combined[[case_col]])
    Main_data_Combined[[case_col]] <- NULL
  }
}

## --- matched_control (+ matched_control_broad): reconciled roster, 0 for
## CaseReport-only patients ---
Main_data_Combined <- Main_data_Combined %>%
  left_join(matched_control_roster, by = "patient_id") %>%
  mutate(
    matched_control = if_else(is.na(matched_control), 0L, as.integer(matched_control)),
    matched_control_broad = if_else(is.na(matched_control_broad), 0L, as.integer(matched_control_broad))
  ) %>%
  relocate(matched_control, matched_control_broad, .after = patient_id)

## --- resistance metrics: 0 (not NA) specifically where patient_id has no
## CaseReport/lab match at all ---
Main_data_Combined <- Main_data_Combined %>%
  mutate(across(all_of(resistance_metric_cols), ~ if_else(is.na(.in_casereport) & is.na(.x), 0L, .x))) %>%
  select(-.in_casereport)

## --- restrict to patients with outcome data (Discharge and/or FollowUp28) ---
# Drops the CaseReport-only patients who have no Discharge/FollowUp28 record
# at all -- no mortality/LOS/complications data exists for them, so they
# can't contribute to any outcome analysis. Verified: keeps 1,302 of 9,912
# (union of Discharge's 1,265 and FollowUp28's 1,253 distinct patient_ids,
# which overlap heavily -- 1,216 in both); removes 8,610, of whom 8,585 are
# CaseReport patients (matched_control == 0) with no outcome form, matching
# the ~8k figure. The other 25 removed are CohortEnrolment (matched_control
# == 1) patients who were enrolled as controls but have no matching
# Discharge/FollowUp28 record either -- dropped for the same reason (no
# outcome data to analyse), not kept just for being a control.
outcome_patient_ids <- union(discharge_dedup$patient_id, followup_dedup$patient_id)
Main_data_Combined <- Main_data_Combined %>%
  filter(patient_id %in% outcome_patient_ids)

rm(discharge_dedup, followup_dedup, cohort_dedup, casereport_dedup, combined_outcomes,
   matched_control_roster, outcome_patient_ids)



save.image("ws_after_sec8.RData")
# =============================================================================
# SECTION 9: Mortality estimates
# =============================================================================
# Estimates 28-day all-cause mortality differences between drug-susceptible
# BSI, drug-resistant BSI, and matched uninfected controls (Main_data_Combined,
# SECTION 8), plus in-hospital mortality with discharge as a competing risk.
# Method: marginal structural Cox models with stabilised inverse-probability
# (of exposure) weights, cluster-robust (hospital) standard errors, pooled
# across multiple-imputed datasets via Rubin's rules; adjusted 28-day risks
# and risk differences from IPW-weighted Kaplan-Meier curves; uncertainty on
# risk differences and national excess deaths via cluster+imputation
# bootstrap. All three pre-specified comparisons throughout:
#   (1) Susceptible BSI    vs Uninfected control
#   (2) Resistant BSI      vs Uninfected control
#   (3) Resistant BSI      vs Susceptible BSI
# Drug resistance = resistant_BSI_positive (resistance to >=1 tested agent,
# organism-agnostic, SECTION 3.6) -- the "at least one tested antibiotic"
# definition specified, not MDR (kept as a secondary/descriptive cut, not a
# separate MSM arm, given only 66 XDR and 7 PDR cases in the whole survey).

# -----------------------------------------------------------------------
# 9.1: Group definition
# -----------------------------------------------------------------------
# Uninfected control = matched_control==1 (BSI_positive is always 0 for
# these, SECTION 8). Susceptible/Resistant BSI = matched_control==0,
# confirmed BSI, with a determinate resistant_BSI_positive. Enrolled
# patients who were never confirmed BSI-positive (matched_control==0 &
# BSI_positive==0, n=186) and BSI-positive patients whose resistance could
# not be determined (AST too incomplete, n=70) are excluded from this
# 3-group comparison -- they remain in Main_data_Combined for other uses.
Main_data_Combined <- Main_data_Combined %>%
  mutate(
    mortality_group = case_when(
      matched_control_broad == 1 ~ "Uninfected control",
      matched_control_broad == 0 & BSI_positive == 1 & resistant_BSI_positive == 1 ~ "Resistant BSI",
      matched_control_broad == 0 & BSI_positive == 1 & resistant_BSI_positive == 0 ~ "Susceptible BSI",
      TRUE ~ NA_character_
    ),
    mortality_group = factor(mortality_group, levels = c("Uninfected control", "Susceptible BSI", "Resistant BSI"))
  )

# -----------------------------------------------------------------------
# 9.2: Survival time construction
# -----------------------------------------------------------------------
# 28-day outcome: time28/event28. event28 = mortality30day (NA for UNABLE,
# excluded). time28 is deliberately NOT computed from the raw "date of
# 28-day review" minus admission date for patients alive at day 28 -- that
# difference is implausible for 63/1106 rows (the same kind of garbled
# source dates already documented in SECTIONS 5-8), whereas the 28-day
# review is BY PROTOCOL DESIGN exactly 28 days after enrolment. Alive
# patients are therefore administratively censored at exactly day 28.
# In-hospital deaths use LOS_days (already cleaned/winsorised, SECTION 5)
# as their event time. Post-discharge deaths within 28 days use the
# post-discharge death date where plausible (0-28 days after admission);
# the 1 patient where that date is missing/implausible is conservatively
# assigned day 28 (a lower-bound approximation on their true hazard, not an
# assumption of higher risk).
followup_raw_death_date_col <- "Status at 28-days review: Date of death (if post-discharge)"
followup_raw_review_date_col <- "Follow-up details: Date of 28-days review"
followup_extra_dates <- tibble(
  patient_id = Main_data_hospitalFollowup28[["Unique survey patient ID #1"]],
  date_of_death_postdischarge = parse_dt(Main_data_hospitalFollowup28[[followup_raw_death_date_col]])
) %>% distinct(patient_id, .keep_all = TRUE)

Main_data_Combined <- Main_data_Combined %>%
  left_join(followup_extra_dates, by = "patient_id") %>%
  mutate(
    event28 = if_else(mortality30day == 1, 1L, if_else(mortality30day == 0, 0L, NA_integer_)),
    .t_death_postdischarge = as.numeric(difftime(date_of_death_postdischarge, Date_of_admission5enro, units = "days")),
    time28 = case_when(
      event28 == 1 & in_hospital_mort == 1 ~ pmin(pmax(LOS_days, 0.5), 28),
      event28 == 1 & (in_hospital_mort == 0 | is.na(in_hospital_mort)) &
        !is.na(.t_death_postdischarge) & .t_death_postdischarge >= 0 & .t_death_postdischarge <= 28 ~ pmax(.t_death_postdischarge, 0.5),
      event28 == 1 ~ 28,
      event28 == 0 ~ 28,
      TRUE ~ NA_real_
    )
  ) %>%
  select(-.t_death_postdischarge)

# In-hospital competing-risk outcome: time_inhosp (= LOS_days) / event_inhosp
# ("death" = died in hospital, "discharge" = left the hospital alive -- the
# competing event that precludes further in-hospital death risk, REGARDLESS
# of what happens after discharge; a post-discharge death is, for this
# specific in-hospital-stay competing-risk framing, preceded by a discharge
# event). NA (in_hospital_mort undetermined) excluded.
Main_data_Combined <- Main_data_Combined %>%
  mutate(
    time_inhosp = pmin(pmax(LOS_days, 0.5), 28),
    event_inhosp = factor(
      case_when(in_hospital_mort == 1 ~ "death", in_hospital_mort == 0 ~ "discharge", TRUE ~ NA_character_),
      levels = c("discharge", "death")
    )
  )

# -----------------------------------------------------------------------
# facility_id backfill for Uninfected/Cohort-3 patients whose CohortEnrolment
# row exists under a DIFFERENT ID suffix than the ID their Discharge/
# FollowUp28 outcome record carries -- the raw join above matches on the
# literal patient_id string, so it misses these even though CohortEnrolment
# itself is 100% complete for facility_id (verified directly against the
# source sheet). Same base-ID normalisation already used by
# is_true_matched_control()/is_cohort3_enrolled() (Cyrillic С vs Latin C,
# "/S1" vs "/C3"), applied here to recover facility_id specifically, since
# those two functions only ever recovered the matched-control FLAG, never
# backfilled this column. Plus one further, single hand-verified correction:
# "HF009/W07/22" (Discharge/FollowUp28) is a zero-padding typo for
# CohortEnrolment's "HF009/W07/022" -- corrected explicitly as a one-off
# remap, not a general zero-padding rule (no other ID in this dataset is
# known to have this issue, and a general rule risks silently merging
# unrelated patients whose IDs happen to differ only by padding).
# Confirmed after this fix: 13 of 29 previously-NA facility_id patients
# recovered via suffix normalisation, 1 more via the padding correction,
# leaving 15 genuinely absent from CohortEnrolment under any ID variant.
normalize_base_id <- function(id) {
  id <- str_replace(id, "/С3$", "/C3")
  str_remove(id, "/(C3|S1)$")
}
cohort_facility_lookup <- Main_data_hospitalCohortEnrolment_clean %>%
  filter(!is.na(facility_id)) %>%
  mutate(base_id = normalize_base_id(patient_id)) %>%
  mutate(base_id = if_else(base_id == "HF009/W07/022", "HF009/W07/22", base_id)) %>%
  distinct(base_id, .keep_all = TRUE) %>%
  select(base_id, facility_id_backfill = facility_id)

Main_data_Combined <- Main_data_Combined %>%
  mutate(.base_id = normalize_base_id(patient_id)) %>%
  left_join(cohort_facility_lookup, by = c(".base_id" = "base_id")) %>%
  mutate(facility_id = coalesce(facility_id, facility_id_backfill)) %>%
  select(-.base_id, -facility_id_backfill)

mortality_analytic <- Main_data_Combined %>% filter(!is.na(mortality_group))

# =========================================================================
# SUPPLEMENTARY: calibrated-design-weighted crude mortality + KM curves,
# mirroring the main-text Figure 3 (crude bar chart) and Figure 4 (KM
# curves) but weighting every estimate by the survey's calibrated design
# weight (SECTION 3.3) instead of reporting raw/unweighted counts. This is
# a fast, MICE-free descriptive companion -- calibrated_weight is a
# facility-level property (valid for cases and controls alike, SECTION
# 9.13's own comment), so every group frame below just needs mortality_
# group/event28/facility_id plus this one joined column. Exits immediately
# after saving, skipping the (expensive) MICE/bootstrap sections below.
cat("=== CALIBRATED-WEIGHTED CRUDE MORTALITY + KM (fast, no MICE) ===\n")
mortality_analytic_cw <- mortality_analytic %>%
  select(-any_of("calibrated_weight")) %>%
  left_join(facility_weights %>% select(UniqueID, calibrated_weight), by = c(facility_id = "UniqueID"))

mortality_cw_mdr <- mortality_analytic_cw %>%
  mutate(mortality_group = case_when(
    matched_control_broad == 1 ~ "Uninfected control",
    matched_control_broad == 0 & mdr_BSI_bacterial == 1 ~ "MDR BSI",
    matched_control_broad == 0 & mdr_BSI_bacterial == 0 ~ "Non-MDR BSI",
    TRUE ~ NA_character_
  ), mortality_group = factor(mortality_group, levels = c("Uninfected control", "Non-MDR BSI", "MDR BSI"))) %>%
  filter(!is.na(mortality_group))

mortality_cw_xdr <- mortality_analytic_cw %>%
  mutate(mortality_group = case_when(
    matched_control_broad == 1 ~ "Uninfected control",
    matched_control_broad == 0 & xdr_BSI_bacterial == 1 ~ "XDR BSI",
    matched_control_broad == 0 & xdr_BSI_bacterial == 0 ~ "Non-XDR BSI",
    TRUE ~ NA_character_
  ), mortality_group = factor(mortality_group, levels = c("Uninfected control", "Non-XDR BSI", "XDR BSI"))) %>%
  filter(!is.na(mortality_group))

mortality_cw_carb <- mortality_analytic_cw %>%
  mutate(mortality_group = case_when(
    matched_control_broad == 1 ~ "Uninfected control",
    matched_control_broad == 0 & is_cre_pa_ab(major_bacterial_group, pathogen1_clean) & carb_BSI_positive == 1 ~ "Carbapenem-R Gram-negative BSI",
    matched_control_broad == 0 & is_cre_pa_ab(major_bacterial_group, pathogen1_clean) & carb_BSI_positive == 0 ~ "Carbapenem-S Gram-negative BSI",
    TRUE ~ NA_character_
  ), mortality_group = factor(mortality_group, levels = c("Uninfected control", "Carbapenem-S Gram-negative BSI", "Carbapenem-R Gram-negative BSI"))) %>%
  filter(!is.na(mortality_group))

mortality_cw_ceph <- mortality_analytic_cw %>%
  mutate(mortality_group = case_when(
    matched_control_broad == 1 ~ "Uninfected control",
    matched_control_broad == 0 & is_cre_pa_ab(major_bacterial_group, pathogen1_clean) & cephalosporin_BSI_positive == 1 ~ "Cephalosporin-R Gram-negative BSI",
    matched_control_broad == 0 & is_cre_pa_ab(major_bacterial_group, pathogen1_clean) & cephalosporin_BSI_positive == 0 ~ "Cephalosporin-S Gram-negative BSI",
    TRUE ~ NA_character_
  ), mortality_group = factor(mortality_group, levels = c("Uninfected control", "Cephalosporin-S Gram-negative BSI", "Cephalosporin-R Gram-negative BSI"))) %>%
  filter(!is.na(mortality_group))

mortality_cw_3gcr <- mortality_analytic_cw %>%
  mutate(mortality_group = case_when(
    matched_control_broad == 1 ~ "Uninfected control",
    matched_control_broad == 0 & is_cre_pa_ab(major_bacterial_group, pathogen1_clean) & c3GC_BSI_positive == 1 ~ "3GCR Gram-negative BSI",
    matched_control_broad == 0 & is_cre_pa_ab(major_bacterial_group, pathogen1_clean) & c3GC_BSI_positive == 0 ~ "3GC-susceptible Gram-negative BSI",
    TRUE ~ NA_character_
  ), mortality_group = factor(mortality_group, levels = c("Uninfected control", "3GC-susceptible Gram-negative BSI", "3GCR Gram-negative BSI"))) %>%
  filter(!is.na(mortality_group))

mortality_cw_mrsa <- mortality_analytic_cw %>%
  mutate(mortality_group = case_when(
    matched_control_broad == 1 ~ "Uninfected control",
    matched_control_broad == 0 & pathogen1_clean == "Staphylococcus aureus" & methicillin_BSI_positive == 1 ~ "MRSA BSI",
    TRUE ~ NA_character_
  ), mortality_group = factor(mortality_group, levels = c("Uninfected control", "MRSA BSI"))) %>%
  filter(!is.na(mortality_group))

mortality_cw_access <- mortality_analytic_cw %>%
  mutate(mortality_group = case_when(
    matched_control_broad == 1 ~ "Uninfected control",
    matched_control_broad == 0 & access_BSI_positive == 1 ~ "Access-R BSI",
    matched_control_broad == 0 & access_BSI_positive == 0 ~ "Access-S BSI",
    TRUE ~ NA_character_
  ), mortality_group = factor(mortality_group, levels = c("Uninfected control", "Access-S BSI", "Access-R BSI"))) %>%
  filter(!is.na(mortality_group))

mortality_cw_watch <- mortality_analytic_cw %>%
  mutate(mortality_group = case_when(
    matched_control_broad == 1 ~ "Uninfected control",
    matched_control_broad == 0 & watch_BSI_positive == 1 ~ "Watch-R BSI",
    matched_control_broad == 0 & watch_BSI_positive == 0 ~ "Watch-S BSI",
    TRUE ~ NA_character_
  ), mortality_group = factor(mortality_group, levels = c("Uninfected control", "Watch-S BSI", "Watch-R BSI"))) %>%
  filter(!is.na(mortality_group))

mortality_cw_reserve <- mortality_analytic_cw %>%
  mutate(mortality_group = case_when(
    matched_control_broad == 1 ~ "Uninfected control",
    matched_control_broad == 0 & reserve_BSI_positive == 1 ~ "Reserve-R BSI",
    matched_control_broad == 0 & reserve_BSI_positive == 0 ~ "Reserve-S BSI",
    TRUE ~ NA_character_
  ), mortality_group = factor(mortality_group, levels = c("Uninfected control", "Reserve-S BSI", "Reserve-R BSI"))) %>%
  filter(!is.na(mortality_group))

mssa_cw <- mortality_analytic_cw %>%
  filter(matched_control == 0, pathogen1_clean == "Staphylococcus aureus", methicillin_BSI_positive == 0) %>%
  mutate(mortality_group = "MSSA")
vre_cw <- mortality_analytic_cw %>%
  filter(matched_control == 0, !is.na(genus1), genus1 == "Enterococcus", !is.na(vancomycin_BSI_positive)) %>%
  mutate(mortality_group = if_else(vancomycin_BSI_positive == 1, "VRE", "VSE"))

weighted_row <- function(data, group_label, display_label) {
  sub <- data %>% dplyr::filter(mortality_group == group_label, !is.na(event28), !is.na(calibrated_weight))
  n <- nrow(sub)
  deaths <- sum(sub$event28 == 1, na.rm = TRUE)
  out <- tryCatch({
    des <- svydesign(ids = ~facility_id, weights = ~calibrated_weight, data = sub)
    cp <- svyciprop(~I(event28 == 1), des, method = "logit")
    ci <- confint(cp)
    tibble(calibrated_pct = as.numeric(cp) * 100, calibrated_lower = ci[1] * 100, calibrated_upper = ci[2] * 100)
  }, error = function(e) tibble(calibrated_pct = 100 * weighted.mean(sub$event28 == 1, sub$calibrated_weight),
                                 calibrated_lower = NA_real_, calibrated_upper = NA_real_))
  tibble(group = display_label, n = n, deaths_28day = deaths) %>%
    bind_cols(round(out, 1))
}

weighted_crude_all <- bind_rows(
  weighted_row(mortality_analytic_cw, "Uninfected control", "Uninfected"),
  weighted_row(mortality_analytic_cw, "Susceptible BSI", "Susceptible (any)"),
  weighted_row(mortality_analytic_cw, "Resistant BSI", "Resistant (any)"),
  weighted_row(mortality_cw_mdr, "Non-MDR BSI", "Non-MDR"),
  weighted_row(mortality_cw_mdr, "MDR BSI", "MDR"),
  weighted_row(mortality_cw_xdr, "Non-XDR BSI", "Non-XDR"),
  weighted_row(mortality_cw_xdr, "XDR BSI", "XDR"),
  weighted_row(mortality_cw_carb, "Carbapenem-S Gram-negative BSI", "Carbapenem-S"),
  weighted_row(mortality_cw_carb, "Carbapenem-R Gram-negative BSI", "Carbapenem-R"),
  weighted_row(mortality_cw_ceph, "Cephalosporin-S Gram-negative BSI", "Cephalosporin-S"),
  weighted_row(mortality_cw_ceph, "Cephalosporin-R Gram-negative BSI", "Cephalosporin-R"),
  weighted_row(mortality_cw_3gcr, "3GC-susceptible Gram-negative BSI", "3GC-S"),
  weighted_row(mortality_cw_3gcr, "3GCR Gram-negative BSI", "3GCR"),
  weighted_row(mssa_cw, "MSSA", "MSSA"),
  weighted_row(mortality_cw_mrsa, "MRSA BSI", "MRSA"),
  weighted_row(vre_cw, "VSE", "VSE"),
  weighted_row(vre_cw, "VRE", "VRE"),
  weighted_row(mortality_cw_access, "Access-S BSI", "Access-S"),
  weighted_row(mortality_cw_access, "Access-R BSI", "Access-R"),
  weighted_row(mortality_cw_watch, "Watch-S BSI", "Watch-S"),
  weighted_row(mortality_cw_watch, "Watch-R BSI", "Watch-R"),
  weighted_row(mortality_cw_reserve, "Reserve-S BSI", "Reserve-S"),
  weighted_row(mortality_cw_reserve, "Reserve-R BSI", "Reserve-R")
)
print(as.data.frame(weighted_crude_all), row.names = FALSE)
write_xlsx(weighted_crude_all, "Mortality_crude_calibrated_weighted_all_mechanisms.xlsx")

# Design-based (Rao-Scott) chi-square for each resistant-vs-susceptible pair,
# calibrated-weighted analogue of the main Figure 3's unweighted chi-square brackets.
cw_sig_pairs <- list(
  c("Susceptible (any)", "Resistant (any)"), c("Non-MDR", "MDR"), c("Non-XDR", "XDR"),
  c("Carbapenem-S", "Carbapenem-R"), c("Cephalosporin-S", "Cephalosporin-R"), c("3GC-S", "3GCR"),
  c("MSSA", "MRSA"), c("VSE", "VRE"), c("Access-S", "Access-R"), c("Watch-S", "Watch-R"), c("Reserve-S", "Reserve-R")
)
cw_pair_frame <- list(
  "Susceptible (any)" = mortality_analytic_cw, "Resistant (any)" = mortality_analytic_cw,
  "Non-MDR" = mortality_cw_mdr, "MDR" = mortality_cw_mdr, "Non-XDR" = mortality_cw_xdr, "XDR" = mortality_cw_xdr,
  "Carbapenem-S" = mortality_cw_carb, "Carbapenem-R" = mortality_cw_carb,
  "Cephalosporin-S" = mortality_cw_ceph, "Cephalosporin-R" = mortality_cw_ceph,
  "3GC-S" = mortality_cw_3gcr, "3GCR" = mortality_cw_3gcr, "MSSA" = mssa_cw, "MRSA" = mortality_cw_mrsa,
  "VSE" = vre_cw, "VRE" = vre_cw, "Access-S" = mortality_cw_access, "Access-R" = mortality_cw_access,
  "Watch-S" = mortality_cw_watch, "Watch-R" = mortality_cw_watch, "Reserve-S" = mortality_cw_reserve, "Reserve-R" = mortality_cw_reserve
)
cw_grp_label <- list("Susceptible (any)" = "Susceptible BSI", "Resistant (any)" = "Resistant BSI",
  "Non-MDR" = "Non-MDR BSI", "MDR" = "MDR BSI", "Non-XDR" = "Non-XDR BSI", "XDR" = "XDR BSI",
  "Carbapenem-S" = "Carbapenem-S Gram-negative BSI", "Carbapenem-R" = "Carbapenem-R Gram-negative BSI",
  "Cephalosporin-S" = "Cephalosporin-S Gram-negative BSI", "Cephalosporin-R" = "Cephalosporin-R Gram-negative BSI",
  "3GC-S" = "3GC-susceptible Gram-negative BSI", "3GCR" = "3GCR Gram-negative BSI",
  "MSSA" = "MSSA", "MRSA" = "MRSA BSI", "VSE" = "VSE", "VRE" = "VRE",
  "Access-S" = "Access-S BSI", "Access-R" = "Access-R BSI", "Watch-S" = "Watch-S BSI", "Watch-R" = "Watch-R BSI",
  "Reserve-S" = "Reserve-S BSI", "Reserve-R" = "Reserve-R BSI")

cw_sig_data <- purrr::map_dfr(cw_sig_pairs, function(p) {
  d1 <- cw_pair_frame[[p[1]]] %>% dplyr::filter(mortality_group == cw_grp_label[[p[1]]], !is.na(event28), !is.na(calibrated_weight))
  d2 <- cw_pair_frame[[p[2]]] %>% dplyr::filter(mortality_group == cw_grp_label[[p[2]]], !is.na(event28), !is.na(calibrated_weight))
  both <- bind_rows(d1 %>% mutate(.pairgrp = p[1]), d2 %>% mutate(.pairgrp = p[2]))
  both$pairgrp_f <- factor(both$.pairgrp, levels = p)
  both$died_f <- factor(both$event28 == 1, levels = c(FALSE, TRUE))
  res <- tryCatch({
    des <- svydesign(ids = ~facility_id, weights = ~calibrated_weight, data = both)
    t <- svychisq(~pairgrp_f + died_f, des, statistic = "Chisq")
    tibble(stat = unname(t$statistic), p = unname(t$p.value))
  }, error = function(e) { cat("SVYCHISQ ERROR for", p[1], "vs", p[2], ":", conditionMessage(e), "\n"); tibble(stat = NA_real_, p = NA_real_) })
  tibble(group1 = p[1], group2 = p[2], chisq_stat = round(res$stat, 2), chisq_p = round(res$p, 4))
})
print(as.data.frame(cw_sig_data), row.names = FALSE)
write_xlsx(cw_sig_data, "Mortality_crude_calibrated_weighted_chisq_pairs.xlsx")

# Calibrated-weighted KM (primary 3-group + MDR mechanism), same weighted-
# survfit() pattern already used elsewhere in this script (SECTIONS 9.9/9.12)
# for IPW-weighted curves, substituting calibrated_weight for ipw.
build_cw_km <- function(data) {
  dd <- data %>% dplyr::filter(!is.na(calibrated_weight))
  fit <- survfit(Surv(time28, event28) ~ mortality_group, weights = calibrated_weight, data = dd)
  tidy_ <- broom::tidy(fit) %>% mutate(mortality_group = sub("^mortality_group=", "", strata))
  tidy0 <- tidy_ %>% distinct(mortality_group) %>% mutate(time = 0, estimate = 1, conf.low = 1, conf.high = 1)
  bind_rows(tidy0, tidy_) %>% arrange(mortality_group, time)
}

km_cw_primary_data <- build_cw_km(mortality_analytic_cw)
km_cw_mdr_data <- build_cw_km(mortality_cw_mdr)

KM_INK_SOFT_CW <- "#4D616A"; KM_LINE_CW <- "#D3DCDD"
km_colors_primary <- c("Uninfected control" = "#17766E", "Susceptible BSI" = "#A8541D", "Resistant BSI" = "#7A1F2B")
km_colors_mdr <- c("Uninfected control" = "#17766E", "Non-MDR BSI" = "#A8541D", "MDR BSI" = "#7A1F2B")

build_cw_km_plot <- function(km_data, colors) {
  ggplot(km_data, aes(x = time, y = 1 - estimate, color = mortality_group, fill = mortality_group)) +
    geom_ribbon(aes(ymin = 1 - conf.low, ymax = 1 - conf.high), alpha = 0.12, color = NA) +
    geom_step(linewidth = 1.1) +
    scale_color_manual(values = colors, name = NULL) +
    scale_fill_manual(values = colors, guide = "none") +
    scale_x_continuous(breaks = seq(0, 28, 7), limits = c(0, 28), expand = c(0, 0)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
    labs(x = "Days since enrolment", y = "Calibrated cumulative mortality") +
    theme_minimal(base_size = 24) +
    theme(
      panel.grid.minor = element_blank(), panel.grid.major = element_line(color = KM_LINE_CW, linewidth = 0.5),
      legend.position = "top", legend.text = element_text(size = 21),
      axis.title = element_text(color = KM_INK_SOFT_CW, size = 23), axis.text = element_text(size = 20)
    )
}

ggsave("Mortality_KM_28day_calibrated_weighted_primary.png", build_cw_km_plot(km_cw_primary_data, km_colors_primary), width = 11, height = 7, dpi = 500)
ggsave("Mortality_KM_28day_calibrated_weighted_primary.tiff", build_cw_km_plot(km_cw_primary_data, km_colors_primary), width = 11, height = 7, dpi = 500, compression = "lzw")
ggsave("Mortality_KM_28day_calibrated_weighted_MDR.png", build_cw_km_plot(km_cw_mdr_data, km_colors_mdr), width = 11, height = 7, dpi = 500)
ggsave("Mortality_KM_28day_calibrated_weighted_MDR.tiff", build_cw_km_plot(km_cw_mdr_data, km_colors_mdr), width = 11, height = 7, dpi = 500, compression = "lzw")

cat("CALIBRATED-WEIGHTED CRUDE MORTALITY + KM DONE\n")

# =========================================================================
# SUPPLEMENTARY: IPW-only LOS model (no calibrated design weight), every
# mechanism -- the main-text LOS model (fit_los_glm, SECTION 11) is IPW x
# calibrated-weight COMBINED throughout; this strips out the calibrated
# weight to give the IPW-only companion requested for the Supplementary.
# Reuses the mortality_cw_* mechanism frames already built above (fast,
# no MICE) -- these all predate SECTION 9.3's Level_of_care/charlson_score_
# combined attachment (only ever joined onto mortality_analytic itself), so
# charlson_score_combined (los_covars' one comorbidity term) is recomputed
# fresh here via add_charlson(), same formula as SECTION 9.3.
cat("=== IPW-ONLY LOS MODEL, ALL MECHANISMS (fast) ===\n")
add_charlson <- function(df) {
  d_pts <- case_when(df$comorb_diabetes_end_organ == 1L ~ 2L, df$comorb_diabetes == 1L ~ 1L, TRUE ~ 0L)
  l_pts <- case_when(df$comorb_liver_severe == 1L ~ 3L, df$comorb_liver_mild == 1L ~ 1L, TRUE ~ 0L)
  m_pts <- case_when(df$comorb_metastatic_tumour == 1L ~ 6L, df$comorb_malignancy_active == 1L ~ 2L, TRUE ~ 0L)
  h_pts <- if_else(df$comorb_hiv_on_art == 1L | df$comorb_hiv_no_art == 1L | df$comorb_hiv_cd4_low == 1L, 6L, 0L, missing = 0L)
  w_pts <- df %>% select(comorb_chf, comorb_cerebrovascular, comorb_dementia, comorb_copd, comorb_connective_tissue, comorb_peptic_ulcer) %>%
    mutate(across(everything(), ~ if_else(.x == 1L, 1L, 0L, missing = 0L))) %>% rowSums()
  he_pts <- if_else(df$comorb_hemi_paraplegia == 1L, 2L, 0L, missing = 0L)
  r_pts <- if_else(df$comorb_renal_disease == 1L, 2L, 0L, missing = 0L)
  df$charlson_score_combined <- w_pts + d_pts + l_pts + m_pts + h_pts + he_pts + r_pts
  df
}

los_covars_ipw <- c("age_years", "sex_at_birth", "charlson_score_combined")

los_analytic_mrsa_cw <- mortality_analytic_cw %>%
  mutate(mortality_group = case_when(
    matched_control_broad == 1 ~ "Uninfected control",
    matched_control_broad == 0 & pathogen1_clean == "Staphylococcus aureus" & methicillin_BSI_positive == 1 ~ "MRSA BSI",
    matched_control_broad == 0 & pathogen1_clean == "Staphylococcus aureus" & methicillin_BSI_positive == 0 ~ "MSSA BSI",
    TRUE ~ NA_character_
  ), mortality_group = factor(mortality_group, levels = c("Uninfected control", "MSSA BSI", "MRSA BSI"))) %>%
  filter(!is.na(mortality_group)) %>% add_charlson()
los_mrsa_comparisons_ipw <- list(
  c(ref = "Uninfected control", comp = "MSSA BSI"),
  c(ref = "Uninfected control", comp = "MRSA BSI"),
  c(ref = "MSSA BSI", comp = "MRSA BSI")
)

fit_los_glm_ipw <- function(data, group_var, ref_level, comp_level, covars) {
  di <- data %>% dplyr::filter(.data[[group_var]] %in% c(ref_level, comp_level)) %>%
    mutate(group_bin = as.integer(.data[[group_var]] == comp_level)) %>%
    filter(if_all(all_of(covars), ~ !is.na(.)))
  form_ps <- as.formula(paste("group_bin ~", paste(covars, collapse = " + ")))
  ps_mod <- glm(form_ps, data = di, family = binomial())
  ps <- pmin(pmax(predict(ps_mod, type = "response"), 0.01), 0.99)
  p_marg <- mean(di$group_bin)
  ipw <- if_else(di$group_bin == 1, p_marg / ps, (1 - p_marg) / (1 - ps))
  q <- quantile(ipw, c(0.01, 0.99), na.rm = TRUE)
  di$combined_w <- pmin(pmax(ipw, q[1]), q[2])
  form <- as.formula(paste("time_inhosp ~ group_bin +", paste(covars, collapse = " + ")))
  glm(form, data = di, family = Gamma(link = "log"), weights = combined_w)
}

standardised_mean_diff_ipw <- function(mod, di) {
  d0 <- di; d0$group_bin <- 0L
  d1 <- di; d1$group_bin <- 1L
  mean_ref <- mean(predict(mod, newdata = d0, type = "response"), na.rm = TRUE)
  mean_comp <- mean(predict(mod, newdata = d1, type = "response"), na.rm = TRUE)
  c(mean_ref = mean_ref, mean_comp = mean_comp, diff = mean_comp - mean_ref)
}

bootstrap_los_diff_ipw <- function(data, group_var, ref_level, comp_level, covars, n_boot = 200) {
  facs <- unique(data$facility_id)
  out <- matrix(NA_real_, nrow = n_boot, ncol = 3, dimnames = list(NULL, c("mean_ref", "mean_comp", "diff")))
  for (b in seq_len(n_boot)) {
    boot_facs <- sample(facs, length(facs), replace = TRUE)
    boot_rows <- purrr::map_dfr(boot_facs, function(f) data[data$facility_id == f, , drop = FALSE])
    r <- tryCatch({
      di <- boot_rows %>% dplyr::filter(.data[[group_var]] %in% c(ref_level, comp_level)) %>%
        mutate(group_bin = as.integer(.data[[group_var]] == comp_level)) %>%
        filter(if_all(all_of(covars), ~ !is.na(.)))
      form_ps <- as.formula(paste("group_bin ~", paste(covars, collapse = " + ")))
      ps_mod <- glm(form_ps, data = di, family = binomial())
      ps <- pmin(pmax(predict(ps_mod, type = "response"), 0.01), 0.99)
      p_marg <- mean(di$group_bin)
      ipw <- if_else(di$group_bin == 1, p_marg / ps, (1 - p_marg) / (1 - ps))
      q <- quantile(ipw, c(0.01, 0.99), na.rm = TRUE)
      di$combined_w <- pmin(pmax(ipw, q[1]), q[2])
      mod <- glm(as.formula(paste("time_inhosp ~ group_bin +", paste(covars, collapse = " + "))),
                 data = di, family = Gamma(link = "log"), weights = combined_w)
      standardised_mean_diff_ipw(mod, di)
    }, error = function(e) c(mean_ref = NA_real_, mean_comp = NA_real_, diff = NA_real_))
    out[b, ] <- r
  }
  out
}

run_los_subgroup_ipw <- function(mechanism, analytic_data, comparisons, covars = los_covars_ipw, n_boot = 200, seed_boot = 45) {
  di_all <- analytic_data %>% filter(!is.na(time_inhosp), !is.na(facility_id))
  purrr::map_dfr(comparisons, function(cmp) {
    mod <- tryCatch(fit_los_glm_ipw(di_all, "mortality_group", cmp["ref"], cmp["comp"], covars), error = function(e) NULL)
    if (is.null(mod)) {
      return(tibble(mechanism = mechanism, comparison = paste(cmp["comp"], "vs", cmp["ref"]), n = NA_integer_,
                     adjusted_mean_los_ref_days = NA_real_, adjusted_mean_los_comp_days = NA_real_,
                     mean_difference_days = NA_real_, ci_lower_days = NA_real_, ci_upper_days = NA_real_, n_bootstrap_ok = 0L))
    }
    di <- di_all %>% dplyr::filter(mortality_group %in% c(cmp["ref"], cmp["comp"])) %>%
      mutate(group_bin = as.integer(mortality_group == cmp["comp"]))
    point <- standardised_mean_diff_ipw(mod, di)
    set.seed(seed_boot)
    boot_mat <- bootstrap_los_diff_ipw(di_all, "mortality_group", cmp["ref"], cmp["comp"], covars, n_boot = n_boot)
    boot_ok <- boot_mat[complete.cases(boot_mat), , drop = FALSE]
    tibble(
      mechanism = mechanism, comparison = paste(cmp["comp"], "vs", cmp["ref"]), n = nrow(di),
      adjusted_mean_los_ref_days = point["mean_ref"], adjusted_mean_los_comp_days = point["mean_comp"],
      mean_difference_days = point["diff"],
      ci_lower_days = if (nrow(boot_ok) > 0) quantile(boot_ok[, "diff"], 0.025, na.rm = TRUE) else NA_real_,
      ci_upper_days = if (nrow(boot_ok) > 0) quantile(boot_ok[, "diff"], 0.975, na.rm = TRUE) else NA_real_,
      n_bootstrap_ok = nrow(boot_ok)
    )
  }) %>% mutate(across(where(is.numeric) & !c(n, n_bootstrap_ok), ~ round(.x, 3)))
}

los_mechanism_list_ipw <- list(
  list(label = "Any resistance", data = mortality_analytic_cw %>% add_charlson(), comparisons = list(
    c(ref = "Uninfected control", comp = "Susceptible BSI"), c(ref = "Uninfected control", comp = "Resistant BSI"),
    c(ref = "Susceptible BSI", comp = "Resistant BSI"))),
  list(label = "MDR (bacterial isolates)", data = mortality_cw_mdr %>% add_charlson(), comparisons = list(
    c(ref = "Uninfected control", comp = "Non-MDR BSI"), c(ref = "Uninfected control", comp = "MDR BSI"),
    c(ref = "Non-MDR BSI", comp = "MDR BSI"))),
  list(label = "XDR (bacterial isolates)", data = mortality_cw_xdr %>% add_charlson(), comparisons = list(
    c(ref = "Uninfected control", comp = "Non-XDR BSI"), c(ref = "Uninfected control", comp = "XDR BSI"),
    c(ref = "Non-XDR BSI", comp = "XDR BSI"))),
  list(label = "Carbapenem-R (Gram-negative)", data = mortality_cw_carb %>% add_charlson(), comparisons = list(
    c(ref = "Uninfected control", comp = "Carbapenem-S Gram-negative BSI"), c(ref = "Uninfected control", comp = "Carbapenem-R Gram-negative BSI"),
    c(ref = "Carbapenem-S Gram-negative BSI", comp = "Carbapenem-R Gram-negative BSI"))),
  list(label = "Cephalosporin-R any-gen (Gram-negative)", data = mortality_cw_ceph %>% add_charlson(), comparisons = list(
    c(ref = "Uninfected control", comp = "Cephalosporin-S Gram-negative BSI"), c(ref = "Uninfected control", comp = "Cephalosporin-R Gram-negative BSI"),
    c(ref = "Cephalosporin-S Gram-negative BSI", comp = "Cephalosporin-R Gram-negative BSI"))),
  list(label = "3GCR (Gram-negative)", data = mortality_cw_3gcr %>% add_charlson(), comparisons = list(
    c(ref = "Uninfected control", comp = "3GC-susceptible Gram-negative BSI"), c(ref = "Uninfected control", comp = "3GCR Gram-negative BSI"),
    c(ref = "3GC-susceptible Gram-negative BSI", comp = "3GCR Gram-negative BSI"))),
  list(label = "Methicillin-R/MRSA (S. aureus)", data = los_analytic_mrsa_cw, comparisons = los_mrsa_comparisons_ipw),
  list(label = "Access-resistant (AWaRe)", data = mortality_cw_access %>% add_charlson(), comparisons = list(
    c(ref = "Uninfected control", comp = "Access-S BSI"), c(ref = "Uninfected control", comp = "Access-R BSI"),
    c(ref = "Access-S BSI", comp = "Access-R BSI"))),
  list(label = "Watch-resistant (AWaRe)", data = mortality_cw_watch %>% add_charlson(), comparisons = list(
    c(ref = "Uninfected control", comp = "Watch-S BSI"), c(ref = "Uninfected control", comp = "Watch-R BSI"),
    c(ref = "Watch-S BSI", comp = "Watch-R BSI"))),
  list(label = "Reserve-resistant (AWaRe)", data = mortality_cw_reserve %>% add_charlson(), comparisons = list(
    c(ref = "Uninfected control", comp = "Reserve-S BSI"), c(ref = "Uninfected control", comp = "Reserve-R BSI"),
    c(ref = "Reserve-S BSI", comp = "Reserve-R BSI")))
)
los_all_results_ipw <- purrr::map_dfr(los_mechanism_list_ipw, ~ run_los_subgroup_ipw(.x$label, .x$data, .x$comparisons))
print(as.data.frame(los_all_results_ipw), row.names = FALSE)
write_xlsx(los_all_results_ipw, "LOS_excess_length_of_stay_IPW_only_all_mechanisms.xlsx")

los_forest_data_ipw <- los_all_results_ipw %>%
  group_by(mechanism) %>%
  mutate(comparison_short = factor(
    c("Susceptible-type vs Uninfected", "Resistant-type vs Uninfected", "Resistant-type vs Susceptible-type")[row_number()],
    levels = rev(c("Susceptible-type vs Uninfected", "Resistant-type vs Uninfected", "Resistant-type vs Susceptible-type"))
  )) %>%
  ungroup() %>%
  mutate(mechanism = factor(mechanism, levels = rev(unique(mechanism))))

CLIP_AT_IPW <- 15
los_forest_data_ipw <- los_forest_data_ipw %>%
  mutate(ci_truncated = ci_upper_days > CLIP_AT_IPW, ci_upper_plot = pmin(ci_upper_days, CLIP_AT_IPW))

los_forest_ipw_plot <- ggplot(los_forest_data_ipw, aes(x = mean_difference_days, y = comparison_short)) +
  geom_vline(xintercept = 0, color = "#4D616A", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = ci_lower_days, xmax = ci_upper_plot), height = 0.25, color = "#4D616A", linewidth = 0.6, na.rm = TRUE) +
  geom_point(size = 3.2, color = "#7A1F2B", na.rm = TRUE) +
  geom_text(data = . %>% dplyr::filter(ci_truncated), aes(x = ci_upper_plot, label = "*"),
            hjust = -0.3, vjust = 0.78, size = 8, color = "#7A1F2B", fontface = "bold") +
  facet_grid(mechanism ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(title = "Adjusted excess length of stay by resistance mechanism (IPW only, no calibrated weight)",
       subtitle = paste0("Every mechanism benchmarked vs Susceptible AND vs Uninfected.\n",
                          "Points = adjusted mean difference (days); bars = 95% bootstrap CI (facility-cluster).\n",
                          "* = CI truncated at ", CLIP_AT_IPW, " days for display."),
       x = "Adjusted mean difference in LOS (days)", y = NULL) +
  theme_minimal(base_size = 17) +
  theme(panel.grid.minor = element_blank(), strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 0, size = 15),
        strip.placement = "outside", plot.subtitle = element_text(color = "#4D616A", size = 13),
        plot.title = element_text(size = 17, face = "bold"),
        axis.text = element_text(size = 14), axis.title = element_text(size = 15))
ggsave("LOS_forest_all_mechanisms_IPWonly.png", los_forest_ipw_plot, width = 17, height = 13, dpi = 500)
ggsave("LOS_forest_all_mechanisms_IPWonly.tiff", los_forest_ipw_plot, width = 17, height = 13, dpi = 500, compression = "lzw")
cat("IPW-ONLY LOS SAVED\n")

# =========================================================================
# SUPPLEMENTARY: sex- and age-group-stratified crude outcomes (prevalence,
# mortality, LOS) across all resistance groups -- 3 new Lancet-style
# figures. Unadjusted/crude throughout (patient-level counts and means, no
# IPW/calibration), per instruction. Sex: 0=Male, 1=Female (SECTION 1's own
# coding). Age group collapsed to two categories only (adult/paediatric,
# adult = age_years>=18), per instruction -- neonates are not broken out
# separately here (they are elsewhere, e.g. SECTION 11's LOS age panels).
cat("=== SEX/AGE-STRATIFIED CRUDE OUTCOMES (fast) ===\n")

crude_row_by_sex <- function(data, group_label, display_label) {
  sub <- data %>% dplyr::filter(mortality_group == group_label, !is.na(sex_at_birth))
  purrr::map_dfr(c(0L, 1L), function(sx) {
    s <- sub %>% dplyr::filter(sex_at_birth == sx)
    n <- nrow(s)
    n_known <- sum(!is.na(s$event28))
    deaths <- sum(s$event28 == 1, na.rm = TRUE)
    mort_pct <- if (n_known > 0) round(100 * deaths / n_known, 1) else NA_real_
    los_vals <- s$time_inhosp[!is.na(s$time_inhosp)]
    los_mean <- if (length(los_vals) > 0) round(mean(los_vals), 2) else NA_real_
    los_sd <- if (length(los_vals) > 1) round(sd(los_vals), 2) else NA_real_
    los_se <- if (length(los_vals) > 1) round(sd(los_vals) / sqrt(length(los_vals)), 2) else NA_real_
    n_adult <- sum(s$age_years >= 18, na.rm = TRUE)
    n_age_known <- sum(!is.na(s$age_years))
    adult_pct <- if (n_age_known > 0) round(100 * n_adult / n_age_known, 1) else NA_real_
    tibble(group = display_label, sex = if_else(sx == 1L, "Female", "Male"), n = n,
           n_outcome_known = n_known, deaths_28day = deaths, mortality_pct = mort_pct,
           n_los_known = length(los_vals), los_mean_days = los_mean, los_sd_days = los_sd, los_se_days = los_se,
           n_adult = n_adult, n_age_known = n_age_known, adult_pct = adult_pct)
  })
}

crude_group_order <- c(
  "Uninfected", "Susceptible (any)", "Resistant (any)",
  "Non-MDR", "MDR", "Non-XDR", "XDR",
  "Carbapenem-S", "Carbapenem-R", "Cephalosporin-S", "Cephalosporin-R",
  "3GC-S", "3GCR", "MSSA", "MRSA", "VSE", "VRE",
  "Access-S", "Access-R", "Watch-S", "Watch-R", "Reserve-S", "Reserve-R"
)
crude_sex_all <- bind_rows(
  crude_row_by_sex(mortality_analytic_cw, "Uninfected control", "Uninfected"),
  crude_row_by_sex(mortality_analytic_cw, "Susceptible BSI", "Susceptible (any)"),
  crude_row_by_sex(mortality_analytic_cw, "Resistant BSI", "Resistant (any)"),
  crude_row_by_sex(mortality_cw_mdr, "Non-MDR BSI", "Non-MDR"),
  crude_row_by_sex(mortality_cw_mdr, "MDR BSI", "MDR"),
  crude_row_by_sex(mortality_cw_xdr, "Non-XDR BSI", "Non-XDR"),
  crude_row_by_sex(mortality_cw_xdr, "XDR BSI", "XDR"),
  crude_row_by_sex(mortality_cw_carb, "Carbapenem-S Gram-negative BSI", "Carbapenem-S"),
  crude_row_by_sex(mortality_cw_carb, "Carbapenem-R Gram-negative BSI", "Carbapenem-R"),
  crude_row_by_sex(mortality_cw_ceph, "Cephalosporin-S Gram-negative BSI", "Cephalosporin-S"),
  crude_row_by_sex(mortality_cw_ceph, "Cephalosporin-R Gram-negative BSI", "Cephalosporin-R"),
  crude_row_by_sex(mortality_cw_3gcr, "3GC-susceptible Gram-negative BSI", "3GC-S"),
  crude_row_by_sex(mortality_cw_3gcr, "3GCR Gram-negative BSI", "3GCR"),
  crude_row_by_sex(mssa_cw %>% mutate(mortality_group = "MSSA"), "MSSA", "MSSA"),
  crude_row_by_sex(mortality_cw_mrsa, "MRSA BSI", "MRSA"),
  crude_row_by_sex(vre_cw, "VSE", "VSE"),
  crude_row_by_sex(vre_cw, "VRE", "VRE"),
  crude_row_by_sex(mortality_cw_access, "Access-S BSI", "Access-S"),
  crude_row_by_sex(mortality_cw_access, "Access-R BSI", "Access-R"),
  crude_row_by_sex(mortality_cw_watch, "Watch-S BSI", "Watch-S"),
  crude_row_by_sex(mortality_cw_watch, "Watch-R BSI", "Watch-R"),
  crude_row_by_sex(mortality_cw_reserve, "Reserve-S BSI", "Reserve-S"),
  crude_row_by_sex(mortality_cw_reserve, "Reserve-R BSI", "Reserve-R")
) %>%
  mutate(group = factor(group, levels = rev(crude_group_order)))
print(as.data.frame(crude_sex_all), row.names = FALSE)
write_xlsx(crude_sex_all, "SexAge_stratified_mortality_LOS_all_mechanisms.xlsx")

prevalence_row_by_sex <- function(data, susceptible_label, resistant_label, display_label) {
  sub <- data %>% dplyr::filter(mortality_group %in% c(susceptible_label, resistant_label), !is.na(sex_at_birth))
  purrr::map_dfr(c(0L, 1L), function(sx) {
    s <- sub %>% dplyr::filter(sex_at_birth == sx)
    n_tested <- nrow(s)
    n_resistant <- sum(s$mortality_group == resistant_label)
    prev_pct <- if (n_tested > 0) round(100 * n_resistant / n_tested, 1) else NA_real_
    n_adult <- sum(s$age_years >= 18, na.rm = TRUE)
    n_age_known <- sum(!is.na(s$age_years))
    adult_pct <- if (n_age_known > 0) round(100 * n_adult / n_age_known, 1) else NA_real_
    tibble(group = display_label, sex = if_else(sx == 1L, "Female", "Male"), n_tested = n_tested,
           n_resistant = n_resistant, prevalence_pct = prev_pct,
           n_adult = n_adult, n_age_known = n_age_known, adult_pct = adult_pct)
  })
}

mrsa_prevalence_frame <- bind_rows(
  mssa_cw %>% mutate(mortality_group = "MSSA"),
  mortality_cw_mrsa %>% dplyr::filter(mortality_group == "MRSA BSI")
)

prevalence_group_order <- c("Any resistance", "MDR", "XDR", "Carbapenem-R", "Cephalosporin-R", "3GCR",
                             "MRSA", "Access-R", "Watch-R", "Reserve-R")
prevalence_sex_all <- bind_rows(
  prevalence_row_by_sex(mortality_analytic_cw, "Susceptible BSI", "Resistant BSI", "Any resistance"),
  prevalence_row_by_sex(mortality_cw_mdr, "Non-MDR BSI", "MDR BSI", "MDR"),
  prevalence_row_by_sex(mortality_cw_xdr, "Non-XDR BSI", "XDR BSI", "XDR"),
  prevalence_row_by_sex(mortality_cw_carb, "Carbapenem-S Gram-negative BSI", "Carbapenem-R Gram-negative BSI", "Carbapenem-R"),
  prevalence_row_by_sex(mortality_cw_ceph, "Cephalosporin-S Gram-negative BSI", "Cephalosporin-R Gram-negative BSI", "Cephalosporin-R"),
  prevalence_row_by_sex(mortality_cw_3gcr, "3GC-susceptible Gram-negative BSI", "3GCR Gram-negative BSI", "3GCR"),
  prevalence_row_by_sex(mrsa_prevalence_frame, "MSSA", "MRSA BSI", "MRSA"),
  prevalence_row_by_sex(mortality_cw_access, "Access-S BSI", "Access-R BSI", "Access-R"),
  prevalence_row_by_sex(mortality_cw_watch, "Watch-S BSI", "Watch-R BSI", "Watch-R"),
  prevalence_row_by_sex(mortality_cw_reserve, "Reserve-S BSI", "Reserve-R BSI", "Reserve-R")
) %>%
  mutate(group = factor(group, levels = rev(prevalence_group_order)))
print(as.data.frame(prevalence_sex_all), row.names = FALSE)
write_xlsx(prevalence_sex_all, "SexAge_stratified_prevalence_all_mechanisms.xlsx")

# --- 3 Lancet-style figures: resistance group on y, dodged bars by sex,
# red adult-% label on top of each bar. ---
sex_colors <- c("Male" = "#2255A4", "Female" = "#87CEEB")
ADULT_RED <- "#C0392B"

build_dodge_bar <- function(data, value_col, value_label, is_pct = TRUE, err_lower = NULL, err_upper = NULL) {
  data$value_ <- data[[value_col]]
  p <- ggplot(data, aes(x = group, y = value_, fill = sex)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7, colour = "black", linewidth = 0.25)
  if (!is.null(err_lower) && !is.null(err_upper)) {
    data$err_lower_ <- data[[err_lower]]; data$err_upper_ <- data[[err_upper]]
    p <- p + geom_errorbar(data = data, aes(ymin = err_lower_, ymax = err_upper_),
                            position = position_dodge(width = 0.75), width = 0.3, linewidth = 0.35, colour = "grey20", na.rm = TRUE)
  }
  p +
    geom_text(aes(label = ifelse(is.na(value_), "", if (is_pct) sprintf("%.1f%%", value_) else sprintf("%.1f", value_))),
              position = position_dodge(width = 0.75), hjust = -0.15, size = 2.6, colour = "black") +
    geom_text(aes(y = 0, label = ifelse(is.na(adult_pct), "", sprintf("%.0f%%A", adult_pct))),
              position = position_dodge(width = 0.75), hjust = 1.05, size = 2.35, colour = ADULT_RED, fontface = "bold") +
    scale_fill_manual(values = sex_colors, name = NULL) +
    coord_flip(clip = "off") +
    labs(x = NULL, y = value_label) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "top", legend.text = element_text(size = 9.5),
      panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.3),
      axis.text.y = element_text(size = 8.8, colour = "black"),
      axis.text.x = element_text(size = 8.5, colour = "black"),
      axis.title.x = element_text(size = 9.5, colour = "black"),
      plot.margin = margin(6, 30, 4, 40)
    )
}

fig_prevalence <- build_dodge_bar(prevalence_sex_all, "prevalence_pct", "Crude prevalence (%)") +
  scale_y_continuous(limits = c(0, max(prevalence_sex_all$prevalence_pct, na.rm = TRUE) * 1.35), expand = expansion(mult = c(0, 0))) +
  labs(caption = "Red label (left of each bar) = % of that bar's patients who are adults (≥ 18y); \"A\" marks the adult-proportion label. Bars: crude (unadjusted) prevalence, dodged by sex (blue = male, sky blue = female).")
ggsave("SexAge_prevalence_by_mechanism.png", fig_prevalence, width = 9.5, height = 7.5, dpi = 500)
ggsave("SexAge_prevalence_by_mechanism.tiff", fig_prevalence, width = 9.5, height = 7.5, dpi = 500, compression = "lzw")

fig_mortality <- build_dodge_bar(crude_sex_all, "mortality_pct", "Crude 28-day mortality (%)") +
  scale_y_continuous(limits = c(0, max(crude_sex_all$mortality_pct, na.rm = TRUE) * 1.35), expand = expansion(mult = c(0, 0))) +
  labs(caption = "Red label (left of each bar) = % of that bar's patients who are adults (≥ 18y); \"A\" marks the adult-proportion label. Bars: crude (unadjusted) 28-day mortality, dodged by sex (blue = male, sky blue = female).")
ggsave("SexAge_mortality_by_mechanism.png", fig_mortality, width = 9.5, height = 10, dpi = 500)
ggsave("SexAge_mortality_by_mechanism.tiff", fig_mortality, width = 9.5, height = 10, dpi = 500, compression = "lzw")

los_plot_data <- crude_sex_all %>%
  mutate(los_ci_lower = los_mean_days - 1.96 * los_se_days, los_ci_upper = los_mean_days + 1.96 * los_se_days)
fig_los <- build_dodge_bar(los_plot_data, "los_mean_days", "Crude mean LOS (days)", is_pct = FALSE,
                            err_lower = "los_ci_lower", err_upper = "los_ci_upper") +
  scale_y_continuous(limits = c(0, max(los_plot_data$los_ci_upper, na.rm = TRUE) * 1.3), expand = expansion(mult = c(0, 0))) +
  labs(caption = "Red label (left of each bar) = % of that bar's patients who are adults (≥ 18y); \"A\" marks the adult-proportion label. Bars: crude (unadjusted) mean LOS (days), 95% CI, dodged by sex (blue = male, sky blue = female).")
ggsave("SexAge_LOS_by_mechanism.png", fig_los, width = 9.5, height = 10, dpi = 500)
ggsave("SexAge_LOS_by_mechanism.tiff", fig_los, width = 9.5, height = 10, dpi = 500, compression = "lzw")

cat("SEX/AGE-STRATIFIED FIGURES SAVED\n")

# =========================================================================
# SUPPLEMENTARY: expanded descriptive ("Table 1"-style) statistics, every
# resistance-mechanism dichotomy -- extends SECTION 9.5's primary
# Resistant/Susceptible/Uninfected table1 to the same 11 mechanism-specific
# comparisons used throughout the mortality/LOS analyses (MDR, XDR,
# carbapenem, cephalosporin, 3GCR, MRSA/MSSA, VRE/VSE, Access/Watch/
# Reserve), plus an overall Confirmed-BSI-vs-Uninfected dichotomy. Reuses
# the fast mortality_cw_* frames (SECTION "CALIBRATED-WEIGHTED CRUDE
# MORTALITY", above) -- these predate charlson_score_combined's attachment
# to mortality_analytic (SECTION 9.5), so it is recomputed here via
# add_charlson() (already defined for the IPW-only LOS block, above).
cat("=== EXPANDED DESCRIPTIVE (TABLE 1) STATISTICS, ALL MECHANISMS (fast) ===\n")

fmt_mean_sd <- function(x) sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
fmt_pct <- function(n, d) if_else(d == 0, "-", sprintf("%d (%.1f%%)", n, 100 * n / d))
# Level_of_care is a facility-level attribute only joined onto mortality_analytic
# at SECTION 9.5/9.6 (below) -- these fast mortality_cw_* frames predate that
# join, so it must be attached fresh here, same pattern as add_charlson().
add_loc <- function(df) {
  df %>% select(-any_of("Level_of_care")) %>%
    left_join(Hospital_characteristics %>% select(`Hospital ID`, Level_of_care), by = c(facility_id = "Hospital ID"))
}
desc_cont_vars <- c("age_years", "charlson_score_combined", "LOS_days")
desc_cat_vars  <- c("sex_at_birth", "ward_type_icu", "specialty_group", "Level_of_care")

desc_group_table <- function(df, labA, labB, nameA, nameB) {
  df <- df %>% mutate(mortality_group = as.character(mortality_group)) %>%
    dplyr::filter(mortality_group %in% c(labA, labB))
  grp <- df$mortality_group
  rows <- list()
  for (v in desc_cont_vars) {
    x <- df[[v]]
    p <- tryCatch(t.test(x[grp == labA], x[grp == labB])$p.value, error = function(e) NA_real_)
    rows[[length(rows) + 1]] <- tibble(
      variable = v, level = "mean (SD)",
      A = fmt_mean_sd(x[grp == labA]), B = fmt_mean_sd(x[grp == labB]),
      n_A = sum(!is.na(x[grp == labA])), n_B = sum(!is.na(x[grp == labB])), p_value = round(p, 4)
    )
  }
  for (v in desc_cat_vars) {
    x <- factor(df[[v]])
    tab <- table(x, grp, useNA = "no")
    if (!all(c(labA, labB) %in% colnames(tab))) next
    sub <- tab[, c(labA, labB), drop = FALSE]
    sub <- sub[rowSums(sub) > 0, , drop = FALSE]
    p <- tryCatch({
      exp_min <- suppressWarnings(min(chisq.test(sub)$expected))
      if (is.na(exp_min) || exp_min < 5) fisher.test(sub, workspace = 2e7)$p.value else chisq.test(sub)$p.value
    }, error = function(e) NA_real_)
    lvls <- rownames(sub)
    for (l in lvls) {
      rows[[length(rows) + 1]] <- tibble(
        variable = v, level = l,
        A = fmt_pct(sub[l, labA], sum(sub[, labA])), B = fmt_pct(sub[l, labB], sum(sub[, labB])),
        n_A = sum(sub[, labA]), n_B = sum(sub[, labB]),
        p_value = if_else(l == lvls[1], round(p, 4), NA_real_)
      )
    }
  }
  # 28-day mortality outcome row (event28), same test-selection rule
  tab_m <- table(df$event28, grp, useNA = "no")
  p_m <- tryCatch({
    exp_min <- suppressWarnings(min(chisq.test(tab_m)$expected))
    if (is.na(exp_min) || exp_min < 5) fisher.test(tab_m)$p.value else chisq.test(tab_m)$p.value
  }, error = function(e) NA_real_)
  rows[[length(rows) + 1]] <- tibble(
    variable = "event28", level = "28-day mortality, n (%)",
    A = fmt_pct(sum(df$event28 == 1 & grp == labA, na.rm = TRUE), sum(grp == labA)),
    B = fmt_pct(sum(df$event28 == 1 & grp == labB, na.rm = TRUE), sum(grp == labB)),
    n_A = sum(grp == labA), n_B = sum(grp == labB), p_value = round(p_m, 4)
  )
  out <- bind_rows(rows)
  names(out)[names(out) == "A"] <- nameA
  names(out)[names(out) == "B"] <- nameB
  names(out)[names(out) == "n_A"] <- paste0("n_", gsub("[^A-Za-z0-9]", "_", nameA))
  names(out)[names(out) == "n_B"] <- paste0("n_", gsub("[^A-Za-z0-9]", "_", nameB))
  out
}

mortality_analytic_cw_bsi <- mortality_analytic_cw %>%
  add_charlson() %>% add_loc() %>%
  mutate(mortality_group = if_else(mortality_group == "Uninfected control", "Uninfected control",
                                    "Confirmed BSI (susceptible + resistant)"))

mrsa_mssa_combined <- bind_rows(
  mssa_cw %>% add_charlson() %>% add_loc(),
  mortality_cw_mrsa %>% add_charlson() %>% add_loc() %>% dplyr::filter(mortality_group == "MRSA BSI")
)

descriptive_all_mechanisms <- list(
  "Confirmed_BSI_vs_Uninfected" = desc_group_table(mortality_analytic_cw_bsi, "Uninfected control", "Confirmed BSI (susceptible + resistant)", "Uninfected", "Confirmed BSI"),
  "Resistant_vs_Susceptible"    = desc_group_table(mortality_analytic_cw %>% add_charlson() %>% add_loc(), "Susceptible BSI", "Resistant BSI", "Susceptible", "Resistant"),
  "NonMDR_vs_MDR"               = desc_group_table(mortality_cw_mdr %>% add_charlson() %>% add_loc(), "Non-MDR BSI", "MDR BSI", "Non-MDR", "MDR"),
  "NonXDR_vs_XDR"               = desc_group_table(mortality_cw_xdr %>% add_charlson() %>% add_loc(), "Non-XDR BSI", "XDR BSI", "Non-XDR", "XDR"),
  "CarbapenemS_vs_CarbapenemR"  = desc_group_table(mortality_cw_carb %>% add_charlson() %>% add_loc(), "Carbapenem-S Gram-negative BSI", "Carbapenem-R Gram-negative BSI", "Carbapenem-S", "Carbapenem-R"),
  "CephalosporinS_vs_CephalosporinR" = desc_group_table(mortality_cw_ceph %>% add_charlson() %>% add_loc(), "Cephalosporin-S Gram-negative BSI", "Cephalosporin-R Gram-negative BSI", "Cephalosporin-S", "Cephalosporin-R"),
  "3GCS_vs_3GCR"                = desc_group_table(mortality_cw_3gcr %>% add_charlson() %>% add_loc(), "3GC-susceptible Gram-negative BSI", "3GCR Gram-negative BSI", "3GC-S", "3GCR"),
  "MSSA_vs_MRSA"                = desc_group_table(mrsa_mssa_combined, "MSSA", "MRSA BSI", "MSSA", "MRSA"),
  "VSE_vs_VRE"                  = desc_group_table(vre_cw %>% add_charlson() %>% add_loc(), "VSE", "VRE", "VSE", "VRE"),
  "AccessS_vs_AccessR"          = desc_group_table(mortality_cw_access %>% add_charlson() %>% add_loc(), "Access-S BSI", "Access-R BSI", "Access-S", "Access-R"),
  "WatchS_vs_WatchR"            = desc_group_table(mortality_cw_watch %>% add_charlson() %>% add_loc(), "Watch-S BSI", "Watch-R BSI", "Watch-S", "Watch-R"),
  "ReserveS_vs_ReserveR"        = desc_group_table(mortality_cw_reserve %>% add_charlson() %>% add_loc(), "Reserve-S BSI", "Reserve-R BSI", "Reserve-S", "Reserve-R")
)
write_xlsx(descriptive_all_mechanisms, "Table1_descriptive_all_mechanisms.xlsx")
cat("EXPANDED DESCRIPTIVE TABLES SAVED\n")

# -----------------------------------------------------------------------
# 9.2b: CoNS-treatment sensitivity check (crude, unadjusted) -- the primary
# mortality analysis throughout this report always uses the unchanged
# mortality_group (above); this is a COMPACT sensitivity check only, not a
# parallel IPW/MSM re-fit under each alternate CoNS population (which would
# require re-deriving the full imputation+bootstrap machinery twice more).
# Reports crude 28-day mortality by group under each of the two alternate
# CoNS-treatment populations (SECTION 2.6), alongside the unchanged primary,
# so a reader can see whether the adjusted headline finding is sensitive to
# how CoNS is treated -- if the crude pattern doesn't move, the adjusted
# estimate is very unlikely to either.
# BSI_positive_CoNS_sensitive only ever DOWNGRADES a primary-confirmed case
# to 0 (never upgrades), so any record still ==1 here has an untouched,
# valid resistant_BSI_positive from the primary flag -- no parallel
# resistance flag is needed for this direction (unlike CoNS_included below).
mortality_group_CoNS_excluded <- case_when(
  Main_data_Combined$matched_control_broad == 1 ~ "Uninfected control",
  Main_data_Combined$BSI_positive_CoNS_sensitive == 1 & Main_data_Combined$resistant_BSI_positive == 1 ~ "Resistant BSI",
  Main_data_Combined$BSI_positive_CoNS_sensitive == 1 & Main_data_Combined$resistant_BSI_positive == 0 ~ "Susceptible BSI",
  TRUE ~ NA_character_
)
mortality_group_CoNS_included <- case_when(
  Main_data_Combined$matched_control_broad == 1 ~ "Uninfected control",
  Main_data_Combined$BSI_positive_CoNS_included == 1 & Main_data_Combined$resistant_BSI_positive_CoNS_included == 1 ~ "Resistant BSI",
  Main_data_Combined$BSI_positive_CoNS_included == 1 & Main_data_Combined$resistant_BSI_positive_CoNS_included == 0 ~ "Susceptible BSI",
  TRUE ~ NA_character_
)
crude_mortality_by_group <- function(group_vec, label) {
  tibble(group = group_vec, event28 = Main_data_Combined$event28) %>%
    filter(!is.na(group), !is.na(event28)) %>%
    group_by(group) %>%
    summarise(n = n(), n_deaths = sum(event28 == 1), crude_pct_mortality = round(100 * mean(event28 == 1), 1), .groups = "drop") %>%
    mutate(population = label, .before = 1)
}
mortality_cons_sensitivity <- bind_rows(
  crude_mortality_by_group(Main_data_Combined$mortality_group, "Primary (unchanged)"),
  crude_mortality_by_group(mortality_group_CoNS_excluded, "CoNS-excluded sensitivity"),
  crude_mortality_by_group(mortality_group_CoNS_included, "CoNS-included sensitivity")
)
write_xlsx(mortality_cons_sensitivity, "Mortality_CoNS_sensitivity_crude.xlsx")

# -----------------------------------------------------------------------
# 9.3: Candidate covariates + missingness report (supplementary table)
# -----------------------------------------------------------------------
# Per the methods spec's candidate list (age, sex, infection origin, ward
# type, specialty, hospital type, pre-culture LOS, severity, comorbidity
# burden, prior antibiotic exposure), mapped to the closest available
# variable in Main_data_Combined, with hospital type from
# Hospital_characteristics and an adapted Charlson score recomputed on
# Main_data_Combined's own (coalesced, SECTION 8) comorb_* columns so it is
# defined for controls too (SECTION 3.1's own charlson_score is
# CaseReport-only). Infection origin, pre-culture LOS and Pitt score are
# structurally BSI-specific (100% missing for the Uninfected control arm by
# definition, not just by chance) -- kept for descriptive comparison
# between the two BSI groups only, never as adjustment covariates.
#
# ward_type_icu and specialty_group are DELIBERATELY EXCLUDED from
# adjustment_covars (previously included, removed on review): both describe
# the patient's ward AT/AROUND enrolment, which for a BSI patient can be a
# CONSEQUENCE of the (possibly resistant) infection's severity -- e.g. ICU
# transfer triggered by septic shock -- not a pre-existing exposure. Using a
# potential post-exposure/collider variable as a propensity-model covariate
# risks biasing the causal estimate rather than correcting it. The
# unambiguously pre-culture versions of these concepts do exist on
# Main_data_CaseReport -- icu_48h_prior ("admitted to ICU >48h prior to
# onset of clinical signs of suspected BSI") and surgery_prior_to_bsi
# ("had surgery...prior to onset of clinical signs of suspected BSI") -- but
# BOTH are 100% missing for the Uninfected control arm (SECTION 7's Cohort
# Enrolment form does not collect them; they are not among the 57 columns
# SECTION 8 could coalesce across sources), so neither can serve as a shared
# adjustment covariate across all 3 mortality_group arms. They are reported
# descriptively instead (below, alongside the other structurally-controls-
# missing candidates), not silently dropped.
mortality_analytic <- mortality_analytic %>%
  left_join(Hospital_characteristics %>% select(`Hospital ID`, Level_of_care), by = c(facility_id = "Hospital ID"))

.charlson_diabetes_pts <- case_when(mortality_analytic$comorb_diabetes_end_organ == 1L ~ 2L, mortality_analytic$comorb_diabetes == 1L ~ 1L, TRUE ~ 0L)
.charlson_liver_pts <- case_when(mortality_analytic$comorb_liver_severe == 1L ~ 3L, mortality_analytic$comorb_liver_mild == 1L ~ 1L, TRUE ~ 0L)
.charlson_malignancy_pts <- case_when(mortality_analytic$comorb_metastatic_tumour == 1L ~ 6L, mortality_analytic$comorb_malignancy_active == 1L ~ 2L, TRUE ~ 0L)
.charlson_hiv_pts <- if_else(mortality_analytic$comorb_hiv_on_art == 1L | mortality_analytic$comorb_hiv_no_art == 1L | mortality_analytic$comorb_hiv_cd4_low == 1L, 6L, 0L, missing = 0L)
.charlson_weight1_pts <- mortality_analytic %>%
  select(comorb_chf, comorb_cerebrovascular, comorb_dementia, comorb_copd, comorb_connective_tissue, comorb_peptic_ulcer) %>%
  mutate(across(everything(), ~ if_else(.x == 1L, 1L, 0L, missing = 0L))) %>%
  rowSums()
.charlson_hemi_pts <- if_else(mortality_analytic$comorb_hemi_paraplegia == 1L, 2L, 0L, missing = 0L)
.charlson_renal_pts <- if_else(mortality_analytic$comorb_renal_disease == 1L, 2L, 0L, missing = 0L)
mortality_analytic$charlson_score_combined <- .charlson_weight1_pts + .charlson_diabetes_pts + .charlson_liver_pts +
  .charlson_malignancy_pts + .charlson_hiv_pts + .charlson_hemi_pts + .charlson_renal_pts
rm(.charlson_diabetes_pts, .charlson_liver_pts, .charlson_malignancy_pts, .charlson_hiv_pts,
   .charlson_weight1_pts, .charlson_hemi_pts, .charlson_renal_pts)

adjustment_covars <- c("age_years", "sex_at_birth", "Level_of_care", "abx_current_therapy", "charlson_score_combined")
descriptive_only_covars <- c("infection_origin_composite", "hours_admission_to_sample1", "pitt_score_partial",
                              "ward_type_icu", "specialty_group", "icu_48h_prior", "surgery_prior_to_bsi")

missingness_report <- purrr::map_dfr(c(adjustment_covars, descriptive_only_covars), function(v) {
  x <- mortality_analytic[[v]]
  tibble(variable = v, n = length(x), n_missing = sum(is.na(x)), pct_missing = round(100 * mean(is.na(x)), 1))
}) %>%
  mutate(
    decision = case_when(
      variable %in% adjustment_covars ~ "Multiply imputed (MICE), used for adjustment",
      pct_missing >= 40 ~ "Missingness >=40% -- descriptive comparison only, NOT used for adjustment",
      TRUE ~ "Excluded from adjustment for reasons other than missingness (e.g. possible post-culture/collider timing) -- descriptive comparison only"
    )
  )
write_xlsx(missingness_report, "Mortality_missingness_supplementary_table.xlsx")

# -----------------------------------------------------------------------
# 9.4: Multiple imputation (MICE) of the adjustment covariates
# -----------------------------------------------------------------------
# Predictive mean matching for continuous variables, logistic/polytomous
# regression for categorical (mice defaults); m=20 imputations given the
# modest missingness (<=33% per SECTION 9.3). Outcome (event28, time28) and
# group are included as imputation-model predictors (recommended practice --
# omitting the outcome from the imputation model biases estimates toward
# the null) but are never themselves imputed (0% missing) -- and must be
# carried through the completed data for SECTION 9.7's Cox models, which
# need them.
set.seed(2026)
mice_input <- mortality_analytic %>%
  mutate(
    sex_at_birth = factor(sex_at_birth), Level_of_care = factor(Level_of_care),
    abx_current_therapy = factor(abx_current_therapy)
  ) %>%
  select(mortality_group, event28, time28, all_of(adjustment_covars))
mortality_mice <- mice(mice_input, m = 20, method = "pmm", printFlag = FALSE,
                        defaultMethod = c("pmm", "logreg", "polyreg", "polr"))

# -----------------------------------------------------------------------
# 9.5: Table 1 -- descriptive statistics by group, Lancet-style long table
# -----------------------------------------------------------------------
# Continuous: mean (SD), Welch t-test. Categorical: n (%), chi-squared (or
# Fisher's exact where any expected cell count <5). Two pre-specified
# comparisons per the request: Resistant vs Susceptible, and Resistant vs
# Uninfected (Susceptible vs Uninfected is the third pre-specified
# comparison for the mortality models themselves, SECTION 9.7, but is not
# duplicated as a third Table 1 p-value column, per what was asked here).
table1_vars_cont <- c("age_years", "charlson_score_combined", "LOS_days")
table1_vars_cat <- c("sex_at_birth", "ward_type_icu", "specialty_group", "Level_of_care", "abx_current_therapy",
                      "infection_origin_composite", "pitt_score_partial")

fmt_pct <- function(n, d) if_else(d == 0, "-", sprintf("%d (%.1f%%)", n, 100 * n / d))

table1_cont_row <- function(varname) {
  x <- mortality_analytic[[varname]]
  grp <- mortality_analytic$mortality_group
  gstat <- function(lvl) sprintf("%.1f (%.1f)", mean(x[grp == lvl], na.rm = TRUE), sd(x[grp == lvl], na.rm = TRUE))
  p_rs <- tryCatch(t.test(x[grp == "Resistant BSI"], x[grp == "Susceptible BSI"])$p.value, error = function(e) NA_real_)
  p_ru <- tryCatch(t.test(x[grp == "Resistant BSI"], x[grp == "Uninfected control"])$p.value, error = function(e) NA_real_)
  tibble(variable = varname, level = "mean (SD)", n_missing = sum(is.na(x)),
         `Uninfected control` = gstat("Uninfected control"), `Susceptible BSI` = gstat("Susceptible BSI"),
         `Resistant BSI` = gstat("Resistant BSI"),
         p_resistant_vs_susceptible = p_rs, p_resistant_vs_uninfected = p_ru)
}

table1_cat_rows <- function(varname) {
  x <- factor(mortality_analytic[[varname]])
  grp <- mortality_analytic$mortality_group
  lvls <- levels(x)
  tab <- table(x, grp, useNA = "no")
  test_p <- function(gA, gB) {
    sub <- tab[, c(gA, gB), drop = FALSE]
    sub <- sub[rowSums(sub) > 0, , drop = FALSE]
    if (nrow(sub) < 2) return(NA_real_)
    exp_min <- suppressWarnings(min(chisq.test(sub)$expected))
    if (is.na(exp_min) || exp_min < 5) {
      tryCatch(fisher.test(sub, workspace = 2e7)$p.value, error = function(e) NA_real_)
    } else {
      chisq.test(sub)$p.value
    }
  }
  p_rs <- test_p("Resistant BSI", "Susceptible BSI")
  p_ru <- test_p("Resistant BSI", "Uninfected control")
  purrr::map_dfr(lvls, function(l) {
    tibble(variable = varname, level = l, n_missing = sum(is.na(x)),
           `Uninfected control` = fmt_pct(tab[l, "Uninfected control"], sum(tab[, "Uninfected control"])),
           `Susceptible BSI` = fmt_pct(tab[l, "Susceptible BSI"], sum(tab[, "Susceptible BSI"])),
           `Resistant BSI` = fmt_pct(tab[l, "Resistant BSI"], sum(tab[, "Resistant BSI"])),
           p_resistant_vs_susceptible = if_else(l == lvls[1], p_rs, NA_real_),
           p_resistant_vs_uninfected = if_else(l == lvls[1], p_ru, NA_real_))
  })
}

mortality_table1 <- bind_rows(
  purrr::map_dfr(table1_vars_cont, table1_cont_row),
  purrr::map_dfr(table1_vars_cat, table1_cat_rows)
) %>%
  mutate(across(c(p_resistant_vs_susceptible, p_resistant_vs_uninfected), ~ round(.x, 4)))
write_xlsx(mortality_table1, "Mortality_table1_descriptive.xlsx")

# -----------------------------------------------------------------------
# SECTION 1 (manuscript): Context-overview descriptive tables. Placed here,
# physically, because the outcome-linked block reuses mortality_analytic
# (built in SECTION 9, above) -- conceptually these tables belong to the
# manuscript's overview/context section, not the mortality section.
# -----------------------------------------------------------------------
fmt_mean_sd <- function(x) sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))

# Level_of_care lives on Hospital_characteristics (facility-level), not
# Main_data_CaseReport directly -- joined here via facility_id_clean, same
# pattern used at SECTIONS 9.3/11/13. id_outcome_conclusion (for the
# contamination/culture-negative rows below) is already present on
# Main_data_CaseReport (added alongside has_CoNS_species, SECTION 2.6) --
# NOT re-joined fresh here, since a redundant join on the same column name
# would get silently suffixed to id_outcome_conclusion.x/.y by left_join(),
# and the two lines below would then silently read NULL via $ rather than
# erroring.
case_with_loc <- Main_data_CaseReport %>%
  left_join(Hospital_characteristics %>% select(`Hospital ID`, Level_of_care), by = c(facility_id_clean = "Hospital ID"))

# --- (A) Overall enrolled population (9k), simple/compact ---
overview_overall <- tibble(
  Characteristic = c(
    "Enrolled patients, n", "Age, years, mean (SD)", "Female sex, n (%)",
    "ICU ward at enrolment, n (%)", "Surgical specialty, n (%)",
    "Tertiary-level facility, n (%)", "Any comorbidity (Charlson>0), n (%)",
    "Culture-confirmed BSI, n (%)", "No culture-confirmed BSI, n (%)"
  ),
  Value = c(
    as.character(nrow(case_with_loc)),
    fmt_mean_sd(case_with_loc$age_years),
    fmt_pct(sum(case_with_loc$sex_at_birth == 1, na.rm = TRUE), sum(!is.na(case_with_loc$sex_at_birth))),
    fmt_pct(sum(case_with_loc$ward_type_icu_num == 1, na.rm = TRUE), sum(!is.na(case_with_loc$ward_type_icu_num))),
    fmt_pct(sum(case_with_loc$specialty_group == "Surgical", na.rm = TRUE), sum(!is.na(case_with_loc$specialty_group))),
    fmt_pct(sum(toupper(case_with_loc$Level_of_care) == "TERTIARY", na.rm = TRUE), sum(!is.na(case_with_loc$Level_of_care))),
    fmt_pct(sum(case_with_loc$charlson_score > 0, na.rm = TRUE), sum(!is.na(case_with_loc$charlson_score))),
    fmt_pct(sum(case_with_loc$BSI_positive == 1, na.rm = TRUE), nrow(case_with_loc)),
    fmt_pct(sum(case_with_loc$BSI_positive == 0, na.rm = TRUE), nrow(case_with_loc))
  )
)

# --- (B) Confirmed BSI cases: susceptible vs resistant, two-column, Lancet-
# style. Demographic/clinical rows use the FULL confirmed-BSI population in
# Main_data_CaseReport (n=487); LOS and 28-day mortality rows require
# outcome-form linkage and use mortality_analytic (n=400 linked), each row's
# own denominator stated explicitly per the report's established
# convention. Continuous: mean (SD), Welch t-test. Categorical: n (%),
# chi-squared (Fisher's exact where any expected cell <5), matching SECTION
# 9.5's method exactly.
bsi_confirmed <- case_with_loc %>% filter(BSI_positive == 1, !is.na(resistant_BSI_positive))
bsi_confirmed$resistance_group <- if_else(bsi_confirmed$resistant_BSI_positive == 1, "Resistant BSI", "Susceptible BSI")

sec1_cont_row_confirmed <- function(varname) {
  x <- bsi_confirmed[[varname]]; grp <- bsi_confirmed$resistance_group
  p <- tryCatch(t.test(x[grp == "Resistant BSI"], x[grp == "Susceptible BSI"])$p.value, error = function(e) NA_real_)
  tibble(variable = varname, n_susceptible = sum(!is.na(x[grp == "Susceptible BSI"])), `Susceptible BSI` = fmt_mean_sd(x[grp == "Susceptible BSI"]),
         n_resistant = sum(!is.na(x[grp == "Resistant BSI"])), `Resistant BSI` = fmt_mean_sd(x[grp == "Resistant BSI"]), p_value = round(p, 4))
}
sec1_cat_rows_confirmed <- function(varname) {
  x <- factor(bsi_confirmed[[varname]]); grp <- bsi_confirmed$resistance_group
  tab <- table(x, grp, useNA = "no")
  p <- tryCatch({
    exp_min <- suppressWarnings(min(chisq.test(tab)$expected))
    if (is.na(exp_min) || exp_min < 5) fisher.test(tab, workspace = 2e7)$p.value else chisq.test(tab)$p.value
  }, error = function(e) NA_real_)
  purrr::map_dfr(levels(x), function(l) {
    tibble(variable = varname, level = l,
           n_susceptible = sum(tab[, "Susceptible BSI"]), `Susceptible BSI` = fmt_pct(tab[l, "Susceptible BSI"], sum(tab[, "Susceptible BSI"])),
           n_resistant = sum(tab[, "Resistant BSI"]), `Resistant BSI` = fmt_pct(tab[l, "Resistant BSI"], sum(tab[, "Resistant BSI"])),
           p_value = if_else(l == levels(x)[1], round(p, 4), NA_real_))
  })
}
sec1_demographic_table <- bind_rows(
  purrr::map_dfr(c("age_years", "charlson_score"), sec1_cont_row_confirmed),
  purrr::map_dfr(c("sex_at_birth", "ward_type_icu", "specialty_group", "Level_of_care", "infection_origin_composite"), sec1_cat_rows_confirmed)
)

sec1_outcome_pop <- mortality_analytic %>% filter(mortality_group %in% c("Susceptible BSI", "Resistant BSI"))
sec1_outcome_pop$resistance_group <- sec1_outcome_pop$mortality_group
sec1_los_row <- {
  x <- sec1_outcome_pop$LOS_days; grp <- sec1_outcome_pop$resistance_group
  p <- tryCatch(t.test(x[grp == "Resistant BSI"], x[grp == "Susceptible BSI"])$p.value, error = function(e) NA_real_)
  tibble(variable = "LOS_days (enrolment to discharge/death, truncated 28d)",
         n_susceptible = sum(!is.na(x[grp == "Susceptible BSI"])), `Susceptible BSI` = fmt_mean_sd(x[grp == "Susceptible BSI"]),
         n_resistant = sum(!is.na(x[grp == "Resistant BSI"])), `Resistant BSI` = fmt_mean_sd(x[grp == "Resistant BSI"]), p_value = round(p, 4))
}
sec1_mortality_row <- {
  tab <- table(sec1_outcome_pop$event28, sec1_outcome_pop$resistance_group)
  p <- tryCatch({
    exp_min <- suppressWarnings(min(chisq.test(tab)$expected))
    if (is.na(exp_min) || exp_min < 5) fisher.test(tab)$p.value else chisq.test(tab)$p.value
  }, error = function(e) NA_real_)
  tibble(variable = "28-day mortality, n (%)",
         n_susceptible = sum(sec1_outcome_pop$resistance_group == "Susceptible BSI"),
         `Susceptible BSI` = fmt_pct(sum(sec1_outcome_pop$event28 == 1 & sec1_outcome_pop$resistance_group == "Susceptible BSI", na.rm = TRUE), sum(sec1_outcome_pop$resistance_group == "Susceptible BSI")),
         n_resistant = sum(sec1_outcome_pop$resistance_group == "Resistant BSI"),
         `Resistant BSI` = fmt_pct(sum(sec1_outcome_pop$event28 == 1 & sec1_outcome_pop$resistance_group == "Resistant BSI", na.rm = TRUE), sum(sec1_outcome_pop$resistance_group == "Resistant BSI")),
         p_value = round(p, 4))
}
sec1_outcome_table <- bind_rows(sec1_los_row, sec1_mortality_row)

sec1_bsi_susceptible_vs_resistant <- bind_rows(
  sec1_demographic_table %>% mutate(denominator = paste0("Confirmed BSI, n=", nrow(bsi_confirmed), " (Main_data_CaseReport)")),
  sec1_outcome_table %>% mutate(denominator = paste0("Confirmed BSI with linked 28-day outcome, n=", nrow(sec1_outcome_pop), " (mortality_analytic)"))
)

# --- (C) Matched (uninfected) control cohort: brief description ---
matched_control_pool <- Main_data_hospitalCohortEnrolment_clean %>% filter(matched_control == 1)
matched_control_description <- tibble(
  Characteristic = c("Pooled matched controls, n", "Age, years, mean (SD)", "Female sex, n (%)", "Facilities represented, n"),
  Value = c(
    as.character(nrow(matched_control_pool)),
    fmt_mean_sd(matched_control_pool$age_years),
    fmt_pct(sum(matched_control_pool$sex_at_birth == 1, na.rm = TRUE), sum(!is.na(matched_control_pool$sex_at_birth))),
    as.character(n_distinct(matched_control_pool$facility_id))
  )
)

write_xlsx(
  list(
    overall_population = overview_overall,
    confirmed_BSI_susceptible_vs_resistant = sec1_bsi_susceptible_vs_resistant,
    matched_control_description = matched_control_description
  ),
  "Section1_context_overview_tables.xlsx"
)

# -----------------------------------------------------------------------
# 9.6: Kaplan-Meier curves, 28-day survival by group
# -----------------------------------------------------------------------
km_fit <- survfit(Surv(time28, event28) ~ mortality_group, data = mortality_analytic)
km_logrank <- survdiff(Surv(time28, event28) ~ mortality_group, data = mortality_analytic)
km_logrank_p <- pchisq(km_logrank$chisq, length(km_logrank$n) - 1, lower.tail = FALSE)

km_tidy <- broom::tidy(km_fit) %>%
  mutate(mortality_group = sub("^mortality_group=", "", strata))
km_tidy0 <- km_tidy %>% distinct(mortality_group) %>%
  mutate(time = 0, estimate = 1, conf.low = 1, conf.high = 1)
km_plot_data <- bind_rows(km_tidy0, km_tidy) %>% arrange(mortality_group, time)

KM_INK <- "#142127"; KM_INK_SOFT <- "#4D616A"; KM_LINE <- "#D3DCDD"
km_colors <- c("Uninfected control" = "#17766E", "Susceptible BSI" = "#A8541D", "Resistant BSI" = "#7A1F2B")

km_plot <- ggplot(km_plot_data, aes(x = time, y = 1 - estimate, color = mortality_group, fill = mortality_group)) +
  geom_ribbon(aes(ymin = 1 - conf.low, ymax = 1 - conf.high), alpha = 0.12, color = NA) +
  geom_step(linewidth = 0.9) +
  scale_color_manual(values = km_colors, name = NULL) +
  scale_fill_manual(values = km_colors, guide = "none") +
  scale_x_continuous(breaks = seq(0, 28, 7), limits = c(0, 28), expand = c(0, 0)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  labs(
    title = "28-day all-cause mortality by infection and resistance status",
    subtitle = sprintf("Log-rank p %s", if_else(km_logrank_p < 0.0001, "<0.0001", sprintf("= %.4f", km_logrank_p))),
    x = "Days since enrolment", y = "Cumulative mortality"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(), panel.grid.major = element_line(color = KM_LINE, linewidth = 0.4),
    legend.position = "top", plot.title = element_text(face = "bold"),
    axis.title = element_text(color = KM_INK_SOFT)
  )
ggsave("Mortality_KM_28day_by_group.png", km_plot, width = 8, height = 6, dpi = 500)

# Large-font, title-free, 500dpi TIFF variant for the Figure 4 multi-panel layout
# (titles are handled by the figure caption there; PNG readability was too small).
km_plot_big <- ggplot(km_plot_data, aes(x = time, y = 1 - estimate, color = mortality_group, fill = mortality_group)) +
  geom_ribbon(aes(ymin = 1 - conf.low, ymax = 1 - conf.high), alpha = 0.12, color = NA) +
  geom_step(linewidth = 1.1) +
  scale_color_manual(values = km_colors, name = NULL) +
  scale_fill_manual(values = km_colors, guide = "none") +
  scale_x_continuous(breaks = seq(0, 28, 7), limits = c(0, 28), expand = c(0, 0)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  labs(x = "Days since enrolment", y = "Cumulative mortality") +
  theme_minimal(base_size = 24) +
  theme(
    panel.grid.minor = element_blank(), panel.grid.major = element_line(color = KM_LINE, linewidth = 0.5),
    legend.position = "top", legend.text = element_text(size = 21),
    axis.title = element_text(color = KM_INK_SOFT, size = 23),
    axis.text = element_text(size = 20)
  )
ggsave("Mortality_KM_28day_by_group_bigfont.tiff", km_plot_big, width = 11, height = 7, dpi = 500, compression = "lzw")
ggsave("Mortality_KM_28day_by_group_bigfont.png", km_plot_big, width = 11, height = 7, dpi = 500)

# -----------------------------------------------------------------------
# 9.7-9.8: IPW-weighted marginal structural Cox models + adjusted risks
# -----------------------------------------------------------------------
# Stabilised IPW from a pairwise logistic propensity model (group membership
# ~ adjustment covariates) per comparison -- pairwise rather than a single
# 3-arm multinomial model, both simpler and more stable given the
# Susceptible BSI arm is small (n=28). Propensities truncated to [0.01,0.99]
# and weights to their [1st,99th] percentile before use, standard practice
# against extreme-weight instability. Fit on each of the 20 imputed
# datasets, pooled via Rubin's rules (between + within-imputation variance)
# on top of hospital cluster-robust (sandwich) variance within each fit --
# combining SECTION 3's design-based clustering logic with MI, as the
# survey calibrated_weight itself was not combined in (see note below).
#
# NOTE on survey weights: the methods spec allows combining calibrated
# survey weights with IPW "where possible", deferring to confounding
# adjustment alone "if combined weighting leads to instability". Given
# SECTIONS 3.3/4's own finding of single- or few-hospital oblast strata,
# multiplying calibrated_weight into these already-small group-specific
# models would concentrate extreme combined weights on a handful of
# clusters -- exactly the instability the spec anticipates. IPW + cluster-
# robust SE (hospital as cluster) is used as the primary analysis; national
# extrapolation (SECTION 9.10) still uses the full survey-weighted estimates
# from SECTIONS 3-4.
mortality_comparisons <- list(
  c(ref = "Uninfected control", comp = "Susceptible BSI"),
  c(ref = "Uninfected control", comp = "Resistant BSI"),
  c(ref = "Susceptible BSI", comp = "Resistant BSI")
)
facility_lookup_mortality <- mortality_analytic$facility_id

fit_ipw_cox_pooled <- function(imp_obj, group_ref, group_comp, covars, facility_ids) {
  form_ps <- as.formula(paste("group_bin ~", paste(covars, collapse = " + ")))
  ests <- numeric(0); vars <- numeric(0)
  for (i in seq_len(imp_obj$m)) {
    fit_result <- tryCatch({
      di <- complete(imp_obj, i)
      di$facility_id <- facility_ids
      di <- di %>% dplyr::filter(mortality_group %in% c(group_ref, group_comp)) %>%
        mutate(group_bin = as.integer(mortality_group == group_comp))
      ps_mod <- glm(form_ps, data = di, family = binomial())
      ps <- pmin(pmax(predict(ps_mod, type = "response"), 0.01), 0.99)
      p_marg <- mean(di$group_bin)
      w <- if_else(di$group_bin == 1, p_marg / ps, (1 - p_marg) / (1 - ps))
      q <- quantile(w, c(0.01, 0.99), na.rm = TRUE)
      di$ipw <- pmin(pmax(w, q[1]), q[2])
      cx <- coxph(Surv(time28, event28) ~ group_bin + cluster(facility_id), data = di, weights = ipw)
      c(est = unname(as.numeric(coef(cx)["group_bin"])), var = unname(as.numeric(vcov(cx)["group_bin", "group_bin"])))
    }, error = function(e) NULL)
    if (is.null(fit_result) || !is.finite(fit_result["est"]) || !is.finite(fit_result["var"])) next
    ests <- c(ests, fit_result["est"]); vars <- c(vars, fit_result["var"])
  }
  m <- length(ests)
  if (m == 0) {
    return(tibble(comparison = paste(group_comp, "vs", group_ref), log_hr = NA_real_, se_log_hr = NA_real_,
                  hazard_ratio = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_, m_used = 0L))
  }
  qbar <- mean(ests); ubar <- mean(vars); b <- if (m > 1) var(ests) else 0
  tvar <- ubar + (1 + 1 / m) * b
  se <- sqrt(tvar)
  tibble(comparison = paste(group_comp, "vs", group_ref), log_hr = qbar, se_log_hr = se,
         hazard_ratio = exp(qbar), hr_lower = exp(qbar - 1.96 * se), hr_upper = exp(qbar + 1.96 * se), m_used = m)
}

compute_ipw_risk_rd <- function(di_full, group_ref, group_comp, covars, facility_ids) {
  na_result <- c(risk_ref = NA_real_, risk_comp = NA_real_, rd = NA_real_)
  result <- tryCatch({
    form_ps <- as.formula(paste("group_bin ~", paste(covars, collapse = " + ")))
    di <- di_full %>% mutate(facility_id = facility_ids) %>%
      dplyr::filter(mortality_group %in% c(group_ref, group_comp)) %>%
      mutate(group_bin = as.integer(mortality_group == group_comp))
    ps_mod <- glm(form_ps, data = di, family = binomial())
    ps <- pmin(pmax(predict(ps_mod, type = "response"), 0.01), 0.99)
    p_marg <- mean(di$group_bin)
    w <- if_else(di$group_bin == 1, p_marg / ps, (1 - p_marg) / (1 - ps))
    q <- quantile(w, c(0.01, 0.99), na.rm = TRUE)
    di$ipw <- pmin(pmax(w, q[1]), q[2])
    sf <- survfit(Surv(time28, event28) ~ group_bin, weights = ipw, data = di)
    s28 <- summary(sf, times = 28, extend = TRUE)
    strata_lab <- as.character(s28$strata)
    risk0 <- 1 - s28$surv[grepl("group_bin=0", strata_lab)]
    risk1 <- 1 - s28$surv[grepl("group_bin=1", strata_lab)]
    if (length(risk0) == 0 || length(risk1) == 0) na_result else c(risk_ref = risk0[1], risk_comp = risk1[1], rd = risk1[1] - risk0[1])
  }, error = function(e) na_result)
  if (any(!is.finite(result))) na_result else result
}

set.seed(42)
N_BOOT <- 200
facs_all <- unique(facility_lookup_mortality)
# Generalised beyond SECTIONS 9.7-9.10 (defaults reproduce those exact primary
# results unchanged) so SECTION 9.15's resistance-mechanism subgroups can
# reuse the same bootstrap engine with their own imputation object/facility
# vector/covariate set/replicate count instead of duplicating this loop.
bootstrap_rd <- function(group_ref, group_comp, imp_obj = mortality_mice, facility_ids = facility_lookup_mortality,
                          covars = adjustment_covars, n_boot = N_BOOT) {
  facs <- unique(facility_ids)
  out <- matrix(NA_real_, nrow = n_boot, ncol = 3, dimnames = list(NULL, c("risk_ref", "risk_comp", "rd")))
  for (b in seq_len(n_boot)) {
    di <- complete(imp_obj, sample.int(imp_obj$m, 1))
    di$facility_id <- facility_ids
    boot_facs <- sample(facs, length(facs), replace = TRUE)
    boot_rows <- purrr::map_dfr(boot_facs, function(f) di[di$facility_id == f, , drop = FALSE])
    r <- tryCatch(compute_ipw_risk_rd(boot_rows, group_ref, group_comp, covars, boot_rows$facility_id),
                  error = function(e) c(risk_ref = NA_real_, risk_comp = NA_real_, rd = NA_real_))
    out[b, ] <- r
  }
  out
}

mortality_msm_results <- purrr::map_dfr(mortality_comparisons, function(cmp) {
  hr_row <- fit_ipw_cox_pooled(mortality_mice, cmp["ref"], cmp["comp"], adjustment_covars, facility_lookup_mortality)
  point <- compute_ipw_risk_rd(complete(mortality_mice, 1), cmp["ref"], cmp["comp"], adjustment_covars, facility_lookup_mortality)
  boot_mat <- bootstrap_rd(cmp["ref"], cmp["comp"])
  boot_ok <- boot_mat[complete.cases(boot_mat), , drop = FALSE]
  hr_row %>% mutate(
    adjusted_risk_ref_pct = point["risk_ref"] * 100, adjusted_risk_comp_pct = point["risk_comp"] * 100,
    risk_difference_pct = point["rd"] * 100,
    rd_lower_pct = quantile(boot_ok[, "rd"], 0.025, na.rm = TRUE) * 100,
    rd_upper_pct = quantile(boot_ok[, "rd"], 0.975, na.rm = TRUE) * 100,
    n_bootstrap_ok = nrow(boot_ok)
  )
}) %>%
  mutate(across(where(is.numeric) & !c(m_used, n_bootstrap_ok), ~ round(.x, 4)))

write_xlsx(mortality_msm_results, "Mortality_MSM_causal_HR_and_risk_differences.xlsx")

# -----------------------------------------------------------------------
# 9.9: Competing-risk analysis -- in-hospital mortality, discharge competing
# -----------------------------------------------------------------------
# Fine-Gray subdistribution hazard (survival::finegray + weighted coxph),
# same 3 pairwise comparisons, unweighted by IPW (descriptive competing-risk
# contrast, not a causal MSM estimate -- the primary causal estimates are
# SECTION 9.7-9.8's 28-day MSM) but still cluster-robust by hospital.
cr_data <- mortality_analytic %>% filter(!is.na(event_inhosp), !is.na(time_inhosp))

# cr_dat defaults to the primary 3-arm cr_data; SECTION 9.15 passes its own
# subgroup-restricted competing-risk dataset instead.
fit_finegray_pairwise <- function(group_ref, group_comp, cr_dat = cr_data) {
  di <- cr_dat %>% filter(mortality_group %in% c(group_ref, group_comp)) %>%
    mutate(mortality_group = droplevels(factor(mortality_group, levels = c(group_ref, group_comp))))
  fg <- finegray(Surv(time_inhosp, event_inhosp) ~ .,
                 data = di %>% select(time_inhosp, event_inhosp, mortality_group, facility_id),
                 etype = "death")
  mod <- coxph(Surv(fgstart, fgstop, fgstatus) ~ mortality_group + cluster(facility_id), data = fg, weights = fgwt)
  s <- summary(mod)
  tibble(
    comparison = paste(group_comp, "vs", group_ref),
    n = s$n, n_events = s$nevent,
    subdistribution_hr = s$conf.int[1, 1], sHR_lower = s$conf.int[1, 3], sHR_upper = s$conf.int[1, 4],
    p_value = s$coefficients[1, ncol(s$coefficients)]
  )
}
mortality_competing_risk_results <- purrr::map_dfr(mortality_comparisons, function(cmp) {
  tryCatch(fit_finegray_pairwise(cmp["ref"], cmp["comp"]), error = function(e) tibble(comparison = paste(cmp["comp"], "vs", cmp["ref"]), n = NA, n_events = NA, subdistribution_hr = NA, sHR_lower = NA, sHR_upper = NA, p_value = NA)) %>%
    mutate(across(where(is.numeric) & !c(n, n_events), ~ round(.x, 4)))
})
write_xlsx(mortality_competing_risk_results, "Mortality_competing_risk_finegray.xlsx")

# -----------------------------------------------------------------------
# 9.10: National excess deaths
# -----------------------------------------------------------------------
# Adjusted risk differences (SECTION 9.8) applied to nationally-estimated
# annual BSI/resistant-BSI case counts, from SECTION 4's survey-weighted
# incidence (per 1,000 admissions) x the national admissions total (ApCL
# sum across the 119-facility frame, SECTION 3.2 -- matches the externally
# reported national inpatient admissions figure, 819,368). Uncertainty
# combines only the risk-difference bootstrap distribution (SECTION 9.8)
# with the point-estimate national case counts -- NOT a full joint
# bootstrap across the prevalence/incidence and mortality estimation
# stages, which is beyond this pipeline's scope; treat the resulting
# interval as a lower bound on total uncertainty.
national_admissions_total <- sum(sampling_frame$ApCL)
rate_bsi_national <- incidence_primary_table %>%
  filter(measure == "Any confirmed BSI (overall)", exposure_type == "Total_inpatient_admissions") %>%
  pull(calibrated_rate_per1000)
rate_resistant_national <- incidence_primary_table %>%
  filter(measure == "Any drug-resistant BSI", exposure_type == "Total_inpatient_admissions") %>%
  pull(calibrated_rate_per1000)
national_bsi_cases <- rate_bsi_national / 1000 * national_admissions_total
national_resistant_cases <- rate_resistant_national / 1000 * national_admissions_total
national_susceptible_cases <- national_bsi_cases - national_resistant_cases

national_case_counts <- tibble(
  quantity = c("National admissions (annual)", "National confirmed BSI cases", "National resistant BSI cases",
               "National susceptible/untested BSI cases"),
  value = round(c(national_admissions_total, national_bsi_cases, national_resistant_cases, national_susceptible_cases))
)

# imp_obj/facility_ids/covars default to the primary analysis's globals;
# SECTION 9.15 passes subgroup-specific versions through to bootstrap_rd().
excess_deaths_from_boot <- function(group_ref, group_comp, national_n, imp_obj = mortality_mice,
                                     facility_ids = facility_lookup_mortality, covars = adjustment_covars) {
  boot_mat <- bootstrap_rd(group_ref, group_comp, imp_obj = imp_obj, facility_ids = facility_ids, covars = covars)
  rd <- boot_mat[, "rd"]; rd <- rd[!is.na(rd)]
  tibble(
    comparison = paste(group_comp, "vs", group_ref), national_cases_applied = round(national_n),
    excess_deaths = round(median(rd) * national_n),
    excess_deaths_lower = round(quantile(rd, 0.025) * national_n),
    excess_deaths_upper = round(quantile(rd, 0.975) * national_n)
  )
}
set.seed(43)
national_excess_deaths <- bind_rows(
  excess_deaths_from_boot("Uninfected control", "Susceptible BSI", national_susceptible_cases),
  excess_deaths_from_boot("Uninfected control", "Resistant BSI", national_resistant_cases),
  excess_deaths_from_boot("Susceptible BSI", "Resistant BSI", national_resistant_cases)
)
write_xlsx(list(national_case_counts = national_case_counts, national_excess_deaths = national_excess_deaths),
           "Mortality_national_excess_deaths.xlsx")

# -----------------------------------------------------------------------
# 9.11: Model assumptions
# -----------------------------------------------------------------------
# Every assumption the SECTION 9.7-9.10 estimates rest on, checked where
# checkable and stated plainly where not:
#
#  1. No unmeasured confounding (exchangeability) -- fundamentally
#     UNTESTABLE from these data. IPW adjusts only for
#     adjustment_covars (age, sex, ward type, specialty, hospital level,
#     prior antibiotic use, comorbidity burden); anything else that
#     differs systematically between groups (e.g. clinical severity at
#     the moment of blood culture, which controls have no equivalent of)
#     is NOT adjusted for and could bias the causal HR/RD in either
#     direction.
#  2. Positivity/overlap -- every patient must have a non-zero, non-one
#     probability of being in EITHER compared group given their
#     covariates. Checked via the propensity score range per comparison
#     below; scores were truncated to [0.01, 0.99] and the resulting
#     weights to their [1st, 99th] percentile specifically because this
#     assumption is imperfect in a comparison this small (Susceptible
#     BSI, n=28) -- truncation trades a small amount of bias for a large
#     reduction in variance, a standard, transparent trade-off, not a
#     free fix.
#  3. Correct propensity-model specification -- a linear-in-covariates
#     logistic model (no interactions/non-linear terms) is assumed
#     sufficient to capture confounding; not formally tested.
#  4. Proportional hazards -- tested via cox.zph() (Schoenfeld residuals)
#     on each comparison's UNWEIGHTED-covariate, IPW-weighted model.
#     Where this rejects (see table below), the reported HR is best read
#     as a weighted-average hazard ratio over the 28-day window, not
#     evidence the relative hazard is constant throughout it.
#  5. Non-informative censoring -- the 28-day outcome uses fixed
#     administrative censoring at day 28 for every patient alive at
#     review (SECTION 9.2), which is non-informative by construction;
#     the alternative (a noisy computed review-date difference) was
#     rejected specifically because it was NOT plausibly non-informative
#     (SECTION 9.2's own diagnostic).
#  6. Multiple imputation (MICE) assumes missingness is Missing At Random
#     (MAR) given the variables in the imputation model -- untestable
#     against the true (unobserved) values; the outcome was included as
#     an imputation-model predictor specifically to make this more
#     plausible (SECTION 9.4).
#  7. Cluster-robust variance (hospital as cluster) requires only that
#     clusters are correctly specified, not independence within a
#     cluster -- satisfied by construction (facility_id is the true PPS
#     sampling unit, SECTION 3.2).
#  8. Fine-Gray subdistribution hazards (SECTION 9.9) do NOT require
#     independence between the competing events (death, discharge) --
#     unlike a cause-specific-hazard approach, this is an assumption the
#     method avoids needing, not one it relies on.
model_assumptions <- tribble(
  ~assumption, ~status,
  "No unmeasured confounding", "Untestable -- adjusts only for measured baseline covariates (SECTION 9.3); interpret causally with caution",
  "Positivity / overlap", "Propensity scores and IPW truncated (1st/99th percentile) where extreme -- see diagnostics table",
  "Propensity model correctly specified", "Assumed (linear logistic, no interactions) -- not formally tested",
  "Proportional hazards (Cox)", "Tested via cox.zph() per comparison -- see diagnostics table; violated for Resistant vs Uninfected",
  "Non-informative censoring", "By construction -- fixed 28-day administrative censoring (SECTION 9.2)",
  "MICE: Missing At Random", "Untestable -- outcome included as an imputation-model predictor to support plausibility",
  "Cluster-robust SE: correct cluster specification", "Satisfied -- facility_id is the true PPS sampling unit",
  "Fine-Gray: independence of competing risks", "NOT required by this method (its advantage over cause-specific hazards)"
)

# imp_obj/facility_ids/covars default to the primary analysis's globals;
# SECTION 9.15 passes subgroup-specific versions.
check_ph_and_positivity <- function(group_ref, group_comp, imp_obj = mortality_mice,
                                     facility_ids = facility_lookup_mortality, covars = adjustment_covars) {
  di <- complete(imp_obj, 1)
  di$facility_id <- facility_ids
  di <- di %>% dplyr::filter(mortality_group %in% c(group_ref, group_comp)) %>%
    mutate(group_bin = as.integer(mortality_group == group_comp))
  form_ps <- as.formula(paste("group_bin ~", paste(covars, collapse = " + ")))
  ps_mod <- glm(form_ps, data = di, family = binomial())
  ps_raw <- predict(ps_mod, type = "response")
  ps <- pmin(pmax(ps_raw, 0.01), 0.99)
  p_marg <- mean(di$group_bin)
  w_raw <- if_else(di$group_bin == 1, p_marg / ps, (1 - p_marg) / (1 - ps))
  cx <- coxph(Surv(time28, event28) ~ group_bin, data = di, weights = pmin(pmax(w_raw, quantile(w_raw, 0.01)), quantile(w_raw, 0.99)))
  zph <- tryCatch(cox.zph(cx), error = function(e) NULL)
  tibble(
    comparison = paste(group_comp, "vs", group_ref),
    n = nrow(di),
    propensity_score_min = round(min(ps_raw), 4), propensity_score_max = round(max(ps_raw), 4),
    pct_propensity_truncated = round(100 * mean(ps_raw < 0.01 | ps_raw > 0.99), 1),
    ipw_min_pretruncation = round(min(w_raw), 2), ipw_max_pretruncation = round(max(w_raw), 2),
    cox_zph_chisq = if (is.null(zph)) NA_real_ else round(zph$table["group_bin", "chisq"], 3),
    cox_zph_p = if (is.null(zph)) NA_real_ else round(zph$table["group_bin", "p"], 4)
  )
}
model_diagnostics <- purrr::map_dfr(mortality_comparisons, ~ check_ph_and_positivity(.x["ref"], .x["comp"]))

write_xlsx(list(assumptions = model_assumptions, diagnostics = model_diagnostics), "Mortality_model_assumptions.xlsx")

# -----------------------------------------------------------------------
# 9.12: Adjusted (model-predicted) survival curves -- overall and by subgroup
# -----------------------------------------------------------------------
# Unlike SECTIONS 9.7-9.8's pairwise IPW (built for the 3 formal causal
# contrasts individually), a single COMBINED weight is needed here to plot
# all three groups together on one standardised curve: a stabilised IPW
# from a multinomial logistic propensity model (group ~ adjustment_covars),
# fit on the first imputed dataset (a representative single completion,
# rather than pooling 20 separate curves -- curves are for illustration of
# the modelled pattern, the formal inferential estimates remain SECTION
# 9.7-9.8's pooled pairwise HRs/risk differences). Subgroup panels use these
# SAME weights (not separately refit per subgroup, which would be unstable
# for a 28-patient arm split further by covariate) -- i.e. the overall
# model's standardisation, viewed within each subgroup.
di_curves <- complete(mortality_mice, 1)
di_curves$facility_id <- facility_lookup_mortality
# ward_type_icu re-attached from the raw (non-imputed, non-modelled) data
# purely for the subgroup-panel SPLIT below -- it is NOT part of the
# propensity/adjustment model (SECTION 9.3 excludes it there as a possible
# post-culture/collider variable). Slicing the already-standardised curves
# by current ward status is a descriptive display choice, not a causal
# adjustment, so the same concern does not apply here.
di_curves$ward_type_icu <- mortality_analytic$ward_type_icu
form_multi <- as.formula(paste("mortality_group ~", paste(adjustment_covars, collapse = " + ")))
multi_ps_mod <- multinom(form_multi, data = di_curves, trace = FALSE)
multi_probs <- predict(multi_ps_mod, type = "probs")
multi_marg <- prop.table(table(di_curves$mortality_group))
own_prob <- multi_probs[cbind(seq_len(nrow(multi_probs)), match(as.character(di_curves$mortality_group), colnames(multi_probs)))]
ipw_combined3 <- as.numeric(multi_marg[as.character(di_curves$mortality_group)]) / own_prob
ipw_combined3 <- pmin(pmax(ipw_combined3, quantile(ipw_combined3, 0.01)), quantile(ipw_combined3, 0.99))
di_curves$ipw_combined3 <- ipw_combined3

MORT_INK <- "#142127"; MORT_INK_SOFT <- "#4D616A"; MORT_LINE <- "#D3DCDD"
mort_colors <- c("Uninfected control" = "#17766E", "Susceptible BSI" = "#A8541D", "Resistant BSI" = "#7A1F2B")

make_adjusted_curve_data <- function(data) {
  # droplevels: a subgroup can genuinely have 0 patients in one of the 3
  # groups (e.g. very few Susceptible BSI neonates) -- fit survfit only on
  # the levels actually present rather than erroring on an empty stratum.
  # Requires >=2 groups present to show anything comparative; a subgroup
  # with only 1 group represented is silently skipped (NULL), not faked.
  data <- data %>% mutate(mortality_group = droplevels(mortality_group))
  if (nlevels(data$mortality_group) < 2) return(NULL)
  sf <- tryCatch(survfit(Surv(time28, event28) ~ mortality_group, weights = ipw_combined3, data = data), error = function(e) NULL)
  if (is.null(sf)) return(NULL)
  td <- broom::tidy(sf) %>% mutate(mortality_group = sub("^mortality_group=", "", strata))
  t0 <- td %>% distinct(mortality_group) %>% mutate(time = 0, estimate = 1, conf.low = 1, conf.high = 1)
  bind_rows(t0, td) %>% arrange(mortality_group, time)
}

plot_adjusted_curve <- function(plot_data, subtitle_txt) {
  ggplot(plot_data, aes(x = time, y = 1 - estimate, color = mortality_group, fill = mortality_group)) +
    geom_ribbon(aes(ymin = 1 - conf.low, ymax = 1 - conf.high), alpha = 0.12, color = NA) +
    geom_step(linewidth = 0.9) +
    scale_color_manual(values = mort_colors, name = NULL) +
    scale_fill_manual(values = mort_colors, guide = "none") +
    scale_x_continuous(breaks = seq(0, 28, 7), limits = c(0, 28), expand = c(0, 0)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
    labs(subtitle = subtitle_txt, x = "Days since enrolment", y = "Adjusted cumulative mortality") +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = MORT_LINE, linewidth = 0.4),
          legend.position = "top", plot.subtitle = element_text(color = MORT_INK_SOFT),
          axis.title = element_text(color = MORT_INK_SOFT))
}

adjusted_overall_data <- make_adjusted_curve_data(di_curves)
adjusted_overall_plot <- plot_adjusted_curve(adjusted_overall_data, "Overall (IPW-standardised)") +
  labs(title = "Model-predicted 28-day mortality by infection and resistance status")
ggsave("Mortality_adjusted_curves_overall.png", adjusted_overall_plot, width = 8, height = 6, dpi = 500)

# =========================================================================
# SUPPLEMENTARY: big-font, title-free 500dpi adjusted-curve figures --
# (1) Cox MSM IPW (reuses adjusted_overall_data above, no recomputation),
# and (2) a new IPW-adjusted Fine-Gray cumulative-incidence curve, primary
# comparison. Fine-Gray's own analytic dataset (cr_data) is complete-case
# restricted to the adjustment covariates and re-weighted by a propensity-
# score IPW (same multinomial-propensity construction as ipw_combined3
# above, refit on cr_data specifically to avoid any row-count mismatch
# between di_curves/mortality_analytic and cr_data's own event_inhosp/
# time_inhosp non-missingness restriction) combined with Fine-Gray's own
# fgwt redistribution weight -- an IPW-adjusted subdistribution/CIF curve,
# not the unadjusted single-covariate fit SECTION 9.9's fit_finegray_
# pairwise() uses for the primary sHR point estimate. Illustrative only
# (single complete-case fit, not pooled across imputations), consistent
# with adjusted_overall_data's own "single completed dataset" convention.
build_bigfont_adjusted_plot <- function(plot_data, colors, ylab = "Adjusted cumulative mortality") {
  ggplot(plot_data, aes(x = time, y = 1 - estimate, color = mortality_group, fill = mortality_group)) +
    geom_ribbon(aes(ymin = 1 - conf.low, ymax = 1 - conf.high), alpha = 0.12, color = NA) +
    geom_step(linewidth = 1.1) +
    scale_color_manual(values = colors, name = NULL) +
    scale_fill_manual(values = colors, guide = "none") +
    scale_x_continuous(breaks = seq(0, 28, 7), limits = c(0, 28), expand = c(0, 0)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
    labs(x = "Days since enrolment", y = ylab) +
    theme_minimal(base_size = 24) +
    theme(
      panel.grid.minor = element_blank(), panel.grid.major = element_line(color = MORT_LINE, linewidth = 0.5),
      legend.position = "top", legend.text = element_text(size = 21),
      axis.title = element_text(color = MORT_INK_SOFT, size = 23), axis.text = element_text(size = 20)
    )
}

ggsave("Mortality_adjusted_curves_CoxMSM_IPW_primary_bigfont.png", build_bigfont_adjusted_plot(adjusted_overall_data, mort_colors), width = 11, height = 7, dpi = 500)
ggsave("Mortality_adjusted_curves_CoxMSM_IPW_primary_bigfont.tiff", build_bigfont_adjusted_plot(adjusted_overall_data, mort_colors), width = 11, height = 7, dpi = 500, compression = "lzw")
cat("Cox MSM IPW primary bigfont curve saved.\n")

fg_ipw_curve_data <- function(cr_dat_full, covars) {
  cc <- cr_dat_full %>%
    mutate(mortality_group = droplevels(mortality_group)) %>%
    filter(if_all(all_of(covars), ~ !is.na(.)))
  form_multi <- as.formula(paste("mortality_group ~", paste(covars, collapse = " + ")))
  mps <- nnet::multinom(form_multi, data = cc, trace = FALSE)
  probs <- predict(mps, type = "probs")
  marg <- prop.table(table(cc$mortality_group))
  own_p <- probs[cbind(seq_len(nrow(probs)), match(as.character(cc$mortality_group), colnames(probs)))]
  ipw3 <- as.numeric(marg[as.character(cc$mortality_group)]) / own_p
  ipw3 <- pmin(pmax(ipw3, quantile(ipw3, 0.01)), quantile(ipw3, 0.99))
  cc$ipw3 <- ipw3
  fg <- finegray(Surv(time_inhosp, event_inhosp) ~ .,
                 data = cc %>% select(time_inhosp, event_inhosp, mortality_group, facility_id, ipw3),
                 etype = "death")
  fg$combined_w <- fg$ipw3 * fg$fgwt
  sf <- survfit(Surv(fgstart, fgstop, fgstatus) ~ mortality_group, weights = combined_w, data = fg)
  td <- broom::tidy(sf) %>% mutate(mortality_group = sub("^mortality_group=", "", strata))
  t0 <- td %>% distinct(mortality_group) %>% mutate(time = 0, estimate = 1, conf.low = 1, conf.high = 1)
  bind_rows(t0, td) %>% arrange(mortality_group, time)
}

fg_primary_data <- fg_ipw_curve_data(cr_data, adjustment_covars)
write_xlsx(fg_primary_data, "Mortality_adjusted_curves_FineGray_primary_data.xlsx")
fg_ylab <- "Adjusted cumulative incidence (death)"
ggsave("Mortality_adjusted_curves_FineGray_primary_bigfont.png", build_bigfont_adjusted_plot(fg_primary_data, mort_colors, fg_ylab), width = 11, height = 7, dpi = 500)
ggsave("Mortality_adjusted_curves_FineGray_primary_bigfont.tiff", build_bigfont_adjusted_plot(fg_primary_data, mort_colors, fg_ylab), width = 11, height = 7, dpi = 500, compression = "lzw")
cat("Fine-Gray IPW-adjusted primary bigfont curve saved.\n")

# subgroup panels: ward type and age group (control variables from the
# adjustment set / core demographics)
subgroup_ward <- purrr::map_dfr(c("ICU", "Non-ICU"), function(lvl) {
  d <- make_adjusted_curve_data(di_curves %>% dplyr::filter(ward_type_icu == lvl))
  if (!is.null(d)) d %>% mutate(panel = lvl)
})
plot_ward <- ggplot(subgroup_ward, aes(x = time, y = 1 - estimate, color = mortality_group, fill = mortality_group)) +
  geom_ribbon(aes(ymin = 1 - conf.low, ymax = 1 - conf.high), alpha = 0.12, color = NA) +
  geom_step(linewidth = 0.9) +
  facet_wrap(~panel) +
  scale_color_manual(values = mort_colors, name = NULL) +
  scale_fill_manual(values = mort_colors, guide = "none") +
  scale_x_continuous(breaks = seq(0, 28, 7), limits = c(0, 28), expand = c(0, 0)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  labs(title = "Model-predicted 28-day mortality, by ward type", x = "Days since enrolment", y = "Adjusted cumulative mortality") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = MORT_LINE, linewidth = 0.4),
        legend.position = "top", strip.text = element_text(face = "bold"), axis.title = element_text(color = MORT_INK_SOFT))
ggsave("Mortality_adjusted_curves_by_wardtype.png", plot_ward, width = 10, height = 5.5, dpi = 500)

# age_group itself was never part of the MICE model (age_years, continuous,
# was used instead) so it isn't in di_curves -- derived here directly from
# the already-imputed age_years using the same NEONATE(<=28 days)/
# PAEDIATRIC(<18y)/ADULT(>=18y) cutoffs as SECTION 7's own age_group.
di_curves <- di_curves %>%
  mutate(age_group_derived = case_when(
    age_years <= 28 / 365.25 ~ "Neonate", age_years < 18 ~ "Paediatric", TRUE ~ "Adult"
  ))
age_group_levels <- c("Neonate", "Paediatric", "Adult")
subgroup_age <- purrr::map_dfr(age_group_levels, function(lvl) {
  d <- make_adjusted_curve_data(di_curves %>% dplyr::filter(age_group_derived == lvl))
  if (!is.null(d)) d %>% mutate(panel = lvl)
})
plot_age <- ggplot(subgroup_age, aes(x = time, y = 1 - estimate, color = mortality_group, fill = mortality_group)) +
  geom_ribbon(aes(ymin = 1 - conf.low, ymax = 1 - conf.high), alpha = 0.12, color = NA) +
  geom_step(linewidth = 0.9) +
  facet_wrap(~ factor(panel, levels = age_group_levels)) +
  scale_color_manual(values = mort_colors, name = NULL) +
  scale_fill_manual(values = mort_colors, guide = "none") +
  scale_x_continuous(breaks = seq(0, 28, 7), limits = c(0, 28), expand = c(0, 0)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  labs(title = "Model-predicted 28-day mortality, by age group", x = "Days since enrolment", y = "Adjusted cumulative mortality") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = MORT_LINE, linewidth = 0.4),
        legend.position = "top", strip.text = element_text(face = "bold"), axis.title = element_text(color = MORT_INK_SOFT))
ggsave("Mortality_adjusted_curves_by_agegroup.png", plot_age, width = 12, height = 5.5, dpi = 500)

# -----------------------------------------------------------------------
# 9.13: Sensitivity model -- IPW combined with the calibrated survey weight
# -----------------------------------------------------------------------
# calibrated_weight (SECTION 3.3) is joined onto EVERY patient here by
# their own facility_id -- not inherited from the SECTION 8 coalesce, which
# only carries calibrated_weight for patients matched to Main_data_CaseReport
# (controls would otherwise be NA). The sampling design weight is a
# property of which hospital a patient is in, so this join is valid for
# cases and controls alike. combined_weight = IPW x calibrated_weight,
# renormalised to mean 1 within each comparison (weight SCALE does not
# affect the Cox partial-likelihood point estimate, but rescaling keeps the
# cluster-robust variance numerically stable). This is the sensitivity
# model the SECTION 9.7-9.8 note anticipated deferring -- run in full here,
# not just described.
calibrated_weight_by_facility <- mortality_analytic %>%
  select(-any_of("calibrated_weight")) %>%
  left_join(facility_weights %>% select(UniqueID, calibrated_weight), by = c(facility_id = "UniqueID")) %>%
  pull(calibrated_weight)

fit_combined_weight_pooled <- function(imp_obj, group_ref, group_comp, covars, facility_ids, cw) {
  form_ps <- as.formula(paste("group_bin ~", paste(covars, collapse = " + ")))
  ests <- numeric(0); vars <- numeric(0)
  for (i in seq_len(imp_obj$m)) {
    fr <- tryCatch({
      di <- complete(imp_obj, i)
      di$facility_id <- facility_ids
      di$cw <- cw
      di <- di %>% dplyr::filter(mortality_group %in% c(group_ref, group_comp), !is.na(cw)) %>%
        mutate(group_bin = as.integer(mortality_group == group_comp))
      ps_mod <- glm(form_ps, data = di, family = binomial())
      ps <- pmin(pmax(predict(ps_mod, type = "response"), 0.01), 0.99)
      p_marg <- mean(di$group_bin)
      ipw <- if_else(di$group_bin == 1, p_marg / ps, (1 - p_marg) / (1 - ps))
      q <- quantile(ipw, c(0.01, 0.99), na.rm = TRUE)
      ipw <- pmin(pmax(ipw, q[1]), q[2])
      combined_w <- ipw * di$cw
      di$combined_w <- combined_w / mean(combined_w)
      cx <- coxph(Surv(time28, event28) ~ group_bin + cluster(facility_id), data = di, weights = combined_w)
      c(est = unname(as.numeric(coef(cx)["group_bin"])), var = unname(as.numeric(vcov(cx)["group_bin", "group_bin"])),
        n = nrow(di))
    }, error = function(e) NULL)
    if (is.null(fr) || !is.finite(fr["est"]) || !is.finite(fr["var"])) next
    ests <- c(ests, fr["est"]); vars <- c(vars, fr["var"])
  }
  m <- length(ests)
  if (m == 0) {
    return(tibble(comparison = paste(group_comp, "vs", group_ref), log_hr = NA_real_, se_log_hr = NA_real_,
                  hazard_ratio = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_, m_used = 0L))
  }
  qbar <- mean(ests); ubar <- mean(vars); b <- if (m > 1) var(ests) else 0
  se <- sqrt(ubar + (1 + 1 / m) * b)
  tibble(comparison = paste(group_comp, "vs", group_ref), log_hr = qbar, se_log_hr = se,
         hazard_ratio = exp(qbar), hr_lower = exp(qbar - 1.96 * se), hr_upper = exp(qbar + 1.96 * se), m_used = m)
}

mortality_msm_combined_weight_results <- purrr::map_dfr(mortality_comparisons, function(cmp) {
  fit_combined_weight_pooled(mortality_mice, cmp["ref"], cmp["comp"], adjustment_covars, facility_lookup_mortality, calibrated_weight_by_facility)
}) %>%
  mutate(across(where(is.numeric) & !m_used, ~ round(.x, 4))) %>%
  mutate(n_missing_calibrated_weight = sum(is.na(calibrated_weight_by_facility)))

write_xlsx(mortality_msm_combined_weight_results, "Mortality_MSM_combined_IPW_survey_weight.xlsx")

# -----------------------------------------------------------------------
# 9.14: IPW propensity-model diagnostics -- coefficients, overlap, balance
# -----------------------------------------------------------------------
# Two things the pairwise IPW pipeline (SECTIONS 9.7-9.8) computed and used
# internally but never surfaced: the propensity (exposure) model itself,
# and whether the resulting weights actually achieved covariate balance.
# Both are reported here on the first imputed dataset (a representative
# single completion, matching SECTION 9.12's convention for the curves --
# the formal HR/RD estimates remain the 20-imputation-pooled figures from
# SECTION 9.7-9.8/9.13).
ps_coef_table <- purrr::map_dfr(mortality_comparisons, function(cmp) {
  di <- complete(mortality_mice, 1) %>%
    dplyr::filter(mortality_group %in% c(cmp["ref"], cmp["comp"])) %>%
    mutate(group_bin = as.integer(mortality_group == cmp["comp"]))
  form_ps <- as.formula(paste("group_bin ~", paste(adjustment_covars, collapse = " + ")))
  ps_mod <- glm(form_ps, data = di, family = binomial())
  broom::tidy(ps_mod, conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(comparison = paste(cmp["comp"], "vs", cmp["ref"]), .before = 1) %>%
    rename(odds_ratio = estimate, ci_lower = conf.low, ci_upper = conf.high)
}) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

# Standardised mean differences (SMD), unweighted vs IPW-weighted, on the
# SAME dummy-coded design matrix the propensity model itself used --
# |SMD| < 0.1 is the conventional threshold for "adequately balanced".
compute_balance <- function(group_ref, group_comp) {
  di <- complete(mortality_mice, 1) %>%
    dplyr::filter(mortality_group %in% c(group_ref, group_comp)) %>%
    mutate(group_bin = as.integer(mortality_group == group_comp))
  form_ps <- as.formula(paste("group_bin ~", paste(adjustment_covars, collapse = " + ")))
  ps_mod <- glm(form_ps, data = di, family = binomial())
  ps <- pmin(pmax(predict(ps_mod, type = "response"), 0.01), 0.99)
  p_marg <- mean(di$group_bin)
  w <- if_else(di$group_bin == 1, p_marg / ps, (1 - p_marg) / (1 - ps))
  q <- quantile(w, c(0.01, 0.99), na.rm = TRUE)
  w <- pmin(pmax(w, q[1]), q[2])
  X <- model.matrix(as.formula(paste("~", paste(adjustment_covars, collapse = "+"))), data = di)[, -1, drop = FALSE]
  wtd_mean <- function(x, g, wt) sum(x[g] * wt[g]) / sum(wt[g])
  wtd_var <- function(x, g, wt) { m <- wtd_mean(x, g, wt); sum(wt[g] * (x[g] - m)^2) / sum(wt[g]) }
  purrr::map_dfr(colnames(X), function(cn) {
    x <- X[, cn]; g1 <- di$group_bin == 1; g0 <- di$group_bin == 0
    sd_u <- sqrt((var(x[g1]) + var(x[g0])) / 2)
    smd_u <- (mean(x[g1]) - mean(x[g0])) / sd_u
    sd_w <- sqrt((wtd_var(x, g1, w) + wtd_var(x, g0, w)) / 2)
    smd_w <- (wtd_mean(x, g1, w) - wtd_mean(x, g0, w)) / sd_w
    tibble(comparison = paste(group_comp, "vs", group_ref), covariate = cn, smd_unweighted = smd_u, smd_weighted = smd_w)
  })
}
mortality_balance_table <- purrr::map_dfr(mortality_comparisons, ~ compute_balance(.x["ref"], .x["comp"])) %>%
  mutate(across(c(smd_unweighted, smd_weighted), ~ round(.x, 4)),
         balanced_after_weighting = abs(smd_weighted) < 0.1)

write_xlsx(list(propensity_model_coefficients = ps_coef_table, covariate_balance = mortality_balance_table),
           "Mortality_IPW_propensity_model_and_balance.xlsx")

## --- balance ("Love") plot ---
# Readable labels for the model.matrix dummy-column names in `covariate`
# (e.g. "Level_of_careTERTIARY") -- kept as a lookup tied to the CURRENT
# adjustment_covars (SECTION 9.3) rather than hardcoded into the plot call,
# so it stays correct if that covariate set changes again.
covariate_label_lookup <- c(
  age_years = "Age (years)",
  sex_at_birth1 = "Sex (female = 1)",
  Level_of_careSECONDARY = "Hospital type: Secondary (vs Primary)",
  Level_of_careTERTIARY = "Hospital type: Tertiary (vs Primary)",
  abx_current_therapy1 = "Currently on antibiotic therapy",
  charlson_score_combined = "Charlson comorbidity score"
)
balance_plot_data <- mortality_balance_table %>%
  mutate(covariate = unname(covariate_label_lookup[covariate])) %>%
  select(comparison, covariate, Unweighted = smd_unweighted, Weighted = smd_weighted) %>%
  tidyr::pivot_longer(c(Unweighted, Weighted), names_to = "stage", values_to = "smd")
balance_plot <- ggplot(balance_plot_data, aes(x = smd, y = covariate, color = stage)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = MORT_LINE, linewidth = 0.5) +
  geom_vline(xintercept = 0, color = MORT_INK_SOFT, linewidth = 0.3) +
  geom_point(size = 2.2) +
  facet_wrap(~comparison, ncol = 1) +
  scale_color_manual(values = c(Unweighted = "#A8541D", Weighted = "#17766E"), name = NULL) +
  labs(title = "Covariate balance before and after IPW", subtitle = "Dashed lines = conventional |SMD| < 0.1 balance threshold",
       x = "Standardised mean difference", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        strip.text = element_text(face = "bold"), plot.subtitle = element_text(color = MORT_INK_SOFT))
ggsave("Mortality_IPW_balance_loveplot.png", balance_plot, width = 8, height = 9, dpi = 500)

## --- propensity score overlap plot ---
ps_overlap_data <- purrr::map_dfr(mortality_comparisons, function(cmp) {
  di <- complete(mortality_mice, 1) %>%
    dplyr::filter(mortality_group %in% c(cmp["ref"], cmp["comp"])) %>%
    mutate(group_bin = as.integer(mortality_group == cmp["comp"]))
  form_ps <- as.formula(paste("group_bin ~", paste(adjustment_covars, collapse = " + ")))
  ps_mod <- glm(form_ps, data = di, family = binomial())
  tibble(comparison = paste(cmp["comp"], "vs", cmp["ref"]), group = as.character(di$mortality_group),
         propensity_score = predict(ps_mod, type = "response"))
})
ps_overlap_plot <- ggplot(ps_overlap_data, aes(x = propensity_score, fill = group, color = group)) +
  geom_density(alpha = 0.35, linewidth = 0.7) +
  facet_wrap(~comparison, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = mort_colors, name = NULL) +
  scale_color_manual(values = mort_colors, guide = "none") +
  scale_x_continuous(limits = c(0, 1)) +
  labs(title = "Propensity score overlap (positivity check)", subtitle = "Kernel density estimate per group",
       x = "Predicted probability of the comparator group", y = "Density") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top", strip.text = element_text(face = "bold"),
        plot.subtitle = element_text(color = MORT_INK_SOFT))
ggsave("Mortality_IPW_propensity_overlap.png", ps_overlap_plot, width = 8, height = 9, dpi = 500)

# -----------------------------------------------------------------------
# 9.15: Resistance-mechanism subgroup analyses (MDR; carbapenem-R among
# Gram-negative BSI; methicillin-R among Staphylococcus BSI)
# -----------------------------------------------------------------------
# SECTIONS 9.1-9.14 pool every resistance mechanism into a single "any
# resistance" Resistant BSI arm (resistant_BSI_positive = resistance to
# >=1 tested agent, organism-agnostic). This section reruns the same
# MICE -> pairwise-IPW-MSM-Cox -> Fine-Gray -> national-excess-deaths ->
# PH/positivity-diagnostics machinery (SECTIONS 9.4/9.7-9.11, generalised
# above to take an explicit imputation object/facility vector/covariate
# set rather than hardcoded globals) on three specific mechanisms, each
# restricted to the organism group the mechanism is clinically defined
# for: MDR (bacterial isolates, any organism), carbapenem resistance
# (Gram-negative isolates only), methicillin resistance (Staphylococcus
# isolates only). Vancomycin resistance (Enterococcus isolates) has only
# 3 VRE patients with 0 deaths among them (SECTION 9.15.4) -- far below
# what any propensity/Cox model can support -- and is reported as a raw
# descriptive count only, per instruction to model a mechanism "only if
# numbers are enough".
#
# run_mortality_subgroup(): one call reproduces the MICE/MSM/Fine-Gray/
# diagnostics/KM-plot bundle for an arbitrary mortality_group factor,
# so it is not hand-duplicated three times. national_cases is an optional
# named list (list("Arm label" = national_case_count)) -- when supplied,
# national excess deaths are computed for any comparison whose `comp` arm
# has an entry; when absent for a given arm (no defensible organism-
# matched national rate exists, e.g. "Carbapenem-S Gram-negative BSI"),
# that comparison is simply skipped from the excess-deaths output, not
# silently estimated from a mismatched denominator.
run_mortality_subgroup <- function(label, analytic_data, comparisons, out_prefix,
                                    national_cases = NULL, seed_mice = 2026, seed_boot = 42, n_boot = 200) {
  set.seed(seed_mice)
  mi <- analytic_data %>%
    mutate(sex_at_birth = factor(sex_at_birth), Level_of_care = factor(Level_of_care),
           abx_current_therapy = factor(abx_current_therapy)) %>%
    select(mortality_group, event28, time28, all_of(adjustment_covars))
  mice_obj <- mice(mi, m = 20, method = "pmm", printFlag = FALSE,
                    defaultMethod = c("pmm", "logreg", "polyreg", "polr"))
  fac_ids <- analytic_data$facility_id

  msm <- purrr::map_dfr(comparisons, function(cmp) {
    hr_row <- fit_ipw_cox_pooled(mice_obj, cmp["ref"], cmp["comp"], adjustment_covars, fac_ids)
    point <- compute_ipw_risk_rd(complete(mice_obj, 1), cmp["ref"], cmp["comp"], adjustment_covars, fac_ids)
    set.seed(seed_boot)
    boot_mat <- bootstrap_rd(cmp["ref"], cmp["comp"], imp_obj = mice_obj, facility_ids = fac_ids,
                              covars = adjustment_covars, n_boot = n_boot)
    boot_ok <- boot_mat[complete.cases(boot_mat), , drop = FALSE]
    hr_row %>% mutate(
      adjusted_risk_ref_pct = point["risk_ref"] * 100, adjusted_risk_comp_pct = point["risk_comp"] * 100,
      risk_difference_pct = point["rd"] * 100,
      rd_lower_pct = if (nrow(boot_ok) > 0) quantile(boot_ok[, "rd"], 0.025, na.rm = TRUE) * 100 else NA_real_,
      rd_upper_pct = if (nrow(boot_ok) > 0) quantile(boot_ok[, "rd"], 0.975, na.rm = TRUE) * 100 else NA_real_,
      n_bootstrap_ok = nrow(boot_ok)
    )
  }) %>% mutate(across(where(is.numeric) & !c(m_used, n_bootstrap_ok), ~ round(.x, 4)))

  cr_sub <- analytic_data %>% filter(!is.na(event_inhosp), !is.na(time_inhosp))
  fg <- purrr::map_dfr(comparisons, function(cmp) {
    tryCatch(
      fit_finegray_pairwise(cmp["ref"], cmp["comp"], cr_dat = cr_sub),
      error = function(e) tibble(comparison = paste(cmp["comp"], "vs", cmp["ref"]), n = NA, n_events = NA,
                                  subdistribution_hr = NA, sHR_lower = NA, sHR_upper = NA, p_value = NA)
    ) %>% mutate(across(where(is.numeric) & !c(n, n_events), ~ round(.x, 4)))
  })

  diagnostics <- purrr::map_dfr(comparisons, ~ check_ph_and_positivity(
    .x["ref"], .x["comp"], imp_obj = mice_obj, facility_ids = fac_ids, covars = adjustment_covars
  ))

  excess <- NULL
  if (!is.null(national_cases)) {
    excess <- purrr::map_dfr(comparisons, function(cmp) {
      nn <- national_cases[[cmp["comp"]]]
      if (is.null(nn)) return(NULL)
      excess_deaths_from_boot(cmp["ref"], cmp["comp"], nn, imp_obj = mice_obj, facility_ids = fac_ids, covars = adjustment_covars)
    })
  }

  km_fit_s <- survfit(Surv(time28, event28) ~ mortality_group, data = analytic_data)
  km_logrank_s <- survdiff(Surv(time28, event28) ~ mortality_group, data = analytic_data)
  km_p_s <- pchisq(km_logrank_s$chisq, length(km_logrank_s$n) - 1, lower.tail = FALSE)
  km_tidy_s <- broom::tidy(km_fit_s) %>% mutate(mortality_group = sub("^mortality_group=", "", strata))
  km_tidy0_s <- km_tidy_s %>% distinct(mortality_group) %>% mutate(time = 0, estimate = 1, conf.low = 1, conf.high = 1)
  km_data_s <- bind_rows(km_tidy0_s, km_tidy_s) %>% arrange(mortality_group, time)
  pal_s <- setNames(c("#17766E", "#A8541D", "#7A1F2B")[seq_along(levels(analytic_data$mortality_group))],
                     levels(analytic_data$mortality_group))
  p <- ggplot(km_data_s, aes(x = time, y = 1 - estimate, color = mortality_group, fill = mortality_group)) +
    geom_ribbon(aes(ymin = 1 - conf.low, ymax = 1 - conf.high), alpha = 0.12, color = NA) +
    geom_step(linewidth = 0.9) +
    scale_color_manual(values = pal_s, name = NULL) + scale_fill_manual(values = pal_s, guide = "none") +
    scale_x_continuous(breaks = seq(0, 28, 7), limits = c(0, 28), expand = c(0, 0)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
    labs(title = paste0("28-day all-cause mortality: ", label),
         subtitle = sprintf("Log-rank p %s", if_else(km_p_s < 0.0001, "<0.0001", sprintf("= %.4f", km_p_s))),
         x = "Days since enrolment", y = "Cumulative mortality") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "top", plot.title = element_text(face = "bold"))
  ggsave(paste0(out_prefix, "_KM.png"), p, width = 8, height = 6, dpi = 500)

  sheets <- list(msm_causal_HR_and_risk_differences = msm, competing_risk_finegray = fg, model_diagnostics = diagnostics)
  if (!is.null(excess) && nrow(excess) > 0) sheets$national_excess_deaths <- excess
  write_xlsx(sheets, paste0(out_prefix, "_results.xlsx"))

  list(msm = msm, finegray = fg, diagnostics = diagnostics, national_excess = excess, km_logrank_p = km_p_s,
       mice_obj = mice_obj, analytic_data = analytic_data, label = label)
}

## --- 9.15.1: MDR (bacterial isolates), any organism ---
# MDR (Magiorakos-style category-count proxy, SECTION 2) replaces the plain
# resistant/susceptible split with a 3rd, more severe tier: Non-MDR BSI
# (susceptible + resistant-but-non-MDR) vs MDR BSI vs Uninfected control.
# n=1027 of mortality_analytic's 1046 (19 excluded: bacterial-domain/AST too
# incomplete to evaluate mdr_resistance, SECTION 9.1's comment).
mortality_analytic_mdr <- mortality_analytic %>%
  mutate(
    mortality_group = case_when(
      matched_control_broad == 1 ~ "Uninfected control",
      matched_control_broad == 0 & mdr_BSI_bacterial == 1 ~ "MDR BSI",
      matched_control_broad == 0 & mdr_BSI_bacterial == 0 ~ "Non-MDR BSI",
      TRUE ~ NA_character_
    ),
    mortality_group = factor(mortality_group, levels = c("Uninfected control", "Non-MDR BSI", "MDR BSI"))
  ) %>%
  filter(!is.na(mortality_group))

# Same large-font, title-free, 500dpi TIFF treatment for the MDR-mechanism KM
# panel (Figure 4 panel B) -- fast standalone KM computation, no MICE/bootstrap
# needed, so exit immediately after to skip the rest of the (expensive)
# mortality pipeline entirely.
km_fit_mdr <- survfit(Surv(time28, event28) ~ mortality_group, data = mortality_analytic_mdr)
km_logrank_mdr <- survdiff(Surv(time28, event28) ~ mortality_group, data = mortality_analytic_mdr)
km_p_mdr <- pchisq(km_logrank_mdr$chisq, length(km_logrank_mdr$n) - 1, lower.tail = FALSE)
km_tidy_mdr <- broom::tidy(km_fit_mdr) %>% mutate(mortality_group = sub("^mortality_group=", "", strata))
km_tidy0_mdr <- km_tidy_mdr %>% distinct(mortality_group) %>% mutate(time = 0, estimate = 1, conf.low = 1, conf.high = 1)
km_data_mdr <- bind_rows(km_tidy0_mdr, km_tidy_mdr) %>% arrange(mortality_group, time)
pal_mdr <- setNames(c("#17766E", "#A8541D", "#7A1F2B")[seq_along(levels(mortality_analytic_mdr$mortality_group))],
                     levels(mortality_analytic_mdr$mortality_group))

km_plot_mdr_big <- ggplot(km_data_mdr, aes(x = time, y = 1 - estimate, color = mortality_group, fill = mortality_group)) +
  geom_ribbon(aes(ymin = 1 - conf.low, ymax = 1 - conf.high), alpha = 0.12, color = NA) +
  geom_step(linewidth = 1.1) +
  scale_color_manual(values = pal_mdr, name = NULL) + scale_fill_manual(values = pal_mdr, guide = "none") +
  scale_x_continuous(breaks = seq(0, 28, 7), limits = c(0, 28), expand = c(0, 0)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  labs(x = "Days since enrolment", y = "Cumulative mortality") +
  theme_minimal(base_size = 24) +
  theme(
    panel.grid.minor = element_blank(), legend.position = "top", legend.text = element_text(size = 21),
    axis.title = element_text(size = 23), axis.text = element_text(size = 20)
  )
ggsave("Mortality_MDR_subgroup_KM_bigfont.tiff", km_plot_mdr_big, width = 11, height = 7, dpi = 500, compression = "lzw")
ggsave("Mortality_MDR_subgroup_KM_bigfont.png", km_plot_mdr_big, width = 11, height = 7, dpi = 500)
cat("BOTH BIG-FONT KM PLOTS SAVED\n")

# MDR-specific Cox MSM IPW adjusted curve, same multinomial-propensity + IPW
# construction as the primary curve (ipw_combined3, above), complete-case
# restricted to the adjustment covariates (illustrative single-fit curve,
# same convention as adjusted_overall_data/fg_primary_data).
mdr_cc <- mortality_analytic_mdr %>% filter(if_all(all_of(adjustment_covars), ~ !is.na(.)))
form_multi_mdr <- as.formula(paste("mortality_group ~", paste(adjustment_covars, collapse = " + ")))
mps_mdr <- nnet::multinom(form_multi_mdr, data = mdr_cc, trace = FALSE)
probs_mdr <- predict(mps_mdr, type = "probs")
marg_mdr <- prop.table(table(mdr_cc$mortality_group))
own_p_mdr <- probs_mdr[cbind(seq_len(nrow(probs_mdr)), match(as.character(mdr_cc$mortality_group), colnames(probs_mdr)))]
ipw_combined3_mdr <- as.numeric(marg_mdr[as.character(mdr_cc$mortality_group)]) / own_p_mdr
ipw_combined3_mdr <- pmin(pmax(ipw_combined3_mdr, quantile(ipw_combined3_mdr, 0.01)), quantile(ipw_combined3_mdr, 0.99))
mdr_cc$ipw_combined3 <- ipw_combined3_mdr
adjusted_mdr_data <- make_adjusted_curve_data(mdr_cc)

mdr_curve_colors <- c("Uninfected control" = "#17766E", "Non-MDR BSI" = "#A8541D", "MDR BSI" = "#7A1F2B")
ggsave("Mortality_adjusted_curves_CoxMSM_IPW_MDR_bigfont.png", build_bigfont_adjusted_plot(adjusted_mdr_data, mdr_curve_colors), width = 11, height = 7, dpi = 500)
ggsave("Mortality_adjusted_curves_CoxMSM_IPW_MDR_bigfont.tiff", build_bigfont_adjusted_plot(adjusted_mdr_data, mdr_curve_colors), width = 11, height = 7, dpi = 500, compression = "lzw")

# MDR-specific Fine-Gray IPW-adjusted cumulative-incidence curve (reuses the
# generic fg_ipw_curve_data() helper defined at the primary-comparison point).
cr_data_mdr <- mortality_analytic_mdr %>% filter(!is.na(event_inhosp), !is.na(time_inhosp))
fg_mdr_data <- fg_ipw_curve_data(cr_data_mdr, adjustment_covars)
write_xlsx(fg_mdr_data, "Mortality_adjusted_curves_FineGray_MDR_data.xlsx")
ggsave("Mortality_adjusted_curves_FineGray_MDR_bigfont.png", build_bigfont_adjusted_plot(fg_mdr_data, mdr_curve_colors, fg_ylab), width = 11, height = 7, dpi = 500)
ggsave("Mortality_adjusted_curves_FineGray_MDR_bigfont.tiff", build_bigfont_adjusted_plot(fg_mdr_data, mdr_curve_colors, fg_ylab), width = 11, height = 7, dpi = 500, compression = "lzw")

# --- Combine into 2-panel Figure-4-style images: one for Cox MSM IPW, one for Fine-Gray ---
suppressMessages(library(patchwork))
tag_theme <- theme(
  plot.tag = element_text(size = 18, face = "bold", margin = margin(b = -8, r = 2)),
  plot.tag.position = c(0.01, 0.99),
  plot.margin = margin(0, 4, 0, 2)
)

coxmsm_combined <- (build_bigfont_adjusted_plot(adjusted_overall_data, mort_colors) |
                     build_bigfont_adjusted_plot(adjusted_mdr_data, mdr_curve_colors)) +
  plot_layout(widths = c(1, 1)) + plot_annotation(tag_levels = "A") & tag_theme
ggsave("SuppFigure_mortality_adjusted_CoxMSM_IPW.png", coxmsm_combined, width = 20, height = 7.2, dpi = 500)
ggsave("SuppFigure_mortality_adjusted_CoxMSM_IPW.tiff", coxmsm_combined, width = 20, height = 7.2, dpi = 500, compression = "lzw")

finegray_combined <- (build_bigfont_adjusted_plot(fg_primary_data, mort_colors, fg_ylab) |
                       build_bigfont_adjusted_plot(fg_mdr_data, mdr_curve_colors, fg_ylab)) +
  plot_layout(widths = c(1, 1)) + plot_annotation(tag_levels = "A") & tag_theme
ggsave("SuppFigure_mortality_adjusted_FineGray_IPW.png", finegray_combined, width = 20, height = 7.2, dpi = 500)
ggsave("SuppFigure_mortality_adjusted_FineGray_IPW.tiff", finegray_combined, width = 20, height = 7.2, dpi = 500, compression = "lzw")

cat("ADJUSTED COX-MSM-IPW AND FINE-GRAY SUPPLEMENTARY FIGURES SAVED\n")

mdr_comparisons <- list(
  c(ref = "Uninfected control", comp = "Non-MDR BSI"),
  c(ref = "Uninfected control", comp = "MDR BSI"),
  c(ref = "Non-MDR BSI", comp = "MDR BSI")
)

rate_mdr_national <- incidence_primary_table %>%
  filter(measure == "MDR BSI (bacterial isolates)", exposure_type == "Total_inpatient_admissions") %>%
  pull(calibrated_rate_per1000)
national_mdr_cases <- rate_mdr_national / 1000 * national_admissions_total
national_nonmdr_cases <- national_bsi_cases - national_mdr_cases
national_cases_mdr <- list(`Non-MDR BSI` = national_nonmdr_cases, `MDR BSI` = national_mdr_cases)

mdr_subgroup_results <- run_mortality_subgroup(
  "MDR vs non-MDR BSI vs uninfected controls", mortality_analytic_mdr, mdr_comparisons,
  "Mortality_MDR_subgroup", national_cases = national_cases_mdr
)

## --- 9.15.1b: XDR (bacterial isolates), any organism -- mirrors 9.15.1's MDR
## construction exactly, substituting xdr_BSI_bacterial for mdr_BSI_bacterial. ---
mortality_analytic_xdr <- mortality_analytic %>%
  mutate(
    mortality_group = case_when(
      matched_control_broad == 1 ~ "Uninfected control",
      matched_control_broad == 0 & xdr_BSI_bacterial == 1 ~ "XDR BSI",
      matched_control_broad == 0 & xdr_BSI_bacterial == 0 ~ "Non-XDR BSI",
      TRUE ~ NA_character_
    ),
    mortality_group = factor(mortality_group, levels = c("Uninfected control", "Non-XDR BSI", "XDR BSI"))
  ) %>%
  filter(!is.na(mortality_group))

xdr_comparisons <- list(
  c(ref = "Uninfected control", comp = "Non-XDR BSI"),
  c(ref = "Uninfected control", comp = "XDR BSI"),
  c(ref = "Non-XDR BSI", comp = "XDR BSI")
)

rate_xdr_national <- incidence_primary_table %>%
  filter(measure == "XDR BSI (bacterial isolates)", exposure_type == "Total_inpatient_admissions") %>%
  pull(calibrated_rate_per1000)
national_xdr_cases <- rate_xdr_national / 1000 * national_admissions_total
national_nonxdr_cases <- national_bsi_cases - national_xdr_cases
national_cases_xdr <- list(`Non-XDR BSI` = national_nonxdr_cases, `XDR BSI` = national_xdr_cases)

xdr_subgroup_results <- run_mortality_subgroup(
  "XDR vs non-XDR BSI vs uninfected controls", mortality_analytic_xdr, xdr_comparisons,
  "Mortality_XDR_subgroup", national_cases = national_cases_xdr
)

## --- 9.15.2: Carbapenem resistance, Gram-negative BSI only ---
# Restricted to Enterobacterales, Pseudomonas aeruginosa, and Acinetobacter
# baumannii specifically (is_cre_pa_ab(), SECTION 3.6) -- the WHO GLASS
# priority-pathogen group carbapenem resistance is surveilled against, NOT
# the broader major_bacterial_group "Non-fermenting Gram-negative bacilli"/
# "Other Gram-negative" buckets, which also sweep in species such as
# Stenotrophomonas maltophilia with well-documented INTRINSIC (not acquired)
# carbapenem resistance -- pooling those in would inflate this rate with
# resistance that reflects species identity, not a true acquired-resistance
# signal. Previously used the broader bucket (major_bacterial_group %in%
# gneg_groups); corrected throughout this section (mortality, risk factors,
# LOS) to match the narrower, WHO-GLASS-consistent restriction SECTION 3.6's
# own prevalence table already used.
mortality_analytic_carb <- mortality_analytic %>%
  mutate(
    mortality_group = case_when(
      matched_control_broad == 1 ~ "Uninfected control",
      matched_control_broad == 0 & is_cre_pa_ab(major_bacterial_group, pathogen1_clean) & carb_BSI_positive == 1 ~ "Carbapenem-R Gram-negative BSI",
      matched_control_broad == 0 & is_cre_pa_ab(major_bacterial_group, pathogen1_clean) & carb_BSI_positive == 0 ~ "Carbapenem-S Gram-negative BSI",
      TRUE ~ NA_character_
    ),
    mortality_group = factor(mortality_group, levels = c("Uninfected control", "Carbapenem-S Gram-negative BSI", "Carbapenem-R Gram-negative BSI"))
  ) %>%
  filter(!is.na(mortality_group))

carb_comparisons <- list(
  c(ref = "Uninfected control", comp = "Carbapenem-S Gram-negative BSI"),
  c(ref = "Uninfected control", comp = "Carbapenem-R Gram-negative BSI"),
  c(ref = "Carbapenem-S Gram-negative BSI", comp = "Carbapenem-R Gram-negative BSI")
)

# National case count: uses SECTION 3.6's own organism-restricted rate-per-
# 100-admissions row ("Carbapenem resistance (Enterobacterales/P.
# aeruginosa/A. baumannii isolates)") rather than the population-wide
# "Carbapenem-resistant BSI" incidence measure -- same fix as MRSA (SECTION
# 9.15.3) and cephalosporin any-generation (SECTION 9.15.5). No comparably
# defensible national rate exists for "carbapenem-susceptible Gram-negative
# BSI" specifically, so only comparisons involving the resistant arm get an
# excess-deaths estimate; the susceptible-vs-uninfected comparison is
# skipped rather than guessed at.
rate_carbR_national_per100 <- amr_rate_per100_table %>%
  filter(measure == "Carbapenem resistance (Enterobacterales/P. aeruginosa/A. baumannii isolates)") %>%
  pull(calibrated_rate_per100)
national_carbR_cases <- rate_carbR_national_per100 / 100 * national_admissions_total
national_cases_carb <- list(`Carbapenem-R Gram-negative BSI` = national_carbR_cases)

carb_subgroup_results <- run_mortality_subgroup(
  "carbapenem resistance, Gram-negative BSI", mortality_analytic_carb, carb_comparisons,
  "Mortality_carbapenemR_gramneg_subgroup", national_cases = national_cases_carb
)

## --- 9.15.3: Methicillin resistance (MRSA), Staphylococcus AUREUS BSI only ---
# Restricted to pathogen1_clean == "Staphylococcus aureus" specifically --
# true MRSA, not the wider "Staphylococcus spp." major bacterial group,
# which also includes coagulase-negative Staphylococcus species (MRSE/MRSH-
# type resistance is a different clinical entity). Only MRSA (n=112) vs
# Uninfected control (n=646) is modelled -- MSSA (n=16) has ZERO 28-day
# deaths among its 14 patients with a determinate outcome, which cannot
# support a Cox HR (quasi-complete separation) in either an MSSA-vs-
# Uninfected or MRSA-vs-MSSA contrast; MSSA is reported descriptively only.
mortality_analytic_mrsa <- mortality_analytic %>%
  mutate(
    mortality_group = case_when(
      matched_control_broad == 1 ~ "Uninfected control",
      matched_control_broad == 0 & pathogen1_clean == "Staphylococcus aureus" & methicillin_BSI_positive == 1 ~ "MRSA BSI",
      TRUE ~ NA_character_
    ),
    mortality_group = factor(mortality_group, levels = c("Uninfected control", "MRSA BSI"))
  ) %>%
  filter(!is.na(mortality_group))

mrsa_comparisons <- list(c(ref = "Uninfected control", comp = "MRSA BSI"))

# National case count: uses Table 2's own S.-aureus-restricted rate-per-100-
# admissions (amr_rate_per100_table's "Methicillin resistance (S. aureus
# isolates)" row, SECTION 3.6.4) rather than incidence_by_amr_class's
# population-wide "Methicillin-resistant BSI" measure (which is not
# organism-restricted at all and would silently fold in coagulase-negative
# Staphylococcus methicillin resistance -- a materially different, and
# typically far more common, resistance pattern than true MRSA). Unlike the
# population-wide-as-proxy approximation already accepted for carbapenem/
# cephalosporin resistance (SECTION 9.15.2's own comment: resistance there
# is "almost exclusively" a Gram-negative finding in practice), methicillin
# resistance is tested and reported across ALL Staphylococcus species, so
# that same approximation would not be defensible here.
rate_mrsa_national_per100 <- amr_rate_per100_table %>%
  filter(measure == "Methicillin resistance (S. aureus isolates)") %>%
  pull(calibrated_rate_per100)
national_mrsa_cases <- rate_mrsa_national_per100 / 100 * national_admissions_total
national_cases_mrsa <- list(`MRSA BSI` = national_mrsa_cases)

mrsa_subgroup_results <- run_mortality_subgroup(
  "MRSA vs uninfected controls (Staphylococcus BSI)", mortality_analytic_mrsa, mrsa_comparisons,
  "Mortality_MRSA_subgroup", national_cases = national_cases_mrsa
)

mssa_descriptive <- mortality_analytic %>%
  filter(matched_control == 0, pathogen1_clean == "Staphylococcus aureus", methicillin_BSI_positive == 0) %>%
  summarise(
    group = "MSSA BSI (methicillin-susceptible)", n = n(),
    n_outcome_known = sum(!is.na(event28)), deaths_28day = sum(event28 == 1, na.rm = TRUE)
  ) %>%
  mutate(mortality_pct_28day = round(100 * deaths_28day / n_outcome_known, 1),
         note = sprintf("Descriptive only -- %d deaths among %d evaluable patients, not modelled (SECTION 9.15.3)", deaths_28day, n_outcome_known))
write_xlsx(mssa_descriptive, "Mortality_MSSA_descriptive_only.xlsx")

## --- 9.15.4: Vancomycin resistance, Enterococcus BSI -- descriptive only ---
# genus1 == "Enterococcus": n=24 total, only 3 VRE (0 deaths among them) and
# 17 VSE (3 deaths) -- an order of magnitude below the smallest arm used
# anywhere else in SECTION 9 (Susceptible BSI, n=28) and with 0 events in
# the exposed arm. No propensity model, Cox model, or bootstrap is fit;
# raw counts only, per instruction to model a mechanism only when numbers
# support it.
vre_descriptive <- mortality_analytic %>%
  filter(matched_control == 0, !is.na(genus1), genus1 == "Enterococcus", !is.na(vancomycin_BSI_positive)) %>%
  mutate(group = if_else(vancomycin_BSI_positive == 1, "VRE (vancomycin-resistant)", "VSE (vancomycin-susceptible)")) %>%
  group_by(group) %>%
  summarise(n = n(), n_outcome_known = sum(!is.na(event28)), deaths_28day = sum(event28 == 1, na.rm = TRUE), .groups = "drop") %>%
  mutate(mortality_pct_28day = round(100 * deaths_28day / n_outcome_known, 1),
         note = "Descriptive only -- insufficient events for IPW/Cox modelling (SECTION 9.15.4)")
write_xlsx(vre_descriptive, "Mortality_VRE_descriptive_only.xlsx")

## --- 9.15.5: Cephalosporin resistance (any generation), Gram-negative BSI only ---
# Same organism restriction as SECTION 9.15.2's corrected carbapenem
# definition: Enterobacterales, Pseudomonas aeruginosa, and Acinetobacter
# baumannii specifically (is_cre_pa_ab()), not the broader major_bacterial_
# group bucket that also sweeps in intrinsically-resistant non-fermenters
# such as Stenotrophomonas maltophilia. Previously used major_bacterial_
# group %in% gneg_groups; corrected for the same reason as 9.15.2.
mortality_analytic_ceph <- mortality_analytic %>%
  mutate(
    mortality_group = case_when(
      matched_control_broad == 1 ~ "Uninfected control",
      matched_control_broad == 0 & is_cre_pa_ab(major_bacterial_group, pathogen1_clean) & cephalosporin_BSI_positive == 1 ~ "Cephalosporin-R Gram-negative BSI",
      matched_control_broad == 0 & is_cre_pa_ab(major_bacterial_group, pathogen1_clean) & cephalosporin_BSI_positive == 0 ~ "Cephalosporin-S Gram-negative BSI",
      TRUE ~ NA_character_
    ),
    mortality_group = factor(mortality_group, levels = c("Uninfected control", "Cephalosporin-S Gram-negative BSI", "Cephalosporin-R Gram-negative BSI"))
  ) %>%
  filter(!is.na(mortality_group))

ceph_comparisons <- list(
  c(ref = "Uninfected control", comp = "Cephalosporin-S Gram-negative BSI"),
  c(ref = "Uninfected control", comp = "Cephalosporin-R Gram-negative BSI"),
  c(ref = "Cephalosporin-S Gram-negative BSI", comp = "Cephalosporin-R Gram-negative BSI")
)

# National case count: uses SECTION 3.6's own organism-restricted rate-per-
# 100-admissions row ("Cephalosporin resistance, any generation
# (Enterobacterales/P. aeruginosa/A. baumannii isolates)") rather than the
# population-wide "Cephalosporin-resistant BSI (any generation)" incidence
# measure -- the same fix applied to MRSA's national rate (SECTION 9.15.3)
# and now possible here because SECTION 3.6 computes this organism-
# restricted rate directly, rather than relying on the "almost exclusively
# a Gram-negative finding in practice, not formally enforced upstream"
# approximation this section previously depended on.
rate_cephR_national_per100 <- amr_rate_per100_table %>%
  filter(measure == "Cephalosporin resistance, any generation (Enterobacterales/P. aeruginosa/A. baumannii isolates)") %>%
  pull(calibrated_rate_per100)
national_cephR_cases <- rate_cephR_national_per100 / 100 * national_admissions_total
national_cases_ceph <- list(`Cephalosporin-R Gram-negative BSI` = national_cephR_cases)

ceph_subgroup_results <- run_mortality_subgroup(
  "cephalosporin resistance (any generation), Gram-negative BSI", mortality_analytic_ceph, ceph_comparisons,
  "Mortality_cephalosporinR_gramneg_subgroup", national_cases = national_cases_ceph
)

## --- 9.15.6: 3rd-generation cephalosporin resistance (3GCR), Gram-negative BSI only ---
# Same organism restriction (is_cre_pa_ab(): Enterobacterales, P. aeruginosa,
# A. baumannii) as SECTIONS 9.15.2/9.15.5, and now consistent with SECTION
# 3.6's own AMR prevalence table, which always used this narrower cre_pa_ab_
# group denominator for 3GCR. Previously this section used the broader
# major_bacterial_group %in% gneg_groups bucket "for internal consistency
# across this section" rather than matching SECTION 3.6 -- corrected so all
# three Gram-negative mechanisms (carbapenem, cephalosporin any-generation,
# 3GCR) now share one restriction, defined once, everywhere.
mortality_analytic_3gcr <- mortality_analytic %>%
  mutate(
    mortality_group = case_when(
      matched_control_broad == 1 ~ "Uninfected control",
      matched_control_broad == 0 & is_cre_pa_ab(major_bacterial_group, pathogen1_clean) & c3GC_BSI_positive == 1 ~ "3GCR Gram-negative BSI",
      matched_control_broad == 0 & is_cre_pa_ab(major_bacterial_group, pathogen1_clean) & c3GC_BSI_positive == 0 ~ "3GC-susceptible Gram-negative BSI",
      TRUE ~ NA_character_
    ),
    mortality_group = factor(mortality_group, levels = c("Uninfected control", "3GC-susceptible Gram-negative BSI", "3GCR Gram-negative BSI"))
  ) %>%
  filter(!is.na(mortality_group))

gcr3_comparisons <- list(
  c(ref = "Uninfected control", comp = "3GC-susceptible Gram-negative BSI"),
  c(ref = "Uninfected control", comp = "3GCR Gram-negative BSI"),
  c(ref = "3GC-susceptible Gram-negative BSI", comp = "3GCR Gram-negative BSI")
)

# National case count: uses SECTION 3.6's own organism-restricted rate-per-
# 100-admissions row directly (same fix as 9.15.2/9.15.3/9.15.5), rather than
# the population-wide "3rd-generation cephalosporin-resistant BSI" incidence
# measure.
rate_3gcr_national_per100 <- amr_rate_per100_table %>%
  filter(measure == "3rd-gen cephalosporin resistance (Enterobacterales/P. aeruginosa/A. baumannii isolates)") %>%
  pull(calibrated_rate_per100)
national_3gcr_cases <- rate_3gcr_national_per100 / 100 * national_admissions_total
national_cases_3gcr <- list(`3GCR Gram-negative BSI` = national_3gcr_cases)

gcr3_subgroup_results <- run_mortality_subgroup(
  "3rd-generation cephalosporin resistance (3GCR), Gram-negative BSI", mortality_analytic_3gcr, gcr3_comparisons,
  "Mortality_3GCR_gramneg_subgroup", national_cases = national_cases_3gcr
)

## --- 9.15.7/8/9: Access/Watch/Reserve (AWaRe), any organism -- mirrors
## 9.15.1's MDR construction exactly (no organism restriction: AWaRe
## categories are antibiotic-class groupings applied BSI-wide, the same
## population any_resistance/MDR/XDR use), substituting access_BSI_positive/
## watch_BSI_positive/reserve_BSI_positive for mdr_BSI_bacterial. Added so
## every mechanism modelled in Risk Factors (SECTION 9.20-adjacent) also has
## the full IPW/Cox + exact-1:1-matched mortality treatment, not just a
## risk-factor logistic model.
build_aware_mortality_subgroup <- function(aware_col, resistant_label, susceptible_label, aware_measure_name, out_prefix, run_label) {
  analytic <- mortality_analytic %>%
    mutate(
      mortality_group = case_when(
        matched_control_broad == 1 ~ "Uninfected control",
        matched_control_broad == 0 & .data[[aware_col]] == 1 ~ resistant_label,
        matched_control_broad == 0 & .data[[aware_col]] == 0 ~ susceptible_label,
        TRUE ~ NA_character_
      ),
      mortality_group = factor(mortality_group, levels = c("Uninfected control", susceptible_label, resistant_label))
    ) %>%
    filter(!is.na(mortality_group))

  comparisons <- list(
    c(ref = "Uninfected control", comp = susceptible_label),
    c(ref = "Uninfected control", comp = resistant_label),
    c(ref = susceptible_label, comp = resistant_label)
  )

  rate_national_per100 <- amr_rate_per100_table %>%
    filter(measure == aware_measure_name) %>%
    pull(calibrated_rate_per100)
  national_resistant_cases <- rate_national_per100 / 100 * national_admissions_total
  national_cases <- setNames(list(national_resistant_cases), resistant_label)

  list(
    analytic = analytic,
    comparisons = comparisons,
    results = run_mortality_subgroup(run_label, analytic, comparisons, out_prefix, national_cases = national_cases)
  )
}

access_mortality <- build_aware_mortality_subgroup(
  "access_BSI_positive", "Access-R BSI", "Access-S BSI",
  "Access-antibiotic resistant BSI (AWaRe)", "Mortality_AccessR_subgroup",
  "Access-antibiotic resistant vs susceptible BSI vs uninfected controls"
)
mortality_analytic_access <- access_mortality$analytic
access_comparisons <- access_mortality$comparisons
access_subgroup_results <- access_mortality$results

watch_mortality <- build_aware_mortality_subgroup(
  "watch_BSI_positive", "Watch-R BSI", "Watch-S BSI",
  "Watch-antibiotic resistant BSI (AWaRe)", "Mortality_WatchR_subgroup",
  "Watch-antibiotic resistant vs susceptible BSI vs uninfected controls"
)
mortality_analytic_watch <- watch_mortality$analytic
watch_comparisons <- watch_mortality$comparisons
watch_subgroup_results <- watch_mortality$results

reserve_mortality <- build_aware_mortality_subgroup(
  "reserve_BSI_positive", "Reserve-R BSI", "Reserve-S BSI",
  "Reserve-antibiotic resistant BSI (AWaRe)", "Mortality_ReserveR_subgroup",
  "Reserve-antibiotic resistant vs susceptible BSI vs uninfected controls"
)
mortality_analytic_reserve <- reserve_mortality$analytic
reserve_comparisons <- reserve_mortality$comparisons
reserve_subgroup_results <- reserve_mortality$results

# =========================================================================
# SUPPLEMENTARY: combined IPW x calibrated-design-weight sensitivity Cox
# model (SECTION 9.13's construction, previously fit only for the primary
# 3-comparison "any resistance" set), extended to every resistance
# mechanism. Reuses each mechanism's own MICE object (mice_obj, already fit
# above inside its run_mortality_subgroup() call) and its own analytic_
# data's facility_id; calibrated_weight is joined per mechanism the same
# way SECTION 9.13 does for the primary comparison (facility-level property,
# valid for cases and controls alike).
cat("=== COMBINED IPW x CALIBRATED WEIGHT SENSITIVITY MODEL, ALL MECHANISMS ===\n")
mechanism_combined_weight_specs <- list(
  list(label = "MDR (bacterial isolates)", results = mdr_subgroup_results, comparisons = mdr_comparisons),
  list(label = "XDR (bacterial isolates)", results = xdr_subgroup_results, comparisons = xdr_comparisons),
  list(label = "Carbapenem-R (Gram-negative)", results = carb_subgroup_results, comparisons = carb_comparisons),
  list(label = "Methicillin-R/MRSA (S. aureus)", results = mrsa_subgroup_results, comparisons = mrsa_comparisons),
  list(label = "Cephalosporin-R any-gen (Gram-negative)", results = ceph_subgroup_results, comparisons = ceph_comparisons),
  list(label = "3GCR (Gram-negative)", results = gcr3_subgroup_results, comparisons = gcr3_comparisons),
  list(label = "Access-resistant (AWaRe)", results = access_subgroup_results, comparisons = access_comparisons),
  list(label = "Watch-resistant (AWaRe)", results = watch_subgroup_results, comparisons = watch_comparisons),
  list(label = "Reserve-resistant (AWaRe)", results = reserve_subgroup_results, comparisons = reserve_comparisons)
)

combined_weight_all_mechanisms <- purrr::map_dfr(mechanism_combined_weight_specs, function(spec) {
  ad <- spec$results$analytic_data
  fac_ids <- ad$facility_id
  cw_by_fac <- ad %>%
    select(-any_of("calibrated_weight")) %>%
    left_join(facility_weights %>% select(UniqueID, calibrated_weight), by = c(facility_id = "UniqueID")) %>%
    pull(calibrated_weight)
  n_missing <- sum(is.na(cw_by_fac))
  purrr::map_dfr(spec$comparisons, function(cmp) {
    tryCatch(
      fit_combined_weight_pooled(spec$results$mice_obj, cmp["ref"], cmp["comp"], adjustment_covars, fac_ids, cw_by_fac),
      error = function(e) tibble(comparison = paste(cmp["comp"], "vs", cmp["ref"]), log_hr = NA_real_, se_log_hr = NA_real_,
                                  hazard_ratio = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_, m_used = 0L)
    ) %>% mutate(mechanism = spec$label, n_missing_calibrated_weight = n_missing, .before = 1)
  })
}) %>%
  mutate(across(where(is.numeric) & !c(m_used, n_missing_calibrated_weight), ~ round(.x, 4)))

print(as.data.frame(combined_weight_all_mechanisms), row.names = FALSE)
write_xlsx(combined_weight_all_mechanisms, "Mortality_MSM_combined_IPW_survey_weight_all_mechanisms.xlsx")

# -----------------------------------------------------------------------
# 9.16: Exact 1:1 matched-pair analysis (counterfactual within-pair design)
# -----------------------------------------------------------------------
# SECTIONS 9.1-9.15 treat "Uninfected control" as one pooled comparator arm
# (is_true_matched_control() already restricts it to genuine matches, but
# still analysed unpaired, via IPW). This section instead exploits the
# survey's own 1:1 matching directly: every matched-control patient_id is
# the SAME base ID as the case they were matched to, with "/C3" appended
# (confirmed directly against the raw sheets -- e.g. case "HF017/W02/005"
# pairs with control "HF017/W02/005/C3"). Stripping "/C3" recovers that
# shared base ID, letting each case be compared against ITS OWN matched
# control specifically, rather than against the pooled control arm -- a
# genuine within-pair counterfactual, additional to (not a replacement for)
# SECTIONS 9.1-9.15's primary IPW analysis.
mortality_analytic$base_id <- str_remove(mortality_analytic$patient_id, "/C3$")
# matched_control == 0 alone is no longer sufficient to mean "is a case":
# mortality_group's "Uninfected control" category now also covers broadly-
# matched (matched_control_broad==1) Cohort-3 patients whose OWN
# matched_control (strict) is 0 -- their nominal paired case was never
# confirmed BSI-positive, so is_true_matched_control() correctly returns 0
# for them, but they are still controls, not cases. Excluding
# mortality_group == "Uninfected control" explicitly (in addition to
# matched_control == 0) keeps pair_cases restricted to genuine BSI cases
# only, as intended -- this 1:1 pairing is deliberately built on the
# STRICT matched_control flag throughout (both sides), unlike SECTIONS
# 9.1-9.15's IPW models.
pair_cases <- mortality_analytic %>% dplyr::filter(matched_control == 0, mortality_group != "Uninfected control")
pair_controls <- mortality_analytic %>% dplyr::filter(matched_control == 1) %>%
  dplyr::select(base_id, control_time28 = time28, control_event28 = event28,
                control_time_inhosp = time_inhosp, control_facility_id = facility_id)

mortality_paired <- pair_cases %>%
  inner_join(pair_controls, by = "base_id") %>%
  dplyr::filter(!is.na(time28), !is.na(event28), !is.na(control_time28), !is.na(control_event28))
cat(sprintf("SECTION 9.16: %d exact 1:1 case-control pairs with complete 28-day outcome data on both sides\n", nrow(mortality_paired)))

# Long format for stratified Cox: each pair (base_id) is its own stratum, so
# the model compares a case only against its own matched control, never
# across pairs -- this is what "strata()" does in coxph, and is the
# standard analysis for a genuinely 1:1 matched cohort. No IPW is needed
# here: matching by design already balances whatever the ward-code/facility
# matching criteria controlled for; this is deliberately a simpler,
# assumption-light complement to SECTION 9.7's covariate-adjusted IPW model,
# not a replacement for it.
mortality_paired_long <- bind_rows(
  mortality_paired %>% transmute(base_id, mortality_group, time28, event28),
  mortality_paired %>% transmute(base_id, mortality_group = "Uninfected control",
                                  time28 = control_time28, event28 = control_event28)
)

# A 2-person-stratified Cox model (one case + its own matched control per
# stratum) is mathematically equivalent to a conditional analysis on
# discordant pairs only -- concordant pairs (both died or both survived)
# contribute no information either way. For a genuinely 1:1 matched BINARY
# outcome this is the classic McNemar's-test setting, and it is the
# textbook-correct tool here, not a stratified Cox: when nearly every
# discordant pair points the same direction (as happens with a rare event
# in a small pair-set), the Cox partial likelihood has no finite maximum
# ("monotone likelihood"/complete separation) and returns a nonsensical,
# numerically-unbounded hazard ratio (this was tried first and produced
# HR > 1 billion for the Susceptible-BSI comparison -- a real small-sample
# statistical limitation, not a coding bug: among 27 Susceptible-BSI pairs
# there were 3 discordant pairs, case-died/control-alive in all 3 and
# control-died/case-alive in 0, i.e. zero pairs favouring survival of the
# case -- see check_susc_pairs.R diagnostic). We report the exact McNemar
# discordant-pair counts and exact binomial p-value (valid regardless of
# how sparse or one-sided the discordant pairs are) plus a
# Haldane-Anscombe-corrected (add 0.5 to each discordant cell) odds ratio,
# which stays finite and appropriately wide even under this near-complete
# separation, instead of reporting the uninterpretable Cox HR.
## Generalised to an arbitrary case-membership vector (not just the primary
## mortality_group), so the SAME 1:1-pair McNemar logic can be applied
## mechanism-by-mechanism (MDR, XDR, organism-restricted carbapenem/
## methicillin/cephalosporin/3GC) -- every mechanism's "case" side is still
## drawn from pair_cases (matched_control==0, strict), so each case is still
## compared only to its own matched control. NA rows (e.g. carbapenem
## resistance never tested, or organism restriction not met) are dropped,
## not coerced to FALSE, mirroring the same convention used for every other
## organism-restricted flag in this script.
mcnemar_paired_generic <- function(case_filter_vec, label) {
  di <- mortality_paired[case_filter_vec & !is.na(case_filter_vec), ]
  n10 <- sum(di$event28 == 1 & di$control_event28 == 0)  # case died, control alive
  n01 <- sum(di$event28 == 0 & di$control_event28 == 1)  # control died, case alive
  n_discordant <- n10 + n01
  exact_p <- if (n_discordant > 0) binom.test(n10, n_discordant, p = 0.5)$p.value else NA_real_
  or_hc <- (n10 + 0.5) / (n01 + 0.5)
  se_log_or <- sqrt(1 / (n10 + 0.5) + 1 / (n01 + 0.5))
  log_or <- log(or_hc)
  tibble(
    comparison = label,
    n_pairs = nrow(di),
    case_died_control_alive = n10,
    control_died_case_alive = n01,
    concordant_pairs = nrow(di) - n_discordant,
    odds_ratio_haldane_corrected = round(or_hc, 3),
    ci_lower = round(exp(log_or - 1.96 * se_log_or), 3),
    ci_upper = round(exp(log_or + 1.96 * se_log_or), 3),
    mcnemar_exact_p = round(exact_p, 4)
  )
}
mcnemar_paired_mortality <- function(comp_level) {
  mcnemar_paired_generic(mortality_paired$mortality_group == comp_level,
                          paste(comp_level, "vs Uninfected control (exact 1:1 pair, McNemar)"))
}
mortality_paired_results <- bind_rows(
  mcnemar_paired_mortality("Susceptible BSI"),
  mcnemar_paired_mortality("Resistant BSI"),
  mcnemar_paired_generic(mortality_paired$mdr_BSI_bacterial == 1, "MDR BSI vs Uninfected control (exact 1:1 pair, McNemar)"),
  mcnemar_paired_generic(mortality_paired$mdr_BSI_bacterial == 0, "Non-MDR BSI vs Uninfected control (exact 1:1 pair, McNemar)"),
  mcnemar_paired_generic(mortality_paired$xdr_BSI_bacterial == 1, "XDR BSI vs Uninfected control (exact 1:1 pair, McNemar)"),
  mcnemar_paired_generic(mortality_paired$xdr_BSI_bacterial == 0, "Non-XDR BSI vs Uninfected control (exact 1:1 pair, McNemar)"),
  mcnemar_paired_generic(is_cre_pa_ab(mortality_paired$major_bacterial_group, mortality_paired$pathogen1_clean) & mortality_paired$carb_BSI_positive == 1, "Carbapenem-R Gram-negative BSI vs Uninfected control (exact 1:1 pair, McNemar)"),
  mcnemar_paired_generic(is_cre_pa_ab(mortality_paired$major_bacterial_group, mortality_paired$pathogen1_clean) & mortality_paired$carb_BSI_positive == 0, "Carbapenem-S Gram-negative BSI vs Uninfected control (exact 1:1 pair, McNemar)"),
  mcnemar_paired_generic(mortality_paired$pathogen1_clean == "Staphylococcus aureus" & mortality_paired$methicillin_BSI_positive == 1, "MRSA BSI vs Uninfected control (exact 1:1 pair, McNemar)"),
  mcnemar_paired_generic(mortality_paired$pathogen1_clean == "Staphylococcus aureus" & mortality_paired$methicillin_BSI_positive == 0, "MSSA BSI vs Uninfected control (exact 1:1 pair, McNemar)"),
  mcnemar_paired_generic(is_cre_pa_ab(mortality_paired$major_bacterial_group, mortality_paired$pathogen1_clean) & mortality_paired$cephalosporin_BSI_positive == 1, "Cephalosporin-R Gram-negative BSI vs Uninfected control (exact 1:1 pair, McNemar)"),
  mcnemar_paired_generic(is_cre_pa_ab(mortality_paired$major_bacterial_group, mortality_paired$pathogen1_clean) & mortality_paired$cephalosporin_BSI_positive == 0, "Cephalosporin-S Gram-negative BSI vs Uninfected control (exact 1:1 pair, McNemar)"),
  mcnemar_paired_generic(is_cre_pa_ab(mortality_paired$major_bacterial_group, mortality_paired$pathogen1_clean) & mortality_paired$c3GC_BSI_positive == 1, "3GCR Gram-negative BSI vs Uninfected control (exact 1:1 pair, McNemar)"),
  mcnemar_paired_generic(is_cre_pa_ab(mortality_paired$major_bacterial_group, mortality_paired$pathogen1_clean) & mortality_paired$c3GC_BSI_positive == 0, "3GC-susceptible Gram-negative BSI vs Uninfected control (exact 1:1 pair, McNemar)"),
  mcnemar_paired_generic(mortality_paired$access_BSI_positive == 1, "Access-R BSI vs Uninfected control (exact 1:1 pair, McNemar)"),
  mcnemar_paired_generic(mortality_paired$access_BSI_positive == 0, "Access-S BSI vs Uninfected control (exact 1:1 pair, McNemar)"),
  mcnemar_paired_generic(mortality_paired$watch_BSI_positive == 1, "Watch-R BSI vs Uninfected control (exact 1:1 pair, McNemar)"),
  mcnemar_paired_generic(mortality_paired$watch_BSI_positive == 0, "Watch-S BSI vs Uninfected control (exact 1:1 pair, McNemar)"),
  mcnemar_paired_generic(mortality_paired$reserve_BSI_positive == 1, "Reserve-R BSI vs Uninfected control (exact 1:1 pair, McNemar)"),
  mcnemar_paired_generic(mortality_paired$reserve_BSI_positive == 0, "Reserve-S BSI vs Uninfected control (exact 1:1 pair, McNemar)")
)

# Paired length-of-stay comparison: per-pair difference (case minus its own
# matched control) in time_inhosp (SECTION 9.2/11's enrolment-to-discharge/
# death LOS, truncated at 28 days), summarised as a mean paired difference
# with a pair-cluster bootstrap CI (resampling PAIRS, not individual
# patients, preserving the matched structure per the pairing itself).
set.seed(47)
paired_los_diff <- function(comp_level, n_boot = 200) {
  di <- mortality_paired %>% dplyr::filter(mortality_group == comp_level, !is.na(time_inhosp), !is.na(control_time_inhosp)) %>%
    mutate(pair_diff = time_inhosp - control_time_inhosp)
  point <- mean(di$pair_diff)
  boot_means <- replicate(n_boot, mean(sample(di$pair_diff, nrow(di), replace = TRUE)))
  tibble(comparison = paste(comp_level, "vs Uninfected control (exact 1:1 pair)"), n_pairs = nrow(di),
         mean_paired_diff_days = point, ci_lower = quantile(boot_means, 0.025), ci_upper = quantile(boot_means, 0.975))
}
los_paired_results <- bind_rows(
  paired_los_diff("Susceptible BSI"),
  paired_los_diff("Resistant BSI")
) %>% mutate(across(c(mean_paired_diff_days, ci_lower, ci_upper), ~ round(.x, 3)))

write_xlsx(
  list(paired_mortality_HR = mortality_paired_results, paired_LOS_mean_diff = los_paired_results),
  "Mortality_exact_1to1_paired_analysis.xlsx"
)

# -----------------------------------------------------------------------
# 9.17: Consolidated mortality summary -- every mechanism, every comparison,
# all four estimators side by side (Cox MSM IPW HR, Fine-Gray subdistribution
# HR, exact 1:1-matched OR, IPW risk difference), plus national extrapolation.
# -----------------------------------------------------------------------
# fit_ipw_cox_pooled()/fit_finegray_pairwise() both build their own
# "comparison" column as paste(group_comp, "vs", group_ref) -- identical
# text in msm and finegray for the same comparison, so they join cleanly on
# that column; mcnemar_paired_generic()'s own comparison text carries an
# extra " (exact 1:1 pair, McNemar)" suffix, stripped here before joining.
# 1:1-matched estimates do not exist for a Resistant-vs-Susceptible-type
# contrast (pairing is always case-vs-its-own-uninfected-control), so that
# row is legitimately blank, not missing data.
consolidate_mechanism <- function(mechanism_label, msm_tbl, fg_tbl, analytic_data) {
  n_by_group_tbl <- analytic_data %>% count(mortality_group)
  n_by_group <- setNames(n_by_group_tbl$n, as.character(n_by_group_tbl$mortality_group))
  msm_tbl %>%
    select(comparison, hazard_ratio, hr_lower, hr_upper, risk_difference_pct, rd_lower_pct, rd_upper_pct) %>%
    left_join(fg_tbl %>% select(comparison, subdistribution_hr, sHR_lower, sHR_upper), by = "comparison") %>%
    left_join(
      mortality_paired_results %>%
        mutate(comparison = str_remove(comparison, " \\(exact 1:1 pair, McNemar\\)")) %>%
        select(comparison, odds_ratio_haldane_corrected, ci_lower, ci_upper) %>%
        rename(matched_or = odds_ratio_haldane_corrected, matched_or_lower = ci_lower, matched_or_upper = ci_upper),
      by = "comparison"
    ) %>%
    mutate(
      mechanism = mechanism_label,
      ref_level = str_remove(comparison, "^.* vs "), comp_level = str_remove(comparison, " vs .*$"),
      n_ref = as.integer(n_by_group[ref_level]), n_comp = as.integer(n_by_group[comp_level]),
      .before = 1
    ) %>%
    select(-ref_level, -comp_level)
}

mortality_summary_all <- bind_rows(
  consolidate_mechanism("Any resistance", mortality_msm_results, mortality_competing_risk_results, mortality_analytic),
  consolidate_mechanism("MDR (bacterial isolates)", mdr_subgroup_results$msm, mdr_subgroup_results$finegray, mortality_analytic_mdr),
  consolidate_mechanism("XDR (bacterial isolates)", xdr_subgroup_results$msm, xdr_subgroup_results$finegray, mortality_analytic_xdr),
  consolidate_mechanism("Carbapenem-R (Gram-negative)", carb_subgroup_results$msm, carb_subgroup_results$finegray, mortality_analytic_carb),
  consolidate_mechanism("Methicillin-R/MRSA (S. aureus)", mrsa_subgroup_results$msm, mrsa_subgroup_results$finegray, mortality_analytic_mrsa),
  consolidate_mechanism("Cephalosporin-R any-gen (Gram-negative)", ceph_subgroup_results$msm, ceph_subgroup_results$finegray, mortality_analytic_ceph),
  consolidate_mechanism("3GCR (Gram-negative)", gcr3_subgroup_results$msm, gcr3_subgroup_results$finegray, mortality_analytic_3gcr),
  consolidate_mechanism("Access-resistant (AWaRe)", access_subgroup_results$msm, access_subgroup_results$finegray, mortality_analytic_access),
  consolidate_mechanism("Watch-resistant (AWaRe)", watch_subgroup_results$msm, watch_subgroup_results$finegray, mortality_analytic_watch),
  consolidate_mechanism("Reserve-resistant (AWaRe)", reserve_subgroup_results$msm, reserve_subgroup_results$finegray, mortality_analytic_reserve)
) %>%
  mutate(across(where(is.numeric) & !c(n_ref, n_comp), ~ round(.x, 4)))

# Vancomycin-R (Enterococcus): descriptive only throughout this script (3 VRE
# isolates, 0 events -- SECTION 9.15.4) -- appended as a labelled descriptive
# row, not a modelled one, rather than silently absent from the summary.
mortality_summary_all <- bind_rows(
  mortality_summary_all,
  tibble(mechanism = "Vancomycin-R (Enterococcus)", comparison = "VRE vs VSE (descriptive only, 3 VRE isolates)",
         n_ref = 17L, n_comp = 3L)
)
write_xlsx(mortality_summary_all, "Mortality_summary_all_mechanisms.xlsx")

# -----------------------------------------------------------------------
# 9.18: National extrapolation, every mechanism -- risk difference (Cox MSM
# IPW, the same adjusted/bias-corrected estimator already used for SECTION
# 9.10's primary national excess-death figure) x national mechanism-specific
# case count (SECTION 4's calibrated incidence x national admissions total).
# Consolidates each mechanism's own national_excess (already computed inside
# run_mortality_subgroup, one workbook per mechanism) into one table; the
# primary comparison's own SECTION 9.10 result is included for completeness.
# -----------------------------------------------------------------------
national_excess_all <- bind_rows(
  national_excess_deaths %>% mutate(mechanism = "Any resistance", .before = 1),
  mdr_subgroup_results$national_excess %>% mutate(mechanism = "MDR (bacterial isolates)", .before = 1),
  xdr_subgroup_results$national_excess %>% mutate(mechanism = "XDR (bacterial isolates)", .before = 1),
  carb_subgroup_results$national_excess %>% mutate(mechanism = "Carbapenem-R (Gram-negative)", .before = 1),
  mrsa_subgroup_results$national_excess %>% mutate(mechanism = "Methicillin-R/MRSA (S. aureus)", .before = 1),
  ceph_subgroup_results$national_excess %>% mutate(mechanism = "Cephalosporin-R any-gen (Gram-negative)", .before = 1),
  gcr3_subgroup_results$national_excess %>% mutate(mechanism = "3GCR (Gram-negative)", .before = 1),
  access_subgroup_results$national_excess %>% mutate(mechanism = "Access-resistant (AWaRe)", .before = 1),
  watch_subgroup_results$national_excess %>% mutate(mechanism = "Watch-resistant (AWaRe)", .before = 1),
  reserve_subgroup_results$national_excess %>% mutate(mechanism = "Reserve-resistant (AWaRe)", .before = 1)
)
write_xlsx(national_excess_all, "Mortality_national_excess_all_mechanisms.xlsx")

# -----------------------------------------------------------------------
# 9.19: Model-predicted (IPW-standardised) 28-day mortality curves, every
# mechanism -- generalises SECTION 9.12's overall combined-IPW curve
# (multinomial group ~ adjustment_covars propensity, standardised weight)
# to each mechanism's own 3-arm mortality_group and its own MICE object
# (reused from SECTION 9.15's run_mortality_subgroup(), not refit), then
# facets all mechanisms into one multi-panel figure.
# -----------------------------------------------------------------------
build_adjusted_curve_for_mechanism <- function(mice_obj_mech, analytic_data_mech, label, resistant_level, susceptible_level = NULL) {
  di <- complete(mice_obj_mech, 1)
  di$facility_id <- analytic_data_mech$facility_id
  form_multi <- as.formula(paste("mortality_group ~", paste(adjustment_covars, collapse = " + ")))
  di$mortality_group <- droplevels(di$mortality_group)
  multi_ps_mod <- tryCatch(multinom(form_multi, data = di, trace = FALSE), error = function(e) NULL)
  if (is.null(multi_ps_mod)) return(NULL)
  multi_probs <- predict(multi_ps_mod, type = "probs")
  if (is.null(dim(multi_probs))) {
    # Exactly 2 mortality_group levels (e.g. MRSA: Uninfected control / MRSA
    # BSI, MSSA excluded) -- multinom()/predict() drop to plain binomial
    # logit and return P(second level) as a bare vector, not a matrix.
    # Reconstruct the 2-column probability matrix multi_probs would give
    # for a >=3-level fit, so the rest of this function doesn't need to
    # branch on level count.
    lv <- levels(di$mortality_group)
    multi_probs <- cbind(1 - multi_probs, multi_probs)
    colnames(multi_probs) <- lv
  }
  multi_marg <- prop.table(table(di$mortality_group))
  own_prob <- multi_probs[cbind(seq_len(nrow(multi_probs)), match(as.character(di$mortality_group), colnames(multi_probs)))]
  ipw <- as.numeric(multi_marg[as.character(di$mortality_group)]) / own_prob
  ipw <- pmin(pmax(ipw, quantile(ipw, 0.01, na.rm = TRUE)), quantile(ipw, 0.99, na.rm = TRUE))
  di$ipw_combined3 <- ipw
  cdata <- make_adjusted_curve_data(di)
  if (!is.null(cdata)) {
    cdata$mechanism <- label
    # Relabel this mechanism's own group names to the 3 generic categories
    # (Resistant BSI / Susceptible BSI / Uninfected control) shared by every
    # panel, purely for a common colour scale across the faceted figure --
    # the mechanism-specific label (e.g. "MDR BSI") is preserved in the
    # panel title itself.
    relabel <- c("Uninfected control" = "Uninfected control")
    relabel[resistant_level] <- "Resistant BSI"
    if (!is.null(susceptible_level)) relabel[susceptible_level] <- "Susceptible BSI"
    cdata$mortality_group <- factor(
      unname(relabel[as.character(cdata$mortality_group)]),
      levels = c("Uninfected control", "Susceptible BSI", "Resistant BSI")
    )
  }
  cdata
}

mechanism_curve_data <- bind_rows(
  adjusted_overall_data %>% mutate(mechanism = "Any resistance"),
  build_adjusted_curve_for_mechanism(mdr_subgroup_results$mice_obj, mortality_analytic_mdr, "MDR", "MDR BSI", "Non-MDR BSI"),
  build_adjusted_curve_for_mechanism(xdr_subgroup_results$mice_obj, mortality_analytic_xdr, "XDR", "XDR BSI", "Non-XDR BSI"),
  build_adjusted_curve_for_mechanism(carb_subgroup_results$mice_obj, mortality_analytic_carb, "Carbapenem-R (Gram-negative)", "Carbapenem-R Gram-negative BSI", "Carbapenem-S Gram-negative BSI"),
  build_adjusted_curve_for_mechanism(mrsa_subgroup_results$mice_obj, mortality_analytic_mrsa, "MRSA (S. aureus)", "MRSA BSI"),
  build_adjusted_curve_for_mechanism(ceph_subgroup_results$mice_obj, mortality_analytic_ceph, "Cephalosporin-R any-gen (Gram-negative)", "Cephalosporin-R Gram-negative BSI", "Cephalosporin-S Gram-negative BSI"),
  build_adjusted_curve_for_mechanism(gcr3_subgroup_results$mice_obj, mortality_analytic_3gcr, "3GCR (Gram-negative)", "3GCR Gram-negative BSI", "3GC-susceptible Gram-negative BSI"),
  build_adjusted_curve_for_mechanism(access_subgroup_results$mice_obj, mortality_analytic_access, "Access-resistant (AWaRe)", "Access-R BSI", "Access-S BSI"),
  build_adjusted_curve_for_mechanism(watch_subgroup_results$mice_obj, mortality_analytic_watch, "Watch-resistant (AWaRe)", "Watch-R BSI", "Watch-S BSI"),
  build_adjusted_curve_for_mechanism(reserve_subgroup_results$mice_obj, mortality_analytic_reserve, "Reserve-resistant (AWaRe)", "Reserve-R BSI", "Reserve-S BSI")
)
mechanism_curve_data$mechanism <- factor(
  mechanism_curve_data$mechanism,
  levels = c("Any resistance", "MDR", "XDR", "Carbapenem-R (Gram-negative)", "MRSA (S. aureus)",
             "Cephalosporin-R any-gen (Gram-negative)", "3GCR (Gram-negative)",
             "Access-resistant (AWaRe)", "Watch-resistant (AWaRe)", "Reserve-resistant (AWaRe)")
)
saveRDS(mechanism_curve_data, "mechanism_curve_data.rds")

mechanism_curves_plot <- ggplot(mechanism_curve_data, aes(x = time, y = 1 - estimate, color = mortality_group, fill = mortality_group)) +
  geom_ribbon(aes(ymin = 1 - conf.low, ymax = 1 - conf.high), alpha = 0.12, color = NA) +
  geom_step(linewidth = 0.7) +
  scale_color_manual(values = mort_colors, name = NULL) +
  scale_fill_manual(values = mort_colors, guide = "none") +
  scale_x_continuous(breaks = seq(0, 28, 14), limits = c(0, 28), expand = c(0, 0)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  facet_wrap(~mechanism, ncol = 3) +
  labs(
    title = "Model-predicted (IPW-standardised) 28-day mortality, by resistance mechanism",
    subtitle = "Each panel: mechanism-specific Resistant / Susceptible-or-Non-mechanism / Uninfected control, own combined-IPW standardisation",
    x = "Days since enrolment", y = "Adjusted cumulative mortality"
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = MORT_LINE, linewidth = 0.4),
        legend.position = "top", strip.text = element_text(face = "bold"),
        plot.subtitle = element_text(color = MORT_INK_SOFT), axis.title = element_text(color = MORT_INK_SOFT))
ggsave("Mortality_adjusted_curves_by_mechanism.png", mechanism_curves_plot, width = 12, height = 12, dpi = 500)

# -----------------------------------------------------------------------
# 9.20: The propensity models themselves -- every mechanism, every pairwise
# comparison. SECTION 9.7-9.19's IPW machinery fits a fresh propensity model
# (group_bin ~ adjustment_covars) inside fit_ipw_cox_pooled()'s own loop
# over all 20 MICE-imputed datasets, but never previously exposed or
# reported those models' own coefficients -- only the resulting weights and
# pooled outcome HR. Reported here on the FIRST imputed dataset only (the
# same "a representative single completion" convention SECTION 9.12 already
# uses for the predicted-curve figures, not a full 20-way Rubin's-rules
# pool), since this table's purpose is transparency about what predicts
# group membership in each model, not a primary inferential estimate --
# the causal HR itself remains SECTION 9.17's pooled Cox MSM estimate.
## fit_propensity_model() fits the SAME single glm() as before, but now
## returns both the tidy coefficients AND the per-patient fitted propensity
## scores (with a generic "Reference group"/"Comparator group" role tag) in
## one call, so the overlap-density figure below and the coefficient table
## come from the identical model fit rather than refitting it twice.
fit_propensity_model <- function(mice_obj, group_ref, group_comp, mechanism_label, comparison_type) {
  di <- complete(mice_obj, 1) %>%
    dplyr::filter(mortality_group %in% c(group_ref, group_comp)) %>%
    mutate(group_bin = as.integer(mortality_group == group_comp))
  form_ps <- as.formula(paste("group_bin ~", paste(adjustment_covars, collapse = " + ")))
  mod <- tryCatch(glm(form_ps, data = di, family = binomial()), error = function(e) NULL)
  comparison_label <- paste(group_comp, "vs", group_ref)
  if (is.null(mod)) {
    return(list(
      coefs = tibble(mechanism = mechanism_label, comparison = comparison_label,
                     term = NA_character_, odds_ratio = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, p_value = NA_real_),
      scores = tibble(mechanism = mechanism_label, comparison = comparison_label, comparison_type = comparison_type,
                       group = character(0), role = character(0), propensity_score = numeric(0))
    ))
  }
  coefs <- broom::tidy(mod, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::filter(term != "(Intercept)") %>%
    transmute(mechanism = mechanism_label, comparison = comparison_label,
              term, odds_ratio = estimate, ci_lower = conf.low, ci_upper = conf.high, p_value = p.value) %>%
    mutate(across(c(odds_ratio, ci_lower, ci_upper, p_value), ~ round(.x, 4)))
  scores <- tibble(
    mechanism = mechanism_label, comparison = comparison_label, comparison_type = comparison_type,
    group = if_else(di$group_bin == 1, group_comp, group_ref),
    role = if_else(di$group_bin == 1, "Comparator group", "Reference group"),
    propensity_score = as.numeric(predict(mod, type = "response"))
  )
  list(coefs = coefs, scores = scores)
}

# comparison_type is passed explicitly (not inferred from list position)
# because MRSA's own comparisons list has only 1 element (Resistant vs
# Uninfected -- MSSA excluded, SECTION 9.15.3) where every other mechanism
# has 3, so position-based inference would mislabel it. These 3 generic
# types are what let the overlap figure below use a mechanism x
# comparison-type GRID (facet_grid) instead of one long, unreadable strip
# of 19 individually-labelled panels.
TYPE_SUS <- "Susceptible / non-mechanism\nvs Uninfected"
TYPE_RES <- "Resistant / mechanism\nvs Uninfected"
TYPE_HEAD <- "Resistant\nvs Susceptible / non-mechanism"

propensity_fits_all <- list(
  fit_propensity_model(mortality_mice, "Uninfected control", "Susceptible BSI", "Any resistance", TYPE_SUS),
  fit_propensity_model(mortality_mice, "Uninfected control", "Resistant BSI", "Any resistance", TYPE_RES),
  fit_propensity_model(mortality_mice, "Susceptible BSI", "Resistant BSI", "Any resistance", TYPE_HEAD),
  fit_propensity_model(mdr_subgroup_results$mice_obj, "Uninfected control", "Non-MDR BSI", "MDR", TYPE_SUS),
  fit_propensity_model(mdr_subgroup_results$mice_obj, "Uninfected control", "MDR BSI", "MDR", TYPE_RES),
  fit_propensity_model(mdr_subgroup_results$mice_obj, "Non-MDR BSI", "MDR BSI", "MDR", TYPE_HEAD),
  fit_propensity_model(xdr_subgroup_results$mice_obj, "Uninfected control", "Non-XDR BSI", "XDR", TYPE_SUS),
  fit_propensity_model(xdr_subgroup_results$mice_obj, "Uninfected control", "XDR BSI", "XDR", TYPE_RES),
  fit_propensity_model(xdr_subgroup_results$mice_obj, "Non-XDR BSI", "XDR BSI", "XDR", TYPE_HEAD),
  fit_propensity_model(carb_subgroup_results$mice_obj, "Uninfected control", "Carbapenem-S Gram-negative BSI", "Carbapenem-R (Gram-neg.)", TYPE_SUS),
  fit_propensity_model(carb_subgroup_results$mice_obj, "Uninfected control", "Carbapenem-R Gram-negative BSI", "Carbapenem-R (Gram-neg.)", TYPE_RES),
  fit_propensity_model(carb_subgroup_results$mice_obj, "Carbapenem-S Gram-negative BSI", "Carbapenem-R Gram-negative BSI", "Carbapenem-R (Gram-neg.)", TYPE_HEAD),
  fit_propensity_model(mrsa_subgroup_results$mice_obj, "Uninfected control", "MRSA BSI", "MRSA (S. aureus)", TYPE_RES),
  fit_propensity_model(ceph_subgroup_results$mice_obj, "Uninfected control", "Cephalosporin-S Gram-negative BSI", "Cephalosporin-R any-gen\n(Gram-neg.)", TYPE_SUS),
  fit_propensity_model(ceph_subgroup_results$mice_obj, "Uninfected control", "Cephalosporin-R Gram-negative BSI", "Cephalosporin-R any-gen\n(Gram-neg.)", TYPE_RES),
  fit_propensity_model(ceph_subgroup_results$mice_obj, "Cephalosporin-S Gram-negative BSI", "Cephalosporin-R Gram-negative BSI", "Cephalosporin-R any-gen\n(Gram-neg.)", TYPE_HEAD),
  fit_propensity_model(gcr3_subgroup_results$mice_obj, "Uninfected control", "3GC-susceptible Gram-negative BSI", "3GCR (Gram-neg.)", TYPE_SUS),
  fit_propensity_model(gcr3_subgroup_results$mice_obj, "Uninfected control", "3GCR Gram-negative BSI", "3GCR (Gram-neg.)", TYPE_RES),
  fit_propensity_model(gcr3_subgroup_results$mice_obj, "3GC-susceptible Gram-negative BSI", "3GCR Gram-negative BSI", "3GCR (Gram-neg.)", TYPE_HEAD),
  fit_propensity_model(access_subgroup_results$mice_obj, "Uninfected control", "Access-S BSI", "Access-resistant\n(AWaRe)", TYPE_SUS),
  fit_propensity_model(access_subgroup_results$mice_obj, "Uninfected control", "Access-R BSI", "Access-resistant\n(AWaRe)", TYPE_RES),
  fit_propensity_model(access_subgroup_results$mice_obj, "Access-S BSI", "Access-R BSI", "Access-resistant\n(AWaRe)", TYPE_HEAD),
  fit_propensity_model(watch_subgroup_results$mice_obj, "Uninfected control", "Watch-S BSI", "Watch-resistant\n(AWaRe)", TYPE_SUS),
  fit_propensity_model(watch_subgroup_results$mice_obj, "Uninfected control", "Watch-R BSI", "Watch-resistant\n(AWaRe)", TYPE_RES),
  fit_propensity_model(watch_subgroup_results$mice_obj, "Watch-S BSI", "Watch-R BSI", "Watch-resistant\n(AWaRe)", TYPE_HEAD),
  fit_propensity_model(reserve_subgroup_results$mice_obj, "Uninfected control", "Reserve-S BSI", "Reserve-resistant\n(AWaRe)", TYPE_SUS),
  fit_propensity_model(reserve_subgroup_results$mice_obj, "Uninfected control", "Reserve-R BSI", "Reserve-resistant\n(AWaRe)", TYPE_RES),
  fit_propensity_model(reserve_subgroup_results$mice_obj, "Reserve-S BSI", "Reserve-R BSI", "Reserve-resistant\n(AWaRe)", TYPE_HEAD)
)
propensity_models_all <- purrr::map_dfr(propensity_fits_all, "coefs")
propensity_scores_all <- purrr::map_dfr(propensity_fits_all, "scores")
write_xlsx(propensity_models_all, "Mortality_propensity_models_all_mechanisms.xlsx")

# -----------------------------------------------------------------------
# 9.20b: Propensity-score overlap, every mechanism and comparison -- kernel
# density of the fitted propensity score, reference vs comparator group, in
# a mechanism (row) x comparison-type (column) GRID -- organised
# "correspondingly with resistance type" per instruction, and legible: an
# earlier facet_wrap version with one long mechanism+comparison string per
# panel produced 19 unreadably truncated strip titles. MRSA has no
# Susceptible/head-to-head panel (MSSA excluded, SECTION 9.15.3) and 3GCR
# has no matched-pair panel dependency here, so those grid cells are simply
# blank -- not missing data, the same convention already used for blank
# cells in Table 7 (main text). This is the standard IPW positivity/common-
# support diagnostic (SECTION 9.11 already reports it numerically --
# propensity_score_min/max -- for the 3 primary comparisons only; this
# figure extends it visually to every mechanism).
mech_row_order <- c("Any resistance", "MDR", "XDR", "Carbapenem-R (Gram-neg.)",
                     "MRSA (S. aureus)", "Cephalosporin-R any-gen\n(Gram-neg.)", "3GCR (Gram-neg.)",
                     "Access-resistant\n(AWaRe)", "Watch-resistant\n(AWaRe)", "Reserve-resistant\n(AWaRe)")
propensity_scores_all$mechanism <- factor(propensity_scores_all$mechanism, levels = mech_row_order)
propensity_scores_all$comparison_type <- factor(propensity_scores_all$comparison_type, levels = c(TYPE_SUS, TYPE_RES, TYPE_HEAD))
PROP_COLORS <- c("Reference group" = "#17766E", "Comparator group" = "#7A1F2B")
propensity_overlap_plot <- ggplot(propensity_scores_all, aes(x = propensity_score, fill = role, color = role)) +
  geom_density(alpha = 0.35, linewidth = 0.55, na.rm = TRUE) +
  scale_fill_manual(values = PROP_COLORS, name = NULL) +
  scale_color_manual(values = PROP_COLORS, name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  facet_grid(mechanism ~ comparison_type, scales = "free_y", switch = "y") +
  labs(
    title = "Propensity-score overlap, every mortality mechanism and comparison",
    subtitle = "Kernel density of the fitted probability of comparator-group membership (Methods); rows = resistance mechanism, columns = comparison type",
    x = "Propensity score", y = "Density"
  ) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        strip.text.y = element_text(face = "bold", size = 7.6, angle = 0),
        strip.text.x = element_text(face = "bold", size = 7.6),
        plot.subtitle = element_text(color = MORT_INK_SOFT))
ggsave("Mortality_propensity_overlap_all_mechanisms.png", propensity_overlap_plot, width = 11, height = 21, dpi = 500)


# =============================================================================
# SECTION 10: Risk factors for drug resistance
# =============================================================================
# Corresponds to the methods spec's "5. Risk Factors for Drug Resistance".
# Objective: characterise patient- and healthcare-related factors associated
# with drug-resistant BSI among CONFIRMED BSI CASES -- this is an
# association/descriptive analysis, NOT the causal mortality estimation of
# SECTION 9 (no IPW, no MICE, no MSM). Design-based survey-weighted
# multivariable logistic regression, PPS facility cluster (facility_id_clean,
# SECTIONS 3.2-3.3's design), each mechanism run against BOTH a "Susceptible"
# comparator (same-mechanism, same-organism-restriction, confirmed BSI) AND
# an "Uninfected" comparator, each as its own pairwise svyglm (not a 3-arm
# multinomial), per the request to do this "separately".
#
# "Uninfected" HERE = Main_data_CaseReport's own BSI-negative patients
# (suspected-BSI, blood culture negative) -- these ARE part of the original
# PPS survey sample with a resolvable design_weight/calibrated_weight
# (SECTION 3.2), unlike SECTION 9's "Uninfected control" arm, which is a
# SEPARATELY RECRUITED matched cohort (SECTION 7) carrying no survey design
# weight at all. A design-based survey-weighted model is only defined for a
# population that is actually part of the sampling frame, so CaseReport's
# own BSI-negative patients -- not SECTION 9's matched controls -- are the
# only valid "Uninfected" comparator for this specific method (design-based
# survey-weighted regression). Do not conflate the two "Uninfected" groups
# across SECTIONS 9 and 10 -- they are different populations answering
# different questions (causal 28-day mortality vs cross-sectional resistance
# association).
#
# Candidate risk factors per the spec: age, sex, infection origin, ward type,
# grouped clinical specialty, hospital type, pre-culture LOS, comorbidity
# burden, prior antibiotic exposure, broad pathogen group -- RESTRICTED here
# to only those that are logically upstream of the blood culture that
# determines resistance status, i.e. genuinely "risk factors" rather than
# things that could themselves be a CONSEQUENCE of (or concurrent with) an
# already-severe/resistant infection:
#   - ward_type_icu (ICU vs Non-ICU) and specialty_group (Medical/Surgical)
#     are DROPPED. Both describe where the patient is AT/AROUND the time of
#     culture -- ICU admission in particular is plausibly a downstream
#     consequence of a severe (possibly resistant) infection, not a
#     pre-existing exposure that caused it, so using it as a "risk factor"
#     risks reverse causation.
#   - broad pathogen group (major_bacterial_group) is DROPPED from every
#     comparison, including vs-Susceptible (the spec listed it as a
#     candidate, but it is identified FROM the very culture whose resistance
#     status is being modelled -- concurrent, not "previous to BSI culture",
#     so it does not belong in a risk-factor model even though it was a
#     defensible covariate in SECTION 9's outcome/prognosis models).
#   - infection_origin_composite (hospital- vs community-acquired) is KEPT:
#     it reflects admission history/timing BEFORE the culture (hospital-
#     acquired = >48h hospital exposure prior to the culture that identified
#     the BSI), a genuine upstream exposure.
#   - Level_of_care (hospital type) is KEPT: a facility-level property fixed
#     before the patient's admission, not affected by their clinical course.
#   - abx_current_therapy is KEPT as the best available proxy for prior
#     antibiotic exposure (collected under "Antimicrobial use history", i.e.
#     exposure leading into this episode, not treatment initiated because of
#     this culture's result).
#   - age_years, sex_at_birth, charlson_score (comorbidity burden) are KEPT:
#     all fixed patient characteristics that predate this infection episode.
#   - riskfx_hospitalisation_3m ("Hospitalisation in the 3 months prior to
#     this admission") is ADDED (not in the spec's own candidate list, but a
#     textbook AMR risk factor and unambiguously pre-admission -- it is
#     collected under "Risk factors AT ADMISSION", i.e. before this episode's
#     culture even exists). 6.7-8.6% missing, well within the 40% ceiling.
#   - surgery_prior_to_bsi ("Did the patient have surgery on this admission
#     prior to onset of clinical signs of suspected BSI?") is ADDED -- also
#     not in the spec's list, but explicitly labelled prior-to-onset in its
#     own source question, i.e. unambiguously upstream of both the clinical
#     suspicion of BSI and the culture. 1.8-2.3% missing.
# Missingness by BSI status (checked against the full CaseReport population):
# pre-culture LOS (hours_admission_to_sample1) -- itself a legitimate
# upstream exposure (longer pre-culture hospital stay = more opportunity for
# nosocomial acquisition) -- is 66.7% missing among BSI-positive and 76.7%
# among BSI-negative patients, far above SECTION 9.3's 40% ceiling, and is
# dropped for that reason (not multiply imputed; this section, unlike
# SECTION 9, does not use MICE, per the spec's own "composite indicators"
# instruction rather than a full MI pipeline). Every retained candidate is
# <=34% missing.
#
# EVENTS-PER-VARIABLE / DESIGN CHANGE (per instruction): "Uninfected" here is
# now SECTION 9's own separately-recruited matched-cohort control population
# (matched_control_broad == 1), NOT CaseReport's own BSI-negative patients as
# in an earlier version of this section -- specifically so that the
# Resistant-vs-Uninfected model IS the same specification (same population,
# same covariates: adjustment_covars, SECTION 9.4) as the propensity-score
# model SECTION 9 fits per MICE-imputed dataset for IPW weight construction
# (fit_ipw_cox_pooled()'s own ps_mod <- glm(form_ps, ...)). This section
# reports a single complete-case fit of that identical specification (cheap,
# interpretable ORs); SECTION 9's causal HR/IPW analysis re-fits the same
# formula across the 20 MICE-imputed datasets pooled via Rubin's rules --
# same model, different purpose (association reporting here vs bias-
# correcting weight construction there), not two different models.
#
# CONSEQUENCE: the matched-cohort Uninfected population carries no PPS
# design/calibrated survey weight at all (SECTION 9's own documented reason
# for using cluster-robust IPW rather than a design-weighted model) -- so
# the Resistant-vs-Uninfected comparison is necessarily UNWEIGHTED
# (facility-cluster-robust svyglm with equal probabilities) only, never
# calibrated. The Resistant-vs-Susceptible comparison stays within
# confirmed-BSI cases (which DO carry calibrated_weight), so it is reported
# both unweighted and calibrated, as everywhere else in this report.
#
# Covariates are the identical adjustment_covars used by SECTION 9's own
# propensity model (age, sex, hospital level of care, prior antibiotic
# exposure, comorbidity burden) for every mechanism and every comparison --
# no separate "reduced"/"full" covariate set by comparison as in an earlier
# version of this section, precisely because the point is now one shared
# specification, not two.
riskfactor_covars <- adjustment_covars

fit_riskfactor_model <- function(mechanism, data, ref_level, comp_level, weight_var) {
  di <- data %>% dplyr::filter(mortality_group_mech %in% c(ref_level, comp_level), !is.na(facility_id))
  if (weight_var != "unweighted") di <- di %>% dplyr::filter(!is.na(.data[[weight_var]]))
  n_ref <- sum(di$mortality_group_mech == ref_level); n_comp <- sum(di$mortality_group_mech == comp_level)
  des <- if (weight_var == "unweighted") {
    svydesign(ids = ~facility_id, data = di)
  } else {
    svydesign(ids = ~facility_id, weights = as.formula(paste0("~", weight_var)), data = di)
  }
  des <- update(des, group_bin = as.integer(mortality_group_mech == comp_level))
  form <- as.formula(paste("group_bin ~", paste(riskfactor_covars, collapse = " + ")))
  result <- tryCatch({
    mod <- suppressWarnings(svyglm(form, design = des, family = quasibinomial()))
    broom::tidy(mod, conf.int = TRUE, exponentiate = TRUE) %>%
      dplyr::filter(term != "(Intercept)") %>%
      transmute(term, odds_ratio = estimate, ci_lower = conf.low, ci_upper = conf.high, p_value = p.value)
  }, error = function(e) tibble(term = NA_character_, odds_ratio = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, p_value = NA_real_))
  result %>% mutate(mechanism = mechanism, comparison = paste(comp_level, "vs", ref_level), weight_type = weight_var,
                     n_ref = n_ref, n_comp = n_comp, .before = 1) %>%
    mutate(across(c(odds_ratio, ci_lower, ci_upper, p_value), ~ round(.x, 4)))
}

## --- mechanism-specific group definitions, built on mortality_analytic
## (SECTION 9's own population: matched-cohort Uninfected + confirmed-BSI
## cases) so "Uninfected" means the same population in every mechanism here
## as it does in SECTION 9's propensity model -- mirrors SECTION 9.15's own
## organism restrictions for the drug-bug-specific mechanisms. ---
riskfactor_mechanisms <- list(
  list(mechanism = "Any resistance", ref = "Susceptible", comp = "Resistant", group_fn = function(d) {
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              d$resistant_BSI_positive == 1 ~ "Resistant",
              d$resistant_BSI_positive == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  list(mechanism = "MDR (bacterial isolates)", ref = "Susceptible", comp = "Resistant", group_fn = function(d) {
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              d$mdr_BSI_bacterial == 1 ~ "Resistant",
              d$mdr_BSI_bacterial == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  list(mechanism = "XDR (bacterial isolates)", ref = "Susceptible", comp = "Resistant", group_fn = function(d) {
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              d$xdr_BSI_bacterial == 1 ~ "Resistant",
              d$xdr_BSI_bacterial == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  # Population-wide ("all confirmed BSI", no organism restriction) risk-
  # factor models for the 4 mechanisms that also have an organism-restricted
  # variant -- mirroring the dual-row pattern SECTION 3.6's own prevalence
  # table (Table 1/2, main text) already uses for these same 4 mechanisms
  # ("(all confirmed BSI)" alongside the organism-restricted "X only" row).
  # Risk Factors previously reported only the organism-restricted variant
  # for each; added so Table 3 is consistent with Table 2's own convention,
  # per instruction.
  list(mechanism = "Carbapenem-R (all confirmed BSI)", ref = "Susceptible", comp = "Resistant", group_fn = function(d) {
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              d$carb_BSI_positive == 1 ~ "Resistant",
              d$carb_BSI_positive == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  list(mechanism = "Carbapenem-R (Gram-negative)", ref = "Susceptible", comp = "Resistant", group_fn = function(d) {
    gneg <- is_cre_pa_ab(d$major_bacterial_group, d$pathogen1_clean)
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              gneg & d$carb_BSI_positive == 1 ~ "Resistant",
              gneg & d$carb_BSI_positive == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  list(mechanism = "Cephalosporin-R any-gen (all confirmed BSI)", ref = "Susceptible", comp = "Resistant", group_fn = function(d) {
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              d$cephalosporin_BSI_positive == 1 ~ "Resistant",
              d$cephalosporin_BSI_positive == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  list(mechanism = "Cephalosporin-R any-gen (Gram-negative)", ref = "Susceptible", comp = "Resistant", group_fn = function(d) {
    gneg <- is_cre_pa_ab(d$major_bacterial_group, d$pathogen1_clean)
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              gneg & d$cephalosporin_BSI_positive == 1 ~ "Resistant",
              gneg & d$cephalosporin_BSI_positive == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  list(mechanism = "3GCR (all confirmed BSI)", ref = "Susceptible", comp = "Resistant", group_fn = function(d) {
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              d$c3GC_BSI_positive == 1 ~ "Resistant",
              d$c3GC_BSI_positive == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  list(mechanism = "3GCR (Gram-negative)", ref = "Susceptible", comp = "Resistant", group_fn = function(d) {
    gneg <- is_cre_pa_ab(d$major_bacterial_group, d$pathogen1_clean)
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              gneg & d$c3GC_BSI_positive == 1 ~ "Resistant",
              gneg & d$c3GC_BSI_positive == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  list(mechanism = "Methicillin-R (all confirmed BSI)", ref = "Susceptible", comp = "Resistant", group_fn = function(d) {
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              d$methicillin_BSI_positive == 1 ~ "Resistant",
              d$methicillin_BSI_positive == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  list(mechanism = "Methicillin-R / MRSA (S. aureus)", ref = "Susceptible", comp = "Resistant", group_fn = function(d) {
    staph <- d$pathogen1_clean == "Staphylococcus aureus"
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              staph & d$methicillin_BSI_positive == 1 ~ "Resistant",
              staph & d$methicillin_BSI_positive == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  list(mechanism = "Access-resistant (AWaRe)", ref = "Susceptible", comp = "Resistant", group_fn = function(d) {
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              d$access_BSI_positive == 1 ~ "Resistant",
              d$access_BSI_positive == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  list(mechanism = "Watch-resistant (AWaRe)", ref = "Susceptible", comp = "Resistant", group_fn = function(d) {
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              d$watch_BSI_positive == 1 ~ "Resistant",
              d$watch_BSI_positive == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  list(mechanism = "Reserve-resistant (AWaRe)", ref = "Susceptible", comp = "Resistant", group_fn = function(d) {
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              d$reserve_BSI_positive == 1 ~ "Resistant",
              d$reserve_BSI_positive == 0 ~ "Susceptible", TRUE ~ NA_character_)
  })
)

## Resistant-vs-Susceptible: both confirmed-BSI arms carry calibrated_weight,
## so reported unweighted AND calibrated. Resistant-vs-Uninfected: the
## matched-cohort Uninfected arm carries no survey weight at all (see
## header comment above), so unweighted only -- this is deliberately the
## exact specification (population + covariates) SECTION 9's own propensity
## model reuses for IPW.
risk_factor_or_table <- purrr::map_dfr(riskfactor_mechanisms, function(mech) {
  d2 <- mortality_analytic %>% mutate(mortality_group_mech = mech$group_fn(.))
  bind_rows(
    fit_riskfactor_model(mech$mechanism, d2, "Susceptible", "Resistant", "unweighted"),
    fit_riskfactor_model(mech$mechanism, d2, "Susceptible", "Resistant", "calibrated_weight"),
    fit_riskfactor_model(mech$mechanism, d2, "Uninfected", "Resistant", "unweighted")
  )
})
write_xlsx(risk_factor_or_table, "RiskFactors_drugresistance_adjusted_OR.xlsx")

## --- CoNS-treatment sensitivity check (genuine adjusted re-fit, "Any
## resistance" mechanism only) -- unlike SECTION 9's mortality models,
## fit_riskfactor_model() is a single fast svyglm() call with no bootstrap
## or imputation, so a real adjusted re-fit under each alternate CoNS
## population (SECTION 2.6) is cheap enough to do properly here, not just
## a crude comparison. The primary risk-factor analysis throughout this
## report always uses the unchanged group definition above; this is an
## additional sensitivity check, never a replacement.
riskfactor_group_fn_CoNS_excluded <- function(d) {
  case_when(d$matched_control_broad == 1 ~ "Uninfected",
            d$BSI_positive_CoNS_sensitive == 1 & d$resistant_BSI_positive == 1 ~ "Resistant",
            d$BSI_positive_CoNS_sensitive == 1 & d$resistant_BSI_positive == 0 ~ "Susceptible", TRUE ~ NA_character_)
}
riskfactor_group_fn_CoNS_included <- function(d) {
  case_when(d$matched_control_broad == 1 ~ "Uninfected",
            d$BSI_positive_CoNS_included == 1 & d$resistant_BSI_positive_CoNS_included == 1 ~ "Resistant",
            d$BSI_positive_CoNS_included == 1 & d$resistant_BSI_positive_CoNS_included == 0 ~ "Susceptible", TRUE ~ NA_character_)
}

# mortality_analytic itself is already filtered to the PRIMARY population
# (!is.na(mortality_group), SECTION 9.1) -- the ~122 patients who are
# confirmed BSI only under the CoNS-included definition (contaminated
# cultures that nonetheless reported a CoNS species) have NA primary
# mortality_group and so are silently ABSENT from mortality_analytic
# already, before riskfactor_group_fn_CoNS_included() ever runs. Every
# CoNS-included group_fn below is therefore built on this broader base
# instead -- every Uninfected patient plus every patient confirmed BSI
# under the CoNS-included definition specifically -- so those ~122 patients
# are actually present to be classified, not invisibly dropped.
mortality_analytic_cons_included_base <- Main_data_Combined %>%
  filter(matched_control_broad == 1 | BSI_positive_CoNS_included == 1)

# Level_of_care and charlson_score_combined are engineered directly onto
# mortality_analytic (SECTION 9.3, above) and never onto Main_data_Combined
# itself -- so this independently-filtered base object needs the same two
# derivations reapplied, or every riskfactor_covars-based model built on it
# fails with "all variables must be in design= argument" (svyglm cannot see
# a formula term that isn't a column of the svydesign's own data).
mortality_analytic_cons_included_base <- mortality_analytic_cons_included_base %>%
  left_join(Hospital_characteristics %>% select(`Hospital ID`, Level_of_care), by = c(facility_id = "Hospital ID"))

.charlson_diabetes_pts2 <- case_when(mortality_analytic_cons_included_base$comorb_diabetes_end_organ == 1L ~ 2L, mortality_analytic_cons_included_base$comorb_diabetes == 1L ~ 1L, TRUE ~ 0L)
.charlson_liver_pts2 <- case_when(mortality_analytic_cons_included_base$comorb_liver_severe == 1L ~ 3L, mortality_analytic_cons_included_base$comorb_liver_mild == 1L ~ 1L, TRUE ~ 0L)
.charlson_malignancy_pts2 <- case_when(mortality_analytic_cons_included_base$comorb_metastatic_tumour == 1L ~ 6L, mortality_analytic_cons_included_base$comorb_malignancy_active == 1L ~ 2L, TRUE ~ 0L)
.charlson_hiv_pts2 <- if_else(mortality_analytic_cons_included_base$comorb_hiv_on_art == 1L | mortality_analytic_cons_included_base$comorb_hiv_no_art == 1L | mortality_analytic_cons_included_base$comorb_hiv_cd4_low == 1L, 6L, 0L, missing = 0L)
.charlson_weight1_pts2 <- mortality_analytic_cons_included_base %>%
  select(comorb_chf, comorb_cerebrovascular, comorb_dementia, comorb_copd, comorb_connective_tissue, comorb_peptic_ulcer) %>%
  mutate(across(everything(), ~ if_else(.x == 1L, 1L, 0L, missing = 0L))) %>%
  rowSums()
.charlson_hemi_pts2 <- if_else(mortality_analytic_cons_included_base$comorb_hemi_paraplegia == 1L, 2L, 0L, missing = 0L)
.charlson_renal_pts2 <- if_else(mortality_analytic_cons_included_base$comorb_renal_disease == 1L, 2L, 0L, missing = 0L)
mortality_analytic_cons_included_base$charlson_score_combined <- .charlson_weight1_pts2 + .charlson_diabetes_pts2 + .charlson_liver_pts2 +
  .charlson_malignancy_pts2 + .charlson_hiv_pts2 + .charlson_hemi_pts2 + .charlson_renal_pts2
rm(.charlson_diabetes_pts2, .charlson_liver_pts2, .charlson_malignancy_pts2, .charlson_hiv_pts2,
   .charlson_weight1_pts2, .charlson_hemi_pts2, .charlson_renal_pts2)

cat("=== RAW CRUDE MORTALITY BY GROUP (for new Results Figure Panel A) ===\n")
crude_mort_row <- function(data, group_label, display_label) {
  sub <- data %>% dplyr::filter(mortality_group == group_label)
  n <- nrow(sub)
  n_known <- sum(!is.na(sub$event28))
  deaths <- sum(sub$event28 == 1, na.rm = TRUE)
  tibble(group = display_label, n = n, n_outcome_known = n_known, deaths_28day = deaths,
         mortality_pct_28day = round(100 * deaths / n_known, 1))
}
vre_split <- mortality_analytic %>%
  dplyr::filter(matched_control == 0, !is.na(genus1), genus1 == "Enterococcus", !is.na(vancomycin_BSI_positive)) %>%
  mutate(grp = if_else(vancomycin_BSI_positive == 1, "VRE", "VSE"))
vre_row <- function(lvl, display_label) {
  sub <- vre_split %>% dplyr::filter(grp == lvl)
  n <- nrow(sub); n_known <- sum(!is.na(sub$event28)); deaths <- sum(sub$event28 == 1, na.rm = TRUE)
  tibble(group = display_label, n = n, n_outcome_known = n_known, deaths_28day = deaths,
         mortality_pct_28day = round(100 * deaths / n_known, 1))
}

crude_mortality_all <- bind_rows(
  crude_mort_row(mortality_analytic, "Uninfected control", "Uninfected"),
  crude_mort_row(mortality_analytic, "Susceptible BSI", "Susceptible (any)"),
  crude_mort_row(mortality_analytic, "Resistant BSI", "Resistant (any)"),
  crude_mort_row(mortality_analytic_mdr, "Non-MDR BSI", "Non-MDR"),
  crude_mort_row(mortality_analytic_mdr, "MDR BSI", "MDR"),
  crude_mort_row(mortality_analytic_xdr, "Non-XDR BSI", "Non-XDR"),
  crude_mort_row(mortality_analytic_xdr, "XDR BSI", "XDR"),
  crude_mort_row(mortality_analytic_carb, "Carbapenem-S Gram-negative BSI", "Carbapenem-S"),
  crude_mort_row(mortality_analytic_carb, "Carbapenem-R Gram-negative BSI", "Carbapenem-R"),
  crude_mort_row(mortality_analytic_ceph, "Cephalosporin-S Gram-negative BSI", "Cephalosporin-S"),
  crude_mort_row(mortality_analytic_ceph, "Cephalosporin-R Gram-negative BSI", "Cephalosporin-R"),
  crude_mort_row(mortality_analytic_3gcr, "3GC-susceptible Gram-negative BSI", "3GC-S"),
  crude_mort_row(mortality_analytic_3gcr, "3GCR Gram-negative BSI", "3GCR"),
  crude_mort_row(mortality_analytic_mrsa, "MRSA BSI", "MRSA"),
  vre_row("VSE", "VSE"),
  vre_row("VRE", "VRE"),
  crude_mort_row(mortality_analytic_access, "Access-S BSI", "Access-S"),
  crude_mort_row(mortality_analytic_access, "Access-R BSI", "Access-R"),
  crude_mort_row(mortality_analytic_watch, "Watch-S BSI", "Watch-S"),
  crude_mort_row(mortality_analytic_watch, "Watch-R BSI", "Watch-R"),
  crude_mort_row(mortality_analytic_reserve, "Reserve-S BSI", "Reserve-S"),
  crude_mort_row(mortality_analytic_reserve, "Reserve-R BSI", "Reserve-R")
)
print(as.data.frame(crude_mortality_all), row.names = FALSE)
write_xlsx(crude_mortality_all, "Mortality_crude_by_group_all_mechanisms.xlsx")

risk_factor_or_cons_sensitivity <- bind_rows(
  {
    d2 <- mortality_analytic %>% mutate(mortality_group_mech = riskfactor_group_fn_CoNS_excluded(.))
    bind_rows(
      fit_riskfactor_model("CoNS-excluded sensitivity", d2, "Susceptible", "Resistant", "unweighted"),
      fit_riskfactor_model("CoNS-excluded sensitivity", d2, "Susceptible", "Resistant", "calibrated_weight"),
      fit_riskfactor_model("CoNS-excluded sensitivity", d2, "Uninfected", "Resistant", "unweighted")
    )
  },
  {
    d2 <- mortality_analytic_cons_included_base %>% mutate(mortality_group_mech = riskfactor_group_fn_CoNS_included(.))
    bind_rows(
      fit_riskfactor_model("CoNS-included sensitivity", d2, "Susceptible", "Resistant", "unweighted"),
      fit_riskfactor_model("CoNS-included sensitivity", d2, "Susceptible", "Resistant", "calibrated_weight"),
      fit_riskfactor_model("CoNS-included sensitivity", d2, "Uninfected", "Resistant", "unweighted")
    )
  }
)
write_xlsx(risk_factor_or_cons_sensitivity, "RiskFactors_CoNS_sensitivity.xlsx")

## --- CoNS-included full mechanism-by-mechanism replication of Table 4 ---
## Same 14 mechanisms, same dual population-wide/organism-restricted
## structure, same 3 comparisons (Susceptible unweighted+calibrated,
## Uninfected unweighted) as riskfactor_mechanisms/risk_factor_or_table
## above, but on the CoNS-included population (contaminated cultures that
## nonetheless reported a CoNS species counted as confirmed BSI, SECTION
## 2.6) instead of the primary group definition -- per instruction, to
## report in the Supplementary Material how the risk-factor estimates shift
## once these ~122 additional patients are included. Uninfected (matched_
## control_broad==1) is unchanged by the CoNS scenario -- only which BSI
## patients count as Resistant/Susceptible changes.
riskfactor_mechanisms_CoNS_included <- list(
  list(mechanism = "Any resistance", group_fn = function(d) {
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              d$BSI_positive_CoNS_included == 1 & d$resistant_BSI_positive_CoNS_included == 1 ~ "Resistant",
              d$BSI_positive_CoNS_included == 1 & d$resistant_BSI_positive_CoNS_included == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  list(mechanism = "MDR (bacterial isolates)", group_fn = function(d) {
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              d$BSI_positive_CoNS_included == 1 & d$mdr_BSI_bacterial_CoNS_included == 1 ~ "Resistant",
              d$BSI_positive_CoNS_included == 1 & d$mdr_BSI_bacterial_CoNS_included == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  list(mechanism = "XDR (bacterial isolates)", group_fn = function(d) {
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              d$BSI_positive_CoNS_included == 1 & d$xdr_BSI_bacterial_CoNS_included == 1 ~ "Resistant",
              d$BSI_positive_CoNS_included == 1 & d$xdr_BSI_bacterial_CoNS_included == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  list(mechanism = "Carbapenem-R (all confirmed BSI)", group_fn = function(d) {
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              d$BSI_positive_CoNS_included == 1 & d$carb_BSI_positive_CoNS_included == 1 ~ "Resistant",
              d$BSI_positive_CoNS_included == 1 & d$carb_BSI_positive_CoNS_included == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  list(mechanism = "Carbapenem-R (Gram-negative)", group_fn = function(d) {
    gneg <- is_cre_pa_ab(d$major_bacterial_group, d$pathogen1_clean)
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              d$BSI_positive_CoNS_included == 1 & gneg & d$carb_resistance_s1 == 1 ~ "Resistant",
              d$BSI_positive_CoNS_included == 1 & gneg & d$carb_resistance_s1 == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  list(mechanism = "Cephalosporin-R any-gen (all confirmed BSI)", group_fn = function(d) {
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              d$BSI_positive_CoNS_included == 1 & d$cephalosporin_BSI_positive_CoNS_included == 1 ~ "Resistant",
              d$BSI_positive_CoNS_included == 1 & d$cephalosporin_BSI_positive_CoNS_included == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  list(mechanism = "Cephalosporin-R any-gen (Gram-negative)", group_fn = function(d) {
    gneg <- is_cre_pa_ab(d$major_bacterial_group, d$pathogen1_clean)
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              d$BSI_positive_CoNS_included == 1 & gneg & d$cephalosporin_resistance_s1 == 1 ~ "Resistant",
              d$BSI_positive_CoNS_included == 1 & gneg & d$cephalosporin_resistance_s1 == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  list(mechanism = "3GCR (all confirmed BSI)", group_fn = function(d) {
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              d$BSI_positive_CoNS_included == 1 & d$c3GC_BSI_positive_CoNS_included == 1 ~ "Resistant",
              d$BSI_positive_CoNS_included == 1 & d$c3GC_BSI_positive_CoNS_included == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  list(mechanism = "3GCR (Gram-negative)", group_fn = function(d) {
    gneg <- is_cre_pa_ab(d$major_bacterial_group, d$pathogen1_clean)
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              d$BSI_positive_CoNS_included == 1 & gneg & d$c3gc_resistance_s1 == 1 ~ "Resistant",
              d$BSI_positive_CoNS_included == 1 & gneg & d$c3gc_resistance_s1 == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  list(mechanism = "Methicillin-R (all confirmed BSI)", group_fn = function(d) {
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              d$BSI_positive_CoNS_included == 1 & d$methicillin_BSI_positive_CoNS_included == 1 ~ "Resistant",
              d$BSI_positive_CoNS_included == 1 & d$methicillin_BSI_positive_CoNS_included == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  list(mechanism = "Methicillin-R / MRSA (S. aureus)", group_fn = function(d) {
    staph <- d$pathogen1_clean == "Staphylococcus aureus"
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              d$BSI_positive_CoNS_included == 1 & staph & d$methicillin_resistance_s1 == 1 ~ "Resistant",
              d$BSI_positive_CoNS_included == 1 & staph & d$methicillin_resistance_s1 == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  list(mechanism = "Access-resistant (AWaRe)", group_fn = function(d) {
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              d$BSI_positive_CoNS_included == 1 & d$access_BSI_positive_CoNS_included == 1 ~ "Resistant",
              d$BSI_positive_CoNS_included == 1 & d$access_BSI_positive_CoNS_included == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  list(mechanism = "Watch-resistant (AWaRe)", group_fn = function(d) {
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              d$BSI_positive_CoNS_included == 1 & d$watch_BSI_positive_CoNS_included == 1 ~ "Resistant",
              d$BSI_positive_CoNS_included == 1 & d$watch_BSI_positive_CoNS_included == 0 ~ "Susceptible", TRUE ~ NA_character_)
  }),
  list(mechanism = "Reserve-resistant (AWaRe)", group_fn = function(d) {
    case_when(d$matched_control_broad == 1 ~ "Uninfected",
              d$BSI_positive_CoNS_included == 1 & d$reserve_BSI_positive_CoNS_included == 1 ~ "Resistant",
              d$BSI_positive_CoNS_included == 1 & d$reserve_BSI_positive_CoNS_included == 0 ~ "Susceptible", TRUE ~ NA_character_)
  })
)

risk_factor_or_table_CoNS_included <- purrr::map_dfr(riskfactor_mechanisms_CoNS_included, function(mech) {
  d2 <- mortality_analytic_cons_included_base %>% mutate(mortality_group_mech = mech$group_fn(.))
  bind_rows(
    fit_riskfactor_model(mech$mechanism, d2, "Susceptible", "Resistant", "unweighted"),
    fit_riskfactor_model(mech$mechanism, d2, "Susceptible", "Resistant", "calibrated_weight"),
    fit_riskfactor_model(mech$mechanism, d2, "Uninfected", "Resistant", "unweighted")
  )
})
write_xlsx(risk_factor_or_table_CoNS_included, "RiskFactors_drugresistance_adjusted_OR_CoNS_included.xlsx")

## --- Vancomycin-R / Enterococcus: too few (VRE n=3) for any regression, descriptive only ---
vre_riskfactor_descriptive <- Main_data_CaseReport %>%
  filter(BSI_positive == 1, !is.na(genus1), genus1 == "Enterococcus", !is.na(vancomycin_BSI_positive)) %>%
  mutate(group = if_else(vancomycin_BSI_positive == 1, "VRE (vancomycin-resistant)", "VSE (vancomycin-susceptible)")) %>%
  count(group, name = "n") %>%
  mutate(note = "Descriptive only -- 3 VRE isolates in the whole survey, not enough for a regression model")
write_xlsx(vre_riskfactor_descriptive, "RiskFactors_VRE_descriptive_only.xlsx")


# =============================================================================
# SECTION 11: Length of hospital stay
# =============================================================================
# Corresponds to the methods spec's "6. Excess Length of Stay". Estimates the
# additional length of hospital stay associated with (1) any confirmed BSI
# and (2) drug-resistant BSI, each vs matched Uninfected controls, using the
# same mortality_analytic population/design as SECTION 9 (same 3 groups,
# same facility clustering) rather than rebuilding it.
#
# Primary outcome: LOS from enrolment to discharge/death, truncated at 28
# days -- this is EXACTLY mortality_analytic$time_inhosp (SECTION 9.2:
# pmin(pmax(LOS_days, 0.5), 28), floored at 0.5 rather than 0 so the log-link
# GLM below is defined for every patient), reused rather than recomputed.
#
# COVARIATES -- the spec's own candidate list for this section is age, sex,
# ward type, infection origin, comorbidity burden. Two of these are DROPPED,
# for reasons that follow directly from decisions already made in SECTION 9
# on your instruction ("same principles as...the mortality analysis"):
#   - ward_type_icu: same reverse-causation/collider concern as SECTION 9 --
#     if anything MORE acute here, since ward transfers happen DURING the
#     very stay whose length is the outcome being modelled.
#   - infection_origin_composite: structurally 100% missing for the
#     Uninfected control arm (SECTION 9.3's finding, unchanged) -- cannot be
#     a shared covariate in a model that compares BSI/Resistant-BSI AGAINST
#     Uninfected controls (unlike SECTION 10's risk-factor models, which
#     never compare against an arm requiring it).
# Adjustment set here is therefore: age, sex, comorbidity burden (Charlson).
# Hospital type (Level_of_care) was not in the spec's candidate list for
# this section and is not added.
#
# -----------------------------------------------------------------------
# 11.1: Outcome + group definitions
# -----------------------------------------------------------------------
los_analytic <- mortality_analytic %>%
  filter(!is.na(time_inhosp), !is.na(facility_id)) %>%
  mutate(
    bsi_vs_uninfected = factor(
      case_when(
        mortality_group == "Uninfected control" ~ "Uninfected control",
        mortality_group %in% c("Susceptible BSI", "Resistant BSI") ~ "Any BSI",
        TRUE ~ NA_character_
      ),
      levels = c("Uninfected control", "Any BSI")
    )
  )
los_comparisons <- list(
  list(group_var = "bsi_vs_uninfected", ref = "Uninfected control", comp = "Any BSI", label = "Any BSI vs Uninfected control"),
  list(group_var = "mortality_group", ref = "Uninfected control", comp = "Resistant BSI", label = "Resistant BSI vs Uninfected control")
)
los_covars <- c("age_years", "sex_at_birth", "charlson_score_combined")

# -----------------------------------------------------------------------
# 11.1b: CoNS-treatment sensitivity check (crude, unadjusted) -- same
# compact-check rationale and same mortality_group_CoNS_included/_excluded
# vectors as SECTION 9.2b, reused directly rather than rebuilt; the
# primary LOS analysis throughout this section always uses the unchanged
# mortality_analytic/los_analytic population above.
crude_los_by_group <- function(group_vec, label) {
  tibble(group = group_vec, time_inhosp = Main_data_Combined$time_inhosp) %>%
    filter(!is.na(group), !is.na(time_inhosp)) %>%
    group_by(group) %>%
    summarise(n = n(), mean_days = round(mean(time_inhosp), 2), median_days = median(time_inhosp), .groups = "drop") %>%
    mutate(population = label, .before = 1)
}
los_cons_sensitivity <- bind_rows(
  crude_los_by_group(Main_data_Combined$mortality_group, "Primary (unchanged)"),
  crude_los_by_group(mortality_group_CoNS_excluded, "CoNS-excluded sensitivity"),
  crude_los_by_group(mortality_group_CoNS_included, "CoNS-included sensitivity")
)
write_xlsx(los_cons_sensitivity, "LOS_CoNS_sensitivity_crude.xlsx")

# -----------------------------------------------------------------------
# 11.2: Adjusted GLM (Gamma, log link) -- coefficients, standardised
# adjusted mean difference, facility-cluster bootstrap CI
# -----------------------------------------------------------------------
# Gamma(link="log") is the standard choice for right-skewed positive LOS
# data (avoids the retransformation-bias correction a log-outcome + OLS/
# log-normal approach would need to get back to the natural day scale).
# The GLM coefficient itself is a log-mean-ratio; the SPEC asks for a mean
# DIFFERENCE in days, so the adjusted mean per arm is obtained by
# standardisation (g-computation: predict every patient's fitted LOS under
# each group value, holding their own covariates fixed, then average) --
# the same logic SECTION 9 used to turn a hazard-scale model into an
# absolute risk difference.
# fit_los_glm: weighted by combined_w = stabilised IPW (from a propensity
# model on `covars`, truncated at the 1st/99th percentile, same construction
# as SECTION 9.7's fit_ipw_cox_pooled) x calibrated_weight (SECTION 3.3's
# PPS design weight, calibrated to oblast population totals), renormalised
# to mean 1 -- identical construction to SECTION 9.13's mortality sensitivity
# model, but used here as the PRIMARY LOS model rather than a sensitivity
# check, since the underlying matched-cohort population and confounding
# structure are the same for LOS as for mortality. Rows with no resolvable
# calibrated_weight (facility not in the sampling frame) are dropped, not
# imputed, matching SECTION 9.13's treatment of the same issue.
fit_los_glm <- function(data, group_var, ref_level, comp_level, covars) {
  di <- data %>% dplyr::filter(.data[[group_var]] %in% c(ref_level, comp_level)) %>%
    mutate(group_bin = as.integer(.data[[group_var]] == comp_level)) %>%
    select(-any_of("calibrated_weight")) %>%
    left_join(facility_weights %>% select(UniqueID, calibrated_weight), by = c(facility_id = "UniqueID")) %>%
    filter(!is.na(calibrated_weight)) %>%
    # Complete-case on the propensity-model covariates -- required so the
    # propensity-model fit (glm()'s own na.action) and `di` stay the same
    # length; without this, a handful of NA age_years/sex_at_birth rows
    # (SECTION 9.3-level missingness, not imputed for the LOS analysis)
    # cause predict(ps_mod) to silently return fewer rows than di has.
    filter(if_all(all_of(covars), ~ !is.na(.)))
  form_ps <- as.formula(paste("group_bin ~", paste(covars, collapse = " + ")))
  ps_mod <- glm(form_ps, data = di, family = binomial())
  ps <- pmin(pmax(predict(ps_mod, type = "response"), 0.01), 0.99)
  p_marg <- mean(di$group_bin)
  ipw <- if_else(di$group_bin == 1, p_marg / ps, (1 - p_marg) / (1 - ps))
  q <- quantile(ipw, c(0.01, 0.99), na.rm = TRUE)
  ipw <- pmin(pmax(ipw, q[1]), q[2])
  combined_w <- ipw * di$calibrated_weight
  di$combined_w <- combined_w / mean(combined_w)
  form <- as.formula(paste("time_inhosp ~ group_bin +", paste(covars, collapse = " + ")))
  glm(form, data = di, family = Gamma(link = "log"), weights = combined_w)
}

standardised_mean_diff <- function(mod, di) {
  d0 <- di; d0$group_bin <- 0L
  d1 <- di; d1$group_bin <- 1L
  # na.rm = TRUE: predict() returns NA for any row with a missing covariate
  # (age_years/sex_at_birth have a few % missing, SECTION 9.3) -- those rows
  # were already excluded from FITTING mod (glm()'s default complete-case
  # na.action), so averaging over only the non-NA predictions here keeps the
  # standardisation population consistent with the fitting population.
  mean_ref <- mean(predict(mod, newdata = d0, type = "response"), na.rm = TRUE)
  mean_comp <- mean(predict(mod, newdata = d1, type = "response"), na.rm = TRUE)
  c(mean_ref = mean_ref, mean_comp = mean_comp, diff = mean_comp - mean_ref)
}

set.seed(44)
N_BOOT_LOS <- 200
bootstrap_los_diff <- function(data, group_var, ref_level, comp_level, covars, n_boot = N_BOOT_LOS) {
  facs <- unique(data$facility_id)
  out <- matrix(NA_real_, nrow = n_boot, ncol = 3, dimnames = list(NULL, c("mean_ref", "mean_comp", "diff")))
  for (b in seq_len(n_boot)) {
    boot_facs <- sample(facs, length(facs), replace = TRUE)
    boot_rows <- purrr::map_dfr(boot_facs, function(f) data[data$facility_id == f, , drop = FALSE])
    r <- tryCatch({
      di <- boot_rows %>% dplyr::filter(.data[[group_var]] %in% c(ref_level, comp_level)) %>%
        mutate(group_bin = as.integer(.data[[group_var]] == comp_level)) %>%
        select(-any_of("calibrated_weight")) %>%
        left_join(facility_weights %>% select(UniqueID, calibrated_weight), by = c(facility_id = "UniqueID")) %>%
        filter(!is.na(calibrated_weight)) %>%
        filter(if_all(all_of(covars), ~ !is.na(.)))
      form_ps <- as.formula(paste("group_bin ~", paste(covars, collapse = " + ")))
      ps_mod <- glm(form_ps, data = di, family = binomial())
      ps <- pmin(pmax(predict(ps_mod, type = "response"), 0.01), 0.99)
      p_marg <- mean(di$group_bin)
      ipw <- if_else(di$group_bin == 1, p_marg / ps, (1 - p_marg) / (1 - ps))
      q <- quantile(ipw, c(0.01, 0.99), na.rm = TRUE)
      ipw <- pmin(pmax(ipw, q[1]), q[2])
      combined_w <- ipw * di$calibrated_weight
      di$combined_w <- combined_w / mean(combined_w)
      mod <- glm(as.formula(paste("time_inhosp ~ group_bin +", paste(covars, collapse = " + "))),
                 data = di, family = Gamma(link = "log"), weights = combined_w)
      standardised_mean_diff(mod, di)
    }, error = function(e) c(mean_ref = NA_real_, mean_comp = NA_real_, diff = NA_real_))
    out[b, ] <- r
  }
  out
}

los_results <- purrr::map_dfr(los_comparisons, function(cmp) {
  mod <- fit_los_glm(los_analytic, cmp$group_var, cmp$ref, cmp$comp, los_covars)
  di <- los_analytic %>% dplyr::filter(.data[[cmp$group_var]] %in% c(cmp$ref, cmp$comp)) %>%
    mutate(group_bin = as.integer(.data[[cmp$group_var]] == cmp$comp))
  point <- standardised_mean_diff(mod, di)
  boot_mat <- bootstrap_los_diff(los_analytic, cmp$group_var, cmp$ref, cmp$comp, los_covars)
  boot_ok <- boot_mat[complete.cases(boot_mat), , drop = FALSE]
  tibble(
    comparison = cmp$label, n = nrow(di),
    adjusted_mean_los_ref_days = point["mean_ref"], adjusted_mean_los_comp_days = point["mean_comp"],
    mean_difference_days = point["diff"],
    ci_lower_days = quantile(boot_ok[, "diff"], 0.025, na.rm = TRUE),
    ci_upper_days = quantile(boot_ok[, "diff"], 0.975, na.rm = TRUE),
    n_bootstrap_ok = nrow(boot_ok)
  )
}) %>%
  mutate(across(where(is.numeric) & !c(n, n_bootstrap_ok), ~ round(.x, 3)))

# Cluster-robust coefficient table (rate-ratio scale, diagnostic/supplementary
# -- the mean-difference table above is the SPEC's requested primary output).
svy_los_coef <- purrr::map_dfr(los_comparisons, function(cmp) {
  di <- los_analytic %>% dplyr::filter(.data[[cmp$group_var]] %in% c(cmp$ref, cmp$comp)) %>%
    mutate(group_bin = as.integer(.data[[cmp$group_var]] == cmp$comp))
  des <- svydesign(ids = ~facility_id, data = di)
  form <- as.formula(paste("time_inhosp ~ group_bin +", paste(los_covars, collapse = " + ")))
  mod <- svyglm(form, design = des, family = Gamma(link = "log"))
  broom::tidy(mod, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::filter(term != "(Intercept)") %>%
    mutate(comparison = cmp$label, .before = 1) %>%
    rename(rate_ratio = estimate, ci_lower = conf.low, ci_upper = conf.high) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))
})

# -----------------------------------------------------------------------
# 11.3: Descriptive statistics -- unadjusted LOS by group, and total LOS
# from TRUE hospital admission (secondary measure)
# -----------------------------------------------------------------------
los_descriptive <- los_analytic %>%
  group_by(mortality_group) %>%
  summarise(n = n(), mean_days = mean(time_inhosp), sd_days = sd(time_inhosp),
            median_days = median(time_inhosp), q25 = quantile(time_inhosp, 0.25), q75 = quantile(time_inhosp, 0.75),
            .groups = "drop") %>%
  mutate(across(where(is.numeric) & !n, ~ round(.x, 2)))

# Secondary measure: total LOS from TRUE hospital admission (distinct from
# enrolment/blood-culture date) to discharge/death, NOT truncated at 28 days
# (winsorised at the 99th percentile instead, same convention SECTION 5 used
# for the enrolment-based LOS_days). The raw admission date is a genuinely
# separate source field from enrolment date ("Follow-up details: Date of
# admission (dd/mm/yyyy)", SECTION 5's own fallback field, dropped from
# Main_data_hospitalFollowup28_clean's final select() there since it was only
# needed to repair a handful of negative-LOS rows) -- re-extracted fresh from
# the raw sheet here, matching SECTION 9.2's convention for pulling
# additional raw Followup28 fields not carried into the _clean table.
raw_admission_dates <- tibble(
  patient_id = Main_data_hospitalFollowup28[["Unique survey patient ID #1"]],
  Date_of_admission_true = parse_dt(Main_data_hospitalFollowup28[["Follow-up details: Date of admission (dd/mm/yyyy)"]])
) %>% distinct(patient_id, .keep_all = TRUE)

los_from_admission_data <- los_analytic %>%
  left_join(raw_admission_dates, by = "patient_id") %>%
  mutate(los_from_admission = as.numeric(difftime(Date_of_discharge, Date_of_admission_true, units = "days"))) %>%
  mutate(los_from_admission = if_else(!is.na(los_from_admission) & los_from_admission >= 0, los_from_admission, NA_real_)) %>%
  mutate(los_from_admission = pmin(los_from_admission, quantile(los_from_admission, 0.99, na.rm = TRUE)))

los_from_admission_descriptive <- los_from_admission_data %>%
  filter(!is.na(los_from_admission)) %>%
  group_by(mortality_group) %>%
  summarise(n = n(), mean_days = mean(los_from_admission), sd_days = sd(los_from_admission),
            median_days = median(los_from_admission), q25 = quantile(los_from_admission, 0.25), q75 = quantile(los_from_admission, 0.75),
            .groups = "drop") %>%
  mutate(across(where(is.numeric) & !n, ~ round(.x, 2)))

write_xlsx(
  list(adjusted_mean_difference = los_results, cluster_robust_rate_ratios = svy_los_coef,
       descriptive_LOS_from_enrolment = los_descriptive, descriptive_LOS_from_admission_secondary = los_from_admission_descriptive),
  "LOS_excess_length_of_stay.xlsx"
)

## --- Lancet-style plot: adjusted mean LOS by group with 95% CI ---
# The reference (Uninfected) bar has no CI plotted -- the bootstrap CI
# (SECTION 11.2) is on the DIFFERENCE, not on either arm's own mean, so it is
# shown on the comparator bar by re-centring it around that arm's point
# estimate: comparator_mean +/- (bootstrap CI bound - point difference).
los_plot_data <- bind_rows(
  los_results %>% transmute(comparison, group = "Uninfected control", mean_days = adjusted_mean_los_ref_days,
                             ci_lower = NA_real_, ci_upper = NA_real_),
  los_results %>% transmute(comparison, group = sub(" vs .*$", "", comparison),
                             mean_days = adjusted_mean_los_comp_days,
                             ci_lower = adjusted_mean_los_comp_days - (mean_difference_days - ci_lower_days),
                             ci_upper = adjusted_mean_los_comp_days + (ci_upper_days - mean_difference_days))
)
los_plot <- ggplot(los_plot_data, aes(x = group, y = mean_days, fill = group)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.15, na.rm = TRUE) +
  facet_wrap(~comparison, scales = "free_x") +
  scale_fill_manual(values = c("Uninfected control" = "#17766E", "Any BSI" = "#A8541D", "Resistant BSI" = "#7A1F2B"), guide = "none") +
  labs(title = "Adjusted mean length of hospital stay (enrolment to discharge/death, truncated at 28 days)",
       subtitle = "95% CI shown on the comparator arm reflects the bootstrap CI on the adjusted mean difference",
       x = NULL, y = "Adjusted mean LOS (days)") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
        plot.subtitle = element_text(color = MORT_INK_SOFT))
ggsave("LOS_adjusted_mean_by_group.png", los_plot, width = 9, height = 5.5, dpi = 500)

# -----------------------------------------------------------------------
# 11.4: Full resistance-mechanism benchmark, each vs Susceptible AND vs
# Uninfected -- reuses SECTION 9.15's already-built mechanism-specific
# populations/comparisons (MDR, carbapenem, cephalosporin, 3GCR, MRSA)
# rather than rebuilding them, matching the same 3 pre-specified
# comparisons per mechanism used throughout SECTION 9 (Susceptible-
# equivalent vs Uninfected, Resistant-equivalent vs Uninfected, Resistant-
# equivalent vs Susceptible-equivalent).
# -----------------------------------------------------------------------
run_los_subgroup <- function(mechanism, analytic_data, comparisons, covars = los_covars, n_boot = N_BOOT_LOS, seed_boot = 45) {
  di_all <- analytic_data %>% filter(!is.na(time_inhosp), !is.na(facility_id))
  purrr::map_dfr(comparisons, function(cmp) {
    mod <- tryCatch(fit_los_glm(di_all, "mortality_group", cmp["ref"], cmp["comp"], covars), error = function(e) NULL)
    if (is.null(mod)) {
      return(tibble(mechanism = mechanism, comparison = paste(cmp["comp"], "vs", cmp["ref"]), n = NA_integer_,
                     adjusted_mean_los_ref_days = NA_real_, adjusted_mean_los_comp_days = NA_real_,
                     mean_difference_days = NA_real_, ci_lower_days = NA_real_, ci_upper_days = NA_real_, n_bootstrap_ok = 0L))
    }
    di <- di_all %>% dplyr::filter(mortality_group %in% c(cmp["ref"], cmp["comp"])) %>%
      mutate(group_bin = as.integer(mortality_group == cmp["comp"]))
    point <- standardised_mean_diff(mod, di)
    set.seed(seed_boot)
    boot_mat <- bootstrap_los_diff(di_all, "mortality_group", cmp["ref"], cmp["comp"], covars, n_boot = n_boot)
    boot_ok <- boot_mat[complete.cases(boot_mat), , drop = FALSE]
    tibble(
      mechanism = mechanism, comparison = paste(cmp["comp"], "vs", cmp["ref"]), n = nrow(di),
      adjusted_mean_los_ref_days = point["mean_ref"], adjusted_mean_los_comp_days = point["mean_comp"],
      mean_difference_days = point["diff"],
      ci_lower_days = if (nrow(boot_ok) > 0) quantile(boot_ok[, "diff"], 0.025, na.rm = TRUE) else NA_real_,
      ci_upper_days = if (nrow(boot_ok) > 0) quantile(boot_ok[, "diff"], 0.975, na.rm = TRUE) else NA_real_,
      n_bootstrap_ok = nrow(boot_ok)
    )
  }) %>% mutate(across(where(is.numeric) & !c(n, n_bootstrap_ok), ~ round(.x, 3)))
}

# MRSA needs its own LOS-specific population that KEEPS MSSA as a 3rd arm --
# SECTION 9.15.3's mortality_analytic_mrsa deliberately excludes MSSA
# entirely (0 deaths among 16-18 patients made the mortality GLM degenerate,
# a binary-outcome separation problem). LOS is a continuous outcome with no
# such degeneracy, so MSSA can be included here.
los_analytic_mrsa <- mortality_analytic %>%
  mutate(
    mortality_group = case_when(
      matched_control_broad == 1 ~ "Uninfected control",
      matched_control_broad == 0 & pathogen1_clean == "Staphylococcus aureus" & methicillin_BSI_positive == 1 ~ "MRSA BSI",
      matched_control_broad == 0 & pathogen1_clean == "Staphylococcus aureus" & methicillin_BSI_positive == 0 ~ "MSSA BSI",
      TRUE ~ NA_character_
    ),
    mortality_group = factor(mortality_group, levels = c("Uninfected control", "MSSA BSI", "MRSA BSI"))
  ) %>%
  filter(!is.na(mortality_group))
los_mrsa_comparisons <- list(
  c(ref = "Uninfected control", comp = "MSSA BSI"),
  c(ref = "Uninfected control", comp = "MRSA BSI"),
  c(ref = "MSSA BSI", comp = "MRSA BSI")
)

los_mechanism_list <- list(
  list(label = "Any resistance", data = mortality_analytic, comparisons = mortality_comparisons),
  list(label = "MDR (bacterial isolates)", data = mortality_analytic_mdr, comparisons = mdr_comparisons),
  list(label = "XDR (bacterial isolates)", data = mortality_analytic_xdr, comparisons = xdr_comparisons),
  list(label = "Carbapenem-R (Gram-negative)", data = mortality_analytic_carb, comparisons = carb_comparisons),
  list(label = "Cephalosporin-R any-gen (Gram-negative)", data = mortality_analytic_ceph, comparisons = ceph_comparisons),
  list(label = "3GCR (Gram-negative)", data = mortality_analytic_3gcr, comparisons = gcr3_comparisons),
  list(label = "Methicillin-R/MRSA (S. aureus)", data = los_analytic_mrsa, comparisons = los_mrsa_comparisons),
  list(label = "Access-resistant (AWaRe)", data = mortality_analytic_access, comparisons = access_comparisons),
  list(label = "Watch-resistant (AWaRe)", data = mortality_analytic_watch, comparisons = watch_comparisons),
  list(label = "Reserve-resistant (AWaRe)", data = mortality_analytic_reserve, comparisons = reserve_comparisons)
)
los_all_results <- purrr::map_dfr(los_mechanism_list, ~ run_los_subgroup(.x$label, .x$data, .x$comparisons))


# Vancomycin-R/Enterococcus: same n=3 VRE problem as SECTIONS 9.15.4/10 --
# too few for any adjusted model. Unadjusted descriptive LOS only.
los_vre_descriptive <- mortality_analytic %>%
  filter(!is.na(time_inhosp), matched_control == 0, !is.na(genus1), genus1 == "Enterococcus", !is.na(vancomycin_BSI_positive)) %>%
  mutate(group = if_else(vancomycin_BSI_positive == 1, "VRE (vancomycin-resistant)", "VSE (vancomycin-susceptible)")) %>%
  group_by(group) %>%
  summarise(n = n(), mean_days = mean(time_inhosp), median_days = median(time_inhosp),
            q25 = quantile(time_inhosp, 0.25), q75 = quantile(time_inhosp, 0.75), .groups = "drop") %>%
  mutate(across(where(is.numeric) & !n, ~ round(.x, 2)),
         note = "Descriptive only (unadjusted) -- insufficient n for a covariate-adjusted model")

write_xlsx(
  # Sheet names kept <=31 chars -- Excel's hard limit silently truncates
  # (and can collide) anything longer rather than erroring.
  list(mean_difference_primary = los_results, cluster_robust_RR_primary = svy_los_coef,
       descriptive_LOS_enrolment = los_descriptive, descriptive_LOS_admission_sec = los_from_admission_descriptive,
       all_mechanisms_mean_diff = los_all_results, VRE_descriptive_only = los_vre_descriptive),
  "LOS_excess_length_of_stay.xlsx"
)

## --- Forest plot: adjusted mean LOS difference, every mechanism x comparison ---
# Short, generic comparator labels by POSITION within each mechanism --
# every comparisons list in SECTION 9.15/11.4 is built in the same fixed
# order (lesser-resistance vs Uninfected, greater-resistance vs Uninfected,
# greater vs lesser), so position alone identifies the comparison type
# without needing the full (verbose, mechanism-specific) arm names.
los_forest_data <- los_all_results %>%
  group_by(mechanism) %>%
  mutate(comparison_short = factor(
    c("Susceptible-type vs Uninfected", "Resistant-type vs Uninfected", "Resistant-type vs Susceptible-type")[row_number()],
    levels = rev(c("Susceptible-type vs Uninfected", "Resistant-type vs Uninfected", "Resistant-type vs Susceptible-type"))
  )) %>%
  ungroup() %>%
  mutate(mechanism = factor(mechanism, levels = rev(unique(mechanism))))
los_forest_plot <- ggplot(los_forest_data, aes(x = mean_difference_days, y = comparison_short)) +
  geom_vline(xintercept = 0, color = MORT_INK_SOFT, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_lower_days, xmax = ci_upper_days), height = 0.2, color = MORT_INK_SOFT, na.rm = TRUE) +
  geom_point(size = 2.2, color = "#7A1F2B", na.rm = TRUE) +
  facet_grid(mechanism ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(title = "Adjusted excess length of stay by resistance mechanism",
       subtitle = "Every mechanism benchmarked vs Susceptible AND vs Uninfected.\nPoints = adjusted mean difference (days); bars = 95% bootstrap CI (facility-cluster).",
       x = "Adjusted mean difference in LOS (days)", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 0),
        strip.placement = "outside", plot.subtitle = element_text(color = MORT_INK_SOFT))
ggsave("LOS_forest_all_mechanisms.png", los_forest_plot, width = 10, height = 8.5, dpi = 500)

# Larger-font, higher-resolution companion (same data/model, IPW x calibrated
# design weight combined -- SECTION 11's primary LOS weighting scheme, per
# fit_los_glm() above -- font size only, no other change).
los_forest_plot_big <- ggplot(los_forest_data, aes(x = mean_difference_days, y = comparison_short)) +
  geom_vline(xintercept = 0, color = MORT_INK_SOFT, linewidth = 0.5) +
  geom_errorbarh(aes(xmin = ci_lower_days, xmax = ci_upper_days), height = 0.25, color = MORT_INK_SOFT, linewidth = 0.6, na.rm = TRUE) +
  geom_point(size = 3.2, color = "#7A1F2B", na.rm = TRUE) +
  facet_grid(mechanism ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(title = "Adjusted excess length of stay by resistance mechanism",
       subtitle = "Every mechanism benchmarked vs Susceptible AND vs Uninfected.\nPoints = adjusted mean difference (days); bars = 95% bootstrap CI (facility-cluster).",
       x = "Adjusted mean difference in LOS (days)", y = NULL) +
  theme_minimal(base_size = 17) +
  theme(panel.grid.minor = element_blank(), strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 0, size = 15),
        strip.placement = "outside", plot.subtitle = element_text(color = MORT_INK_SOFT, size = 13),
        plot.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14), axis.title = element_text(size = 15))
ggsave("LOS_forest_all_mechanisms_bigfont.png", los_forest_plot_big, width = 14, height = 13, dpi = 500)
ggsave("LOS_forest_all_mechanisms_bigfont.tiff", los_forest_plot_big, width = 14, height = 13, dpi = 500, compression = "lzw")

# -----------------------------------------------------------------------
# 11.5: Overall + subgroup panels (ward type, age group) -- adjusted mean
# LOS per arm with 95% CI, for the 3 PRIMARY (any-resistance) comparisons
# -----------------------------------------------------------------------
# fit_los_arm_means(): pairwise-fits Uninfected vs each other level present
# (droplevels-safe, mirrors SECTION 9.12's make_adjusted_curve_data()) and
# returns EVERY arm's own adjusted mean + bootstrap CI (not just the
# difference) by pulling both mean_ref and mean_comp out of each pairwise
# bootstrap -- the Uninfected arm's CI only needs computing once even though
# it appears as `ref` in multiple pairwise fits.
fit_los_arm_means <- function(data, covars, n_boot = 100) {
  data <- data %>% mutate(mortality_group = droplevels(mortality_group)) %>%
    filter(!is.na(time_inhosp), !is.na(facility_id))
  lv <- levels(data$mortality_group)
  if (length(lv) < 2) return(NULL)
  ref_lv <- if ("Uninfected control" %in% lv) "Uninfected control" else lv[1]
  other_lv <- setdiff(lv, ref_lv)
  arm_stats <- list()
  for (ol in other_lv) {
    mod <- tryCatch(fit_los_glm(data, "mortality_group", ref_lv, ol, covars), error = function(e) NULL)
    if (is.null(mod)) next
    di <- data %>% dplyr::filter(mortality_group %in% c(ref_lv, ol)) %>%
      mutate(group_bin = as.integer(mortality_group == ol))
    point <- tryCatch(standardised_mean_diff(mod, di), error = function(e) NULL)
    if (is.null(point)) next
    boot_mat <- bootstrap_los_diff(data, "mortality_group", ref_lv, ol, covars, n_boot = n_boot)
    boot_ok <- boot_mat[complete.cases(boot_mat), , drop = FALSE]
    if (nrow(boot_ok) == 0) next
    if (is.null(arm_stats[[ref_lv]])) {
      arm_stats[[ref_lv]] <- tibble(group = ref_lv, mean_days = point["mean_ref"],
                                     ci_lower = quantile(boot_ok[, "mean_ref"], 0.025), ci_upper = quantile(boot_ok[, "mean_ref"], 0.975))
    }
    arm_stats[[ol]] <- tibble(group = ol, mean_days = point["mean_comp"],
                               ci_lower = quantile(boot_ok[, "mean_comp"], 0.025), ci_upper = quantile(boot_ok[, "mean_comp"], 0.975))
  }
  if (length(arm_stats) == 0) return(NULL)
  bind_rows(arm_stats)
}

los_subgroup_colors <- c("Uninfected control" = "#17766E", "Susceptible BSI" = "#A8541D", "Resistant BSI" = "#7A1F2B")

los_by_ward <- purrr::map_dfr(c("ICU", "Non-ICU"), function(lvl) {
  d <- fit_los_arm_means(los_analytic %>% dplyr::filter(ward_type_icu == lvl), los_covars)
  if (!is.null(d)) d %>% mutate(panel = lvl)
})
plot_los_ward <- ggplot(los_by_ward, aes(x = group, y = mean_days, fill = group)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.15) +
  facet_wrap(~panel) +
  scale_fill_manual(values = los_subgroup_colors, guide = "none") +
  labs(title = "Adjusted mean LOS by ward type (display stratification, not a causal covariate)",
       x = NULL, y = "Adjusted mean LOS (days)") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"), axis.text.x = element_text(angle = 20, hjust = 1))
ggsave("LOS_by_wardtype.png", plot_los_ward, width = 9, height = 5.5, dpi = 500)

los_analytic <- los_analytic %>%
  mutate(age_group_derived = case_when(age_years <= 28 / 365.25 ~ "Neonate", age_years < 18 ~ "Paediatric", TRUE ~ "Adult"))
los_by_age <- purrr::map_dfr(age_group_levels, function(lvl) {
  d <- fit_los_arm_means(los_analytic %>% dplyr::filter(age_group_derived == lvl), los_covars)
  if (!is.null(d)) d %>% mutate(panel = lvl)
})
plot_los_age <- ggplot(los_by_age, aes(x = group, y = mean_days, fill = group)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.15) +
  facet_wrap(~ factor(panel, levels = age_group_levels)) +
  scale_fill_manual(values = los_subgroup_colors, guide = "none") +
  labs(title = "Adjusted mean LOS by age group", x = NULL, y = "Adjusted mean LOS (days)") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"), axis.text.x = element_text(angle = 20, hjust = 1))
ggsave("LOS_by_agegroup.png", plot_los_age, width = 11, height = 5.5, dpi = 500)


# =============================================================================
# SECTION 12: Diagnostic stewardship and assessment of measurement bias
# =============================================================================
# Corresponds to the methods spec's "7. Diagnostic Stewardship and
# Assessment of Measurement Bias". Population = Main_data_CaseReport (all
# 9247 enrolled suspected-BSI patients) -- the original PPS survey sample,
# same design (facility_id_clean cluster, design_weight/calibrated_weight)
# as SECTIONS 3.2-3.3/10.
#
# Several needed fields were dropped from Main_data_CaseReport in SECTION 1's
# own column selection and are re-extracted here from a raw re-read of the
# CaseReport sheet, joined by patient_id -- the same pattern SECTIONS
# 9.2/11.3 already use for raw fields not preserved in a _clean table.
#
# "Prior antibiotic exposure ... within 48 hours prior to blood culture
# collection" (the spec's own definition) is built EXACTLY per the source
# form (Kyrgyzstan Master Protocol VS2, Annex 13.a, Section 4 "Antimicrobial
# use history"): 1 if abx_current_therapy==1 (on antibiotics AT the time of
# sampling), OR abx_current_therapy==0 AND the form's own follow-up question
# ("time between blood sampling and the last antibiotic intake") = "<=48
# hours" (option OP1 of 5: OP1 <=48h, OP2 49h-7d, OP3 8-14d, OP4 >14d, OP5
# Unknown -- option text confirmed from the protocol PDF, not guessed).
diag_raw_cr <- read_excel(DHIS2_XLSX, sheet = "CaseReport", guess_max = 1e6) %>%
  transmute(
    patient_id = `Unique survey patient ID`,
    case_def_adult = `Adult patient (≥ 18 years): Case definition met with ≥ 2 of the following signs (qSOFA)`,
    case_def_paed = `Paediatric patient (≥ 29 days and <18 years): Case definition met with ≥ 1 of the following signs (LqSOFA)`,
    case_def_neo = `Neonate patient (≤ 28 days): Case definition met with ≥ 2 of the following signs`,
    ind_abx_adult = `Adult patient (≥ 18 years): Indication for parenteral antibiotics`,
    ind_abx_paed = `Paediatric patient (≥ 29 days and <18 years): Indication for parenteral antibiotics`,
    ind_abx_neo = `Neonate patient (≤ 28 days): Indication for parenteral antibiotics`,
    abx_period_since_last = `Antimicrobial use history: If no, indicate the period between the sampling and the last antibiotic intake`
  ) %>%
  mutate(
    case_definition_met = as.integer(!is.na(coalesce(case_def_adult, case_def_paed, case_def_neo))),
    indication_parenteral_abx = as.integer(!is.na(coalesce(ind_abx_adult, ind_abx_paed, ind_abx_neo))),
    enrolment_pathway = case_when(
      case_definition_met == 1 ~ "Formal case definition met",
      indication_parenteral_abx == 1 ~ "Clinician judgement only (parenteral abx indicated)",
      TRUE ~ "Neither documented"
    )
    # NOTE: abx_history is NOT re-derived here -- Main_data_CaseReport
    # already carries it (SECTION 1) and diagnostic_base inherits it
    # directly below. Re-deriving it here under the same name previously
    # caused a silent dplyr join-collision (abx_history.x/.y), which is why
    # the model below referencing plain `abx_history` failed.
  ) %>%
  select(patient_id, case_definition_met, indication_parenteral_abx, enrolment_pathway, abx_period_since_last) %>%
  distinct(patient_id, .keep_all = TRUE)

diagnostic_base <- Main_data_CaseReport %>%
  filter(!is.na(calibrated_weight)) %>%
  left_join(diag_raw_cr, by = "patient_id") %>%
  left_join(Hospital_characteristics %>% select(`Hospital ID`, Level_of_care), by = c(facility_id_clean = "Hospital ID")) %>%
  mutate(
    # Unknown/unrecorded timing among patients not currently on antibiotics
    # (period = "Unknown" or blank) is coded 0 (not recently exposed), not
    # NA -- per explicit instruction, a patient who wasn't asked or couldn't
    # recall the timing is treated as lacking the recent exposure rather than
    # dropped from the analysis. abx_current_therapy is never NA (see its own
    # definition, SECTION 1), so the only remaining ambiguity this resolves
    # is the timing question.
    abx_exposure_48h = case_when(
      abx_current_therapy == 1 ~ 1L,
      abx_current_therapy == 0 & abx_period_since_last == "OP1" ~ 1L,
      abx_current_therapy == 0 & abx_period_since_last %in% c("OP2", "OP3", "OP4") ~ 0L,
      TRUE ~ 0L
    )
  )
diag_design <- svydesign(ids = ~facility_id_clean, weights = ~calibrated_weight, data = diagnostic_base)

# -----------------------------------------------------------------------
# 12.1: Blood culture utilisation + diagnostic yield
# -----------------------------------------------------------------------
# Testing/utilisation rate reuses SECTION 4's incidence_row() unchanged --
# the "indicator" here is simply "was a blood culture taken" (=1 for every
# enrolled CaseReport patient, by construction of the survey), so the same
# facility-level Poisson-rate machinery gives the blood-culture testing rate
# per 1,000 patient-days/admissions. Diagnostic yield (% positive among
# tested) is SECTION 3.2's prevalence_overall_calibrated -- referenced, not
# recomputed.
bc_utilisation_table <- bind_rows(
  incidence_row(rep(1L, nrow(Main_data_CaseReport)), "Total_patientDays", "Blood culture testing rate"),
  incidence_row(rep(1L, nrow(Main_data_CaseReport)), "Total_inpatient_admissions", "Blood culture testing rate")
) %>% mutate(across(where(is.numeric) & !c(n_facilities, n_cases, total_exposure), ~ round(.x, 4)))

diagnostic_yield_summary <- tibble(
  measure = "Diagnostic yield (BSI-positive among tested)",
  n_tested = nrow(Main_data_CaseReport), n_positive = sum(Main_data_CaseReport$BSI_positive == 1, na.rm = TRUE),
  unweighted_pct = 100 * mean(Main_data_CaseReport$BSI_positive == 1, na.rm = TRUE),
  calibrated_pct = prevalence_overall_calibrated * 100,
  ci_lower = prevalence_overall_calibrated_ci[1] * 100, ci_upper = prevalence_overall_calibrated_ci[2] * 100
)

# -----------------------------------------------------------------------
# 12.2: Enrolment pathway (formal case definition vs clinician judgement)
# -----------------------------------------------------------------------
enrolment_pathway_summary <- diagnostic_base %>%
  count(enrolment_pathway) %>%
  mutate(pct = round(100 * n / sum(n), 1))
# Positivity only here (CaseReport-scoped); the mortality comparison across
# enrolment pathways needs Main_data_Combined's outcome linkage and is built
# in SECTION 13 instead, alongside its own data-availability caveat.
enrolment_pathway_positivity <- diagnostic_base %>%
  filter(!is.na(enrolment_pathway)) %>%
  group_by(enrolment_pathway) %>%
  summarise(n = n(), n_bsi_positive = sum(BSI_positive == 1, na.rm = TRUE),
            unweighted_positivity_pct = round(100 * mean(BSI_positive == 1, na.rm = TRUE), 1), .groups = "drop")

# -----------------------------------------------------------------------
# 12.3: Prior antibiotic exposure (48h) -- prevalence + 2 associations
# -----------------------------------------------------------------------
abx_exposure_summary <- tibble(
  measure = "Antibiotic exposure within 48h prior to blood culture",
  n_determinate = sum(!is.na(diagnostic_base$abx_exposure_48h)),
  n_exposed = sum(diagnostic_base$abx_exposure_48h == 1, na.rm = TRUE),
  unweighted_pct = 100 * mean(diagnostic_base$abx_exposure_48h == 1, na.rm = TRUE),
  calibrated_pct = svyciprop(~I(abx_exposure_48h == 1), subset(diag_design, !is.na(abx_exposure_48h)), method = "logit") * 100
)

# Association 1: prior abx (48h) vs BSI positivity, ALL enrolled patients.
# Covariates per the spec: age, grouped clinical specialty, hospital type,
# infection origin, available severity indicator (Pitt score).
diag_covars_all <- c("age_years", "specialty_group", "Level_of_care", "infection_origin_composite", "pitt_score_partial")
form_abx_vs_positivity <- as.formula(paste("I(BSI_positive == 1) ~ abx_exposure_48h +", paste(diag_covars_all, collapse = " + ")))
mod_abx_vs_positivity <- svyglm(form_abx_vs_positivity, design = subset(diag_design, !is.na(abx_exposure_48h)), family = quasibinomial())
abx_vs_positivity_or <- broom::tidy(mod_abx_vs_positivity, conf.int = TRUE, exponentiate = TRUE) %>%
  dplyr::filter(term != "(Intercept)") %>%
  rename(odds_ratio = estimate, ci_lower = conf.low, ci_upper = conf.high) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

# abx_history vs abx_current_therapy ("current_abx"): same model, swapping
# the exposure term, to decide which prior-antibiotic-use covariate is the
# better one to carry forward. abx_history = 1 if a historical drug name was
# recorded (joined in via diag_raw_cr above); abx_current_therapy = 1 if on
# antibiotics at the moment of sampling. Both are complete (no NA), so both
# fit on the same n as the primary model above.
form_abxhist_vs_positivity <- as.formula(paste("I(BSI_positive == 1) ~ abx_history +", paste(diag_covars_all, collapse = " + ")))
mod_abxhist_vs_positivity <- svyglm(form_abxhist_vs_positivity, design = diag_design, family = quasibinomial())
abxhist_vs_positivity_or <- broom::tidy(mod_abxhist_vs_positivity, conf.int = TRUE, exponentiate = TRUE) %>%
  dplyr::filter(term != "(Intercept)") %>%
  rename(odds_ratio = estimate, ci_lower = conf.low, ci_upper = conf.high) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

form_current_vs_positivity <- as.formula(paste("I(BSI_positive == 1) ~ abx_current_therapy +", paste(diag_covars_all, collapse = " + ")))
mod_current_vs_positivity <- svyglm(form_current_vs_positivity, design = diag_design, family = quasibinomial())
current_vs_positivity_or <- broom::tidy(mod_current_vs_positivity, conf.int = TRUE, exponentiate = TRUE) %>%
  dplyr::filter(term != "(Intercept)") %>%
  rename(odds_ratio = estimate, ci_lower = conf.low, ci_upper = conf.high) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

abx_variable_comparison <- tibble(
  variable = c("abx_history (historical drug name recorded)", "abx_current_therapy (on antibiotics at sampling)"),
  n_exposed = c(sum(diagnostic_base$abx_history == 1, na.rm = TRUE), sum(diagnostic_base$abx_current_therapy == 1, na.rm = TRUE)),
  pct_exposed = round(100 * c(mean(diagnostic_base$abx_history == 1, na.rm = TRUE), mean(diagnostic_base$abx_current_therapy == 1, na.rm = TRUE)), 1),
  aic = c(AIC(mod_abxhist_vs_positivity)[1], AIC(mod_current_vs_positivity)[1]),
  exposure_or = c(abxhist_vs_positivity_or$odds_ratio[1], current_vs_positivity_or$odds_ratio[1]),
  exposure_p = c(abxhist_vs_positivity_or$p.value[1], current_vs_positivity_or$p.value[1])
)
cat("=== abx_history vs abx_current_therapy: which covariate to use ===\n")
print(as.data.frame(abx_variable_comparison))

# CoNS-treatment sensitivity check (genuine re-fit, headline model only) --
# same cheap-svyglm() rationale as SECTION 10's risk-factor sensitivity
# check; the primary model above always uses the unchanged BSI_positive.
fit_abx_vs_positivity_sensitivity <- function(bsi_col, label) {
  form <- as.formula(paste0("I(", bsi_col, " == 1) ~ abx_exposure_48h +", paste(diag_covars_all, collapse = " + ")))
  mod <- tryCatch(svyglm(form, design = subset(diag_design, !is.na(abx_exposure_48h)), family = quasibinomial()), error = function(e) NULL)
  if (is.null(mod)) return(tibble(population = label, term = NA_character_, odds_ratio = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_))
  broom::tidy(mod, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::filter(term != "(Intercept)") %>%
    rename(odds_ratio = estimate, ci_lower = conf.low, ci_upper = conf.high) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)), population = label, .before = 1)
}
abx_vs_positivity_cons_sensitivity <- bind_rows(
  fit_abx_vs_positivity_sensitivity("BSI_positive_CoNS_sensitive", "CoNS-excluded sensitivity"),
  fit_abx_vs_positivity_sensitivity("BSI_positive_CoNS_included", "CoNS-included sensitivity")
)
write_xlsx(abx_vs_positivity_cons_sensitivity, "Diagnostic_stewardship_CoNS_sensitivity.xlsx")

# Association 2: prior abx (48h) vs drug-resistant BSI, CONFIRMED BSI only.
# NOT MODELLED: abx_exposure_48h==0 & resistant_BSI_positive==0 has only 1
# patient (checked directly -- 2x2 cross-tab: {0,0}=1, {0,1}=14, {1,0}=19,
# {1,1}=286). A near-empty cell like this makes even an UNADJUSTED odds
# ratio unstable (Haldane-Anscombe territory), let alone a 4-covariate
# adjusted model -- attempting the adjusted model produces exactly the kind
# of degenerate estimate (OR collapsing to 0 with a huge SE) already seen
# and rejected earlier in this script for other near-empty cells (VRE,
# MSSA). Reported as a crude cross-tab only, consistent with that
# precedent, not as a modelled association.
diag_bsi_only <- diagnostic_base %>% filter(BSI_positive == 1, !is.na(resistant_BSI_positive), !is.na(abx_exposure_48h))
abx_vs_resistance_crosstab <- diag_bsi_only %>%
  count(abx_exposure_48h, resistant_BSI_positive) %>%
  tidyr::pivot_wider(names_from = resistant_BSI_positive, values_from = n, values_fill = 0,
                      names_prefix = "resistant_BSI_positive_") %>%
  mutate(note = "Descriptive only -- {abx_exposure_48h=0, resistant=0} cell has n=1, too sparse for any regression (crude or adjusted)")

# -----------------------------------------------------------------------
# 12.4: Laboratory turnaround times + association with positivity
# -----------------------------------------------------------------------
# SampleShipment (facility-level processing chain, one row per isolate slot,
# not currently loaded elsewhere in this script) + CentralRefLab's own
# result-communication fields (already loaded as Main_data_hospitalLab).
sample_shipment_raw <- read_excel(DHIS2_XLSX, sheet = "SampleShipment", guess_max = 1e6) %>%
  transmute(
    patient_id = `Unique survey patient ID`,
    collection_dt = parse_dt(`Collection point 1: Sample collection date, time (24h) - S1`),
    referral_dt = parse_dt(`Collection point 1: Referral/shipment date, time (24h) - S1`),
    incub_arrival_dt = parse_dt(`Incubation point without gram stain test 1: Sample arrival date, time (24h) - S1`),
    incub_start_dt = parse_dt(`Incubation point without gram stain test 1: Incubation start date, time (24h) - S1`),
    incub_end_dt = parse_dt(`Incubation point without gram stain test 1: Incubation end date, time (24h) - S1`)
  ) %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  mutate(
    hrs_collection_to_incubation_start = as.numeric(difftime(incub_start_dt, collection_dt, units = "hours")),
    hrs_incubation_duration = as.numeric(difftime(incub_end_dt, incub_start_dt, units = "hours"))
  ) %>%
  mutate(across(starts_with("hrs_"), ~ if_else(.x >= 0 & .x <= 24 * 30, .x, NA_real_)))  # drop implausible negatives/multi-month gaps

lab_result_timing_raw <- Main_data_hospitalLab %>%
  transmute(
    patient_id = `Unique survey patient ID`,
    result_date = parse_dt(`ID outcome: Date of test result`),
    final_idast_date = parse_dt(`Responsible officer: Date final ID/AST results`),
    sent_to_clinical_date = parse_dt(`Responsible officer: Date final ID/AST results sent to clinical`)
  ) %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  mutate(hrs_result_to_communicated = as.numeric(difftime(sent_to_clinical_date, final_idast_date, units = "hours"))) %>%
  mutate(hrs_result_to_communicated = if_else(hrs_result_to_communicated >= 0 & hrs_result_to_communicated <= 24 * 30, hrs_result_to_communicated, NA_real_))

diagnostic_base <- diagnostic_base %>%
  left_join(sample_shipment_raw, by = "patient_id") %>%
  left_join(lab_result_timing_raw, by = "patient_id") %>%
  # hrs_incubation_to_final_idast: incubation END to final ID/AST result --
  # the ToR's own "time from incubation to pathogen identification and AST",
  # distinct from hrs_incubation_duration (the incubation START-to-END
  # period itself, i.e. how long the isolate sat incubating, not how long
  # identification/AST took afterwards). Previously conflated: only
  # incubation duration was computed, not this post-incubation ID/AST step.
  mutate(
    hrs_incubation_to_final_idast = as.numeric(difftime(final_idast_date, incub_end_dt, units = "hours")),
    hrs_incubation_to_final_idast = if_else(hrs_incubation_to_final_idast >= 0 & hrs_incubation_to_final_idast <= 24 * 30, hrs_incubation_to_final_idast, NA_real_)
  )

lab_turnaround_descriptive <- purrr::map_dfr(
  c("hrs_collection_to_incubation_start", "hrs_incubation_duration", "hrs_incubation_to_final_idast", "hrs_result_to_communicated"),
  function(v) {
    x <- diagnostic_base[[v]]
    tibble(measure = v, n = sum(!is.na(x)), median_hours = median(x, na.rm = TRUE),
           q25_hours = quantile(x, 0.25, na.rm = TRUE), q75_hours = quantile(x, 0.75, na.rm = TRUE))
  }
) %>% mutate(across(where(is.numeric) & !c(n), ~ round(.x, 1)))

# Association: laboratory-processing delay vs BSI positivity, run for each
# of the ToR's three named turnaround intervals (collection-to-incubation
# start; incubation duration; incubation-to-final-ID/AST) -- previously only
# the first of these three was tested. hrs_result_to_communicated is
# deliberately NOT tested here: that interval is only observed for already-
# positive cultures (communication of a result presupposes a result), so
# regressing it against positive/negative status would be circular. Same
# adjustment covariates as 12.3's Association 1.
turnaround_vs_positivity_or <- purrr::map_dfr(
  c("hrs_collection_to_incubation_start", "hrs_incubation_duration", "hrs_incubation_to_final_idast"),
  function(v) {
    di <- diagnostic_base %>% dplyr::filter(!is.na(.data[[v]]))
    des <- svydesign(ids = ~facility_id_clean, weights = ~calibrated_weight, data = di)
    form <- as.formula(paste("I(BSI_positive == 1) ~", v, "+", paste(diag_covars_all, collapse = " + ")))
    mod <- tryCatch(svyglm(form, design = des, family = quasibinomial()), error = function(e) NULL)
    if (is.null(mod)) return(tibble(turnaround_measure = v, term = v, odds_ratio = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, p.value = NA_real_, n = nrow(di)))
    broom::tidy(mod, conf.int = TRUE, exponentiate = TRUE) %>%
      dplyr::filter(term == v) %>%
      rename(odds_ratio = estimate, ci_lower = conf.low, ci_upper = conf.high) %>%
      mutate(turnaround_measure = v, n = nrow(di), .before = 1) %>%
      mutate(across(where(is.numeric) & !n, ~ round(.x, 4)))
  }
)

# -----------------------------------------------------------------------
# 12.5b: Blood culture utilisation and diagnostic yield, by oblast (region)
# -----------------------------------------------------------------------
# Complements 12.5's hospital-type breakdown -- the ToR asks for variation
# "across hospitals AND regions"; Table 1 already reports weighted
# prevalence by oblast, but not testing/utilisation counts specifically.
diag_variation_by_oblast <- diagnostic_base %>%
  filter(!is.na(NameAU)) %>%
  group_by(NameAU) %>%
  summarise(n_enrolled = n(), n_bsi_positive = sum(BSI_positive == 1, na.rm = TRUE),
            unweighted_positivity_pct = round(100 * mean(BSI_positive == 1, na.rm = TRUE), 1),
            n_facilities = n_distinct(facility_id_clean), .groups = "drop") %>%
  arrange(desc(n_enrolled))

# -----------------------------------------------------------------------
# 12.5: Variation in testing rate and yield across hospital type
# -----------------------------------------------------------------------
# Facility-level testing-rate/yield variation is summarised by hospital type
# (Level_of_care) rather than per-individual-facility (119 facilities would
# make an individual-facility table unwieldy and mostly single-patient-driven
# noise) -- bounded scope, consistent with SECTION 3's own subgroup pattern.
diag_variation_by_hospitaltype <- diagnostic_base %>%
  filter(!is.na(Level_of_care)) %>%
  group_by(Level_of_care) %>%
  summarise(n_enrolled = n(), n_bsi_positive = sum(BSI_positive == 1, na.rm = TRUE),
            unweighted_positivity_pct = round(100 * mean(BSI_positive == 1, na.rm = TRUE), 1),
            n_facilities = n_distinct(facility_id_clean), .groups = "drop")

write_xlsx(
  list(bc_utilisation = bc_utilisation_table, diagnostic_yield = diagnostic_yield_summary,
       enrolment_pathway_counts = enrolment_pathway_summary, enrolment_pathway_positivity = enrolment_pathway_positivity,
       abx_exposure_48h_summary = abx_exposure_summary, abx_vs_BSI_positivity_OR = abx_vs_positivity_or,
       abx_history_vs_current_comparison = abx_variable_comparison,
       abx_vs_resistance_crosstab = abx_vs_resistance_crosstab, lab_turnaround_descriptive = lab_turnaround_descriptive,
       turnaround_vs_positivity_OR = turnaround_vs_positivity_or, variation_by_hospital_type = diag_variation_by_hospitaltype,
       variation_by_oblast = diag_variation_by_oblast),
  "Diagnostic_stewardship_measurement_bias.xlsx"
)


# =============================================================================
# SECTION 13: Value of collecting CRF data for culture-negative/contaminated cases
# =============================================================================
# Corresponds to the methods spec's "8. Value of Collecting CRF Data for
# Blood-Culture-Negative and Contaminated Cases".
#
# -----------------------------------------------------------------------
# 13.1: Microbiological outcome classification (all 9247 enrolled patients)
# -----------------------------------------------------------------------
# CentralRefLab's "ID outcome: Conclusion" already distinguishes
# CONTAMINATION from a true negative -- SECTION 2's existing BSI_positive/
# lab_match flags only ever check for =="MICROBIOLOGICALLY" and fold
# everything else (contamination, unconclusive, no CentralRefLab referral at
# all) into "0", which is exactly the distinction this section needs to
# undo. Main_data_CaseReport already carries id_outcome_conclusion (added
# alongside has_CoNS_species/has_CoNS_species_if_contaminated, SECTION 2.6 --
# same source field, same join key), so no fresh re-join is needed here;
# a redundant one would silently get suffixed to id_outcome_conclusion.x/.y
# by left_join() and break the bare-name case_when() below.
crf_value_base <- Main_data_CaseReport %>%
  mutate(
    microbio_outcome = factor(case_when(
      BSI_positive == 1 ~ "Confirmed BSI",
      TRUE ~ "No confirmed BSI"
    ), levels = c("Confirmed BSI", "No confirmed BSI"))
  )
microbio_outcome_counts <- crf_value_base %>%
  count(microbio_outcome) %>%
  mutate(pct = round(100 * n / sum(n), 1))

# -----------------------------------------------------------------------
# 13.1b: CoNS-treatment sensitivity check on the microbiological-outcome
# classification itself -- directly relevant here, since this section's own
# subject is exactly the culture-negative/contaminated boundary that CoNS
# treatment (SECTION 2.6) moves patients across. Two alternate
# classifications, each moving patients in the OPPOSITE direction from the
# unchanged primary above: CoNS-excluded moves some "Confirmed BSI" patients
# into "Contaminated culture"; CoNS-included moves some "Contaminated
# culture" patients (with a CoNS species reported) into "Confirmed BSI".
microbio_outcome_CoNS_excluded <- factor(case_when(
  crf_value_base$BSI_positive_CoNS_sensitive == 1 ~ "Confirmed BSI",
  crf_value_base$BSI_positive == 1 & crf_value_base$BSI_positive_CoNS_sensitive == 0 ~ "Contaminated culture",
  crf_value_base$id_outcome_conclusion == "CONTAMINATION" ~ "Contaminated culture",
  crf_value_base$id_outcome_conclusion == "UNCONCLUSIVE" ~ "Unconclusive",
  TRUE ~ "Culture-negative"
), levels = c("Confirmed BSI", "Contaminated culture", "Culture-negative", "Unconclusive"))
microbio_outcome_CoNS_included <- factor(case_when(
  crf_value_base$BSI_positive_CoNS_included == 1 ~ "Confirmed BSI",
  crf_value_base$id_outcome_conclusion == "CONTAMINATION" ~ "Contaminated culture",
  crf_value_base$id_outcome_conclusion == "UNCONCLUSIVE" ~ "Unconclusive",
  TRUE ~ "Culture-negative"
), levels = c("Confirmed BSI", "Contaminated culture", "Culture-negative", "Unconclusive"))
microbio_outcome_cons_sensitivity <- bind_rows(
  microbio_outcome_counts %>% mutate(population = "Primary (unchanged)", .before = 1),
  tibble(microbio_outcome = microbio_outcome_CoNS_excluded) %>% count(microbio_outcome) %>%
    mutate(pct = round(100 * n / sum(n), 1), population = "CoNS-excluded sensitivity", .before = 1),
  tibble(microbio_outcome = microbio_outcome_CoNS_included) %>% count(microbio_outcome) %>%
    mutate(pct = round(100 * n / sum(n), 1), population = "CoNS-included sensitivity", .before = 1)
)
write_xlsx(microbio_outcome_cons_sensitivity, "CRF_value_microbio_outcome_CoNS_sensitivity.xlsx")

# -----------------------------------------------------------------------
# 13.2: Outcome-form linkage availability -- a finding in its own right
# -----------------------------------------------------------------------
# SECTION 8 kept only CaseReport patients ALSO present in
# Discharge/Followup28 (matched_control==1, or a genuine outcome-form match)
# when building Main_data_Combined -- 8,610 of 9,912 rows were dropped
# there specifically because they were matched_control==0 with no outcome
# form at all. That gap lands very unevenly across the 4 microbiological
# outcome groups checked here, which is itself directly relevant to this
# section's question: if culture-negative patients essentially never have
# follow-up data available in the first place, that is a strong signal
# about (lack of) current operational value, independent of what any
# mortality COMPARISON on the tiny linked subset can show.
microbio_outcome_linkage <- crf_value_base %>%
  mutate(has_outcome_linkage = patient_id %in% Main_data_Combined$patient_id) %>%
  count(microbio_outcome, has_outcome_linkage) %>%
  tidyr::pivot_wider(names_from = has_outcome_linkage, values_from = n, values_fill = 0, names_prefix = "linked_") %>%
  rename(n_no_outcome_form = linked_FALSE, n_has_outcome_form = linked_TRUE) %>%
  mutate(pct_linked = round(100 * n_has_outcome_form / (n_has_outcome_form + n_no_outcome_form), 1))

# -----------------------------------------------------------------------
# 13.3: Demographic/clinical profile across microbiological outcome groups
# -----------------------------------------------------------------------
# Uses CaseReport-native fields (available for all 9247 regardless of
# outcome-form linkage) -- age, sex, comorbidity burden (SECTION 3.1's own
# charlson_score), prior antibiotic exposure (SECTION 12's abx_exposure_48h,
# re-joined), ward type.
crf_value_base <- crf_value_base %>%
  left_join(diag_raw_cr %>% select(patient_id, abx_period_since_last), by = "patient_id") %>%
  mutate(
    # Same 0-not-NA convention as SECTION 12's abx_exposure_48h -- see that
    # definition's comment for rationale.
    abx_exposure_48h = case_when(
      abx_current_therapy == 1 ~ 1L,
      abx_current_therapy == 0 & abx_period_since_last == "OP1" ~ 1L,
      abx_current_therapy == 0 & abx_period_since_last %in% c("OP2", "OP3", "OP4") ~ 0L,
      TRUE ~ 0L
    )
  )
microbio_profile_cont <- purrr::map_dfr(c("age_years", "charlson_score", "pitt_score_partial"), function(v) {
  crf_value_base %>% group_by(microbio_outcome) %>%
    summarise(mean_val = round(mean(.data[[v]], na.rm = TRUE), 2), sd_val = round(sd(.data[[v]], na.rm = TRUE), 2), .groups = "drop") %>%
    mutate(variable = v, .before = 1)
})
microbio_profile_cat <- purrr::map_dfr(c("sex_at_birth", "ward_type_icu_num", "abx_exposure_48h"), function(v) {
  crf_value_base %>% filter(!is.na(.data[[v]])) %>% group_by(microbio_outcome) %>%
    summarise(pct_value1 = round(100 * mean(.data[[v]] == 1, na.rm = TRUE), 1), n = n(), .groups = "drop") %>%
    mutate(variable = v, .before = 1)
})

# -----------------------------------------------------------------------
# 13.3b: Survey-weighted (calibrated) version of the profile comparison
# -----------------------------------------------------------------------
# The ToR explicitly requires this comparison to "account for the PPS
# cluster sampling design and hospital-level clustering through survey-
# weighted estimation" -- unlike the SECTION 9 mortality/LOS matched-cohort
# analyses (whose separately-recruited uninfected arm has no design weight
# and legitimately uses cluster-robust SEs instead, per Methods), this
# comparison uses only the full CaseReport survey population, for which
# calibrated weights exist for every row. The unweighted table above is kept
# alongside this one because it is easier to read directly; this is the
# ToR-compliant primary estimate. (PPS design-weight-only reporting was
# removed per instruction; the calibrated weight is reported here instead.)
crf_design_data <- crf_value_base %>% filter(!is.na(calibrated_weight))
# Per-(variable, group) estimate built from a FRESH svydesign() on an
# already dplyr-filtered data frame, deliberately NOT subset() on a single
# shared svydesign object: subset.survey.design()'s non-standard evaluation
# of a compound closure-variable condition (`microbio_outcome == lvl &
# !is.na(get(v))`, both lvl and v supplied by purrr::map_dfr) produced
# silently wrong subsets here during development -- e.g. abx_exposure_48h
# among culture-negative patients came out as 3.4% weighted against a
# well-established, independently cross-checked 96.4% weighted estimate for
# the same variable overall (SECTION 12's abx_exposure_summary). Filtering
# the plain data frame with ordinary dplyr first removes that ambiguity
# entirely.
microbio_outcome_levels <- levels(crf_design_data$microbio_outcome)
microbio_profile_cont_weighted <- purrr::map_dfr(c("age_years", "charlson_score", "pitt_score_partial"), function(v) {
  purrr::map_dfr(microbio_outcome_levels, function(lvl) {
    di <- crf_design_data %>% dplyr::filter(as.character(microbio_outcome) == lvl, !is.na(.data[[v]]))
    m <- if (nrow(di) > 0) tryCatch({
      des <- svydesign(ids = ~facility_id_clean, weights = ~calibrated_weight, data = di)
      svymean(as.formula(paste0("~", v)), des, na.rm = TRUE)
    }, error = function(e) NULL) else NULL
    tibble(variable = v, microbio_outcome = lvl, n = nrow(di),
           calibrated_mean = if (is.null(m)) NA_real_ else round(as.numeric(coef(m)[1]), 2),
           calibrated_se = if (is.null(m)) NA_real_ else round(as.numeric(SE(m)[1]), 3))
  })
})
microbio_profile_cat_weighted <- purrr::map_dfr(c("sex_at_birth", "ward_type_icu_num", "abx_exposure_48h"), function(v) {
  purrr::map_dfr(microbio_outcome_levels, function(lvl) {
    di <- crf_design_data %>% dplyr::filter(as.character(microbio_outcome) == lvl, !is.na(.data[[v]]))
    # svyciprop(), not svymean()+coef() -- the SAME construction already
    # proven correct elsewhere in this script for this exact variable
    # (SECTION 12's abx_exposure_summary gives 96.4% weighted overall via
    # svyciprop(~I(abx_exposure_48h==1), ...); the svymean()+coef(m)[1]
    # version tried first here gave a wildly different, implausible 3.4%
    # for the culture-negative subgroup -- svymean()'s coefficient-naming
    # for an I(...)-wrapped logical term is evidently not as reliable as
    # assumed. svyciprop is also the more appropriate estimator for a
    # bounded [0,1] proportion in the first place.
    p <- if (nrow(di) > 0) tryCatch({
      des <- svydesign(ids = ~facility_id_clean, weights = ~calibrated_weight, data = di)
      svyciprop(as.formula(paste0("~I(", v, "==1)")), des, method = "logit")
    }, error = function(e) NULL) else NULL
    tibble(variable = v, microbio_outcome = lvl, n = nrow(di),
           calibrated_pct = if (is.null(p)) NA_real_ else round(as.numeric(p[1]) * 100, 1),
           calibrated_se_pct = if (is.null(p)) NA_real_ else round(as.numeric(SE(p)[1]) * 100, 2))
  })
})

# -----------------------------------------------------------------------
# 13.4: 28-day mortality across microbiological outcome groups AND vs
# matched Uninfected controls
# -----------------------------------------------------------------------
# CAVEAT (directly follows from 13.2): only 18 of 8,452 Culture-negative
# patients (0.2%) and 3 of 4 Unconclusive patients have any linked outcome
# data at all -- these are NOT a random subset (patients who happened to
# also get a Discharge/Followup28 form completed likely differ systematically,
# e.g. longer/more complicated admissions), so the Culture-negative and
# Unconclusive mortality rows below are exploratory at best, not a
# defensible population estimate. Contaminated culture has much better
# linkage (146 of 302, 48%) and Confirmed BSI/Uninfected control match
# SECTION 9's own established populations.
mortality_by_microbio <- Main_data_Combined %>%
  left_join(crf_value_base %>% select(patient_id, microbio_outcome), by = "patient_id") %>%
  mutate(microbio_outcome = if_else(matched_control == 1, "Uninfected control", as.character(microbio_outcome))) %>%
  filter(!is.na(microbio_outcome), !is.na(mortality30day)) %>%
  group_by(microbio_outcome) %>%
  summarise(n = n(), n_deaths = sum(mortality30day == 1), mortality30day_pct = round(100 * mean(mortality30day == 1), 1), .groups = "drop") %>%
  mutate(microbio_outcome = factor(microbio_outcome, levels = c("Uninfected control", "Confirmed BSI", "No confirmed BSI"))) %>%
  arrange(microbio_outcome) %>%
  mutate(reliability = case_when(
    microbio_outcome %in% c("Uninfected control", "Confirmed BSI") ~ "Established population (SECTION 9)",
    microbio_outcome == "Contaminated culture" ~ "Moderate linkage (48%) -- interpret with some caution",
    TRUE ~ "SEVERE linkage gap (<1%) -- exploratory only, not a population estimate"
  ))

# -----------------------------------------------------------------------
# 13.5: Enrolment pathway (formal case definition vs clinician judgement) --
# mortality comparison (positivity already reported in SECTION 12.2)
# -----------------------------------------------------------------------
enrolment_pathway_mortality <- Main_data_Combined %>%
  left_join(diag_raw_cr %>% select(patient_id, enrolment_pathway), by = "patient_id") %>%
  filter(!is.na(enrolment_pathway), !is.na(mortality30day)) %>%
  group_by(enrolment_pathway) %>%
  summarise(n = n(), n_deaths = sum(mortality30day == 1), mortality30day_pct = round(100 * mean(mortality30day == 1), 1), .groups = "drop")

# -----------------------------------------------------------------------
# 13.5b: Clinical characteristics by enrolment pathway (previously only
# positivity, SECTION 12.2, and mortality, above, were compared; the ToR's
# 8.4 output list also asks for demographic/clinical-profile comparison
# between the formal-case-definition and clinician-judgement-only groups).
# -----------------------------------------------------------------------
enrolment_pathway_profile_cont <- purrr::map_dfr(c("age_years", "charlson_score", "pitt_score_partial"), function(v) {
  crf_value_base %>%
    left_join(diag_raw_cr %>% select(patient_id, enrolment_pathway), by = "patient_id") %>%
    filter(!is.na(enrolment_pathway)) %>%
    group_by(enrolment_pathway) %>%
    summarise(mean_val = round(mean(.data[[v]], na.rm = TRUE), 2), sd_val = round(sd(.data[[v]], na.rm = TRUE), 2), n = n(), .groups = "drop") %>%
    mutate(variable = v, .before = 1)
})
enrolment_pathway_profile_cat <- purrr::map_dfr(c("sex_at_birth", "ward_type_icu_num", "abx_exposure_48h"), function(v) {
  crf_value_base %>%
    left_join(diag_raw_cr %>% select(patient_id, enrolment_pathway), by = "patient_id") %>%
    filter(!is.na(enrolment_pathway), !is.na(.data[[v]])) %>%
    group_by(enrolment_pathway) %>%
    summarise(pct_value1 = round(100 * mean(.data[[v]] == 1, na.rm = TRUE), 1), n = n(), .groups = "drop") %>%
    mutate(variable = v, .before = 1)
})

# -----------------------------------------------------------------------
# 13.7: Sensitivity analysis -- potential effect of microbiological
# under-ascertainment on prevalence and incidence estimates
# -----------------------------------------------------------------------
# Per the ToR's SECTION 8.2 closing paragraph: "Where findings suggest
# substantial under-ascertainment of true BSI among culture-negative
# patients, sensitivity analyses may be considered to explore the potential
# effect of missed infections on prevalence, incidence, and mortality
# estimates under alternative assumptions." SECTION 12 found near-universal
# (96%) prior antibiotic exposure before culture and an established
# mechanism (reduced culture sensitivity after antibiotic exposure) by
# which this could cause false-negative cultures; this analysis quantifies
# how much the headline prevalence/incidence estimates would shift if a
# fraction of antibiotic-exposed culture-negative patients were, in truth,
# missed BSIs. A mortality-side sensitivity analysis is NOT attempted:
# SECTION 13.2 already established that only 0.3% of culture-negative
# patients have any linked 28-day outcome data at all, so there is no
# usable outcome data to reclassify or reweight for that group -- reported
# explicitly as a limitation of this sensitivity exercise itself, not
# silently omitted.
underascertainment_scenarios <- c(0.05, 0.10, 0.20)
cn_abx_exposed_pool <- crf_value_base %>%
  filter(microbio_outcome == "No confirmed BSI", abx_exposure_48h == 1, !is.na(calibrated_weight))
n_cn_abx_exposed <- nrow(cn_abx_exposed_pool)

underascertainment_sensitivity <- purrr::map_dfr(underascertainment_scenarios, function(p_missed) {
  n_reclassified <- round(n_cn_abx_exposed * p_missed)
  sens_data <- Main_data_CaseReport
  if (n_reclassified > 0) {
    reclass_ids <- cn_abx_exposed_pool$patient_id[seq_len(n_reclassified)]
    sens_data <- sens_data %>% mutate(BSI_positive_sens = if_else(patient_id %in% reclass_ids, 1L, BSI_positive))
  } else {
    sens_data <- sens_data %>% mutate(BSI_positive_sens = BSI_positive)
  }
  sens_design <- svydesign(ids = ~facility_id_clean, weights = ~calibrated_weight, data = sens_data %>% filter(!is.na(calibrated_weight)))
  prev <- svyciprop(~I(BSI_positive_sens == 1), sens_design, method = "logit")
  ci <- confint(prev)
  prev_num <- as.numeric(prev[1]); ci_lo <- as.numeric(ci[1]); ci_hi <- as.numeric(ci[2])
  tibble(
    scenario = paste0(scales::percent(p_missed), " of antibiotic-exposed culture-negative patients assumed missed BSI"),
    n_reclassified_to_confirmed = as.integer(n_reclassified),
    n_confirmed_BSI_sensitivity = as.integer(sum(sens_data$BSI_positive_sens == 1, na.rm = TRUE)),
    calibrated_prevalence_pct = round(prev_num * 100, 2),
    calibrated_prevalence_ci_lower = round(ci_lo * 100, 2),
    calibrated_prevalence_ci_upper = round(ci_hi * 100, 2),
    national_cases_extrapolated = round(prev_num * national_admissions_total)
  )
})
underascertainment_sensitivity <- bind_rows(
  tibble(scenario = "Primary estimate (no reclassification)", n_reclassified_to_confirmed = 0L,
         n_confirmed_BSI_sensitivity = sum(Main_data_CaseReport$BSI_positive == 1, na.rm = TRUE),
         calibrated_prevalence_pct = round(prevalence_overall_calibrated[1] * 100, 2),
         calibrated_prevalence_ci_lower = round(prevalence_overall_calibrated_ci[1] * 100, 2),
         calibrated_prevalence_ci_upper = round(prevalence_overall_calibrated_ci[2] * 100, 2),
         national_cases_extrapolated = round(prevalence_overall_calibrated[1] * national_admissions_total)),
  underascertainment_sensitivity
)

write_xlsx(
  list(microbio_outcome_counts = microbio_outcome_counts, outcome_linkage_availability = microbio_outcome_linkage,
       profile_continuous = microbio_profile_cont, profile_categorical = microbio_profile_cat,
       profile_continuous_calibrated = microbio_profile_cont_weighted, profile_categorical_calibrated = microbio_profile_cat_weighted,
       mortality_by_microbio_outcome = mortality_by_microbio, enrolment_pathway_mortality = enrolment_pathway_mortality,
       enrolment_pathway_profile_cont = enrolment_pathway_profile_cont, enrolment_pathway_profile_cat = enrolment_pathway_profile_cat,
       underascertainment_sensitivity = underascertainment_sensitivity),
  "CRF_value_culturenegative_contaminated.xlsx"
)

# -----------------------------------------------------------------------
# 13.6: Structured recommendation (per the spec's 8.3) -- synthesised from
# the actual numbers computed above, not a template.
# -----------------------------------------------------------------------
# 1. Scale: of 9,247 enrolled patients, 489 (5.3%) were confirmed BSI, 302
#    (3.3%) were contaminated cultures, 8,452 (91.4%) were culture-negative,
#    4 unconclusive -- so culture-negative patients are the overwhelming
#    majority of the CRF-collection burden this section is asking about.
# 2. The single clearest finding is the outcome-linkage gap itself (13.2):
#    the CURRENT survey infrastructure links 28-day outcome data to only
#    ~0.2% of culture-negative patients, vs 96% of confirmed-BSI patients
#    and 48% of contaminated-culture patients. Whatever full CRF collection
#    for culture-negative patients is meant to buy analytically, THIS
#    dataset mostly hasn't been able to cash it in -- not because
#    culture-negative patients are unimportant, but because the operational
#    linkage between their CRF and their Discharge/Followup28 forms is
#    almost entirely missing in practice.
# 3. On the (very thin, n=18) mortality signal that IS available:
#    Contaminated-culture patients' mortality (SECTION 13.4 output) sits
#    much closer to Confirmed-BSI than to Uninfected control, on
#    reasonably-linked data (48%) -- suggesting contaminated-culture
#    patients are, on average, clinically similar to true BSI patients
#    (plausible: a positive-but-contaminated culture is still typically
#    drawn because the patient looked septic) and are NOT simply
#    "false alarms" safe to deprioritise. Culture-negative patients' own
#    mortality figure is not a defensible population estimate (point 2)
#    and should not by itself be used to argue either for or against
#    continued full CRF collection.
# 4. Enrolment pathway (13.5, SECTION 12.2): clinician-judgement-only
#    patients (parenteral antibiotics indicated without meeting the formal
#    scored case definition) still have a non-trivial BSI positivity rate
#    (SECTION 12.2's enrolment_pathway_positivity) -- clinician-only
#    indication is not simply noise, supporting continued enrolment via
#    that pathway rather than restricting to the formal case definition
#    alone.
# NET RECOMMENDATION: the analytical case for full CRF collection on
# contaminated-culture patients is reasonably well supported (they resemble
# confirmed BSI clinically and are a small group, ~3% of enrolment -- low
# marginal burden for real signal). The case for culture-negative patients
# (91% of enrolment, the actual operational burden) cannot be assessed on
# THIS dataset because their outcome data was essentially never linked --
# the priority recommendation for future surveys is fixing that linkage
# (ensuring Discharge/Followup28 forms are completed and joinable for a
# representative sample of culture-negative patients), not simplifying
# their CRF, since the current gap is a follow-up/linkage problem rather
# than evidence that the CRF itself lacks value.


# =============================================================================
# SECTION 14: Data collection rationalisation and variable reduction
# =============================================================================
# Corresponds to the methods spec's "9. Data Collection Rationalisation and
# Variable Reduction". Per the spec's own framing, this characterises
# analytical contribution and redundancy -- it does NOT determine causal
# importance or unilaterally recommend removal (SECTION 9.2's own words),
# so no causal-timing restriction is applied here (unlike SECTIONS 9-11):
# every candidate variable, including ward/specialty/organism fields
# excluded from the CAUSAL models elsewhere in this script, is fair game
# for a pure redundancy/prediction-contribution exercise.
#
# SCOPE: the spec names 4 source forms (CRF/CaseReport, Cohort 3, discharge,
# 28-day follow-up). Main_data_CaseReport's own already-curated clinical
# fields remain the primary candidate set -- comorbidities (comorb_*), risk
# factors (riskfx_*), admission reasons (reason_*), demographics, ward/
# hospital type, severity, and prior antibiotic exposure -- since this is
# the single richest, most redundancy-prone form (many of the ~240+ raw
# CaseReport columns alone are IDs, free text, or documented in SECTION 15's
# own codebook as structurally near-empty, e.g. the INN drug-name fields).
# EXTENDED here (previously scoped out) to also include the discharge
# form's own clinical covariate fields -- ICU stay during the hospital
# episode, in-hospital complications (any/count), and treatment-intensity
# markers (invasive ventilation, ECMO, vasopressors, renal replacement
# therapy) -- joined from Main_data_hospitalDischarge_clean and a fresh
# extraction of the 4 treatment-history flags (not previously built
# anywhere in this script). The 28-day follow-up form's own fields are
# overwhelmingly outcome-ascertainment (status/location at day 28, care-
# seeking since discharge) rather than clinical covariates comparable to
# the other three forms' candidates, and are deliberately not added here as
# PREDICTORS of the very BSI/resistance/mortality outcomes this section
# models -- doing so would risk exactly the reverse-causation/collider
# issue SECTIONS 9-11 already avoid for ward type. Cohort 3's own
# comorbidity/risk-factor/admission-reason fields use an IDENTICAL variable
# taxonomy to CaseReport's (SECTION 7's comorb_map_cohort/riskfx_map_cohort/
# reason_map_cohort reuse CaseReport's own `new` names) but cannot share
# this section's outcome models -- the matched-control population is, by
# construction, 100% uninfected, so confirmed-BSI and resistant-BSI have no
# variation to model there. Cohort 3 accordingly gets its own, separate
# redundancy-only analysis (14.1b, below), directly using the same
# cramers_v()/pairwise_association() machinery on its own population.
discharge_treatment_flags <- read_excel(DHIS2_XLSX, sheet = "Discharge", guess_max = 1e6) %>%
  transmute(
    patient_id = `Unique survey patient ID #1`,
    discharge_invasive_ventilation = case_when(
      toupper(trimws(`Treatment history: Invasive ventilation (any)?`)) %in% c("YES", "TRUE") ~ 1L,
      toupper(trimws(`Treatment history: Invasive ventilation (any)?`)) %in% c("NO", "FALSE") ~ 0L, TRUE ~ NA_integer_),
    discharge_ecmo = case_when(
      toupper(trimws(`Treatment history: Extracorporeal (ECMO) support?`)) %in% c("YES", "TRUE") ~ 1L,
      toupper(trimws(`Treatment history: Extracorporeal (ECMO) support?`)) %in% c("NO", "FALSE") ~ 0L, TRUE ~ NA_integer_),
    discharge_vasopressors = case_when(
      toupper(trimws(`Treatment history: Inotropes/ vasopressors?`)) %in% c("YES", "TRUE") ~ 1L,
      toupper(trimws(`Treatment history: Inotropes/ vasopressors?`)) %in% c("NO", "FALSE") ~ 0L, TRUE ~ NA_integer_),
    discharge_rrt = case_when(
      toupper(trimws(`Treatment history: Renal replacement therapy (RRT) or dialysis?`)) %in% c("YES", "TRUE") ~ 1L,
      toupper(trimws(`Treatment history: Renal replacement therapy (RRT) or dialysis?`)) %in% c("NO", "FALSE") ~ 0L, TRUE ~ NA_integer_)
  ) %>%
  filter(!is.na(patient_id) & !grepl("DATABASE|TEST", patient_id, ignore.case = TRUE)) %>%
  distinct(patient_id, .keep_all = TRUE)

# The 4 "Treatment history" fields (invasive ventilation, ECMO,
# vasopressors, RRT) are 100% missing across all 1275 Discharge form rows
# in this extract -- verified directly against the raw sheet, not a join
# artifact. Reported as its own finding (discharge_treatment_history_
# completeness, below) since 100% non-completion of a defined CRF field is
# itself exactly the kind of data-collection-rationalisation evidence this
# section exists to surface -- but excluded from the modelling candidate
# set below: a column that is always NA carries zero information and, when
# included in ANY complete-case model (even alongside just the other 3
# discharge fields, tried first), collapses that model's usable sample to
# zero rather than simply contributing nothing.
discharge_treatment_history_completeness <- tibble(
  field = c("discharge_invasive_ventilation", "discharge_ecmo", "discharge_vasopressors", "discharge_rrt"),
  n_discharge_forms = nrow(Main_data_hospitalDischarge), n_non_missing = 0L,
  note = "100% missing in this extract -- verified against the raw Discharge sheet, not a linkage/join artifact"
)

candidate_vars_flags <- grep("^comorb_|^riskfx_|^reason_", names(Main_data_CaseReport), value = TRUE)
candidate_vars_core <- c("age_years", "sex_at_birth", "ward_type_icu", "specialty_group", "charlson_score",
                          "abx_current_therapy", "infection_origin_composite", "pitt_score_partial")
candidate_vars_discharge <- c("ICU_entry", "any_complication", "complications_count")
# candidate_vars: the CRF-only main set (unchanged scope, used for the
# redundancy matrix/FAMD/LASSO/RF below) -- candidate_vars_discharge is
# deliberately NOT folded in here. Requiring joint complete-case-ness across
# the discharge form's own linkage gap (most patients, especially culture-
# negative ones, never have ANY discharge form) AND the ~20 already-fairly-
# missing CRF fields compounds multiplicatively: tried first, and it
# collapsed the resistant-vs-susceptible and mortality LASSO/RF complete-
# case n to literally zero. Discharge-form variable contribution is instead
# assessed as its own smaller, separate model (14.3b/14.4b, below) --
# exactly the same reasoning already applied to Cohort 3's fields.
candidate_vars <- c(candidate_vars_core, candidate_vars_flags)
# candidate_vars_redundancy: discharge vars ARE included here, since the
# pairwise redundancy matrix (14.1) computes each pair's own complete-case
# overlap independently (pairwise_association()'s own `sum(ok) < 10` guard)
# rather than requiring one joint complete-case set across every candidate
# -- no collapse risk, and it is a directly ToR-requested output.
candidate_vars_redundancy <- c(candidate_vars, candidate_vars_discharge)

varred_base <- Main_data_CaseReport %>%
  select(patient_id, BSI_positive, resistant_BSI_positive, all_of(candidate_vars)) %>%
  left_join(Main_data_hospitalDischarge_clean %>% select(patient_id, ICU_entry, any_complication, complications_count), by = "patient_id") %>%
  left_join(discharge_treatment_flags, by = "patient_id")

# -----------------------------------------------------------------------
# 14.1: Pairwise redundancy matrix -- Cramer's V (categorical/binary pairs),
# Spearman (continuous pairs or continuous-vs-binary)
# -----------------------------------------------------------------------
cramers_v <- function(x, y) {
  tab <- table(x, y)
  if (any(dim(tab) < 2)) return(NA_real_)
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
  n <- sum(tab)
  sqrt(as.numeric(chi) / (n * (min(dim(tab)) - 1)))
}
continuous_vars <- c("age_years", "charlson_score", "pitt_score_partial", "complications_count")
categorical_vars <- setdiff(candidate_vars, continuous_vars)

pairwise_association <- function(v1, v2, data) {
  x <- data[[v1]]; y <- data[[v2]]
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) < 10) return(NA_real_)
  x <- x[ok]; y <- y[ok]
  if (v1 %in% continuous_vars && v2 %in% continuous_vars) {
    suppressWarnings(cor(x, y, method = "spearman"))
  } else if (!(v1 %in% continuous_vars) && !(v2 %in% continuous_vars)) {
    cramers_v(x, y)
  } else {
    # mixed continuous/categorical pair: Spearman on the categorical side
    # numerically coded (as.numeric(factor(.)), not just as-is since some
    # candidates, e.g. ward_type_icu/specialty_group, are character strings
    # cor() can't take directly) -- equivalent to a rank-biserial correlation
    # for the binary case, a pragmatic single-metric choice for a mixed pair
    # (the spec's Cramer's V/Pearson-Spearman split covers same-type pairs;
    # mixed pairs are not separately specified).
    if (!is.numeric(x)) x <- as.numeric(as.factor(x))
    if (!is.numeric(y)) y <- as.numeric(as.factor(y))
    suppressWarnings(cor(x, y, method = "spearman"))
  }
}

var_pairs <- combn(candidate_vars_redundancy, 2, simplify = FALSE)
redundancy_long <- purrr::map_dfr(var_pairs, function(p) {
  tibble(var1 = p[1], var2 = p[2], association = pairwise_association(p[1], p[2], varred_base))
})
redundancy_matrix_wide <- redundancy_long %>%
  bind_rows(redundancy_long %>% rename(var1 = var2, var2 = var1)) %>%
  tidyr::pivot_wider(names_from = var2, values_from = association)

# Highly correlated/redundant pairs: |association| >= 0.6 (conventional
# "strong" threshold for both Cramer's V and Spearman), excluding
# structurally paired variables (e.g. comorb_diabetes/comorb_diabetes_end_organ
# are DESIGNED to be correlated -- severity tiers of the same condition, not
# redundant data-collection burden).
high_redundancy_pairs <- redundancy_long %>%
  filter(!is.na(association), abs(association) >= 0.6) %>%
  arrange(desc(abs(association)))

# -----------------------------------------------------------------------
# 14.1b: Cohort 3 form -- separate redundancy matrix on its own population
# -----------------------------------------------------------------------
# Cohort 3's comorbidity/risk-factor/admission-reason fields reuse
# CaseReport's own variable names (comorb_map_cohort etc., SECTION 7), so
# the SAME candidate_vars_flags/cramers_v()/pairwise_association() machinery
# applies directly -- just to the 464 pooled matched controls instead of
# the 9222 CaseReport patients. No FAMD/LASSO/RF here: those need an
# outcome to model against (confirmed-BSI, resistant-BSI), and the matched-
# control population is, by construction, 100% uninfected -- no variation
# in either outcome to model. Redundancy alone is still a well-defined,
# ToR-requested question independent of any outcome.
cohort3_candidate_vars <- intersect(c(candidate_vars_core[1:2], candidate_vars_flags), names(Main_data_hospitalCohortEnrolment_clean))
cohort3_varred_base <- Main_data_hospitalCohortEnrolment_clean %>%
  filter(matched_control == 1) %>%
  select(patient_id, all_of(cohort3_candidate_vars))
cohort3_var_pairs <- combn(cohort3_candidate_vars, 2, simplify = FALSE)
cohort3_redundancy_long <- purrr::map_dfr(cohort3_var_pairs, function(p) {
  tibble(var1 = p[1], var2 = p[2], association = pairwise_association(p[1], p[2], cohort3_varred_base))
})
cohort3_high_redundancy_pairs <- cohort3_redundancy_long %>%
  filter(!is.na(association), abs(association) >= 0.6) %>%
  arrange(desc(abs(association)))

# -----------------------------------------------------------------------
# 14.2: FAMD (Factor Analysis of Mixed Data) -- overall structure, not a
# removal mechanism (per the spec's own caveat, SECTION 9.2)
# -----------------------------------------------------------------------
# Complete-case only for this exploratory structural pass (a full MICE
# pipeline, as used for the causal models in SECTION 9, would be
# disproportionate for a structure-characterisation step the spec itself
# says should not mechanically drive removal decisions).
famd_data <- varred_base %>%
  select(all_of(candidate_vars)) %>%
  mutate(across(all_of(categorical_vars), as.factor)) %>%
  tidyr::drop_na()
famd_n_complete <- nrow(famd_data)
famd_result <- tryCatch(
  FactoMineR::FAMD(famd_data, ncp = 5, graph = FALSE),
  error = function(e) NULL
)
famd_var_contributions <- if (is.null(famd_result)) {
  tibble(variable = candidate_vars, note = "FAMD did not converge on the complete-case subset")
} else {
  as.data.frame(famd_result$var$contrib) %>%
    tibble::rownames_to_column("variable") %>%
    mutate(across(where(is.numeric), ~ round(.x, 2)))
}

# -----------------------------------------------------------------------
# 14.3: Penalised regression (LASSO, glmnet) -- variable contribution to 3
# key outcomes: confirmed BSI vs culture-negative, resistant vs susceptible
# BSI, 28-day mortality
# -----------------------------------------------------------------------
# Reports coefficients at BOTH lambda.min (the CV-optimal fit -- more
# informative, typically retains a working-sized variable set) and
# lambda.1se (the conventional "most parsimonious model within 1 SE of
# optimal" -- often much more conservative, and can legitimately retain
# ZERO variables when no single predictor's signal clearly beats the noise
# floor, which is itself a valid, reportable finding, not a bug to work
# around by only ever showing the more generous lambda).
run_lasso_importance <- function(data, outcome_var, covars) {
  di <- data %>% select(all_of(c(outcome_var, covars))) %>% tidyr::drop_na()
  # Zero-variance guard: a rare candidate (e.g. discharge_ecmo) can be
  # constant (a single level/value) within a small outcome-specific
  # complete-case subset even though it varies in the full dataset --
  # model.matrix()'s contrast construction errors outright on a
  # single-level factor, so such columns are dropped here (not included,
  # not silently miscoded) and reported as such via the returned tibble's
  # dropped_zero_variance rows.
  n_distinct_by_col <- sapply(di[covars], function(x) length(unique(x)))
  zero_var_covars <- names(n_distinct_by_col)[n_distinct_by_col < 2]
  covars_used <- setdiff(covars, zero_var_covars)
  if (length(covars_used) == 0) {
    return(tibble(term = covars, coefficient_lambda_min = NA_real_, coefficient_lambda_1se = NA_real_,
                  retained_lambda_min = NA, retained_lambda_1se = NA, n_used = nrow(di)))
  }
  x <- model.matrix(as.formula(paste("~", paste(covars_used, collapse = " + "))), data = di)[, -1, drop = FALSE]
  y <- di[[outcome_var]]
  cvfit <- glmnet::cv.glmnet(x, y, family = "binomial", alpha = 1)
  coefs_min <- as.matrix(coef(cvfit, s = "lambda.min"))
  coefs_1se <- as.matrix(coef(cvfit, s = "lambda.1se"))
  out <- tibble(term = rownames(coefs_min), coefficient_lambda_min = coefs_min[, 1], coefficient_lambda_1se = coefs_1se[, 1]) %>%
    dplyr::filter(term != "(Intercept)") %>%
    mutate(retained_lambda_min = coefficient_lambda_min != 0, retained_lambda_1se = coefficient_lambda_1se != 0, n_used = nrow(di)) %>%
    arrange(desc(abs(coefficient_lambda_min)))
  if (length(zero_var_covars) > 0) {
    out <- bind_rows(out, tibble(term = zero_var_covars, coefficient_lambda_min = NA_real_, coefficient_lambda_1se = NA_real_,
                                  retained_lambda_min = NA, retained_lambda_1se = NA, n_used = nrow(di)))
  }
  out
}

# Discharge-form-derived variables (candidate_vars_discharge) are excluded
# from the confirmed-BSI-vs-culture-negative model specifically: culture-
# negative patients essentially never have a linked Discharge form
# (SECTION 13.2 -- 0.3% linkage), so requiring complete cases on these
# variables here would silently collapse this outcome's whole population
# down to the tiny, non-representative linked subset. They are kept for the
# resistant-vs-susceptible and mortality outcomes below, where the
# population already requires outcome-form linkage by construction.
candidate_vars_no_discharge <- setdiff(candidate_vars, candidate_vars_discharge)
lasso_bsi_vs_negative <- run_lasso_importance(
  varred_base %>% mutate(BSI_positive = as.integer(BSI_positive)), "BSI_positive", candidate_vars_no_discharge
) %>% mutate(outcome = "Confirmed BSI vs culture-negative", .before = 1)

lasso_resistant_vs_susceptible <- run_lasso_importance(
  varred_base %>% filter(BSI_positive == 1, !is.na(resistant_BSI_positive)) %>% mutate(resistant_BSI_positive = as.integer(resistant_BSI_positive)),
  "resistant_BSI_positive", candidate_vars
) %>% mutate(outcome = "Resistant vs susceptible BSI", .before = 1)

mortality_varred_data <- mortality_analytic %>%
  mutate(charlson_score = charlson_score_combined) %>%
  left_join(discharge_treatment_flags, by = "patient_id") %>%
  select(event28, all_of(intersect(candidate_vars, names(.))))
lasso_mortality <- run_lasso_importance(
  mortality_varred_data %>% mutate(event28 = as.integer(event28)), "event28",
  intersect(candidate_vars, names(mortality_varred_data))
) %>% mutate(outcome = "28-day mortality", .before = 1)

lasso_all_results <- bind_rows(lasso_bsi_vs_negative, lasso_resistant_vs_susceptible, lasso_mortality) %>%
  mutate(across(starts_with("coefficient_"), ~ round(.x, 4)))

# -----------------------------------------------------------------------
# 14.4: Random forest permutation importance (ranger) -- supportive measure
# -----------------------------------------------------------------------
run_rf_importance <- function(data, outcome_var, covars) {
  di <- data %>% select(all_of(c(outcome_var, covars))) %>% tidyr::drop_na()
  # Same zero-variance guard as run_lasso_importance() -- a constant column
  # within this outcome's complete-case subset carries no information and
  # some ranger code paths mishandle a single-level factor predictor.
  n_distinct_by_col <- sapply(di[covars], function(x) length(unique(x)))
  covars_used <- setdiff(covars, names(n_distinct_by_col)[n_distinct_by_col < 2])
  if (length(covars_used) == 0) return(tibble(term = character(0), importance = numeric(0), n_used = nrow(di)))
  di <- di %>% select(all_of(c(outcome_var, covars_used))) %>%
    mutate(across(where(is.character), as.factor), !!outcome_var := as.factor(.data[[outcome_var]]))
  form <- as.formula(paste(outcome_var, "~ ."))
  rf <- ranger::ranger(form, data = di, importance = "permutation", probability = FALSE, num.trees = 500, seed = 46)
  tibble(term = names(rf$variable.importance), importance = rf$variable.importance, n_used = nrow(di)) %>%
    arrange(desc(importance))
}
rf_bsi_vs_negative <- run_rf_importance(varred_base %>% mutate(BSI_positive = as.character(BSI_positive)), "BSI_positive", candidate_vars_no_discharge) %>%
  mutate(outcome = "Confirmed BSI vs culture-negative", .before = 1)
rf_resistant_vs_susceptible <- run_rf_importance(
  varred_base %>% filter(BSI_positive == 1, !is.na(resistant_BSI_positive)) %>% mutate(resistant_BSI_positive = as.character(resistant_BSI_positive)),
  "resistant_BSI_positive", candidate_vars
) %>% mutate(outcome = "Resistant vs susceptible BSI", .before = 1)
rf_mortality <- run_rf_importance(
  mortality_varred_data %>% mutate(event28 = as.character(event28)), "event28", intersect(candidate_vars, names(mortality_varred_data))
) %>% mutate(outcome = "28-day mortality", .before = 1)
rf_all_results <- bind_rows(rf_bsi_vs_negative, rf_resistant_vs_susceptible, rf_mortality) %>%
  mutate(importance = round(importance, 5))

# -----------------------------------------------------------------------
# 14.3b/14.4b: Discharge-form variable contribution, as its own smaller
# model (not pooled with the main CRF candidate set -- see the comment at
# candidate_vars_discharge's definition, above, for why). Confirmed-BSI-
# vs-culture-negative is not attempted here either: culture-negative
# patients essentially never have a linked Discharge form (SECTION 13.2),
# so that contrast has no meaningful discharge-side population to model.
# -----------------------------------------------------------------------
discharge_varred_resistant <- varred_base %>%
  filter(BSI_positive == 1, !is.na(resistant_BSI_positive)) %>%
  select(resistant_BSI_positive, all_of(candidate_vars_discharge))
discharge_varred_mortality <- mortality_varred_data %>%
  select(event28, all_of(intersect(candidate_vars_discharge, names(mortality_varred_data))))

lasso_discharge_resistant <- tryCatch(
  run_lasso_importance(discharge_varred_resistant %>% mutate(resistant_BSI_positive = as.integer(resistant_BSI_positive)),
                        "resistant_BSI_positive", candidate_vars_discharge),
  error = function(e) tibble(term = candidate_vars_discharge, coefficient_lambda_min = NA_real_, coefficient_lambda_1se = NA_real_,
                              retained_lambda_min = NA, retained_lambda_1se = NA, n_used = nrow(discharge_varred_resistant))
) %>% mutate(outcome = "Resistant vs susceptible BSI (discharge-form variables)", .before = 1)
lasso_discharge_mortality <- tryCatch(
  run_lasso_importance(discharge_varred_mortality %>% mutate(event28 = as.integer(event28)),
                        "event28", intersect(candidate_vars_discharge, names(discharge_varred_mortality))),
  error = function(e) tibble(term = candidate_vars_discharge, coefficient_lambda_min = NA_real_, coefficient_lambda_1se = NA_real_,
                              retained_lambda_min = NA, retained_lambda_1se = NA, n_used = nrow(discharge_varred_mortality))
) %>% mutate(outcome = "28-day mortality (discharge-form variables)", .before = 1)
lasso_discharge_results <- bind_rows(lasso_discharge_resistant, lasso_discharge_mortality) %>%
  mutate(across(starts_with("coefficient_"), ~ round(.x, 4)))

rf_discharge_resistant <- tryCatch(
  run_rf_importance(discharge_varred_resistant %>% mutate(resistant_BSI_positive = as.character(resistant_BSI_positive)),
                     "resistant_BSI_positive", candidate_vars_discharge),
  error = function(e) tibble(term = character(0), importance = numeric(0), n_used = integer(0))
) %>% mutate(outcome = "Resistant vs susceptible BSI (discharge-form variables)", .before = 1)
rf_discharge_mortality <- tryCatch(
  run_rf_importance(discharge_varred_mortality %>% mutate(event28 = as.character(event28)),
                     "event28", intersect(candidate_vars_discharge, names(discharge_varred_mortality))),
  error = function(e) tibble(term = character(0), importance = numeric(0), n_used = integer(0))
) %>% mutate(outcome = "28-day mortality (discharge-form variables)", .before = 1)
rf_discharge_results <- bind_rows(rf_discharge_resistant, rf_discharge_mortality) %>%
  mutate(importance = round(importance, 5))

# -----------------------------------------------------------------------
# 14.5: Classification -- essential / useful / candidate for removal
# -----------------------------------------------------------------------
# Synthesised (not mechanical, per the spec's own instruction). Uses
# lambda.min (not lambda.1se, which retained 0 variables for every outcome
# -- see 14.3's own comment) for LASSO retention, and a RELATIVE RF
# importance threshold (>=50% of that outcome's own top variable) rather
# than a fixed top-N cutoff -- with only 53 candidate variables, a fixed
# "top 10 per outcome x 3 outcomes" cutoff flagged 41/53 (77%) as
# "Essential" on the first pass, which defeats the point of this exercise;
# a relative threshold self-adjusts to how concentrated (or flat) each
# outcome's actual importance signal is.
# "Essential" = LASSO(lambda.min)-retained for >=2 of 3 outcomes, OR RF
# importance >=50% of the top variable for that outcome, for any outcome.
# "Candidate for removal" = retained by no outcome's LASSO, not RF-notable
# for any outcome, AND paired (14.1, |association|>=0.6) with a variable
# that IS essential (i.e. genuinely substitutable, not just unimportant on
# its own). Everything else is "Useful".
lasso_retained_count <- lasso_all_results %>% filter(retained_lambda_min, term %in% candidate_vars) %>% count(term, name = "n_outcomes_retained")
rf_notable <- rf_all_results %>%
  group_by(outcome) %>%
  mutate(rel_importance = importance / max(importance, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(rel_importance >= 0.5) %>%
  distinct(term)
redundant_with_essential <- high_redundancy_pairs %>%
  filter(var1 %in% lasso_retained_count$term | var2 %in% lasso_retained_count$term) %>%
  transmute(candidate = if_else(var1 %in% lasso_retained_count$term, var2, var1)) %>%
  distinct()

variable_classification <- tibble(variable = candidate_vars) %>%
  left_join(lasso_retained_count, by = c(variable = "term")) %>%
  mutate(
    n_outcomes_retained = coalesce(n_outcomes_retained, 0L),
    is_rf_notable = variable %in% rf_notable$term,
    is_redundant_with_essential = variable %in% redundant_with_essential$candidate,
    classification = case_when(
      n_outcomes_retained >= 2 | is_rf_notable ~ "Essential",
      n_outcomes_retained == 0 & !is_rf_notable & is_redundant_with_essential ~ "Candidate for removal (redundant)",
      TRUE ~ "Useful"
    )
  ) %>%
  arrange(match(classification, c("Essential", "Useful", "Candidate for removal (redundant)")))

write_xlsx(
  list(redundancy_matrix = redundancy_matrix_wide, high_redundancy_pairs = high_redundancy_pairs,
       famd_variable_contributions = famd_var_contributions, lasso_variable_contribution = lasso_all_results,
       rf_permutation_importance = rf_all_results, variable_classification = variable_classification,
       cohort3_redundancy_matrix = cohort3_redundancy_long, cohort3_high_redundancy_pairs = cohort3_high_redundancy_pairs,
       discharge_lasso_variable_contribution = lasso_discharge_results, discharge_rf_permutation_importance = rf_discharge_results,
       discharge_treatment_history_completeness = discharge_treatment_history_completeness),
  "VariableReduction_redundancy_and_importance.xlsx"
)


# =============================================================================
# SECTION 15: Automated codebook
# =============================================================================
# variable_registry maps every variable name to a human description, built
# automatically from the same (orig -> new) mapping tables used above so it
# can't drift out of sync with the recoding logic, plus a manual tribble for
# one-off variables. Anything present in Main_data_CaseReport but missing
# from the registry still appears in the codebook with a placeholder
# description instead of breaking the export -- re-run this SECTION any
# time after adding new variables.

manual_registry <- tribble(
  ~variable, ~description,
  "patient_id", "Unique patient identifier for this survey.",
  "facility_id", "Unique facility identifier.",
  "ward_id", "Unique identifier of the ward/department that enrolled the patient.",
  "enrolment_date", "Date the patient was enrolled into the survey (parsed to datetime).",
  "transfer_date_F1", "Date of transfer to inpatient ward, case report form field F1 (parsed to datetime).",
  "days_enrolment_to_transfer", "transfer_date_F1 minus enrolment_date, in days (with decimals).",
  "dob", "Patient date of birth (parsed to datetime).",
  "age_years", "Current chronological age = (script run time - dob) / 365.25, rounded to 2 decimals. NOT age at a fixed clinical event -- changes on every re-run. Some values are implausible due to source data-entry errors (see DATA QUALITY NOTE in SECTION 1).",
  "sex_at_birth", "1 = FEMALE, 0 = MALE, NA = blank/other in source.",
  "age_group", "0 = NEONATE, 1 = PAEDIATRIC, 3 = ADULT (coding as explicitly requested; note value 2 is intentionally unused), NA = blank in source.",
  "speciality_code", "Raw ward speciality code (e.g. AMW-OBG) -- see ward_type_icu/specialty_group for the derived groupings.",
  "ward_type_icu", "Derived from speciality_code via ward_code_lookup: 'ICU' vs 'Non-ICU'.",
  "specialty_group", "Derived from speciality_code via ward_code_lookup: 'Medical' vs 'Surgical'; NA for bare ICU codes that mix both.",
  "bsi_onset_datetime", "Date/time of onset of clinical signs and/or symptoms of BSI (parsed to datetime).",
  "hospital_admissionDate", "Hospital admission date AND time. Used (not hospital_admission_date_only) for hosp_comm_acquis48 since it has time-of-day precision.",
  "hospital_admission_date_only", "Hospital admission date without time-of-day (a separate source field from hospital_admissionDate).",
  "blood_sample1_date", "Date/time of first blood culture sample collection.",
  "blood_sample2_date", "Date/time of second blood culture sample collection.",
  "blood_sample1_ID", "First blood culture sample ID.",
  "blood_sample2_ID", "Second blood culture sample ID.",
  "hours_admission_to_sample1", "Hours from hospital_admissionDate to blood_sample1_date. NA if either date missing or admission is after sampling (implausible sequence).",
  "hosp_comm_acquis48", "1 if hours_admission_to_sample1 > 48 (hospital-acquired), 0 if <=48 (community-acquired), NA if hours_admission_to_sample1 is NA.",
  "infection_origin_hosp", "From 'Likely infection origin: Please indicate likely origin of the infection' (clinician-recorded, distinct from the timing-based hosp_comm_acquis48). 1 = HOSPITAL, 0 = COMMUNITY, NA = blank.",
  "infection_origin_composite", "Hierarchical hospital- vs community-acquired classification combining infection_origin_hosp, hosp_comm_acquis48, and 5 structural proxy variables (transferred_other_ward, icu_48h_prior, surgery_prior_to_bsi, transferred_other_hospital, discharged_48h_prior), in that order of confidence. 1 = hospital-acquired, 0 = community-acquired. Coverage 97.9% (vs 48.5% for infection_origin_hosp alone) -- see infection_origin_source for which tier produced each value, and the SECTION 1 comment for full reasoning/caveats.",
  "infection_origin_source", "Which tier produced infection_origin_composite for this patient: 'clinician' (most direct), 'timing_48h', 'proxy_positive', or 'proxy_all_negative' (weakest tier -- absence of 5 specific exposures inferred as community-acquired, not confirmed). NA if unclassifiable by any tier (2.1% of patients).",
  "infection_origin_composite_strict", "Conservative companion to infection_origin_composite: identical values for tiers 1-2 (clinician + 48h timing), but NA wherever the full composite relied on tier 3 (the structural-proxy inference). Coverage 49.8% (vs 97.9% for the full composite, 48.5% for infection_origin_hosp alone). Use this to sensitivity-check any result that changes materially when the proxy tier is included.",
  "abx_current_therapy", "From 'Antimicrobial use history: Is the patient currently under antibiotic therapy?'. 1 = true, 0 = false, NA = blank.",
  "abx_current_1", "Renamed from 'Antimicrobial use history: Indicate the drug (INN) #1'. 100% EMPTY in this extract -- kept for structural completeness only. No duration/gap-to-culture variables were built for this domain (no start/end dates exist in the source).",
  "abx_current_2", "Renamed from 'Antimicrobial use history: Indicate the drug (INN) #2'. 100% EMPTY in this extract.",
  "abx_current_3", "Renamed from 'Antimicrobial use history: Indicate the drug (INN) #3'. 100% EMPTY in this extract.",
  "abx_current_4", "Renamed from 'Antimicrobial use history: Indicate the drug (INN) #4'. 100% EMPTY in this extract.",
  "pitt_temp_points", "1 point if Pitt-domain temperature <36C or >=40C, else 0 (or 0 if temperature missing). Sub-component of pitt_score_partial.",
  "pitt_score_partial", "PARTIAL Pitt bacteraemia score = pitt_temp_points + 2 if mech_vent_48h + 4 if cardiac_arrest_48h. Hypotension and Mental status axes excluded (100% missing in this extract, see SECTION 1 comment); Respiratory rate/Heart rate/Systolic BP also excluded (100% missing). NOT a complete/validated Pitt score.",
  "pitt_score_partial_n_components", "Count (0-3) of how many of {temperature, mech_vent_48h, cardiac_arrest_48h} were non-missing for this patient -- low counts mean pitt_score_partial is based on sparse data.",
  "lab_match", "1 if this patient ID matched an 'ID outcome: Conclusion' == MICROBIOLOGICALLY record in CentralRefLab, 0 otherwise.",
  "pathogen1", "Species 1 identified by CentralRefLab, only for lab_match == 1 patients; NA otherwise.",
  "pathogen2", "Species 2 identified by CentralRefLab (polymicrobial/second-sample cases only); NA otherwise.",
  "test_result_date", "From CentralRefLab 'ID outcome: Date of test result' (parsed to datetime), only for lab_match == 1 patients; NA otherwise. Note: at least one source row has a garbled year (e.g. '0025-...' instead of '2025-...').",
  "resistance_any", "1 if resistant (ast_s1_*/ast_s2_* == 1) to ANY of the up-to-56 antibiotics tested for this patient, 0 if tested but never resistant, NA if never tested at all.",
  "mdr_resistance", "Multidrug-resistant (MDR): 1 if resistant to >=1 agent in >=3 antimicrobial categories (organism-agnostic category-count proxy, NOT the organism-specific Magiorakos definition -- see SECTION 2 comment), 0 if tested but <3 categories resistant, NA if never tested.",
  "xdr_resistance", "Extensively drug-resistant (XDR): 1 if mdr_resistance == 1 AND susceptible-only in <=2 of the categories actually tested (i.e. non-susceptible to all-but-<=2 tested categories), 0 otherwise (when tested), NA if never tested. Same organism-agnostic caveat as mdr_resistance.",
  "carb_resistance", "1 if resistant to any carbapenem (any antibiotic name containing 'penem': meropenem, imipenem, ertapenem, meropenem-vaborbactam, etc.), 0 if a carbapenem was tested but not resistant, NA if no carbapenem was ever tested.",
  "methicilin_resistance", "1 if resistant to oxacillin, methicillin, or cefoxitin, 0 if tested but not resistant, NA if none of these three were ever tested.",
  "cephalosporin_resis", "1 if resistant to any cephalosporin (any antibiotic name starting with 'Cef', all generations pooled), 0 if a cephalosporin was tested but not resistant, NA if no cephalosporin was ever tested.",
  "3GC_resistance", "Third-generation cephalosporin resistance specifically: 1 if resistant to ceftriaxone, cefotaxime, ceftazidime, or ceftazidime-avibactam (excludes 2nd-gen cefoxitin/cefuroxime, 4th-gen cefepime, 5th-gen ceftaroline -- narrower than cephalosporin_resis). 0 if a 3GC agent was tested but not resistant, NA if none were tested. Accessed via [[\"3GC_resistance\"]], not $3GC_resistance (invalid R symbol, starts with a digit).",
  "vancomycin_resistance", "1 if resistant to vancomycin specifically (does NOT include teicoplanin or other glycopeptides), 0 if vancomycin was tested but not resistant, NA if vancomycin was never tested.",
  "pdr_resistance", "Pandrug-resistant (PDR, Magiorakos-style, same category-count proxy as mdr_resistance/xdr_resistance): 1 if xdr_resistance == 1 AND resistant in ALL antimicrobial categories tested (0 susceptible), 0 otherwise (when tested), NA if never tested. See n_amr_categories_tested before trusting a PDR=1 value -- a narrow tested panel can qualify spuriously under this proxy.",
  "n_amr_categories_tested", "Count of distinct antimicrobial categories (of the ~20 in SECTION 2's antibiotic_category_map) actually tested for this patient's isolate(s), across all up-to-56 AST slots. Context for mdr_resistance/xdr_resistance/pdr_resistance -- a low count means those flags rest on a narrow panel.",
  "n_amr_categories_resistant", "Count of distinct antimicrobial categories with >=1 resistant result, out of n_amr_categories_tested.",
  "pathogen1_clean", "pathogen1 with nominotypical-subspecies name duplication collapsed (e.g. 'Klebsiella pneumoniae pneumoniae' -> 'Klebsiella pneumoniae'); genuine subspecies/serovars (different trailing word, e.g. 'Klebsiella pneumoniae ozaenae') are untouched. See SECTION 2 comment for the 2 organisms ('Micrococcus kristinae' / 'Rothia kristinae') flagged as a likely naming-convention duplicate but deliberately NOT auto-merged.",
  "pathogen2_clean", "pathogen2, cleaned the same way as pathogen1_clean.",
  "genus1", "First word of pathogen1_clean (the isolate's genus). NA if pathogen1_clean is NA or is the literal 'Other' catch-all.",
  "major_bacterial_group", "genus1 mapped to one of 8 groups (Enterobacterales, Non-fermenting Gram-negative bacilli, Other Gram-negative, Staphylococcus spp., Streptococcus/Enterococcus spp., Other Gram-positive, Fungi, Unspecified/Other) -- see SECTION 2's genus_group_map for the full genus assignment. Used for SECTION 3.6's AMR heatmap group-level rows.",
  "pathogen_domain", "'Bacteria' or 'Fungi', from major_bacterial_group. NA for the 'Other' free-text catch-all (genuinely unknown, not assumed bacterial) and for lab_match == 0 patients. Used to restrict MDR/XDR/PDR (bacteria-only concepts) away from the 3 fungal isolates in SECTION 3.6.",
  "BSI_positive", "SECTION 3.1: 1 if patient_id matched a confirmed CentralRefLab record in lab_full, 0 otherwise. Numerically identical to lab_match; kept as its own explicitly-named variable for the prevalence analysis.",
  "facility_id_clean", "SECTION 3.1: facility_id normalized to Hospital_characteristics's 'HFNNN' format (extracts the HF+digits prefix, re-pads to 3 digits). Recovers 24 of ~175 malformed facility_id values; NA for the rest, including ~151 rows referencing facility codes absent from Hospital_characteristics altogether (mostly 'HF034', 149 rows) or with no recoverable code ('W12'). Used only for the facility-level prevalence aggregation, not elsewhere.",
  "charlson_score", "SECTION 3.1: adapted Charlson comorbidity index built from available comorb_* flags. NOT a complete validated Charlson score -- myocardial infarction and peripheral vascular disease have no corresponding source variable and are omitted entirely. Mutually-exclusive severity pairs (diabetes, liver disease, malignancy) resolved to the higher weight only; missing comorbidity flags contribute 0.",
  "charlson_group", "SECTION 3.1: 'Charlson = 0' vs 'Charlson >0'. NOT a median split as such splits would normally be defined -- the actual median charlson_score is 0 (73% of patients), so a literal <median/>=median split would leave one group empty; this 0-vs->0 split is what that collapses to in practice.",
  "design_weight", "SECTION 3.2: Horvitz-Thompson sampling design weight = 1/(85 * p_i), where p_i is the facility's PPS-with-replacement selection probability (its ApCL / sum(ApCL) across the 119-facility national sampling frame). Joined from NewSamplingFrameDec2022.csv via facility_id_clean; NA for the 7 patients with no resolvable facility_id_clean, who are consequently excluded from the SECTION 3.2 weighted analysis.",
  "NameAU", "SECTION 3.1: patient's facility's oblast/region, from the sampling frame's NameAU field (Osh oblast and Osh city are not distinguished in this field -- both are 'Osh'). NA for the ~7 patients with no resolvable facility_id_clean. Used as the 'oblast' subgroup throughout SECTIONS 3.1-3.4 and 4, and as the calibration margin in SECTION 3.3.",
  "raking_factor", "SECTION 3.3: per-region calibration multiplier = (true 2024 population share of that region) / (current design-weighted share of that region), computed from Demographic Yearbook of the Kyrgyz Republic Table I.1.7.1. Applied to design_weight to correct regional imbalance left over from this single PPS draw (e.g. Issyk-Kul and Naryn were both under-represented, Talas over-represented -- see SECTION 3.3 comment for the verified numbers).",
  "calibrated_weight", "SECTION 3.3: design_weight * raking_factor -- the final, region-calibrated sampling weight. Weighted regional shares built from this column exactly match true 2024 population shares by construction (verified in SECTION 3.3).",
  "resistant_BSI_positive", "SECTION 3.6: population-level flag = 1 if this patient has a confirmed BSI (BSI_positive==1) AND resistance_any==1, 0 if BSI-negative OR BSI-positive-but-susceptible, NA if BSI-positive but AST data was too incomplete to determine resistance_any at all. Denominator for this variable is ALL tested patients, same population as BSI_positive itself.",
  "mdr_BSI_bacterial", "SECTION 3.6: population-level flag = 1 if confirmed bacterial BSI (pathogen_domain=='Bacteria') AND mdr_resistance==1, 0 if BSI-negative OR bacterial-BSI-but-not-MDR, NA if BSI-positive but non-bacterial/unknown-domain (fungal isolates and the 'Other' catch-all are not MDR-evaluable) or if mdr_resistance itself is NA (AST too incomplete).",
  "xdr_BSI_bacterial", "SECTION 3.6: same construction as mdr_BSI_bacterial, using xdr_resistance instead of mdr_resistance.",
  "pdr_BSI_bacterial", "SECTION 3.6: same construction as mdr_BSI_bacterial, using pdr_resistance instead of mdr_resistance -- inherits pdr_resistance's narrow-panel caveat (see that variable's description).",
  "carb_BSI_positive", "SECTION 3.6: population-level flag = 1 if confirmed BSI (BSI_positive==1) AND carb_resistance==1, 0 if BSI-negative OR BSI-positive-but-not-carbapenem-resistant, NA if BSI-positive but carb_resistance itself is NA (no carbapenem ever tested). NOT bacteria-restricted like mdr/xdr/pdr_BSI_bacterial -- carb_resistance is already a narrow, named-class flag.",
  "methicillin_BSI_positive", "SECTION 3.6: same construction as carb_BSI_positive, using methicilin_resistance (oxacillin/methicillin/cefoxitin) instead of carb_resistance.",
  "cephalosporin_BSI_positive", "SECTION 3.6: same construction as carb_BSI_positive, using cephalosporin_resis (any generation) instead of carb_resistance.",
  "c3GC_BSI_positive", "SECTION 3.6: same construction as carb_BSI_positive, using [[\"3GC_resistance\"]] (3rd-generation cephalosporins specifically) instead of carb_resistance.",
  "vancomycin_BSI_positive", "SECTION 3.6: same construction as carb_BSI_positive, using vancomycin_resistance instead of carb_resistance."
)

s1_ast_registry <- bind_rows(
  tibble(
    variable = paste0("abx_s1_", 1:28),
    description = paste0(
      "Species 1 antibiotic #", 1:28, " tested, homologated from CentralRefLab \"S1 - Antibiotic ", 1:28,
      ": Specify the antibiotic\" (parenthetical alternate names stripped, e.g. dropdown option ",
      "'Trimethoprim-sulfamethoxazole (Co-trimoxazole)' -> 'Trimethoprim-sulfamethoxazole'). ",
      "NA if that slot was never used for this isolate, or the slot column doesn't exist in this extract at all (slots 20/23/26/27)."
    )
  ),
  tibble(
    variable = paste0("ast_s1_", 1:28),
    description = paste0(
      "AST result for abx_s1_", 1:28, " (from \"S1 - Antibiotic ", 1:28,
      ": AST results\"). 1 = R or I (resistant/intermediate), 0 = S (susceptible), NA = untested/unknown ('U')/blank, ",
      "or the slot doesn't exist in this extract (slots 20/23/26/27)."
    )
  )
)

s2_ast_registry <- bind_rows(
  tibble(
    variable = paste0("abx_s2_", 1:28),
    description = paste0(
      "Species 2 antibiotic #", 1:28, " tested, homologated from CentralRefLab \"S2 - Antibiotic ", 1:28,
      ": Specify the antibiotic\" (parenthetical alternate names stripped). NA if that slot was never used ",
      "for this isolate, or the slot column doesn't exist in this extract at all (slots 20/23/24/26/27 -- ",
      "a different missing set than Species 1)."
    )
  ),
  tibble(
    variable = paste0("ast_s2_", 1:28),
    description = paste0(
      "AST result for abx_s2_", 1:28, " (from \"S2 - Antibiotic ", 1:28,
      ": AST results\"). 1 = R or I (resistant/intermediate), 0 = S (susceptible), NA = untested/unknown ('U')/blank, ",
      "or the slot doesn't exist in this extract."
    )
  )
)

abx_hist_registry <- tibble(
  variable = paste0("abx_hist_", 1:5),
  description = paste0(
    "Antibiotic #", 1:5, " from 'History of antibiotic therapy NO: Indicate the drug (INN) ", 1:5,
    "', name homologated to a canonical English drug name from free-text source; 'UNMAPPED: <value>' where not confidently recognized."
  )
)

integrated_registry <- integrated_map %>%
  distinct(new, group) %>%
  mutate(description = paste0(
    "Integrated flag: 1 if TRUE in ANY of the age-domain source columns for '", group,
    "' (Adult/Paediatric/Neonate, as applicable -- a patient only has one age domain populated), 0 otherwise."
  )) %>%
  select(variable = new, description)

agespecific_registry <- tibble(
  variable = agespecific_map$new,
  description = paste0(
    "Age-specific flag from source column \"", agespecific_map$orig,
    "\". 1 = flagged true in source, 0 = not flagged/blank (source has no explicit false value)."
  )
)

yesno_registry <- tibble(
  variable = yesno_map$new,
  description = paste0("From source column \"", yesno_map$orig, "\". 1 = YES, 0 = NO, NA = UNK/blank in source.")
)

device_registry <- tibble(
  variable = device_map$new,
  description = paste0(
    "From source column \"", device_map$orig,
    "\". 1 = flagged true in source, 0 = not flagged/blank."
  )
)

reason_registry <- tibble(
  variable = reason_map$new,
  description = paste0(
    "Primary reason for admission flag from source column \"", reason_map$orig,
    "\". 1 = flagged true in source, 0 = not flagged/blank."
  )
)

site_registry <- tibble(
  variable = site_map$new,
  description = paste0(
    "Clinically suspected primary site flag from source column \"", site_map$orig,
    "\". 1 = flagged true in source, 0 = not flagged/blank."
  )
)

comorb_registry <- tibble(
  variable = comorb_map$new,
  description = paste0(
    "Comorbidity flag from source column \"", comorb_map$orig,
    "\". 1 = YES, 0 = NO, NA = UNK/blank in source."
  )
)

riskfx_registry <- tibble(
  variable = riskfx_map$new,
  description = paste0(
    "Risk factor flag from source column \"", riskfx_map$orig,
    "\". 1 = YES, 0 = NO, NA = UNK/blank in source."
  )
)

occur48h_registry <- tibble(
  variable = occur48h_map$new,
  description = paste0(
    "Occurrence-in-48h flag from source column \"", occur48h_map$orig,
    "\". 1 = YES, 0 = NO, NA = UNK/blank in source."
  )
)

sample_id_registry <- tibble(
  variable = sample_id_map$new,
  description = paste0("From source column \"", sample_id_map$orig, "\".")
)

date_registry <- tibble(
  variable = date_map$new,
  description = paste0("Renamed and parsed to datetime from source column \"", date_map$orig, "\".")
)

variable_registry <- bind_rows(
  manual_registry, abx_hist_registry, s1_ast_registry, s2_ast_registry, integrated_registry, agespecific_registry,
  yesno_registry, device_registry, reason_registry, site_registry, comorb_registry,
  riskfx_registry, occur48h_registry, sample_id_registry, date_registry
) %>%
  # manual_registry gives a few variables (blood_sample*_ID/date,
  # bsi_onset_datetime, hospital_admissionDate*) a more detailed description
  # than sample_id_registry/date_registry auto-generate for the same names;
  # distinct() keeps the FIRST (manual_registry's, since it's listed first
  # above) and drops the generic duplicate, so the codebook join can't fan out.
  distinct(variable, .keep_all = TRUE)

# Summarize the actual values found in a column, for the codebook's
# "values" column: numeric -> observed range; everything else -> up to 12
# distinct non-missing values (more than that -> range + count).
summarize_values <- function(x, max_show = 12) {
  x_nonNA <- x[!is.na(x)]
  if (length(x_nonNA) == 0) return("(all missing)")
  if (is.numeric(x)) {
    return(paste0(round(min(x_nonNA), 2), " to ", round(max(x_nonNA), 2)))
  }
  vals <- sort(unique(x_nonNA))
  if (length(vals) <= max_show) {
    paste(vals, collapse = "; ")
  } else {
    paste0(paste(vals[seq_len(max_show)], collapse = "; "), "; ... (", length(vals), " distinct values total)")
  }
}

build_codebook <- function(df, registry) {
  tibble(
    variable    = names(df),
    r_class     = map_chr(df, ~ paste(class(.x), collapse = "/")),
    n_distinct  = map_int(df, ~ n_distinct(.x, na.rm = TRUE)),
    pct_missing = round(map_dbl(df, ~ mean(is.na(.x))) * 100, 1),
    values      = map_chr(df, summarize_values)
  ) %>%
    left_join(registry, by = "variable") %>%
    mutate(description = if_else(
      is.na(description),
      "(no description recorded -- add a row to variable_registry in SECTION 7)",
      description
    )) %>%
    select(variable, description, r_class, n_distinct, pct_missing, values)
}

Codebook_Main_data_CaseReport <- build_codebook(Main_data_CaseReport, variable_registry)
write_xlsx(Codebook_Main_data_CaseReport, "Codebook_Main_data_CaseReport.xlsx")

# Re-run this SECTION 7 block any time after adding new variables to
# Main_data_CaseReport -- it will regenerate the codebook automatically.

# =============================================================================
# SECTION 16: DESCRIPTIVE STATISTICS ACROSS ANALYTIC SAMPLES (Supplementary Material)
# =============================================================================
# Plain descriptive characterisation (no test, no model) of the CRF/baseline
# variable universe already used for LASSO/RF variable reduction (SECTION 14's
# candidate_vars_core/candidate_vars_flags/candidate_vars_discharge), split
# across the 5 patient groups repeatedly compared throughout this report --
# with each variable's own missingness rate reported per group, never
# imputed or hidden.
#
# DELIBERATELY NOT built from Main_data_Combined: SECTION 8 restricts that
# object to patients with EITHER a Discharge OR a Followup28 record
# (outcome_patient_ids, SECTION 8's own filter just before its cleanup rm())
# -- fine for the mortality/LOS analyses that need linked outcome data, but
# it would silently shrink the BSI-negative group from ~8,452 to the 18
# (0.2%) with any linked outcome at all (SECTION 13.2's own documented
# linkage-gap finding), which is not what a baseline/CRF descriptive table
# is for. Built instead directly from the two CRF-collection sources
# (Main_data_CaseReport, Main_data_hospitalCohortEnrolment_clean), each
# already carrying the same 57 harmonised column names (SECTION 7/8) and
# each de-duplicated by patient_id the same way SECTION 8 itself does
# (dedupe_by_patient(), still in scope -- the specific discharge_dedup/
# casereport_dedup/cohort_dedup OBJECTS were rm()'d by SECTION 8's own
# cleanup, so re-applied fresh here rather than re-derived as standing
# objects nothing else needs).
#
# Group definitions (mutually exclusive by construction, EXCEPT that a
# confirmed-BSI patient with an indeterminate resistance result -- AST too
# incomplete to call -- falls into BSI-positive but neither the BSI+AMR nor
# BSI+Susceptible subgroup, exactly as everywhere else this script stratifies
# by resistant_BSI_positive):
#   BSI-negative (culture-negative)             : BSI_positive == 0
#   BSI-positive (confirmed BSI)                : BSI_positive == 1
#   BSI-positive, resistant (BSI+AMR)           : BSI_positive == 1 & resistant_BSI_positive == 1
#   BSI-positive, susceptible (BSI+Susceptible) : BSI_positive == 1 & resistant_BSI_positive == 0
#   Uninfected cohort                            : matched_control_broad == 1 (SECTION 9's
#                                                   broad IPW population -- every genuinely
#                                                   enrolled uninfected participant, not just
#                                                   the strict 1:1 pairs)
# These use the unchanged, protocol-defined BSI_positive/resistant_BSI_positive
# throughout -- unaffected by the CoNS-treatment sensitivity analyses
# (SECTION 2.6/throughout this script), consistent with every other primary
# analysis in this report.
#
# ICU_entry/any_complication/complications_count are DISCHARGE-form fields
# (SECTION 14's own varred_base joins them the same way) -- present for
# whichever CaseReport/CohortEnrolment patients happen to also have a linked
# Discharge record (roughly half of each, per SECTION 8), genuinely NA
# (never imputed) for the rest; that partial linkage is exactly what the
# missingness column below is for.
desc_base <- bind_rows(
  dedupe_by_patient(Main_data_CaseReport) %>%
    select(patient_id, BSI_positive, resistant_BSI_positive, any_of(candidate_vars_core), any_of(candidate_vars_flags)) %>%
    mutate(matched_control_broad = 0L),
  dedupe_by_patient(Main_data_hospitalCohortEnrolment_clean) %>%
    select(patient_id, matched_control_broad, any_of(candidate_vars_core), any_of(candidate_vars_flags)) %>%
    mutate(BSI_positive = NA_integer_, resistant_BSI_positive = NA_integer_)
) %>%
  left_join(dedupe_by_patient(Main_data_hospitalDischarge_clean) %>%
              select(patient_id, ICU_entry, any_complication, complications_count),
            by = "patient_id")
stopifnot(all(c(candidate_vars_core, candidate_vars_flags) %in% names(desc_base)))

desc_groups <- list(
  `BSI-negative (culture-negative)` = which(desc_base$BSI_positive == 0),
  `BSI-positive (confirmed BSI)` = which(desc_base$BSI_positive == 1),
  `BSI-positive, resistant (BSI+AMR)` = which(desc_base$BSI_positive == 1 & desc_base$resistant_BSI_positive == 1),
  `BSI-positive, susceptible (BSI+Susceptible)` = which(desc_base$BSI_positive == 1 & desc_base$resistant_BSI_positive == 0),
  `Uninfected cohort` = which(desc_base$matched_control_broad == 1)
)
group_n_summary <- tibble(group = names(desc_groups), n = purrr::map_int(desc_groups, length))

desc_continuous_vars <- intersect(c("age_years", "charlson_score", "pitt_score_partial", "complications_count"), names(desc_base))
desc_categorical_vars <- setdiff(c(candidate_vars_core, candidate_vars_flags, "ICU_entry", "any_complication"), desc_continuous_vars)

desc_cont_row <- function(varname, group_name, idx) {
  x <- desc_base[[varname]][idx]
  n_total <- length(idx); n_nonmiss <- sum(!is.na(x))
  tibble(
    variable = varname, group = group_name, n_total = n_total, n_nonmissing = n_nonmiss,
    pct_missing = round(100 * (n_total - n_nonmiss) / n_total, 1),
    mean_sd = if (n_nonmiss > 0) sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)) else NA_character_,
    median_iqr = if (n_nonmiss > 0) sprintf("%.1f [%.1f-%.1f]", median(x, na.rm = TRUE), quantile(x, 0.25, na.rm = TRUE), quantile(x, 0.75, na.rm = TRUE)) else NA_character_
  )
}

# Binary 0/1 CRF flags (comorb_/riskfx_/reason_/ICU_entry/any_complication/
# abx_current_therapy, etc.) collapse to a single "1 (present)" row per
# variable; genuine multi-level factors (sex_at_birth, ward_type_icu,
# specialty_group, infection_origin_composite) get one row per observed
# level -- mirrors SECTION 1's own sec1_cat_rows_confirmed() convention.
desc_cat_rows <- function(varname, group_name, idx) {
  x <- desc_base[[varname]][idx]
  n_total <- length(idx); n_nonmiss <- sum(!is.na(x)); pct_missing <- round(100 * (n_total - n_nonmiss) / n_total, 1)
  if (n_nonmiss == 0) {
    return(tibble(variable = varname, group = group_name, level = NA_character_,
                  n_total = n_total, n_nonmissing = 0L, pct_missing = pct_missing,
                  n_level = NA_integer_, pct_level = NA_real_))
  }
  vals <- sort(unique(x))
  if (setequal(vals, c(0, 1)) || length(vals) == 1) {
    tibble(variable = varname, group = group_name, level = "1 (present)",
           n_total = n_total, n_nonmissing = n_nonmiss, pct_missing = pct_missing,
           n_level = sum(x == 1, na.rm = TRUE), pct_level = round(100 * mean(x == 1, na.rm = TRUE), 1))
  } else {
    purrr::map_dfr(as.character(sort(unique(x))), function(l) {
      tibble(variable = varname, group = group_name, level = l,
             n_total = n_total, n_nonmissing = n_nonmiss, pct_missing = pct_missing,
             n_level = sum(as.character(x) == l, na.rm = TRUE),
             pct_level = round(100 * mean(as.character(x) == l, na.rm = TRUE), 1))
    })
  }
}

descriptive_stats_continuous <- purrr::map_dfr(names(desc_groups), function(g) {
  purrr::map_dfr(desc_continuous_vars, desc_cont_row, group_name = g, idx = desc_groups[[g]])
})
descriptive_stats_categorical <- purrr::map_dfr(names(desc_groups), function(g) {
  purrr::map_dfr(desc_categorical_vars, desc_cat_rows, group_name = g, idx = desc_groups[[g]])
})

# Wide layout for readability -- one row per variable(/level), one column
# triple (summary + missing%) per group, mirroring the compact wide tables
# used throughout SECTIONS 3-14's manuscript-facing outputs.
descriptive_stats_continuous_wide <- descriptive_stats_continuous %>%
  select(variable, group, mean_sd, median_iqr, pct_missing) %>%
  tidyr::pivot_wider(names_from = group, values_from = c(mean_sd, median_iqr, pct_missing))
descriptive_stats_categorical_wide <- descriptive_stats_categorical %>%
  select(variable, level, group, n_level, pct_level, pct_missing) %>%
  tidyr::pivot_wider(names_from = group, values_from = c(n_level, pct_level, pct_missing))

write_xlsx(
  list(
    group_sizes = group_n_summary,
    continuous_long = descriptive_stats_continuous,
    continuous_wide = descriptive_stats_continuous_wide,
    categorical_long = descriptive_stats_categorical,
    categorical_wide = descriptive_stats_categorical_wide
  ),
  "Descriptive_statistics_by_analytic_sample.xlsx"
)
save.image("ws_final.RData")
