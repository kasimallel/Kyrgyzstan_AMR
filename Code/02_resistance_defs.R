# 02_resistance_defs.R -- AST cleaning, drug/class maps (heatmap template), AWaRe, Magiorakos MDR/XDR/PDR, isolate flags
source(file.path(Sys.getenv("KGZ_PKG"), "R", "config.R"))
x <- readRDS(file.path(OUT, "01_pathogens_ast.rds"))
iso <- x$pathogen_rows %>% filter(!is.na(sample_used)) %>%
  transmute(record_id = `Record ID`, case_rid = `Case report record ID`, facility_raw = `Unique facility ID (from case report)`,
            ward_raw = `Unique ward ID (from case report)`, patient_raw = `Unique survey patient ID`,
            species, pathogen_group, name_source, sample_used, polymicrobial, comments = `Comments: Comments`)
a <- x$ast_long

# ---- 1. drug dictionary: name -> code, AWaRe --------------------------------------------------
drug <- tribble(
 ~pattern, ~code, ~drug, ~aware,
 "^Amikacin","AMK","Amikacin","Access", "^Gentamicin","GEN","Gentamicin","Access", "^Netilmicin","NET","Netilmicin","Access",
 "^Tobramycin","TOB","Tobramycin","Access", "^Streptomycin","STR","Streptomycin","Watch",
 "^Ampicillin-sulbactam","SAM","Ampicillin-sulbactam","Access", "^Ampicillin","AMP","Ampicillin","Access",
 "^Amoxicillin-clav","AMC","Amoxicillin-clavulanic acid","Access", "^Benzylpenicillin","PEN","Benzylpenicillin","Access",
 "^Oxacillin","OXA","Oxacillin","Access", "^Piperacillin","TZP","Piperacillin-tazobactam","Watch",
 "^Cefoxitin","FOX","Cefoxitin","Watch", "^Cefuroxime","CXM","Cefuroxime","Watch", "^Cefotaxime","CTX","Cefotaxime","Watch",
 "^Ceftriaxone","CRO","Ceftriaxone","Watch", "^Ceftazidime-avi","CZA","Ceftazidime-avibactam","Reserve", "^Ceftazidime","CAZ","Ceftazidime","Watch",
 "^Cefepime","FEP","Cefepime","Watch", "^Ceftaroline","CPT","Ceftaroline","Reserve",
 "^Ertapenem","ETP","Ertapenem","Watch", "^Meropenem-vabo","MEV","Meropenem-vaborbactam","Reserve", "^Meropenem","MEM","Meropenem","Watch", "^Imipenem","IPM","Imipenem","Watch",
 "^Ciprofloxacin","CIP","Ciprofloxacin","Watch", "^Levofloxacin","LVX","Levofloxacin","Watch", "^Moxifloxacin","MXF","Moxifloxacin","Watch", "^Nalidixic","NAL","Nalidixic acid","Watch",
 "^Trimethoprim-sulf","SXT","Trimethoprim-sulfamethoxazole","Access", "^Colistin","COL","Colistin","Reserve",
 "^Vancomycin","VAN","Vancomycin","Watch", "^Teicoplanin","TEC","Teicoplanin","Watch", "^Daptomycin","DAP","Daptomycin","Reserve", "^Linezolid","LZD","Linezolid","Reserve",
 "^Tetracycline","TET","Tetracycline","Access", "^Tigecycline","TGC","Tigecycline","Reserve", "^Clindamycin","CLI","Clindamycin","Access",
 "^Erythromycin","ERY","Erythromycin","Watch", "^Rifampicin","RIF","Rifampicin","Watch", "^Chloramphenicol","CHL","Chloramphenicol","Access",
 "^Fusidic","FUS","Fusidic acid","Watch", "^Fosfomycin","FOS","Fosfomycin","Reserve", "^Nitrofurantoin","NIT","Nitrofurantoin","Access")
code_of <- function(nm) { out <- rep(NA_character_, length(nm)); for (i in seq_len(nrow(drug))) { h <- is.na(out) & grepl(drug$pattern[i], nm, ignore.case = TRUE); out[h] <- drug$code[i] }; out }
a$agent <- code_of(a$antibiotic_raw)
stopifnot(!any(is.na(a$agent[!is.na(a$antibiotic_raw)])))

# ---- 2. result vocabulary ---------------------------------------------------------------------
# S: Susceptible, S, S*, SÜ (encoding variants S†)  | I: Intermediate, I*, I (non-PNA) | R: Resistant, R, HLR+
# NA (not interpretable): No BP, IE, UTI only, HLR?, Uninterpretable result, missing; "Check FOX" handled via MRSA inference below
a$res <- case_when(
  a$result_raw == "Susceptible" | grepl("^S($|[^a-z])", a$result_raw) ~ "S",
  a$result_raw %in% c("Intermediate","I*","I (non-PNA)") ~ "R",   # PLAN DECISION: I and I* are counted as resistant (non-susceptible) throughout
  a$result_raw %in% c("Resistant","R","HLR+") ~ "R",
  TRUE ~ NA_character_)
a$mic <- suppressWarnings(as.numeric(sub(",", ".", gsub("[^0-9,.]", "", a$value_raw))))
a$res_note <- case_when(a$result_raw == "Check FOX" ~ "check_fox", a$result_raw %in% c("No BP","IE","UTI only","HLR?","Uninterpretable result") ~ "not_interpretable", is.na(a$result_raw) ~ "missing", TRUE ~ NA_character_)

# ---- 3. organism gram / intrinsic & non-interpretable agents (EUCAST expected-resistance logic) --
gram_neg <- c("Escherichia coli","Klebsiella pneumoniae","Acinetobacter baumannii","Pseudomonas aeruginosa","Pseudomonas stutzeri","Enterobacter cloacae","Salmonella Typhi","Salmonella spp. (non-Typhi)","Proteus mirabilis","Citrobacter youngae","Pantoea agglomerans","Pantoea sp.","Enterobacter sp.","Stenotrophomonas maltophilia","Sphingomonas paucimobilis","Moraxella catarrhalis","Aeromonas salmonicida","Agrobacterium tumefaciens","Chromobacterium violaceum","Rodentibacter pneumotropicus")
gneg_excl <- c("PEN","OXA","VAN","TEC","LZD","DAP","CLI","ERY","FUS","RIF","CPT","NIT")           # NIT: urine-only (Magiorakos criterion ii)
gpos_excl <- c("COL","CAZ","NAL","NIT","AMK","NET","TOB","MEV","CZA")
org_excl <- list(
  "Klebsiella pneumoniae" = c("AMP"),
  "Enterobacter cloacae" = c("AMP","AMC","SAM","FOX","CXM"),
  "Proteus mirabilis" = c("COL","TGC","TET"),
  "Pseudomonas aeruginosa" = c("ETP","CTX","CRO","GEN","SXT","AMP","AMC","SAM","CHL","TET","TGC","FOS"),
  "Pseudomonas stutzeri" = c("ETP","CTX","CRO","GEN","SXT","AMP","AMC","SAM","CHL","TET","TGC","FOS"),
  "Acinetobacter baumannii" = c("AMP","AMC","ETP","CTX","CRO","CAZ","FOX","CXM","TGC","FOS","SXT_none"),   # EUCAST: no BP / intrinsic (template dagger)
  "Stenotrophomonas maltophilia" = c("AMP","AMC","SAM","CTX","CRO","CAZ","FEP","ETP","MEM","IPM","TZP","AMK","GEN","NET","TOB","COL","FOS"),
  "Enterococcus faecalis" = c("CTX","CRO","CAZ","FEP","SXT","CLI","COL","CXM","FOX","GEN_low","AMK","NET","TOB","LVX","CIP","TGC_none"),
  "Enterococcus faecium"  = c("CTX","CRO","CAZ","FEP","SXT","CLI","COL","CXM","FOX","AMK","NET","TOB","LVX","CIP"),
  "Enterococcus gallinarum" = c("VAN","CTX","CRO","CAZ","FEP","SXT","CLI","COL","AMK","NET","TOB"),
  "Viridans group streptococci" = c("GEN","AMK","NET","TOB","COL","SXT_none"),
  "Streptococcus pyogenes" = c("GEN","AMK","NET","TOB","COL"),
  "Streptococcus agalactiae" = c("GEN","AMK","NET","TOB","COL"),
  "Streptococcus pneumoniae" = c("GEN","AMK","NET","TOB","COL"))
interpretable <- function(species, grp, agent) {
  g <- if (species %in% gram_neg) gneg_excl else gpos_excl
  ex <- unique(c(g, org_excl[[species]], org_excl[[grp]]))
  !(agent %in% ex)
}
a <- a %>% left_join(iso %>% select(record_id, species, pathogen_group, comments), by = "record_id")
a$interp <- mapply(interpretable, a$species, a$pathogen_group, a$agent)
# enterococci gentamicin / streptomycin: only high-level (HLR+) is a valid result; plain S/I/R from low-level test not applicable
a <- a %>% mutate(res_final = ifelse(interp & !is.na(res), res, NA_character_))
# duplicate agent on one isolate -> worst result (R > I > S)
rank_res <- c(S = 1, I = 2, R = 3)
ast <- a %>% filter(!is.na(agent)) %>% group_by(record_id, agent) %>%
  summarise(res_final = { v <- res_final[!is.na(res_final)]; if (!length(v)) NA_character_ else names(rank_res)[max(rank_res[v])] },
            check_fox = any(res_note == "check_fox", na.rm = TRUE), .groups = "drop")
cat("Isolate-agent pairs:", nrow(ast), " interpretable:", sum(!is.na(ast$res_final)), "\n")

# ---- 4. MRSA inference (S. aureus) --------------------------------------------------------------
sa <- iso %>% filter(pathogen_group == "Staphylococcus aureus")
# Oxacillin reported as "Check FOX" carries the oxacillin MIC in the value column: MIC > 2 mg/L => MRSA (EUCAST/CLSI S. aureus oxacillin breakpoint)
oxa_mic <- a %>% filter(agent == "OXA", res_note %in% "check_fox") %>% group_by(record_id) %>% summarise(mic = max(mic, na.rm = TRUE), .groups = "drop") %>% mutate(mic = ifelse(is.infinite(mic), NA, mic))
oxa <- ast %>% filter(agent %in% c("OXA","FOX")) %>% group_by(record_id) %>%
  summarise(oxa_fox = { v <- res_final[!is.na(res_final)]; if (!length(v)) NA_character_ else if (any(v=="R")) "R" else "S" },
            check_fox = any(check_fox), .groups = "drop") %>% left_join(oxa_mic, by = "record_id") %>%
  mutate(oxa_fox = ifelse(is.na(oxa_fox) & check_fox & !is.na(mic), ifelse(mic > 2, "R", "S"), oxa_fox),
         mic_call = ifelse(check_fox & !is.na(mic), ifelse(mic > 2, "MRSA", "MSSA"), NA))
cm <- toupper(iso$comments)
iso$cm_mrsa_pos <- grepl("MRS|CEFOXITIN[- ]?POS|FOX[- ]?POS", cm)
iso$cm_mrsa_neg <- grepl("CEFOXITIN[- ]?NEG|SCREENING NEGATIVE|MRSA[- ]?NEG|MSSA", cm)
iso <- iso %>% left_join(oxa, by = "record_id") %>%
  mutate(mrsa = case_when(pathogen_group != "Staphylococcus aureus" ~ NA_character_,
                          oxa_fox == "R" ~ "MRSA", cm_mrsa_pos ~ "MRSA", oxa_fox == "S" ~ "MSSA", cm_mrsa_neg ~ "MSSA", TRUE ~ "unknown"))
cat("\nMRSA status (S. aureus):\n"); print(table(iso$mrsa[iso$pathogen_group=="Staphylococcus aureus"], useNA="ifany"))
cat("MIC-based call vs comment evidence (Check FOX isolates):\n"); print(table(mic_call = iso$mic_call, comment_pos = iso$cm_mrsa_pos, useNA = "ifany"))

# push the inferred methicillin status back into the AST table as an oxacillin result (so heatmap, drug-bug and MDR use one consistent source)
sa_status <- iso %>% filter(mrsa %in% c("MRSA","MSSA")) %>% transmute(record_id, ox_inf = ifelse(mrsa == "MRSA", "R", "S"))
ast <- ast %>% left_join(sa_status, by = "record_id") %>%
  mutate(res_final = ifelse(agent == "OXA" & is.na(res_final) & !is.na(ox_inf), ox_inf, res_final)) %>% select(-ox_inf)
need <- sa_status %>% filter(!record_id %in% ast$record_id[ast$agent == "OXA"]) %>% transmute(record_id, agent = "OXA", res_final = ox_inf, check_fox = FALSE)
ast <- bind_rows(ast, need)
cat("S. aureus oxacillin/methicillin results after inference:", sum(ast$agent == "OXA" & !is.na(ast$res_final)), "\n")

# ---- 5. heatmap template cells (rows x class -> agents) ------------------------------------------
cls <- c("Penicillins","Anti-staph β-lactams (MR screen)","Penicillin + BLI","Cephalosporins 3rd gen","Cephalosporins 4th gen","Carbapenems","Aminoglycosides","Fluoroquinolones","Folate pathway inhibitors","Polymyxins","Glycopeptides","Lipopeptides","Oxazolidinones","Tetracyclines","Lincosamides","Macrolides","Rifamycins","Phenicols")
class_of <- c(AMP="Penicillins",PEN="Penicillins",FOX="Anti-staph β-lactams (MR screen)",OXA="Anti-staph β-lactams (MR screen)",AMC="Penicillin + BLI",SAM="Penicillin + BLI",TZP="Penicillin + BLI",
  CTX="Cephalosporins 3rd gen",CAZ="Cephalosporins 3rd gen",CRO="Cephalosporins 3rd gen",FEP="Cephalosporins 4th gen",ETP="Carbapenems",MEM="Carbapenems",IPM="Carbapenems",
  AMK="Aminoglycosides",GEN="Aminoglycosides",NET="Aminoglycosides",TOB="Aminoglycosides",CIP="Fluoroquinolones",LVX="Fluoroquinolones",MXF="Fluoroquinolones",SXT="Folate pathway inhibitors",
  COL="Polymyxins",VAN="Glycopeptides",TEC="Glycopeptides",DAP="Lipopeptides",LZD="Oxazolidinones",TET="Tetracyclines",CLI="Lincosamides",ERY="Macrolides",RIF="Rifamycins",CHL="Phenicols")
H <- function(...) { l <- list(...); stack(l) }
tpl <- list(
 "Escherichia coli" = list(`Penicillins`="AMP",`Penicillin + BLI`="AMC",`Cephalosporins 3rd gen`=c("CTX","CAZ"),`Cephalosporins 4th gen`="FEP",Carbapenems=c("ETP","MEM"),Aminoglycosides="AMK",Fluoroquinolones="CIP",`Folate pathway inhibitors`="SXT",Polymyxins="COL"),
 "Klebsiella pneumoniae" = list(`Penicillin + BLI`="AMC",`Cephalosporins 3rd gen`=c("CTX","CAZ"),`Cephalosporins 4th gen`="FEP",Carbapenems=c("ETP","MEM"),Aminoglycosides="AMK",Fluoroquinolones="CIP",`Folate pathway inhibitors`="SXT",Polymyxins="COL"),
 "Acinetobacter baumannii" = list(`Cephalosporins 3rd gen`=c("CTX","CAZ"),`Cephalosporins 4th gen`="FEP",Carbapenems="MEM",Aminoglycosides="AMK",Fluoroquinolones="CIP",`Folate pathway inhibitors`="SXT",Polymyxins="COL"),
 "Pseudomonas spp." = list(`Cephalosporins 3rd gen`="CAZ",`Cephalosporins 4th gen`="FEP",Carbapenems="MEM",Aminoglycosides="AMK",Fluoroquinolones="CIP",Polymyxins="COL"),
 "Staphylococcus aureus" = list(Penicillins="PEN",`Anti-staph β-lactams (MR screen)`=c("FOX","OXA"),Aminoglycosides="GEN",Fluoroquinolones="LVX",Glycopeptides=c("VAN","TEC"),Lipopeptides="DAP",Oxazolidinones="LZD",Tetracyclines="TET",Lincosamides="CLI",Macrolides="ERY",Rifamycins="RIF"),
 "Other Staphylococcus spp. (pathogen-confirmed)" = list(`Anti-staph β-lactams (MR screen)`=c("FOX","OXA"),Aminoglycosides="GEN",Fluoroquinolones="LVX",Glycopeptides=c("VAN","TEC"),Lipopeptides="DAP",Oxazolidinones="LZD",Tetracyclines="TET",Lincosamides="CLI",Macrolides="ERY",Rifamycins="RIF"),
 "Enterococcus faecalis" = list(Aminoglycosides="GEN",Fluoroquinolones="LVX",Glycopeptides=c("VAN","TEC"),Oxazolidinones="LZD"),
 "Enterococcus faecium" = list(Aminoglycosides="GEN",Fluoroquinolones="LVX",Glycopeptides=c("VAN","TEC"),Oxazolidinones="LZD"),
 "Viridans group streptococci" = list(Penicillins=c("PEN","AMP"),Aminoglycosides="GEN",Fluoroquinolones="LVX",`Folate pathway inhibitors`="SXT",Lincosamides="CLI",Macrolides="ERY",Rifamycins="RIF"),
 "Streptococcus pyogenes" = list(Penicillins=c("PEN","AMP"),Aminoglycosides="GEN",Fluoroquinolones="LVX",`Folate pathway inhibitors`="SXT",Lincosamides="CLI",Macrolides="ERY",Rifamycins="RIF"))
# template dagger cells (tested per protocol but not interpretable): enterococci FQ (+ GEN low-level), VGS/pyogenes GEN, A. baumannii CTX/CAZ
dagger <- list("Enterococcus faecalis"=c("Fluoroquinolones"), "Enterococcus faecium"=c("Fluoroquinolones"), "Viridans group streptococci"="Aminoglycosides","Streptococcus pyogenes"="Aminoglycosides","Acinetobacter baumannii"="Cephalosporins 3rd gen")
hm_row <- function(sp, grp) case_when(
  grp == "Escherichia coli" ~ "Escherichia coli", grp == "Klebsiella pneumoniae" ~ "Klebsiella pneumoniae", grp == "Acinetobacter baumannii" ~ "Acinetobacter baumannii",
  grp %in% c("Pseudomonas aeruginosa","Pseudomonas stutzeri") ~ "Pseudomonas spp.", grp == "Staphylococcus aureus" ~ "Staphylococcus aureus",
  grp %in% c("Staphylococcus pseudintermedius","Staphylococcus sp.") ~ "Other Staphylococcus spp. (pathogen-confirmed)",
  grp %in% c("Enterococcus faecalis","Enterococcus faecium","Viridans group streptococci","Streptococcus pyogenes") ~ grp, TRUE ~ NA_character_)
iso$hm_row <- hm_row(iso$species, iso$pathogen_group)

# ---- 6. Magiorakos Tables 1-5 (categories -> agents) ------------------------------------------------
mag <- list(
 "Staphylococcus aureus" = list(Aminoglycosides="GEN",Ansamycins="RIF",`Anti-MRSA cephalosporins`="CPT",`Anti-staphylococcal beta-lactams`=c("OXA","FOX"),Fluoroquinolones=c("CIP","MXF"),`Folate pathway inhibitors`="SXT",Fucidanes="FUS",Glycopeptides=c("VAN","TEC"),Glycylcyclines="TGC",Lincosamides="CLI",Lipopeptides="DAP",Macrolides="ERY",Oxazolidinones="LZD",Phenicols="CHL",`Phosphonic acids`="FOS",Tetracyclines="TET"),
 "Enterococcus" = list(`Aminoglycosides (HL gentamicin)`="GEN",`Streptomycin (HL)`="STR",Carbapenems=c("IPM","MEM"),Fluoroquinolones=c("CIP","LVX","MXF"),Glycopeptides=c("VAN","TEC"),Glycylcyclines="TGC",Lipopeptides="DAP",Oxazolidinones="LZD",Penicillins="AMP",Tetracyclines="TET"),
 "Enterobacterales" = list(Aminoglycosides=c("GEN","TOB","AMK","NET"),`Anti-MRSA cephalosporins`="CPT",`Antipseudomonal penicillins + BLI`="TZP",Carbapenems=c("ETP","IPM","MEM"),`1st/2nd gen cephalosporins`="CXM",`Extended-spectrum cephalosporins`=c("CTX","CRO","CAZ","FEP"),Cephamycins="FOX",Fluoroquinolones="CIP",`Folate pathway inhibitors`="SXT",Glycylcyclines="TGC",Penicillins="AMP",`Penicillins + BLI`=c("AMC","SAM"),Phenicols="CHL",`Phosphonic acids`="FOS",Polymyxins="COL",Tetracyclines="TET"),
 "Pseudomonas aeruginosa" = list(Aminoglycosides=c("GEN","TOB","AMK","NET"),`Antipseudomonal carbapenems`=c("IPM","MEM"),`Antipseudomonal cephalosporins`=c("CAZ","FEP"),`Antipseudomonal fluoroquinolones`=c("CIP","LVX"),`Antipseudomonal penicillins + BLI`="TZP",`Phosphonic acids`="FOS",Polymyxins="COL"),
 "Acinetobacter" = list(Aminoglycosides=c("GEN","TOB","AMK","NET"),`Antipseudomonal carbapenems`=c("IPM","MEM"),`Antipseudomonal fluoroquinolones`=c("CIP","LVX"),`Antipseudomonal penicillins + BLI`="TZP",`Extended-spectrum cephalosporins`=c("CTX","CRO","CAZ","FEP"),`Folate pathway inhibitors`="SXT",`Penicillins + BLI`="SAM",Polymyxins="COL",Tetracyclines="TET"))
mag_table <- function(sp, grp) case_when(
  grp == "Staphylococcus aureus" ~ "Staphylococcus aureus", grepl("^Enterococcus", grp) ~ "Enterococcus",
  grp %in% c("Escherichia coli","Klebsiella pneumoniae","Enterobacter cloacae","Proteus mirabilis","Citrobacter youngae","Salmonella Typhi","Salmonella spp. (non-Typhi)","Pantoea agglomerans","Pantoea sp.","Enterobacter sp.") ~ "Enterobacterales",
  grp == "Pseudomonas aeruginosa" ~ "Pseudomonas aeruginosa", grp == "Acinetobacter baumannii" ~ "Acinetobacter", TRUE ~ NA_character_)
iso$mag_table <- mag_table(iso$species, iso$pathogen_group)
# Magiorakos intrinsic categories (Tables 2-3) removed before classification
mag_intr <- list("Klebsiella pneumoniae"=c("AMP"),"Enterobacter cloacae"=c("AMP","AMC","SAM","FOX","CXM"),"Proteus mirabilis"=c("COL","TGC","TET"),
                 "Enterococcus faecium"=c("IPM","MEM"),"Enterococcus faecalis"=character(0))

gen_cat <- c(AMP="Penicillins",PEN="Penicillins",OXA="Anti-staphylococcal beta-lactams",FOX="Anti-staphylococcal beta-lactams",AMC="Penicillins + BLI",SAM="Penicillins + BLI",TZP="Penicillins + BLI",CZA="Penicillins + BLI",
  CXM="Cephalosporins",CTX="Cephalosporins",CRO="Cephalosporins",CAZ="Cephalosporins",FEP="Cephalosporins",CPT="Anti-MRSA cephalosporins",ETP="Carbapenems",MEM="Carbapenems",IPM="Carbapenems",MEV="Carbapenems",
  AMK="Aminoglycosides",GEN="Aminoglycosides",NET="Aminoglycosides",TOB="Aminoglycosides",STR="Aminoglycosides",CIP="Fluoroquinolones",LVX="Fluoroquinolones",MXF="Fluoroquinolones",NAL="Fluoroquinolones",
  SXT="Folate pathway inhibitors",COL="Polymyxins",VAN="Glycopeptides",TEC="Glycopeptides",DAP="Lipopeptides",LZD="Oxazolidinones",TET="Tetracyclines",TGC="Glycylcyclines",CLI="Lincosamides",ERY="Macrolides",
  RIF="Ansamycins",CHL="Phenicols",FUS="Fucidanes",FOS="Phosphonic acids",NIT="Nitrofurans")
# ---- 7. isolate-level flags -----------------------------------------------------------------------------
agents_wide <- ast %>% select(record_id, agent, res_final) %>% filter(!is.na(res_final))
getres <- function(rid) { s <- agents_wide[agents_wide$record_id == rid, ]; setNames(s$res_final, s$agent) }
panel_cov <- function(tab, grp_sel, thr = 0.30) {   # panel categories = categories tested in >=30% of that organism's isolates
  ids <- iso$record_id[iso$mag_table == tab]
  cats <- mag[[tab]]
  cov <- sapply(cats, function(ag) mean(sapply(ids, function(r) any(names(getres(r)) %in% ag))))
  names(cats)[cov >= thr]
}
flag_row <- function(i) {
  r <- iso[i, ]; res <- getres(r$record_id); sp <- r$species
  agents <- names(res)
  out <- list(R_any = NA, R_agents = NA_character_, n_agents = length(res))
  # --- any resistance (R only) among interpretable agents
  out$R_any <- if (length(res)) any(res == "R") else NA
  # --- AWaRe
  aw <- drug$aware[match(agents, drug$code)]
  for (k in c("Access","Watch","Reserve")) { sel <- aw == k; out[[paste0("R_", k)]] <- if (any(sel)) any(res[sel] == "R") else NA }
  # --- MDR/XDR/PDR (Magiorakos; non-susceptible = I or R)
  tab <- r$mag_table
  out$mdr <- out$xdr <- out$pdr <- NA_character_; out$n_cat_tested <- out$n_cat_ns <- NA_integer_
  {
    if (!is.na(tab)) { cats <- mag[[tab]]; ex <- mag_intr[[sp]]; if (!is.null(ex)) cats <- lapply(cats, setdiff, ex); cats <- cats[lengths(cats) > 0]
    } else { gc_ <- gen_cat[names(res)]; cats <- split(names(res), gc_) }   # organisms without a published category list: same criteria on the antimicrobial classes tested
    tested <- vapply(cats, function(ag) any(names(res) %in% ag), logical(1)); ns <- vapply(cats, function(ag) any(res[names(res) %in% ag] %in% c("I","R")), logical(1))
    n_t <- sum(tested); n_ns <- sum(ns)
    out$n_cat_tested <- n_t; out$n_cat_ns <- n_ns
    mrsa <- identical(r$mrsa, "MRSA")
    # conventional Magiorakos criteria applied to the categories actually tested for the isolate (intrinsic categories removed)
    out$mdr <- if (n_t == 0) NA_character_ else if (mrsa || n_ns >= 3) "Yes" else "No"
    out$xdr <- if (n_t == 0) NA_character_ else if (n_ns >= 3 && (n_t - n_ns) <= 2) "Yes" else "No"
    allns <- length(res) > 0 && all(res == "R")
    out$pdr <- if (n_t == 0) NA_character_ else if (n_ns >= 3 && allns) "Yes" else "No"
  }
  # --- key phenotypes
  g <- function(a) if (a %in% names(res)) res[[a]] else NA_character_
  out$ESBL_3GCR <- if (r$mag_table %in% "Enterobacterales") { v <- res[names(res) %in% c("CTX","CRO","CAZ")]; if (length(v)) any(v == "R") else NA } else NA
  out$CR <- if (r$species %in% gram_neg) { v <- res[names(res) %in% c("MEM","IPM","ETP")]; if (length(v)) any(v == "R") else NA } else NA
  out$CephR <- if (r$species %in% gram_neg) { v <- res[names(res) %in% c("CTX","CRO","CAZ","FEP","CXM")]; if (length(v)) any(v == "R") else NA } else NA
  out$MRSA <- if (r$pathogen_group == "Staphylococcus aureus") (if (r$mrsa == "unknown") NA else r$mrsa == "MRSA") else NA
  out$VR <- if (grepl("^Enterococcus", r$pathogen_group)) { v <- res[names(res) %in% c("VAN")]; if (length(v)) any(v == "R") else NA } else NA
  as_tibble(out)
}
flags <- bind_cols(iso %>% select(record_id), bind_rows(lapply(seq_len(nrow(iso)), flag_row)))
iso_f <- iso %>% left_join(flags, by = "record_id")
cat("\nR_any among isolates with >=1 interpretable result:\n"); print(table(iso_f$R_any, useNA = "ifany"))
cat("\nMDR / XDR / PDR:\n"); print(table(MDR = iso_f$mdr, useNA = "ifany")); print(table(XDR = iso_f$xdr, useNA = "ifany")); print(table(PDR = iso_f$pdr, useNA = "ifany"))
cat("\nMDR by organism table:\n"); print(table(iso_f$mag_table, iso_f$mdr, useNA = "ifany"))
saveRDS(list(iso = iso_f, ast = ast, drug = drug, tpl = tpl, dagger = dagger, class_of = class_of, cls = cls, mag = mag, gen_cat = gen_cat), file.path(OUT, "02_resistance.rds"))
