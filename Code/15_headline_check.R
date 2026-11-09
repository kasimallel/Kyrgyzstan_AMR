# 15_headline_check.R -- self-check: the headline numbers of the manuscript must be reproduced exactly.
source(file.path(Sys.getenv("KGZ_PKG"), "R", "config.R"))
x <- file.path(OUT, "04_prevalence_incidence_tables.xlsx")
L1 <- readxl::read_excel(x, "Level1_agnostic"); M <- readxl::read_excel(x, "Prevalence_incidence_merged"); dist <- readxl::read_excel(x, "Pathogen_distribution")
g <- function(pat, col) L1[[col]][grepl(pat, L1$indicator)][1]
checks <- list(
  `culture-confirmed bacterial BSI isolates` = c(sum(dist$n), 294),
  `isolates with an interpretable susceptibility result` = c(L1$n_denominator[grepl("^Any resistance", L1$indicator)], 252),
  `any resistance (n)` = c(g("^Any resistance", "n_resistant"), 229),
  `MDR (n)` = c(g("^MDR", "n_resistant"), 138),
  `XDR (n)` = c(g("^XDR", "n_resistant"), 46),
  `3GC-resistant Enterobacterales (n resistant)` = c(g("^3rd-gen", "n_resistant"), 88),
  `carbapenem-resistant Gram-negatives (n resistant)` = c(g("^Carbapenem", "n_resistant"), 27),
  `carbapenem-resistant Gram-negatives (n evaluable)` = c(g("^Carbapenem", "n_denominator"), 119),
  `MRSA (n resistant)` = c(g("^MRSA", "n_resistant"), 22),
  `incidence, 1 in x admissions (all BSI)` = c(M$one_in[M$indicator == "All culture-confirmed bacterial BSI"], 2629))
ok <- vapply(checks, function(v) isTRUE(all.equal(as.numeric(v[1]), as.numeric(v[2]))), logical(1))
print(data.frame(check = names(checks), reproduced = vapply(checks, `[`, 1, 1), expected = vapply(checks, `[`, 1, 2), ok = ok), row.names = FALSE)
if (!all(ok)) stop("Headline check FAILED: results differ from the published numbers.")
cat("Headline check passed: all", length(ok), "numbers reproduced.\n")
