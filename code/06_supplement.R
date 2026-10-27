#!/usr/bin/env Rscript
# ==============================================================================
# Builds the supplemental dataset for the Mendeley repository.
#
# Included:
#   Table S1  all significant sex x drug x event interactions (full statistics)
#   Table S2  sex composition of reports by dermatologic drug category
#   Table S3  serious outcomes, sex contrast stratified by drug
#   Table S4  sensitivity scenarios for the isotretinoin signals
#   Table S5  every drug-event pair evaluated (full exploratory scan)
#   README    data dictionary for the above
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "C:/Users/lenovo/WorkBuddy/2026-10-17-00-49-45"
PROC <- file.path(ROOT, "data", "processed")
TAB  <- file.path(ROOT, "output", "tables")
SUP  <- file.path(ROOT, "supplement")
dir.create(SUP, recursive = TRUE, showWarnings = FALSE)

CATALOG0 <- fread(file.path(ROOT, "code", "drug_catalog.csv"))
CATLBL <- CATALOG0[, .(drug_group, drug_label, category)][!duplicated(drug_group)]

# cohorts vocabularies below are read back from the analysis script so that the
# supplement can never drift away from what was actually executed
rx_block <- function(txt, name) {
  i <- grep(paste0("^", name, " <- paste0\\("), txt)
  if (!length(i)) return(NA_character_)
  j <- i[1]
  # the block closes on the first later line that ends in a closing paren
  later <- (j + 1):length(txt)
  end <- later[grep("\\)\\s*$", trimws(txt[later]))[1]]
  body <- txt[(j + 1):end]
  # each line is a quoted R string fragment: drop the quote/comma/paren wrapper
  clean <- trimws(body)
  clean <- sub('^"', "", clean)
  clean <- sub('[",)]+$', "", clean)
  paste(clean, collapse = "")
}
ASCRIPT <- readLines(file.path(ROOT, "code", "03_analysis.R"))
DERM_TERMS   <- rx_block(ASCRIPT, "DERM")
NON_AE_TERMS <- rx_block(ASCRIPT, "NON_AE")
stopifnot(!is.na(DERM_TERMS), !is.na(NON_AE_TERMS))

A <- readRDS(file.path(PROC, "analysis_objects.rds"))
demo <- as.data.table(A$demo)
ser  <- as.data.table(A$ser)
t1b  <- as.data.table(A$t1b)

r2 <- function(x) round(as.numeric(x), 2)
# p values below ~1e-308 underflow to 0 in double precision, so the -log10
# transform is computed from the z statistic in log space, where it stays finite
neglog10p <- function(z) round(-(log(2) + pnorm(-abs(z), log.p = TRUE)) / log(10), 1)

# ------------------------------------------------------------------- Table S1
sig <- as.data.table(A$sig)
# keep all statistics on the table BEFORE attaching readable labels:
# a join can reorder rows, and referring to an external sig$column afterwards
# would silently pair the wrong numbers with the wrong event
s1 <- copy(sig)
s1 <- CATLBL[s1, on = "drug_group", nomatch = NULL]
s1 <- s1[, .(category, drug_group, drug_label, pt,
             n_female = a_F, n_male = a_M,
             ROR_female = r2(ror_F),
             ROR_female_lci = r2(lci_F),
             ROR_female_uci = r2(uci_F),
             ROR_male = r2(ror_M),
             ROR_male_lci = r2(lci_M),
             ROR_male_uci = r2(uci_M),
             ROR_ratio = r2(rorr),
             ROR_ratio_lci = r2(rorr_lci),
             ROR_ratio_uci = r2(rorr_uci),
             neglog10_p_interaction = neglog10p(z),
             p_interaction = signif(p_int, 4),
             q_interaction_BH = signif(q_int, 4))]
setorder(s1, -neglog10_p_interaction)
fwrite(s1, file.path(SUP, "TableS1_significant_sex_interactions.csv"))
cat("[S1] significant interactions:", nrow(s1), "\n")

# ------------------------------------------------------------------- Table S2
s2 <- t1b[, .(category,
              reports_female = n_F,
              reports_male   = n_M,
              total_reports  = tot,
              pct_female     = pct_F)]
setorder(s2, -total_reports)
fwrite(s2, file.path(SUP, "TableS2_category_sex_composition.csv"))
cat("[S2] categories:", nrow(s2), "\n")

# ------------------------------------------------------------------- Table S3
evs <- c("DE", "HO", "LT", "DS")
s3 <- ser[, .(n_total = .N,
              n_female = sum(sex == "F"),
              n_male   = sum(sex == "M"),
              deaths_female = sum(DE > 0 & sex == "F"),
              deaths_male   = sum(DE > 0 & sex == "M"),
              hospitalised_female = sum(HO > 0 & sex == "F"),
              hospitalised_male   = sum(HO > 0 & sex == "M"),
              lifethreatening_female = sum(LT > 0 & sex == "F"),
              lifethreatening_male   = sum(LT > 0 & sex == "M"),
              disability_female = sum(DS > 0 & sex == "F"),
              disability_male   = sum(DS > 0 & sex == "M")), by = drug_group]
s3 <- CATLBL[s3, on = "drug_group", nomatch = NULL]
s3 <- s3[n_total >= 200]
setorder(s3, -n_total)
fwrite(s3, file.path(SUP, "TableS3_serious_outcomes_by_drug.csv"))
cat("[S3] drugs with >=200 reports:", nrow(s3), "\n")

# ------------------------------------------------------------------- Table S4
sf <- file.path(TAB, "sensitivity_isotretinoin.csv")
if (file.exists(sf)) {
  s4 <- fread(sf)
  setorder(s4, pt, scenario)
  fwrite(s4, file.path(SUP, "TableS4_sensitivity_isotretinoin.csv"))
  cat("[S4] sensitivity rows:", nrow(s4), "\n")
}

# ------------------------------------------------------------------- Table S5
R <- as.data.table(A$R)
s5 <- R[, .(drug_group, pt, a_F, a_M,
            ROR_female = r2(ror_F), ROR_male = r2(ror_M),
            ROR_ratio = r2(rorr),
            ROR_ratio_lci = r2(rorr_lci), ROR_ratio_uci = r2(rorr_uci),
            neglog10_p_interaction = neglog10p(z),
            p_interaction = signif(p_int, 4),
            q_interaction_BH = signif(q_int, 4))]
s5 <- CATLBL[s5, on = "drug_group", nomatch = NULL]
setorder(s5, -neglog10_p_interaction)
fwrite(s5, file.path(SUP, "TableS5_all_drug_event_pairs.csv.gz"))
cat("[S5] all evaluated pairs:", nrow(s5), "\n")

# ------------------------------------------------------------------- Table S6
# RECORD-PE item 7.1: complete list of codes and algorithms used to define
# exposure. One row per source string per drug, so a reader can reproduce the
# cohort without reading the code.
s6 <- CATALOG0[, {
  # NB: regex split - fixed = TRUE would treat "\\|" as two literal characters
  v <- trimws(unlist(strsplit(as.character(synonyms), "\\|")))
  v <- v[nzchar(v) & !is.na(v)]
  .(category = category[1], drug_label = drug_label[1], search_string = v)
}, by = drug_group]
setorder(s6, category, drug_group, search_string)
fwrite(s6, file.path(SUP, "TableS6_drug_search_strings.csv"))
cat("[S6] drug search strings:", nrow(s6), "\n")

# ------------------------------------------------------------------- Table S7
# RECORD-PE item 7.1: cohort-defining vocabularies
#  - the indication strings that qualified a report as dermatologic
#  - the preferred terms excluded as non-adverse events
ind_list <- trimws(strsplit(DERM_TERMS, "\\|")[[1]])
ind_list <- ind_list[nzchar(ind_list)]
# NOTE: the non-AE list is a grouped regular expression, not a flat term list,
# so it is reproduced verbatim rather than split on the alternation operator
s7 <- rbindlist(list(
  data.table(vocabulary = "dermatologic_indication",
             match_rule = "case-insensitive substring match on INDI.indi_pt",
             term = ind_list),
  data.table(vocabulary = "excluded_non_adverse_event_pt",
             match_rule = "anchored case-insensitive regular expression on REAC.pt",
             term = NON_AE_TERMS)
))
s7 <- s7[term != ""]
fwrite(s7, file.path(SUP, "TableS7_cohort_definition_terms.csv"))
cat("[S7] cohort-definition terms:", nrow(s7), "\n")

# --------------------------------------------------------------------- README
# The README is generated here, in full, and nowhere else. Earlier versions were
# assembled by having 06 write a skeleton and 07 and 09 append to it; the result
# was out-of-order table entries and, on one Windows configuration, a truncated
# file. One writer, one order.
readme <- c(
  "# Supplemental data - Male predominance in depression reporting with isotretinoin after rebuilding the comparator pool",
  "",
  "Companion dataset for the JAAD Brief Report *Male predominance in",
  "depression reporting with isotretinoin after rebuilding the comparator pool*.",
  "",
  "## Source",
  "FDA Adverse Event Reporting System (FAERS) quarterly ASCII extracts,",
  "2012 Q4 to 2026 Q2 (55 quarters). FDA retains no earlier quarters at the",
  "public download endpoint. Reports were retained if they carried at least one",
  "dermatologic indication and listed a dermatologic drug as primary suspect.",
  "",
  "## Cohort",
  "3,464,909 dermatologic drug reports; 1,148,083 after restricting to",
  "dermatologic indications. Female 709,690 (61.8%); male 438,393 (38.2%).",
  "Deduplication: one record per caseid, retaining the highest caseversion.",
  "1,091,621 of these carry at least one primary-suspect dermatologic drug and",
  "one retained adverse-event term and form the analysis population; 13,829 are",
  "isotretinoin reports (7,496 female, 6,333 male).",
  "",
  "## Two comparators",
  "Every ratio of ratios is reported twice, against two reference pools built from",
  "the same analysis population:",
  "",
  "- *all-drug* - all other dermatologic drug reports, the conventional choice,",
  "  retained in the supplement because it shows the size of the artefact;",
  "- *sex-neutral* - the primary analysis, in which comparator drugs whose reports",
  "  come <20% or >80% from one sex are removed and a comparator report is kept",
  "  only if every primary-suspect drug on it is kept.",
  "",
  "The restriction removes 30,311 of 1,077,792 comparator reports - 2.8% - and",
  "keeps 1,047,481, or 97.2% of the pool (TableS34).",
  "",
  "`NOTE_comparator_design.md` documents the two faults that motivated the second",
  "pool: an age-structure mismatch producing Simpson's paradox, and sex-specific",
  "comparator drugs (chiefly finasteride) setting the male baseline.",
  "",
  "## Files",
  "",
  "**TableS1_significant_sex_interactions.csv** - drug-event pairs with a",
  "statistically significant sex difference (p_interaction < 0.05, ratio",
  "confidence interval excluding 1, at least 3 reports in each sex and 20 in",
  "total).",
  "",
  "**TableS2_category_sex_composition.csv** - number and percentage of reports",
  "from female patients, by dermatologic drug category. Shows how far the sex",
  "composition of reporting is driven by who receives each drug class.",
  "",
  "**TableS3_serious_outcomes_by_drug.csv** - counts of death, hospitalisation,",
  "life-threatening and disability outcomes by sex, for each drug with at least",
  "200 reports. Underlies the drug-stratified pooled estimates. Caveat: the",
  "death contrast is within drug and therefore free of drug mixing, but it is a",
  "crude sex contrast - men die more than women at every age in the general",
  "population and FAERS holds no denominator, so nothing here identifies how much",
  "of the difference is the drug. That is why it is not in the Brief Report body.",
  "",
  "**TableS4_sensitivity_isotretinoin.csv** - the isotretinoin signals",
  "re-estimated under three restrictions (single suspect drug, health-",
  "professional reporter, ages 12-45 years) and their combination. Built on the",
  "superseded design; see TableS30 for the same restrictions on the corrected",
  "design.",
  "",
  "**TableS5_all_drug_event_pairs.csv.gz** - every drug-event pair evaluated,",
  "regardless of statistical significance. Provided so that the full",
  "exploratory scan can be inspected and re-analysed.",
  "",
  "**TableS6_drug_search_strings.csv** - every ingredient name, trade name and",
  "synonym used to identify exposure, one row per search string. Matching was",
  "case-insensitive with word boundaries against the concatenation of DRUG.",
  "drugname and DRUG.prod_ai. This is the exposure code list.",
  "",
  "**TableS7_cohort_definition_terms.csv** - the indication terms that qualified",
  "a report as dermatologic, plus the regular expression defining preferred",
  "terms excluded as non-adverse events (medication errors, product-use issues,",
  "lack of efficacy).",
  "",
  "**TableS8_cohort_selection_flow.csv** - stepwise selection of the analysis",
  "cohort, with the count remaining and the count dropped at each step and the",
  "reason for the drop.",
  "",
  "**TableS9_reporter_by_sex_isotretinoin.csv** - who submitted each",
  "isotretinoin report, by patient sex (FAERS occupation code: CN consumer,",
  "MD physician, PH pharmacist, HP other health professional, RN nurse,",
  "LW lawyer, OT other).",
  "",
  "**TableS10_indication_by_sex_isotretinoin.csv** - what the drug was taken",
  "for, by patient sex. Both indication rows and distinct reports are given,",
  "since a report may carry several indication terms.",
  "",
  "**TableS11_top_indication_terms_isotretinoin.csv** - the unchanged literal",
  "indication strings behind TableS10, so the grouping can be inspected.",
  "",
  "**TableS12_reports_by_year.csv** - reports and sex composition by year of",
  "receipt. Four rows carry impossible years (321, 1012, 2029, 5015) that come",
  "from the source extracts and are reproduced rather than cleaned, because",
  "editing a date would change a count. They total five reports.",
  "event, which shows how drug availability changed over the 2012-2026 window.",
  "",
  "**TableS13_lawyer_channel_by_sex.csv** - reports reaching FDA through a",
  "lawyer, by patient sex, across the whole dermatologic cohort.",
  "",
  "**TableS14_lawyer_channel_by_category_sex.csv** - the same, split by",
  "therapeutic class. Shows the litigation channel is concentrated almost",
  "entirely in systemic retinoids rather than spread across the specialty.",
  "",
  "**TableS15_significant_pairs_with_without_lawyers.csv** - every",
  "significant sex interaction re-estimated with lawyer-filed reports removed,",
  "so a reader can see which signals depend on the reporting channel.",
  "",
  "**Column convention.** In Tables S19-S23 and S30-S40, `a_F` and `a_M` are the",
  "index arm's event counts, `n_iF` and `n_iM` its denominators, `c_F` and `c_M`",
  "the comparator's event counts, and `n_cF` and `n_cM` the comparator's",
  "denominators. Every reporting odds ratio in those tables can be recomputed",
  "from these six cells; nothing is left implicit.",
  "",
  "**TableS16_comparator_baseline_by_sex.csv** - an early diagnostic, computed on",
  "the wider pre-restriction cohort (index arm 7,684 women / 6,398 men, TableS9)",
  "and therefore NOT the denominator of the primary panels. Kept for the record;",
  "the denominators actually used are in TableS19 and TableS20.",
  "",
  "**TableS17_active_comparator_tetracycline.csv** - an early version of the",
  "active comparator, computed on a slightly wider cohort and without the acne",
  "indication restriction. Superseded by TableS42, which is the version quoted.",
  "",
  "**TableS18_comparator_baseline_age12_45.csv** - the same early diagnostic",
  "restricted to ages 12-45; it inherits the wider cohort of TableS16.",
  "",
  "## Tables S19-S29: corrected comparator design",
  "",
  "- TableS19 primary ratios of ratios, all ages (sex-neutral comparator)",
  "- TableS20 primary ratios of ratios, ages 12-45",
  "- TableS21 uncorrected ratios of ratios, all ages",
  "- TableS22 uncorrected ratios of ratios, ages 12-45",
  "- TableS23 the same estimates under five comparator trimming rules",
  "- TableS24 sex composition of every comparator drug, with retention flag",
  "- TableS25 share of male comparator depression reports supplied by finasteride",
  "- TableS26 age-band composition of both arms after trimming",
  "- TableS27 formal age-band interaction (12-45 versus 46+), uncorrected pool",
  "- TableS28 pooled estimates restricted to reports with a recorded age",
  "- TableS29 age-band composition of both arms before trimming",
  "",
  "See `NOTE_comparator_design.md` for why the comparator was rebuilt.",
  "",
  "## Tables S30-S34: the corrected design under restriction",
  "",
  "- TableS30 the seven restriction scenarios re-run on the corrected design",
  "- TableS31 age-band interaction (12-45 versus 46+) on the sex-neutral comparator",
  "- TableS32 ratios of ratios by age band on the sex-neutral comparator",
  "- TableS33 how far each term moved when the comparator was rebuilt, by panel",
  "- TableS34 comparator reports before and after the restriction, and the fraction",
  "  retained",
  "",
  "TableS31 supersedes TableS27, which was computed on the uncorrected",
  "comparator; both are reported because the interaction is null in one and",
  "significant in the other, and the difference is the point. TableS33 is the",
  "internal control: the correction moves the psychiatric and sexual terms",
  "(median |log shift| 0.67 and 2.26) and leaves the mucocutaneous and",
  "gastrointestinal terms where they were (0.02).",
  "",
  "## Tables S35-S40: stability over time, the male baseline, monitoring",
  "",
  "- TableS35 ratios of ratios by era, sex-neutral comparator",
  "- TableS36 heterogeneity of the ratio of ratios across eras",
  "- TableS37 finasteride's share of male comparator depression reports, by era",
  "- TableS38 depression reports per 1,000 male reports for the 35 comparator drugs",
  "- TableS39 reports carrying only monitoring or pregnancy terms, by arm and sex",
  "- TableS40 ratios of ratios with and without those reports",
  "- TableS41 the 15-term table formerly called Table I",
  "- TableS42 isotretinoin versus tetracyclines, both arms restricted to acne",
  "- TableS43 leave-one-drug-out refits across the 64 retained comparator drugs",
  "",
  "**TableS41_primary_estimates_all_terms.csv** - the manuscript's Table I, moved",
  "here because a two-panel figure already uses both display slots. Sex-stratified",
  "disproportionality for 15 adverse events with isotretinoin as primary suspect",
  "among 13,829 such reports in 1,091,621 dermatologic-indication reports (7,496",
  "from women, 6,333 from men), FAERS 2012 Q4 to 2026 Q2. \"All-drug\" uses every",
  "other dermatologic drug as comparator; \"sex-neutral\" is the primary analysis.",
  "The last two columns decompose the ratio of ratios: \"comparator OR, men versus",
  "women\" multiplied by \"isotretinoin OR, women versus men\" equals the sex-neutral",
  "ratio of ratios up to rounding of the printed factors, so a comparator OR far",
  "is being set by the pool rather than by the drug. P values are unadjusted;",
  "Benjamini-Hochberg q values are in TableS19.",
  "",
  "Depression is below 1 in all three eras - 0.77 (0.67-0.88), 0.65 (0.54-0.79)",
  "and 0.48 (0.35-0.66), heterogeneity p = 0.018 - but at ages 12-45 only the last",
  "era excludes 1 (0.87, 0.90, 0.37; p = 0.002). The direction is consistent and",
  "the magnitude is not, and terms with no psychiatric mechanism drift as much or",
  "more, so the drift is a property of the comparator rather than of the drug.",
  "See `NOTE_comparator_design.md`.",
  "",
  "## Column definitions",
  "",
  "- `drug_group`: internal identifier for the active ingredient or product.",
  "- `drug_label`: readable name of that ingredient or product.",
  "- `category`: therapeutic class.",
  "- `pt`: MedDRA preferred term as reported.",
  "- `n_female` / `n_male`: number of reports from each sex carrying that",
  "  drug-event combination.",
  "- `ROR_female` / `ROR_male`: reporting odds ratio computed within that sex,",
  "  against all other dermatologic drug reports of the same sex.",
  "- `ROR_ratio`: ROR_female divided by ROR_male. Values above 1 denote stronger",
  "  disproportionality in women, values below 1 stronger in men.",
  "- `p_interaction`: two-sided z test on the difference of log RORs. Values",
  "  below about 1e-308 underflow to 0 in double precision; use the companion",
  "  `neglog10_p_interaction` column for those rows, which is derived from the",
  "  z statistic in log space and remains finite at any magnitude.",
  "- `q_interaction_BH`: Benjamini-Hochberg adjusted p value across all pairs.",
  "- `a_F` / `a_M`: reports with the event, in women and in men.",
  "- `n_iF` / `n_iM`: isotretinoin reports in women and in men.",
  "- `RORR` with `lci` / `uci`: ratio of ratios and its 95% confidence interval,",
  "  interval from the standard error of the sum of the two log odds ratios.",
  "- `index_OR_F_vs_M`: the drug's own odds ratio, women versus men, within the",
  "  index arm. Free of the comparator pool and therefore the quantity the",
  "  correction is meant to leave alone.",
  "- `comp_OR_M_vs_F`: the comparator factor, taken in the opposite direction.",
  "  `RORR = index_OR_F_vs_M x comp_OR_M_vs_F` before rounding; a comparator factor far",
  "  from 1 means the ratio of ratios is mostly an artefact of the reference pool.",
  "- `scenario`: which restriction produced the row (main, single suspect,",
  "  health-professional reporter, ages 12-45, their combination, lawyers",
  "  removed).",
  "- `abs_log_shift`: absolute change in log RORR between the two comparators,",
  "  used as the internal control in TableS33.",
  "",
  "## Reproducibility",
  "",
  "Analysis code is in the code directory of the associated project, run in this",
  "order:",
  "",
  "1. `01b_download_parallel.py` - retrieve the 55 quarterly extracts",
  "2. `02_build_cohort.py` - parse, deduplicate, restrict to dermatologic",
  "   indications, build the report-level tables",
  "3. `03_analysis.R` - exploratory disproportionality scan, sex interaction",
  "4. `03b_serious_outcomes.R` - drug-stratified serious outcomes (TableS3)",
  "5. `12_correct_2x2.R`, `13_primary_numbers.R`, `14_simpson_check.R`,",
  "   `16_pool_mechanism.R`, `17_final_analysis.R` - the corrected comparator",
  "   design and the two faults it addresses",
  "6. `18_robustness.R` - comparator trimming rules (TableS23)",
  "7. `20_primary_sensitivity.R` - restrictions and the age interaction on the",
  "   corrected design (TablesS30-S33)",
  "8. `08_verify_numbers.R` - re-derives every number quoted in the manuscript",
  "   from the output tables and exits non-zero on any mismatch",
  "9. `04_figures.R`, `99_check_pdf.R` - figure and vector-text verification",
  "",
  "This README is written in full by `06_supplement.R`, which is the only script",
  "that writes it; `19_supplement_update.R` copies Tables S19-S34 into this",
  "directory and writes `NOTE_comparator_design.md`.",
  "",
  "`10_diagnose_comparator.R`, `11_active_comparator.R` and",
  "`15_acne_comparator.R` are retained but describe rejected approaches; their",
  "outputs are TableS16-S18 and the note that no active comparator is claimed.",
  "",
  "## Ethics",
  "De-identified public data; no institutional review board approval required."
)
writeLines(readme, file.path(SUP, "README.md"))
cat("[README] written:", length(readme), "lines\n")

cat("\n[done] supplemental dataset in", SUP, "\n")
for (f in list.files(SUP)) {
  cat("  ", f, sprintf("%.1f KB", file.size(file.path(SUP, f)) / 1024), "\n")
}
