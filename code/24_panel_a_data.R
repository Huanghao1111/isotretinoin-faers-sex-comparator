#!/usr/bin/env Rscript
# ==============================================================================
# Panel A data: depression reports per 1,000 male reports, by comparator drug,
# ages 12-45, in the untrimmed pool.
#
# The manuscript states the contamination as a pair of counts (2,572 of 92,976
# reports, 779 of 1,385 depression reports). Those numbers are true and hard to
# feel. Put every comparator drug on the same axis and the male baseline is
# visibly one drug.
#
# Rates need a denominator, so a drug enters the panel only if it carries at
# least 200 male reports at these ages; below that a rate is two or three events
# and noise. The threshold is stated in the legend.
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "C:/Users/lenovo/WorkBuddy/2026-10-17-00-49-45"
PROC <- file.path(ROOT, "data", "processed")
TAB  <- file.path(ROOT, "output", "tables")

DERM <- paste0(
  "psoriasis|psoriatic|dermatitis|eczema|acne|urticaria|angioedema|",
  "hidradenitis|alopecia|hair loss|vitiligo|melanoma|basal cell carcinoma|",
  "squamous cell carcinoma|actinic keratosis|skin cancer|skin neoplasm|",
  "rosacea|pemphigus|pemphigoid|lupus|pruritus|itch|onychomycosis|tinea|",
  "lichen planus|morphea|scleroderma|pyoderma|photo|seborrhoe|",
  "hyperhidrosis|cutane|skin|wart|cellulitis|impetigo|scabies|herpes|",
  "erythema nodosum|granuloma annulare|sarcoid|keratosis")
NON_AE <- paste0(
  "^(product use issue|off label use|unapproved drug use|",
  "product used? for unapproved indication|wrong (patient|drug|dose)|",
  "drug (dispensing|administration) error|medication error|",
  "incorrect dose administered|intercepted (drug|medication)|",
  "intentional product (misuse|use issue)|product (misuse|quality issue)|",
  "no adverse (event|drug effect)|treatment noncompliance|",
  "therapy (interrupted|cessation)|drug (withdrawn|dose omission)|",
  "circumstance or information capable of leading to|",
  "product (substitution|label|packaging) issue|expired product)")

INDEX <- "isotretinoin"
MIN_MALE_N <- 200L

demo <- fread(file.path(PROC, "demo.csv.gz"))[sex %in% c("M", "F"),
                                              .(primaryid, sex, age_yr)]
ps   <- unique(fread(file.path(PROC, "drug.csv.gz"))[, .(primaryid, role_cod, drug_group)])
rc   <- unique(fread(file.path(PROC, "reac.csv.gz"))[, .(primaryid, pt)])
indi <- unique(fread(file.path(PROC, "indi.csv.gz"))[, .(primaryid, indi_pt)])

demo <- demo[primaryid %in% unique(indi[grepl(DERM, indi_pt, ignore.case = TRUE),
                                        primaryid])]
ps <- ps[primaryid %in% demo$primaryid & toupper(role_cod) == "PS"]
rc <- rc[primaryid %in% demo$primaryid & !is.na(pt) & pt != "" &
           !grepl(NON_AE, pt, ignore.case = TRUE)]
pop <- intersect(ps$primaryid, rc$primaryid)
demo <- demo[primaryid %in% pop]; ps <- ps[primaryid %in% pop]
rc <- rc[primaryid %in% pop]

arm <- ps[, .(arm = if (any(drug_group == INDEX)) "index" else "comparator"),
          by = primaryid]
arm <- arm[!duplicated(primaryid)]
demo <- arm[demo, on = "primaryid", nomatch = NULL]
demo <- demo[arm == "comparator" & !is.na(age_yr) & age_yr >= 12 & age_yr <= 45]

dep <- rc[pt == "Depression", unique(primaryid)]
# a comparator report can carry several drugs; count it once per drug
cmp <- unique(ps[demo, on = "primaryid", nomatch = NULL,
                 .(primaryid, drug_group, sex)])
cmp[, dep := primaryid %in% dep]

agg <- cmp[, .(n_reports = .N, n_depression = sum(dep)), by = .(drug_group, sex)]
dcast_wide <- dcast(agg, drug_group ~ sex, value.var = c("n_reports", "n_depression"),
                    fill = 0)
setnames(dcast_wide,
         c("n_reports_F", "n_reports_M", "n_depression_F", "n_depression_M"),
         c("n_F", "n_M", "dep_F", "dep_M"))
dcast_wide[, `:=`(rate_F = 1000 * dep_F / n_F, rate_M = 1000 * dep_M / n_M,
                  pct_female = 100 * n_F / (n_F + n_M))]
dcast_wide <- dcast_wide[n_M >= MIN_MALE_N]

cat <- fread(file.path(ROOT, "code", "drug_catalog.csv"),
             encoding = "UTF-8")[, .(drug_group, drug_label, category)]
dcast_wide <- merge(dcast_wide, cat, by = "drug_group", all.x = TRUE)
dcast_wide[is.na(drug_label), drug_label := drug_group]
dcast_wide <- dcast_wide[order(-rate_M)]

fwrite(dcast_wide, file.path(TAB, "panel_a_drug_rates.csv"))
cat("drugs with >=", MIN_MALE_N, "male reports at ages 12-45:", nrow(dcast_wide), "\n")
print(dcast_wide[1:12, .(drug_label, n_M, dep_M, rate_M, rate_F, pct_female)])

# Pooled rate across the whole (untrimmed) comparator at these ages. This is a
# report-level quantity: counting `cmp` would count a report once per comparator
# drug it names, which inflates both the denominator and the numerator and gave
# 15.1 instead of 14.9 per 1,000.
pool <- demo[, .(n = .N, d = sum(primaryid %in% dep)), by = sex]
pool[, rate_per_1000 := round(1000 * d / n, 2)]
fwrite(pool, file.path(TAB, "panel_a_pooled.csv"))
print(pool)
cat("\npooled male comparator rate per 1,000:",
    pool[sex == "M", rate_per_1000], "\n")
cat("pooled female comparator rate per 1,000:",
    pool[sex == "F", rate_per_1000], "\n")
cat("\n[written] panel_a_drug_rates.csv\n[done]\n")
