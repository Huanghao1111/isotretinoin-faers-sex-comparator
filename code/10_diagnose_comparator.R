#!/usr/bin/env Rscript
# ==============================================================================
# Diagnostic for the one threat the ROR ratio does NOT remove.
#
# Cancellation argument: because each ROR is computed within one sex, the sex
# composition of treated patients and the overall willingness to report cancel
# in the ratio. What does NOT cancel is a sex difference in the BASELINE RATE of
# a given event across the comparator pool.
#
#   ROR = (a / b) / (c / d),  c = reports carrying the event but NOT the drug
#
# If men accumulate many reports of, say, sexual dysfunction from other drugs
# (finasteride, antidepressants), c_male is inflated, ROR_male is deflated, and
# the female-to-male ROR ratio is inflated for reasons that have nothing to do
# with the index drug. This measures that term directly.
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "C:/Users/lenovo/WorkBuddy/2026-10-17-00-49-45"
PROC <- file.path(ROOT, "data", "processed")
TAB  <- file.path(ROOT, "output", "tables")
SUP  <- file.path(ROOT, "supplement")

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

EVENTS <- c("Loss of libido", "Sexual dysfunction", "Libido decreased",
            "Depression", "Anxiety", "Completed suicide", "Suicidal ideation",
            "Major depression", "Anhedonia", "Lip dry", "Dry skin",
            "Epistaxis", "Acne", "Crohn's disease",
            "Inflammatory bowel disease")

demo <- fread(file.path(PROC, "demo.csv.gz"))
drug <- fread(file.path(PROC, "drug.csv.gz"))
reac <- fread(file.path(PROC, "reac.csv.gz"))
indi <- fread(file.path(PROC, "indi.csv.gz"))

demo <- demo[sex %in% c("M", "F"), .(primaryid, sex)]
keep <- intersect(demo$primaryid,
                  unique(indi[grepl(DERM, indi_pt, ignore.case = TRUE), primaryid]))
demo <- demo[primaryid %in% keep]
iso_pid <- drug[toupper(role_cod) == "PS" & drug_group == "isotretinoin" &
                  primaryid %in% keep, unique(primaryid)]
reac <- unique(reac[primaryid %in% keep &
                      !is.na(pt) & pt != "" &
                      !grepl(NON_AE, pt, ignore.case = TRUE), .(primaryid, pt)])
cat("cohort:", nrow(demo), " | isotretinoin PS reports:", length(iso_pid), "\n")

# comparator pool: every report in the cohort NOT flagged as index drug
comp <- demo[!(primaryid %in% iso_pid)]
oniso <- demo[primaryid %in% iso_pid]

tab <- rbindlist(lapply(EVENTS, function(ev) {
  pid <- reac[pt == ev, unique(primaryid)]
  cf <- comp[sex == "F" & primaryid %in% pid, .N]
  cm <- comp[sex == "M" & primaryid %in% pid, .N]
  nf <- comp[sex == "F", .N]
  nm <- comp[sex == "M", .N]
  af <- oniso[sex == "F" & primaryid %in% pid, .N]
  am <- oniso[sex == "M" & primaryid %in% pid, .N]
  data.table(
    pt = ev,
    n_index_female = af, n_index_male = am,
    n_comparator_female = cf, n_comparator_male = cm,
    rate_index_female = round(1000 * af / oniso[sex == "F", .N], 2),
    rate_index_male   = round(1000 * am / oniso[sex == "M", .N], 2),
    rate_comparator_female = round(1000 * cf / nf, 2),
    rate_comparator_male   = round(1000 * cm / nm, 2),
    comparator_rate_ratio_M_vs_F = round((cm / nm) / (cf / nf), 2)
  )
}))

t1 <- fread(file.path(TAB, "table1_isotretinoin.csv"))
tab <- t1[, .(pt = Event, ror_ratio = `ROR ratio (95% CI)`)][tab, on = "pt"]
tab[, srt := -abs(log(comparator_rate_ratio_M_vs_F))]
setorder(tab, srt)
tab[, srt := NULL]
fwrite(tab, file.path(SUP, "TableS16_comparator_baseline_by_sex.csv"))

cat("\n== comparator-pool baseline rate, men versus women ==\n")
cat("ratio > 1 means the event is commoner in the male comparator pool, which\n")
cat("deflates ROR_male and inflates the female-to-male ROR ratio\n\n")
print(tab[, .(pt, ror_ratio, comparator_rate_ratio_M_vs_F,
              rate_comparator_female, rate_comparator_male,
              n_index_female, n_index_male)], nrows = 20)

# ------------------------------------------- same diagnostic at prescribing ages
# The one positive finding in the manuscript is that depression shifts to female
# predominance at ages 12-45. That survives only if the comparator baseline for
# depression is balanced within that age band, so it is checked here rather than
# assumed.
cat("\n== restricted to ages 12-45 years ==\n")
demo_a <- fread(file.path(PROC, "demo.csv.gz"))[sex %in% c("M", "F"),
                .(primaryid, sex, age_yr)]
demo_a <- demo_a[primaryid %in% keep & age_yr >= 12 & age_yr <= 45]
iso_a <- intersect(iso_pid, demo_a$primaryid)
comp_a <- demo_a[!(primaryid %in% iso_a)]
oniso_a <- demo_a[primaryid %in% iso_a]
cat("cohort 12-45:", nrow(demo_a), "| isotretinoin:", length(iso_a), "\n")

tab_a <- rbindlist(lapply(c("Depression", "Anxiety", "Suicidal ideation",
                            "Completed suicide", "Loss of libido", "Lip dry",
                            "Dry skin"), function(ev) {
  pid <- reac[pt == ev, unique(primaryid)]
  cf <- comp_a[sex == "F" & primaryid %in% pid, .N]
  cm <- comp_a[sex == "M" & primaryid %in% pid, .N]
  nf <- comp_a[sex == "F", .N]; nm <- comp_a[sex == "M", .N]
  af <- oniso_a[sex == "F" & primaryid %in% pid, .N]
  am <- oniso_a[sex == "M" & primaryid %in% pid, .N]
  data.table(
    pt = ev,
    n_index_female = af, n_index_male = am,
    n_comparator_female = cf, n_comparator_male = cm,
    rate_index_female = round(1000 * af / oniso_a[sex == "F", .N], 2),
    rate_index_male   = round(1000 * am / oniso_a[sex == "M", .N], 2),
    rate_comparator_female = round(1000 * cf / nf, 2),
    rate_comparator_male   = round(1000 * cm / nm, 2),
    comparator_rate_ratio_M_vs_F = round((cm / nm) / (cf / nf), 2),
    index_rate_ratio_F_vs_M = round((af / oniso_a[sex == "F", .N]) /
                                      (am / oniso_a[sex == "M", .N]), 2)
  )
}))
fwrite(tab_a, file.path(SUP, "TableS18_comparator_baseline_age12_45.csv"))
print(tab_a[, .(pt, n_index_female, n_index_male,
                 rate_index_female, rate_index_male,
                 rate_comparator_female, rate_comparator_male,
                 comparator_rate_ratio_M_vs_F, index_rate_ratio_F_vs_M)],
      nrows = 10)

cat("\n[done]\n")
