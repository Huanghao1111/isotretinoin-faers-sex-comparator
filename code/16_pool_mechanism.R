#!/usr/bin/env Rscript
# ==============================================================================
# Why does the comparator pool behave differently in the two sexes?
#
# The 12-45 female comparator pool is 180,872 reports against 92,976 male. If
# the female pool is dominated by drugs prescribed almost exclusively to women
# and those drugs carry few psychiatric reports, the female comparator rate is
# diluted, the female ROR is inflated, and a female signal appears for the index
# drug without the drug having anything to do with it. That mechanism is testable
# and it is tested here.
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
demo <- demo[primaryid %in% pop]
ps <- ps[primaryid %in% pop]
rc <- rc[primaryid %in% pop]
arm <- ps[, .(arm = if (any(drug_group == INDEX)) "index" else "comparator"),
          by = primaryid]
arm <- arm[!duplicated(primaryid)]
demo <- arm[demo, on = "primaryid", nomatch = NULL]

young <- demo[!is.na(age_yr) & age_yr >= 12 & age_yr <= 45, primaryid]
dep <- rc[pt == "Depression", unique(primaryid)]

cmp <- ps[primaryid %in% young][demo[arm == "comparator" & primaryid %in% young],
                                on = "primaryid", nomatch = NULL]
cmp <- cmp[!duplicated(primaryid)]
cmp[, depr := primaryid %in% dep]

for (sx in c("F", "M")) {
  d <- cmp[demo[cmp, on = "primaryid", sex == sx]]
  tot <- nrow(d)
  cat(sprintf("\n== 12-45 comparator pool, sex = %s : %s reports, %s with depression (%.2f%%) ==\n",
              sx, format(tot, big.mark = ","),
              format(d[, sum(depr)], big.mark = ","),
              100 * d[, mean(depr)]))
  agg <- d[, .(n = .N, depr = sum(depr)), by = drug_group]
  agg[, pct_pool := round(100 * n / tot, 1)]
  agg[, rate := round(1000 * depr / n, 1)]
  setorder(agg, -n)
  print(agg[1:15, .(drug_group, n, pct_pool, depr, rate)])
}

# how much of the female-male comparator gap do the top women-heavy drugs explain
cat("\n== share of the 12-45 comparator pool held by drugs whose use is >80% female ==\n")
sx_mix <- cmp[demo[cmp, on = "primaryid"], nomatch = NULL]
sx_mix <- merge(cmp[, .(primaryid, drug_group)],
                demo[, .(primaryid, sex)], by = "primaryid", all.x = FALSE)
sx_mix <- sx_mix[!duplicated(primaryid)]
mix <- sx_mix[, .(n = .N), by = .(drug_group, sex)]
mix <- dcast(mix, drug_group ~ sex, value.var = "n", fill = 0)
mix <- mix[F + M >= 100]
mix[, pct_female := round(100 * F / (F + M), 1)]
fheavy <- mix[pct_female >= 80 | pct_female <= 20]
print(fheavy[order(-(F + M))][1:20], nrows = 20)
cat("\nshare of female comparator reports on drugs >=80% female:",
    round(100 * sum(mix[pct_female >= 80, F]) / sum(mix$F), 1), "%\n")
cat("share of male comparator reports on drugs <=20% female:",
    round(100 * sum(mix[pct_female <= 20, M]) / sum(mix$M), 1), "%\n")

fwrite(mix, file.path(TAB, "comparator_drug_sex_mix.csv"))
cat("\n[written] comparator_drug_sex_mix.csv\n[done]\n")
