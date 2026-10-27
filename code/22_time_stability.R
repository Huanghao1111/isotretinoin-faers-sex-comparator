#!/usr/bin/env Rscript
# ==============================================================================
# Calendar-time stability of the corrected ratio of ratios.
#
# The correction rests on finasteride sitting in the male comparator pool. Two
# things could make that an artefact of when the reports arrived: finasteride's
# own psychiatric literature appeared during the window (Welk 2017, Nguyen 2021),
# so reporting of it may be concentrated late, and isotretinoin's reporting
# volume moved as iPLEDGE and the acne market changed. If the male predominance
# is present in every era, the finding is not a period effect. If it is present
# only where finasteride reports are dense, the paper has to say so.
#
# The comparator definition is fixed once, on the full period, and then applied
# inside each era: letting the drug mix be re-estimated per era would change the
# definition of the comparator between columns.
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
MIN_DRUG_N <- 200L; LO_SX <- 20; HI_SX <- 80

TERMS <- c("Depression", "Anxiety", "Suicidal ideation", "Completed suicide",
           "Major depression", "Anhedonia",
           # non-neuropsychiatric controls: no sex difference is expected, so a
           # shift here would mean the correction is doing something general
           "Lip dry", "Dry skin", "Acne")
ERAS <- c("2012-2016", "2017-2021", "2022-2026")

demo <- fread(file.path(PROC, "demo.csv.gz"))[sex %in% c("M", "F"),
                                              .(primaryid, sex, age_yr, event_year)]
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
cat("analysis population:", nrow(demo), "\n")

demo[, era := fifelse(event_year <= 2016, ERAS[1],
               fifelse(event_year <= 2021, ERAS[2], ERAS[3]))]
demo[, in_band := !is.na(age_yr) & age_yr >= 12 & age_yr <= 45]
print(demo[, .N, by = era][order(era)])

arm <- ps[, .(arm = if (any(drug_group == INDEX)) "index" else "comparator"),
          by = primaryid]
arm <- arm[!duplicated(primaryid)]
demo <- arm[demo, on = "primaryid", nomatch = NULL]

# ------------------------------------------------- comparator, fixed globally
cmp <- ps[demo[arm == "comparator"], on = "primaryid", nomatch = NULL]
cmp <- cmp[!duplicated(primaryid)][!is.na(drug_group)]
cmp <- merge(cmp[, .(primaryid, drug_group)], demo[, .(primaryid, sex)],
             by = "primaryid", all.x = FALSE)
mix <- dcast(cmp[, .(n = .N), by = .(drug_group, sex)], drug_group ~ sex,
             value.var = "n", fill = 0)
mix[, tot := F + M]
mix[, pct_female := round(100 * F / tot, 1)]
keep_drug <- mix[tot >= MIN_DRUG_N & pct_female >= LO_SX & pct_female <= HI_SX,
                 drug_group]
cat("comparator drugs retained:", length(keep_drug), "\n")

allps <- ps[, .(bad = sum(!(drug_group %in% keep_drug))), by = primaryid]
demo <- merge(demo, allps, by = "primaryid", all.x = TRUE)
demo[, usable := arm == "index" | (is.na(bad) | bad == 0)]

# One pass per band. Subsetting and recomputing both the denominators and the
# event counts is the only way to keep "all ages" meaning all ages: labelling
# every row with a band and grouping on it would make the "all" rows the
# complement of the 12-45 rows rather than the whole population.
run_band <- function(bd) {
  d <- if (bd == "all") demo else demo[in_band == TRUE]
  base <- d[usable == TRUE, .(n = .N), by = .(arm, sex, era)]
  setkey(base, arm, sex, era)
  event <- merge(rc[pt %in% TERMS, .(pt, primaryid)],
                 d[, .(primaryid, arm, sex, era, usable)], by = "primaryid")
  event <- event[usable == TRUE, .N, by = .(pt, arm, sex, era)]
  setkey(event, pt, arm, sex, era)
  cell <- function(p, armx, s, ey) {
    n <- base[.(armx, s, ey), n]
    a <- event[.(p, armx, s, ey), N]
    list(n = n, a = if (is.na(a)) 0L else a)
  }
  rbindlist(lapply(TERMS, function(p) rbindlist(lapply(ERAS, function(ey) {
    f <- cell(p, "index", "F", ey); m <- cell(p, "index", "M", ey)
    cf <- cell(p, "comparator", "F", ey); cm <- cell(p, "comparator", "M", ey)
    if (any(c(f$a, m$a, cf$a, cm$a) == 0)) return(NULL)
    ror_f <- (f$a / (f$n - f$a)) / (cf$a / (cf$n - cf$a))
    ror_m <- (m$a / (m$n - m$a)) / (cm$a / (cm$n - cm$a))
    rorr  <- ror_f / ror_m
    se <- sqrt(sum(1 / c(f$a, f$n - f$a, cf$a, cf$n - cf$a,
                         m$a, m$n - m$a, cm$a, cm$n - cm$a)))
    data.table(pt = p, era = ey, band = bd, a_F = f$a, a_M = m$a, n_iF = f$n,
               n_iM = m$n, n_cF = cf$n, n_cM = cm$n, c_F = cf$a, c_M = cm$a,
               RORR = rorr, lci = exp(log(rorr) - 1.96 * se),
               uci = exp(log(rorr) + 1.96 * se), se_log = se)
  }))))
}

res <- rbindlist(lapply(c("all", "12-45"), run_band))

res[, `:=`(lci = round(lci, 2), uci = round(uci, 2), RORR = round(RORR, 3))]
fwrite(res, file.path(TAB, "time_stability_rorr.csv"))

# ------------------------------------------------------ heterogeneity by term
het <- res[, {
  w <- 1 / se_log^2
  pooled <- sum(w * log(RORR)) / sum(w)
  se_p <- sqrt(1 / sum(w))
  Q <- sum(w * (log(RORR) - pooled)^2)
  list(eras = .N, pooled_RORR = round(exp(pooled), 3),
       pooled_lci = round(exp(pooled - 1.96 * se_p), 3),
       pooled_uci = round(exp(pooled + 1.96 * se_p), 3),
       Q = round(Q, 2), df = .N - 1L, p_heterogeneity = round(pchisq(Q, .N - 1L,
                                                                     lower.tail = FALSE), 3))
}, by = .(pt, band)]
fwrite(het, file.path(TAB, "time_stability_heterogeneity.csv"))

cat("\n== ratio of ratios by era, sex-neutral comparator ==\n")
print(dcast(res[band == "all", .(pt, era, v = sprintf("%.2f (%.2f-%.2f)", RORR, lci, uci))],
             pt ~ era, value.var = "v"), nrows = 20)
cat("\n== ages 12-45 only ==\n")
print(dcast(res[band == "12-45", .(pt, era, v = sprintf("%.2f (%.2f-%.2f)", RORR, lci, uci))],
             pt ~ era, value.var = "v"), nrows = 20)
cat("\n== heterogeneity across eras ==\n"); print(het)

# ------------------------------- finasteride's share of male depression by era
# Measured in the untrimmed pool: after trimming there is no finasteride left to
# count, and the question here is how much contamination each era carried.
fin_ids <- ps[drug_group == "finasteride", unique(primaryid)]
dep_ids <- rc[pt == "Depression", unique(primaryid)]
demo[, band := fifelse(is.na(age_yr), "unknown",
                fifelse(age_yr < 12, "<12", fifelse(age_yr <= 45, "12-45", "46+")))]

fin_era <- rbindlist(lapply(c("all", "12-45"), function(bx) {
  d <- if (bx == "all") demo[arm == "comparator" & sex == "M"]
       else demo[arm == "comparator" & sex == "M" & band == "12-45"]
  d[, {
    ids <- primaryid
    dp <- intersect(ids, dep_ids)
    list(band = bx,
         n_comparator = .N,
         n_depression = length(dp),
         n_depression_finasteride = sum(dp %in% fin_ids),
         n_reports_finasteride = sum(ids %in% fin_ids))
  }, by = era]
}))[order(band, era)]
fin_era[, pct_depression_on_finasteride := round(100 * n_depression_finasteride /
                                                   n_depression, 1)]
fwrite(fin_era, file.path(TAB, "finasteride_share_by_era.csv"))
cat("\n== finasteride share of male comparator depression reports, by era ==\n")
print(fin_era)

cat("\n[written] time_stability_rorr.csv, time_stability_heterogeneity.csv,\n")
cat("          finasteride_share_by_era.csv\n[done]\n")
