#!/usr/bin/env Rscript
# ==============================================================================
# Does the iPLEDGE monitoring channel explain the male predominance?
#
# Isotretinoin is dispensed under a pregnancy-prevention programme, so a share of
# the reports from women are monitoring reports - pregnancy tests, contraception,
# lipids, liver function - with no psychiatric event on them. Those reports enter
# the female denominator and add no depression cases, which lowers the female
# reporting odds and therefore makes the ratio of ratios look more male. The
# direction of that bias favours the paper's conclusion, which is exactly why it
# has to be measured rather than assumed away.
#
# The exclusion rule is applied to both arms symmetrically: a report is dropped if
# every retained adverse-event term on it is a monitoring term.
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
# terms that a pregnancy-prevention programme generates
MON <- paste0("pregnan|contracep|condom|birth control|abortion|foetal|fetal|",
              "menstrua|blood test|liver function|hepatic function|lipid|",
              "cholesterol|triglycerid|laborator|serum|haematolog|hematolog")

INDEX <- "isotretinoin"
MIN_A <- 3L; MIN_PT <- 50L
MIN_DRUG_N <- 200L; LO_SX <- 20; HI_SX <- 80
TERMS <- c("Depression", "Anxiety", "Suicidal ideation")

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
# a boolean, not a label: labelling the complement "all" invites the mistake
# of grouping on it and calling the complement the whole population
demo[, in_band := !is.na(age_yr) & age_yr >= 12 & age_yr <= 45]

# trimming rule, fixed over the whole period as in the primary analysis
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
allps <- ps[, .(bad = sum(!(drug_group %in% keep_drug))), by = primaryid]
demo <- merge(demo, allps, by = "primaryid", all.x = TRUE)
demo[, usable := arm == "index" | (is.na(bad) | bad == 0)]

# ------------------------------------------------- monitoring-only flag
rc2 <- merge(rc[, .(primaryid, pt)], demo[, .(primaryid, arm, sex)],
             by = "primaryid")
rc2[, mon := grepl(MON, pt, ignore.case = TRUE)]
flag <- rc2[, .(n_terms = .N, n_mon = sum(mon)), by = primaryid]
flag[, monitoring_only := n_terms == n_mon]
demo <- merge(demo, flag[, .(primaryid, monitoring_only)], by = "primaryid",
              all.x = TRUE)
demo[is.na(monitoring_only), monitoring_only := FALSE]

cat("\n== monitoring-only reports, by arm and sex ==\n")
cnt <- demo[usable == TRUE, .(n = .N, mon = sum(monitoring_only),
                              pct = round(100 * sum(monitoring_only) / .N, 1)),
            by = .(arm, sex)][order(arm, sex)]
print(cnt)
fwrite(cnt, file.path(TAB, "monitoring_report_counts.csv"))

# ------------------------------------------------------------- estimator
est <- function(d, label) {
  d <- d[usable == TRUE]
  rcc <- rc[primaryid %in% d$primaryid]
  freq <- rcc[, .(tot = .N), by = pt][tot >= MIN_PT, pt]
  rcc <- rcc[pt %in% freq]
  out <- rbindlist(lapply(TERMS, function(ev) {
    pid <- rcc[pt == ev, unique(primaryid)]
    one <- function(sx) {
      ni <- d[arm == "index" & sex == sx, .N]
      ai <- d[arm == "index" & sex == sx & primaryid %in% pid, .N]
      nc <- d[arm == "comparator" & sex == sx, .N]
      ac <- d[arm == "comparator" & sex == sx & primaryid %in% pid, .N]
      if (min(ai, ni - ai, ac, nc - ac) <= 0) return(NULL)
      A <- as.numeric(ai); B <- as.numeric(ni - ai)
      C <- as.numeric(ac); D <- as.numeric(nc - ac)
      list(a = ai, n = ni, c = ac, ln = log((A * D) / (B * C)),
           se = sqrt(1/A + 1/B + 1/C + 1/D))
    }
    f <- one("F"); m <- one("M")
    if (is.null(f) || is.null(m) || min(f$a, m$a) < MIN_A) return(NULL)
    lo <- f$ln - m$ln; se <- sqrt(f$se^2 + m$se^2)
    data.table(pt = ev, scenario = label, a_F = f$a, a_M = m$a,
               n_iF = f$n, n_iM = m$n,
               RORR = round(exp(lo), 3), lci = round(exp(lo - 1.96 * se), 3),
               uci = round(exp(lo + 1.96 * se), 3), c_F = f$c, c_M = m$c)
  }))
  out
}

res_all <- rbind(est(demo, "all reports"),
                 est(demo[monitoring_only == FALSE], "monitoring-only reports excluded"))
res_all[, band := "all"]
res_45 <- rbind(est(demo[in_band == TRUE], "all reports"),
                est(demo[in_band == TRUE & monitoring_only == FALSE],
                    "monitoring-only reports excluded"))
res_45[, band := "12-45"]
out <- rbind(res_all, res_45)
setcolorder(out, c("pt", "band", "scenario", "a_F", "a_M", "n_iF", "n_iM",
                   "RORR", "lci", "uci"))
fwrite(out, file.path(TAB, "monitoring_sensitivity.csv"))

cat("\n== ratio of ratios with and without monitoring-only reports ==\n")
print(dcast(out[, .(pt, band, scenario,
                    v = sprintf("%.2f (%.2f-%.2f)", RORR, lci, uci))],
             pt + band ~ scenario, value.var = "v"), nrows = 20)

cat("\n[written] monitoring_sensitivity.csv, monitoring_report_counts.csv\n[done]\n")
