#!/usr/bin/env Rscript
# ==============================================================================
# Sensitivity analyses for the isotretinoin sex-specific signals.
#
# The headline concern: the extreme female predominance for sexual events
# (RORR 23-64) rests on 10-40 reports, and may arise because men already have a
# high baseline reporting rate for sexual dysfunction in FAERS (finasteride,
# antidepressants, ...), which deflates the MALE ROR and inflates the ratio.
# ROR cancels the sex composition of users, but it does NOT cancel sex
# differences in the BASELINE RATE of a given event across the whole database.
#
# Scenarios probe robustness from different angles:
#   A  main analysis (reference)
#   B  single suspected dermatologic drug only  -> removes comedication mixing
#   C  health-professional reports only         -> removes consumer self-report
#   D  age 12-45 years                          -> acne-typical window
#   E  B + C combined                           -> strictest
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "C:/Users/lenovo/WorkBuddy/2026-10-17-00-49-45"
PROC <- file.path(ROOT, "data", "processed")
TAB  <- file.path(ROOT, "output", "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

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

KEY <- c("Depression", "Anxiety", "Suicidal ideation", "Completed suicide",
         "Major depression", "Anhedonia",
         "Loss of libido", "Sexual dysfunction", "Libido decreased",
         "Lip dry", "Crohn's disease", "Inflammatory bowel disease",
         "Dry skin", "Epistaxis", "Acne")

cat("== loading ==\n")
demo <- fread(file.path(PROC, "demo.csv.gz"))
drug <- fread(file.path(PROC, "drug.csv.gz"))
reac <- fread(file.path(PROC, "reac.csv.gz"))
indi <- fread(file.path(PROC, "indi.csv.gz"))

demo <- demo[sex %in% c("M", "F"),
             .(primaryid, sex, age_yr, event_year, occp_cod)]
derm_pid <- unique(indi[!is.na(indi_pt) & indi_pt != ""][
  grepl(DERM, indi_pt, ignore.case = TRUE), primaryid])
keep <- intersect(demo$primaryid, derm_pid)
demo <- demo[primaryid %in% keep]
drug <- drug[primaryid %in% keep]
reac <- reac[primaryid %in% keep]
cat("base cohort:", nrow(demo), "\n")

drug_ps <- unique(drug[toupper(role_cod) == "PS", .(primaryid, drug_group)])
reac    <- unique(reac[!is.na(pt) & pt != "" &
                        !grepl(NON_AE, pt, ignore.case = TRUE), .(primaryid, pt)])
cat("PS drug rows:", nrow(drug_ps), "| reaction rows:", nrow(reac), "\n")

# number of distinct primary-suspect dermatologic drugs per report
ndrug <- drug_ps[, .(ndrug = uniqueN(drug_group)), by = primaryid]

compute <- function(ids, label) {
  dv <- demo[primaryid %in% ids]
  dv <- dv[sex %in% c("M", "F")]
  if (nrow(dv) < 1000) { cat("[skip] too few cases in", label, "\n"); return(NULL) }
  if (uniqueN(dv$sex) < 2) { cat("[skip] single sex in", label, "\n"); return(NULL) }
  dp <- drug_ps[primaryid %in% ids]
  rc <- reac[primaryid %in% ids]

  de <- dp[, .(primaryid, drug_group)][rc, on = "primaryid",
                                       allow.cartesian = TRUE, nomatch = NULL]
  de <- dv[, .(primaryid, sex)][de, on = "primaryid", nomatch = NULL]
  de <- unique(de[, .(primaryid, sex, drug_group, pt)])

  N_sex <- dv[, .(N = .N), by = sex]
  Nd  <- de[, .(Nd = uniqueN(primaryid)), by = .(drug_group, sex)]
  Ne  <- de[, .(Ne = uniqueN(primaryid)), by = .(pt, sex)]
  Nde <- de[, .(a = uniqueN(primaryid)), by = .(drug_group, pt, sex)]

  build <- function(sx) {
    dd <- Nde[sex == sx][Nd[sex == sx], on = "drug_group", nomatch = NULL][
      Ne[sex == sx], on = "pt", nomatch = NULL]
    dd[, `:=`(b = Nd - a, c = Ne - a, nn = N_sex[sex == sx, N] - Nd - Ne + a)]
    setnames(dd, "nn", "dcell")
    dd <- dd[a >= 3 & b >= 0 & c >= 0 & dcell >= 0]
    dd[, `:=`(A = as.numeric(a), B = as.numeric(b),
              C = as.numeric(c), D = as.numeric(dcell))]
    dd[, `:=`(ror = (A * D) / (B * C), se = sqrt(1/A + 1/B + 1/C + 1/D))]
    dd[, .(drug_group, pt, a, ror, ln = log(ror), se)]
  }
  Fs <- build("F"); Ms <- build("M")
  R  <- merge(Fs, Ms, by = c("drug_group", "pt"), suffixes = c("_F", "_M"))
  R[, `:=`(lorr = ln_F - ln_M, se_r = sqrt(se_F^2 + se_M^2))]
  R[, `:=`(rorr = exp(lorr),
           lci = exp(lorr - 1.96 * se_r),
           uci = exp(lorr + 1.96 * se_r))]
  R[, p_int := 2 * pnorm(-abs(lorr / se_r))]

  out <- R[drug_group == "isotretinoin" & pt %in% KEY,
           .(pt, a_F, a_M, ror_F = round(ror_F, 2), ror_M = round(ror_M, 2),
             rorr = round(rorr, 2), lci = round(lci, 2), uci = round(uci, 2),
             p = signif(p_int, 3))]
  # `sex` does not survive the F/M merge; cohort sizes must come from `dv`
  if (nrow(out)) out[, `:=`(scenario = label, n_cohort = nrow(dv),
                            n_women = dv[sex == "F", .N])]
  cat(sprintf("[%-16s] cohort=%-9s events=%d\n", label,
              format(nrow(dv), big.mark = ","), nrow(out)))
  out
}

single <- ndrug[ndrug == 1, primaryid]
hp     <- demo[toupper(occp_cod) %in% c("MD", "PH", "HP", "RN"), primaryid]
young  <- demo[!is.na(age_yr) & age_yr >= 12 & age_yr <= 45, primaryid]

# Scenario F: drop lawyer-mediated reports. Isotretinoin carries a long US
# litigation history over inflammatory bowel disease, and litigation files at a
# different rate in men and women, so this is not a hypothetical concern.
nolaw  <- demo[toupper(occp_cod) != "LW", primaryid]

cat("\n== running scenarios ==\n")
res <- rbindlist(list(
  compute(demo$primaryid,                          "A_main"),
  compute(single,                                  "B_single_drug"),
  compute(hp,                                      "C_healthprof"),
  compute(young,                                   "D_age12_45"),
  compute(intersect(single, hp),                   "E_single+hp"),
  compute(nolaw,                                   "F_no_lawyer"),
  compute(intersect(nolaw, single),                "G_no_lawyer_single")
))

fwrite(res, file.path(TAB, "sensitivity_isotretinoin.csv"))

cat("\n================ ApJ sensitivity table ================\n")
disp <- dcast(res[, .(scenario, pt, rorr_when = sprintf("%.2f [%.2f-%.2f]", rorr, lci, uci))],
              pt ~ scenario, value.var = "rorr_when")
ord <- c("Depression", "Anxiety", "Suicidal ideation", "Completed suicide",
         "Major depression", "Anhedonia",
         "Loss of libido", "Sexual dysfunction", "Libido decreased",
         "Lip dry", "Crohn's disease", "Inflammatory bowel disease",
         "Dry skin", "Epistaxis", "Acne")
disp[, pt := factor(pt, levels = ord)]
setorder(disp, pt)

# how much of the sex difference lives in litigation-driven reports
cat("\n== share of isotretinoin reports filed by lawyers, by sex ==\n")
iso_p <- demo[drug_ps[drug_group == "isotretinoin"], on = "primaryid", nomatch = NULL]
iso_p[, is_lawyer := toupper(occp_cod) == "LW"]
print(iso_p[, .(n_reports = .N,
                pct_lawyer = round(100 * sum(is_lawyer) / .N, 1)), by = sex])
print(disp, nrows = 30)

cat("\n--- interpretive summary ---\n")
for (ev in KEY) {
  s <- res[pt == ev][order(scenario)]
  if (nrow(s) == 0) next
  rr <- s$rorr
  cat(sprintf("%-28s main=%6.2f  range across scenarios = %6.2f to %6.2f  %s\n",
              ev, s[scenario == "A_main", rorr], min(rr), max(rr),
              ifelse(all(rr > 1) | all(rr < 1), "direction STABLE", "DIRECTION VARIES")))
}
cat("\n[done]\n")
