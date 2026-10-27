#!/usr/bin/env Rscript
# ==============================================================================
# An indication-matched active comparator: tetracyclines in acne.
#
# The primary comparator is every other dermatologic drug, which treats psoriasis,
# eczema, urticaria and so on. If the comparator's psychiatric reporting differs
# because its patients differ, the ratio of ratios inherits that. The cleanest
# answer available in FAERS is to compare isotretinoin with the other systemic
# drugs given for the same disease: the tetracyclines, restricted to reports that
# carry an acne indication.
#
# The comparison is small and nobody should read a precise effect from it. It is
# worth running because it is the one design here in which the indication is
# matched by construction rather than argued.
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
ACTIVE <- c("tetracycline", "doxycycline", "minocycline")
TERMS <- c("Depression", "Anxiety", "Suicidal ideation", "Completed suicide",
           "Major depression", "Anhedonia", "Loss of libido", "Dry skin", "Acne")

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

# acne indication, matched across arms
acne_ids <- unique(indi[primaryid %in% demo$primaryid &
                          grepl("acne", indi_pt, ignore.case = TRUE), primaryid])
cat("reports carrying an acne indication:", length(acne_ids), "\n")

arm <- ps[, .(drugs = list(drug_group)), by = primaryid]
arm[, index := vapply(drugs, function(x) INDEX %in% x, logical(1))]
arm[, active := vapply(drugs, function(x) any(ACTIVE %in% x), logical(1))]
arm <- arm[primaryid %in% acne_ids]
# a report naming both would contaminate either arm; there should be none
cat("reports naming both arms:", arm[index == TRUE & active == TRUE, .N], "\n")
arm <- arm[xor(index, active)]
arm[, grp := fifelse(index, "index", "active")]
d <- merge(demo, arm[, .(primaryid, grp)], by = "primaryid")
print(d[, .N, by = .(grp, sex)])

est <- function(ev) {
  pid <- rc[pt == ev, unique(primaryid)]
  one <- function(sx) {
    ni <- d[grp == "index" & sex == sx, .N]
    ai <- d[grp == "index" & sex == sx & primaryid %in% pid, .N]
    nc <- d[grp == "active" & sex == sx, .N]
    ac <- d[grp == "active" & sex == sx & primaryid %in% pid, .N]
    if (min(ai, ni - ai, ac, nc - ac) <= 0) return(NULL)
    A <- as.numeric(ai); B <- as.numeric(ni - ai)
    C <- as.numeric(ac); D <- as.numeric(nc - ac)
    list(a = ai, n = ni, ca = ac, cn = nc,
         ror = (A * D) / (B * C), ln = log((A * D) / (B * C)),
         se = sqrt(1/A + 1/B + 1/C + 1/D))
  }
  f <- one("F"); m <- one("M")
  if (is.null(f) || is.null(m)) return(NULL)
  lo <- f$ln - m$ln; se <- sqrt(f$se^2 + m$se^2)
  data.table(pt = ev,
             a_F = f$a, c_F = f$ca, n_iF = f$n, n_cF = f$cn, ROR_F = f$ror,
             a_M = m$a, c_M = m$ca, n_iM = m$n, n_cM = m$cn, ROR_M = m$ror,
             RORR = exp(lo), lci = exp(lo - 1.96 * se), uci = exp(lo + 1.96 * se),
             p = 2 * pnorm(-abs(lo / se)))
}

res <- rbindlist(lapply(TERMS, est))
res[, `:=`(ROR_F = round(ROR_F, 2), ROR_M = round(ROR_M, 2),
           RORR = round(RORR, 3), lci = round(lci, 3), uci = round(uci, 3))]
fwrite(res, file.path(TAB, "active_comparator_acne.csv"))
cat("\n== isotretinoin versus tetracyclines, acne reports only ==\n")
print(res[, .(pt, a_F, a_M, c_F, c_M, n_iF, n_iM, RORR, lci, uci, p)], nrows = 20)

cat("\n[written] active_comparator_acne.csv\n[done]\n")
