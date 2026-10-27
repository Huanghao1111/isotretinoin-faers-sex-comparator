#!/usr/bin/env Rscript
# ==============================================================================
# Sex-stratified disproportionality analysis of dermatologic drugs (FAERS)
#
# Design decisions (all pre-specified):
#  1) Cohort is restricted to reports carrying at least one DERMATOLOGIC
#     indication, so TNF-blocker reports from rheumatology/gastroenterology do
#     not contaminate a dermatology-facing analysis.
#  2) Medication-error / product-use preferred terms are removed: they describe
#     handling mistakes, not drug toxicity.
#  3) Core metric: ROR ratio (RORR) = ROR_female / ROR_male. Each ROR is
#     computed WITHIN one sex, so the unknown sex composition of the treated
#     population and overall reporting propensity cancel in the ratio.
#  4) Serious outcomes are compared with a Mantel-Haenszel OR stratified by
#     drug, so the contrast is not driven by which drugs each sex receives.
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

ROOT  <- "C:/Users/lenovo/WorkBuddy/2026-10-17-00-49-45"
PROC  <- file.path(ROOT, "data", "processed")
TAB   <- file.path(ROOT, "output", "tables")
FIG   <- file.path(ROOT, "output", "figures")
for (d in c(TAB, FIG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

MIN_A    <- 3L     # min co-occurrences per sex
MIN_PT   <- 50L    # min PT frequency in the cohort (exploratory scan)

# --------------------------------------------------- indication vocabulary
DERM <- paste0(
  "psoriasis|psoriatic|dermatitis|eczema|acne|urticaria|angioedema|",
  "hidradenitis|alopecia|hair loss|vitiligo|melanoma|basal cell carcinoma|",
  "squamous cell carcinoma|actinic keratosis|skin cancer|skin neoplasm|",
  "rosacea|pemphigus|pemphigoid|lupus|pruritus|itch|onychomycosis|tinea|",
  "lichen planus|morphea|scleroderma|pyoderma|photo|seborrhoe|",
  "hyperhidrosis|cutane|skin|wart|cellulitis|impetigo|scabies|herpes|",
  "erythema nodosum|granuloma annulare|sarcoid|keratosis")

# ----------------------------------------------- non-adverse-event PTs to drop
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

cat("== reading cleaned FAERS tables ==\n")
demo <- fread(file.path(PROC, "demo.csv.gz"))
drug <- fread(file.path(PROC, "drug.csv.gz"))
reac <- fread(file.path(PROC, "reac.csv.gz"))
outc <- fread(file.path(PROC, "outc.csv.gz"))
indi <- fread(file.path(PROC, "indi.csv.gz"))

demo <- demo[, .(primaryid, sex, age_yr, event_year, occp_cod)]
demo <- demo[sex %in% c("M", "F")]
setkey(demo, primaryid)

# ------------------------------------------------- restrict to derm indication
cat("\n== restricting to dermatologic indications ==\n")
indi_u <- unique(indi[!is.na(indi_pt) & indi_pt != "", .(primaryid, indi_pt)])
derm_pid <- indi_u[grepl(DERM, indi_pt, ignore.case = TRUE), unique(primaryid)]
cat("cases with >=1 dermatologic indication:", length(derm_pid), "\n")

keep <- intersect(demo$primaryid, derm_pid)
demo <- demo[primaryid %in% keep]
drug <- drug[primaryid %in% keep]
reac <- reac[primaryid %in% keep]
outc <- outc[primaryid %in% keep]
cat("analysis cohort:", nrow(demo), "\n")

# ------------------------------------------------- tidy drug / reaction
drug_ps <- unique(drug[toupper(role_cod) == "PS", .(primaryid, drug_group, category)])
reac    <- unique(reac[, .(primaryid, pt)])
reac    <- reac[!is.na(pt) & pt != "" & !grepl(NON_AE, pt, ignore.case = TRUE)]
cat("PS drug rows:", nrow(drug_ps), " | reaction rows (non-AE dropped):", nrow(reac), "\n")

cat("\n== Table 1: sex distribution overall ==\n")
t1 <- demo[, .(N = .N), by = sex][, `:=`(pct = round(100 * N / sum(N), 1))]
print(t1)

drug_sex <- drug_ps[demo, on = "primaryid", nomatch = NULL]
t1b <- drug_sex[, .(N = uniqueN(primaryid)), by = .(category, sex)]
t1b <- dcast(t1b, category ~ sex, value.var = "N", fill = 0)
if (!"F" %in% names(t1b)) t1b[, F := 0L]
if (!"M" %in% names(t1b)) t1b[, M := 0L]
setnames(t1b, c("F", "M"), c("n_F", "n_M"))
t1b[, tot := n_F + n_M][, pct_F := round(100 * n_F / tot, 1)]
setorder(t1b, -tot)
print(t1b, nrows = 40)
fwrite(t1b, file.path(TAB, "table1_sex_by_category.csv"))

# ------------------------------------------------- drug x event pairs
cat("\n== building drug x event pairs ==\n")
de <- drug_ps[, .(primaryid, drug_group)][reac, on = "primaryid",
                                          allow.cartesian = TRUE, nomatch = NULL]
de <- demo[, .(primaryid, sex)][de, on = "primaryid", nomatch = NULL]
de <- unique(de[, .(primaryid, sex, drug_group, pt)])
cat("drug-event pairs:", nrow(de), "\n")

N_sex <- demo[, .(N = .N), by = sex]
Nd    <- de[, .(Nd = uniqueN(primaryid)), by = .(drug_group, sex)]
Ne    <- de[, .(Ne = uniqueN(primaryid)), by = .(pt, sex)]
Nde   <- de[, .(a = uniqueN(primaryid)), by = .(drug_group, pt, sex)]

build <- function(sx) {
  dd <- Nde[sex == sx][Nd[sex == sx], on = "drug_group", nomatch = NULL][
    Ne[sex == sx], on = "pt", nomatch = NULL]
  dd[, `:=`(b = Nd - a, c = Ne - a,
            nn = N_sex[sex == sx, N] - Nd - Ne + a)]
  setnames(dd, "nn", "dcell")
  dd <- dd[a >= MIN_A & b >= 0 & c >= 0 & dcell >= 0]
  dd[, `:=`(A = as.numeric(a), B = as.numeric(b),
            C = as.numeric(c), D = as.numeric(dcell))]
  dd[, `:=`(ror = (A * D) / (B * C),
            se  = sqrt(1/A + 1/B + 1/C + 1/D))]
  dd[, `:=`(ln_ror = log(ror),
            lci = exp(log(ror) - 1.96 * se),
            uci = exp(log(ror) + 1.96 * se))]
  dd[, p_ror := 2 * pnorm(-abs(log(ror) / se))]
  dd[, .(drug_group, pt, a, b, c, d = dcell, ror, lci, uci, ln_ror, se, p_ror)]
}

cat("== computing sex-specific RORs ==\n")
F_set <- build("F"); M_set <- build("M")

pt_freq <- de[, .(tot = uniqueN(primaryid)), by = pt]
F_set <- F_set[pt_freq, on = "pt"][tot >= MIN_PT & !is.na(tot)][, tot := NULL]
M_set <- M_set[pt_freq, on = "pt"][tot >= MIN_PT & !is.na(tot)][, tot := NULL]

cat("\n== joining F/M and computing ROR ratio ==\n")
R <- merge(F_set, M_set, by = c("drug_group", "pt"), suffixes = c("_F", "_M"))
R[, `:=`(lorr     = ln_ror_F - ln_ror_M,
         se_r     = sqrt(se_F^2 + se_M^2))]
R[, `:=`(rorr     = exp(lorr),
         rorr_lci = exp(lorr - 1.96 * se_r),
         rorr_uci = exp(lorr + 1.96 * se_r))]
R[, z := lorr / se_r]
R[, p_int := 2 * pnorm(-abs(z))]
R[, q_int := p.adjust(p_int, method = "fdr")]

sig <- R[p_int < 0.05 & (rorr_uci < 1 | rorr_lci > 1) &
           a_F >= MIN_A & a_M >= MIN_A & (a_F + a_M) >= 20]
setorder(sig, p_int)
cat("\nsignificant sex x drug x event interactions:", nrow(sig), "\n")
print(head(sig[, .(drug_group, pt, a_F, a_M, ror_F = round(ror_F, 2),
                   ror_M = round(ror_M, 2), rorr = round(rorr, 2),
                   lci = round(rorr_lci, 2), uci = round(rorr_uci, 2),
                   p = signif(p_int, 3), q = signif(q_int, 3))], 30), nrows = 30)

fwrite(R,   file.path(TAB, "all_ror_sex_stratified.csv.gz"))
fwrite(sig, file.path(TAB, "significant_rorr.csv"))

# ------------------------------------------------- serious outcomes, drug-stratified
cat("\n== serious outcomes: Mantel-Haenszel OR (male vs female), stratified by drug ==\n")
flag <- unique(outc[!is.na(outc_cod), .(primaryid, outc_cod)])
ser  <- dcast(flag, primaryid ~ outc_cod, fun.aggregate = length, value.var = "outc_cod")
ser  <- demo[, .(primaryid, sex)][ser, on = "primaryid", nomatch = NULL]
for (cc in c("DE", "HO", "LT", "DS", "CA", "RI"))
  if (!cc %in% names(ser)) ser[, (cc) := 0L]

# primary drug per case: the most frequently matched PS drug in that report
prim <- drug_ps[, .(n = .N), by = .(primaryid, drug_group)]
setorder(prim, primaryid, -n)
prim <- prim[!duplicated(primaryid), .(primaryid, drug_group)]
ser  <- prim[ser, on = "primaryid", nomatch = NULL]

mh_or <- function(dt, ev) {
  d <- dt[, .(m_n = sum(get(ev) > 0 & sex == "M"),
              m_t = sum(sex == "M"),
              f_n = sum(get(ev) > 0 & sex == "F"),
              f_t = sum(sex == "F")), by = drug_group]
  d <- d[m_n > 0 & f_n > 0 & m_t > m_n & f_t > f_n]
  if (nrow(d) < 3) return(NULL)
  d[, `:=`(a = m_n, b = m_t - m_n, c = f_n, dd = f_t - f_n, n = m_t + f_t)]
  P <- sum((d$a * d$dd) / d$n); Q <- sum((d$b * d$c) / d$n)
  if (Q <= 0) return(NULL)
  or <- P / Q
  # Robins-Breslow-Greenland variance
  s1 <- sum((d$a + d$dd) * (d$a * d$dd) / d$n^2)
  s2 <- sum((d$a * d$dd + d$b * d$c) * (d$a + d$dd) / d$n^2)
  s3 <- sum((d$a * d$dd) * (d$a + d$dd) / d$n^2)
  var <- s1 / (2 * P^2) + s2 / (2 * P * Q) + s3 / (2 * Q^2)
  se <- sqrt(var)
  data.table(outcome = ev, strata = nrow(d),
             OR_M_vs_F = round(or, 3),
             lci = round(exp(log(or) - 1.96 * se), 3),
             uci = round(exp(log(or) + 1.96 * se), 3),
             p = signif(2 * pnorm(-abs(log(or) / se)), 3))
}

cmp <- rbindlist(lapply(c("DE", "HO", "LT", "DS"), function(e) mh_or(ser, e)))
lbl <- c(DE = "Death", HO = "Hospitalisation",
         LT = "Life-threatening", DS = "Disability")
cmp[, outcome := lbl[outcome]]
print(cmp)
fwrite(cmp, file.path(TAB, "serious_outcome_MH_drugstratified.csv"))

# crude (unstratified) numbers for the manuscript text
raw <- ser[, .(n = .N,
               death = round(100 * sum(DE > 0) / .N, 1),
               hosp  = round(100 * sum(HO > 0) / .N, 1),
               lt    = round(100 * sum(LT > 0) / .N, 1),
               disab = round(100 * sum(DS > 0) / .N, 1)), by = sex]
cat("\n-- crude rates (%) by sex --\n"); print(raw)
fwrite(raw, file.path(TAB, "serious_outcome_crude.csv"))

# ------------------------------------------------- top PTs by sex
cat("\n== top PTs by sex ==\n")
for (sx in c("F", "M")) {
  cat("\n--", sx, "--\n")
  tt <- de[sex == sx, .(n = uniqueN(primaryid)), by = pt][order(-n)]
  print(head(tt, 12), nrows = 12)
}

saveRDS(list(R = R, sig = sig, de = de, demo = demo,
             drug_sex = drug_sex, ser = ser, cmp = cmp, t1b = t1b),
        file.path(PROC, "analysis_objects.rds"))
cat("\n[done] analysis complete\n")
