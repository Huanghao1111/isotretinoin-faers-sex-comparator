#!/usr/bin/env Rscript
# ==============================================================================
# Rebuilds the 2x2 tables as a proper contingency table.
#
# The original implementation drew its cells from three differently restricted
# denominators:
#   a  reports with the index drug as primary suspect AND the event
#   b  index-drug reports MINUS a, among reports carrying >=1 retained term
#   c  event reports MINUS a, among reports carrying >=1 primary suspect drug
#   d  everything left over from the full cohort
#
# Cells b and c therefore exclude different sets of reports, and the reports
# they exclude land in d. That is not a partition of one population, and the
# exclusion is not sex-symmetric: 17.0% of depression reports from women carry
# no primary-suspect dermatologic drug, against 2.6% from men, which pushes the
# depression ROR ratio from 0.82 to 0.98 and manufactures the null result.
#
# Corrected: restrict the analysis population first, then partition it.
#   population = reports with >=1 primary suspect dermatologic drug AND >=1
#                retained adverse event term
#   exposed    = index drug is the primary suspect
#   comparator = another dermatologic drug is the primary suspect
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

MIN_A  <- 3L
MIN_PT <- 50L
INDEX  <- "isotretinoin"

demo <- fread(file.path(PROC, "demo.csv.gz"))
drug <- fread(file.path(PROC, "drug.csv.gz"))
reac <- fread(file.path(PROC, "reac.csv.gz"))
indi <- fread(file.path(PROC, "indi.csv.gz"))

demo <- demo[sex %in% c("M", "F"), .(primaryid, sex, age_yr, occp_cod)]
derm_pid <- unique(indi[grepl(DERM, indi_pt, ignore.case = TRUE), primaryid])
demo <- demo[primaryid %in% derm_pid]
cat("dermatologic-indication reports with known sex:", nrow(demo), "\n")

ps <- unique(drug[primaryid %in% demo$primaryid & toupper(role_cod) == "PS",
                  .(primaryid, drug_group)])
rc <- unique(reac[primaryid %in% demo$primaryid &
                    !is.na(pt) & pt != "" &
                    !grepl(NON_AE, pt, ignore.case = TRUE), .(primaryid, pt)])

# ---------------------------------------------------------------- population
pop <- intersect(ps$primaryid, rc$primaryid)   # PS drug AND retained term
dv <- demo[primaryid %in% pop]
ps <- ps[primaryid %in% pop]
rc <- rc[primaryid %in% pop]
cat("analysis population (PS dermatologic drug AND retained term):",
    nrow(dv), "\n")
print(dv[, .(n = .N), by = sex])

# one exposed arm per report: index drug if present, else comparator
arm <- ps[, .(primaryid, drug_group)]
arm <- arm[, .(arm = if (any(drug_group == INDEX)) "index" else "comparator"),
           by = primaryid]
arm <- arm[!duplicated(primaryid)]
dv <- arm[dv, on = "primaryid", nomatch = NULL]
cat("index arm:", dv[arm == "index", .N],
    "| comparator arm:", dv[arm == "comparator", .N], "\n")
print(dv[, .(n = .N), by = .(arm, sex)])

# ------------------------------------------------------- full scan of events
pt_freq <- rc[, .(tot = .N), by = pt][tot >= MIN_PT]
cat("events meeting the", MIN_PT, "-report floor:", nrow(pt_freq), "\n")

compute <- function(dvv, rcc, pss, label) {
  rcc <- rcc[pt %in% pt_freq$pt]
  N_sex <- dvv[, .(N = .N), by = sex]
  Nd <- dvv[arm == "index", .(Nd = .N), by = sex]
  Ne <- rcc[, .(pt, primaryid)][dvv[, .(primaryid, sex)], on = "primaryid",
                                nomatch = NULL][, .(Ne = .N), by = .(pt, sex)]
  a <- rcc[, .(pt, primaryid)][dvv[arm == "index", .(primaryid, sex)],
                               on = "primaryid", nomatch = NULL][
    , .(a = .N), by = .(pt, sex)]

  one <- function(sx) {
    dd <- a[sex == sx][Ne[sex == sx], on = "pt", nomatch = NULL]
    dd[, `:=`(Nd = Nd[sex == sx, Nd], N = N_sex[sex == sx, N])]
    dd[, `:=`(b = Nd - a, c = Ne - a, d = N - Nd - Ne + a)]
    dd <- dd[a >= MIN_A & b > 0 & c > 0 & d > 0]
    dd[, `:=`(A = as.numeric(a), B = as.numeric(b),
              C = as.numeric(c), D = as.numeric(d))]
    dd[, `:=`(ror = (A * D) / (B * C), se = sqrt(1/A + 1/B + 1/C + 1/D))]
    dd[, .(pt, a, ror, ln = log(ror), se)]
  }
  Fs <- one("F"); Ms <- one("M")
  R <- merge(Fs, Ms, by = "pt", suffixes = c("_F", "_M"))
  R[, `:=`(lorr = ln_F - ln_M, se_r = sqrt(se_F^2 + se_M^2))]
  R[, `:=`(rorr = exp(lorr), rorr_lci = exp(lorr - 1.96 * se_r),
           rorr_uci = exp(lorr + 1.96 * se_r))]
  R[, p_int := 2 * pnorm(-abs(lorr / se_r))]
  R[, q_int := p.adjust(p_int, method = "fdr")]
  cat(sprintf("[%s] events=%d  significant=%d\n", label, nrow(R),
              R[p_int < 0.05 & (rorr_uci < 1 | rorr_lci > 1) &
                  a_F >= MIN_A & a_M >= MIN_A, .N]))
  R
}

R <- compute(dv, rc, ps, "corrected, all ages")
young <- dv[!is.na(age_yr) & age_yr >= 12 & age_yr <= 45, primaryid]
Ry <- compute(dv[primaryid %in% young], rc[primaryid %in% young],
              ps[primaryid %in% young], "corrected, ages 12-45")

KEY <- c("Depression", "Anxiety", "Suicidal ideation", "Completed suicide",
         "Major depression", "Anhedonia", "Loss of libido",
         "Sexual dysfunction", "Libido decreased", "Lip dry",
         "Crohn's disease", "Inflammatory bowel disease", "Dry skin",
         "Epistaxis", "Acne")
show <- function(x, ttl) {
  cat("\n==", ttl, "==\n")
  d <- x[pt %in% KEY, .(pt, a_F, a_M, ROR_F = round(ror_F, 2),
                        ROR_M = round(ror_M, 2), RORR = round(rorr, 2),
                        lci = round(rorr_lci, 2), uci = round(rorr_uci, 2),
                        p = signif(p_int, 3))]
  setorder(d, RORR)
  print(d, nrows = 20)
}
show(R, "corrected ROR ratio, all ages")
show(Ry, "corrected ROR ratio, ages 12-45")

# side by side with what the previous implementation produced
old <- fread(file.path(TAB, "all_ror_sex_stratified.csv.gz"))
old <- old[drug_group == "isotretinoin" & pt %in% KEY,
           .(pt, RORR_old = round(rorr, 2), lci_old = round(rorr_lci, 2),
             uci_old = round(rorr_uci, 2))]
new <- R[pt %in% KEY, .(pt, RORR_new = round(rorr, 2),
                        lci_new = round(rorr_lci, 2),
                        uci_new = round(rorr_uci, 2))]
cmp <- merge(old, new, by = "pt")
cmp[, flipped := (RORR_old > 1) != (RORR_new > 1)]
cat("\n== old versus corrected ==\n")
print(cmp[order(-abs(log(RORR_new / RORR_old)))], nrows = 20)
cat("\nsignals whose conclusion about direction changed:", cmp[, sum(flipped)], "\n")

# ------------------------------------------------------------- decomposition
# ROR ratio = (index odds, women / men) x (comparator odds, men / women).
# The first factor is what the question is about; the second is an artefact of
# how frequently each sex reports that event across the comparator pool. Any
# ratio whose second factor is far from 1 cannot be read as a drug effect.
decomp <- function(dvv, rcc, pss) {
  rcc <- rcc[pt %in% pt_freq$pt]
  idx <- dvv[arm == "index", .(primaryid, sex)]
  cmp <- dvv[arm == "comparator", .(primaryid, sex)]
  out <- rbindlist(lapply(KEY, function(ev) {
    pid <- rcc[pt == ev, unique(primaryid)]
    # counts are kept alongside the odds so the figure can draw an interval
    od <- function(x, sx) {
      n <- x[sex == sx, .N]; a <- x[sex == sx & primaryid %in% pid, .N]
      if (n - a <= 0) return(NULL)
      list(n = n, a = a, odds = a / (n - a))
    }
    fi <- od(idx, "F"); mi <- od(idx, "M")
    fc <- od(cmp, "F"); mc <- od(cmp, "M")
    if (any(vapply(list(fi, mi, fc, mc), is.null, TRUE)) || fc$odds == 0)
      return(NULL)
    data.table(pt = ev,
               n_index_F = fi$n, a_index_F = fi$a,
               n_index_M = mi$n, a_index_M = mi$a,
               n_comp_F  = fc$n, a_comp_F  = fc$a,
               n_comp_M  = mc$n, a_comp_M  = mc$a,
               index_odds_female = round(1000 * fi$odds, 3),
               index_odds_male   = round(1000 * mi$odds, 3),
               comparator_odds_female = round(1000 * fc$odds, 3),
               comparator_odds_male   = round(1000 * mc$odds, 3),
               index_OR_F_vs_M   = round(fi$odds / mi$odds, 2),
               comparator_OR_M_vs_F = round(mc$odds / fc$odds, 2))
  }))
  out[, product := round(index_OR_F_vs_M * comparator_OR_M_vs_F, 2)]
  rr <- R[pt %in% KEY, .(pt, rorr = round(rorr, 2))]
  out <- rr[out, on = "pt"]
  out[, drug_specific := abs(log(index_OR_F_vs_M)) >= abs(log(comparator_OR_M_vs_F))]
  setorder(out, -comparator_OR_M_vs_F)
  out
}
dc <- decomp(dv, rc, ps)
cat("\n== decomposition of the ROR ratio (corrected population) ==\n")
cat("product should equal ROR ratio; a comparator factor far from 1 means the\n")
cat("ratio is mostly reporting context rather than drug effect\n\n")
print(dc, nrows = 20)
fwrite(dc, file.path(TAB, "decomposition_rorr.csv"))

fwrite(R,  file.path(TAB, "corrected_rorr_all_ages.csv"))
fwrite(Ry, file.path(TAB, "corrected_rorr_age12_45.csv"))
fwrite(cmp, file.path(TAB, "corrected_vs_previous.csv"))
cat("\n[done]\n")
