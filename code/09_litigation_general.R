#!/usr/bin/env Rscript
# ==============================================================================
# Does the litigation channel contaminate sex contrasts beyond isotretinoin?
#
# The second Threat to a within-sex disproportionality contrast is that the two
# sexes reach the regulator by different routes. TableS9 showed lawyer-filed
# reports are 25.0% of male versus 15.1% of female isotretinoin reports, and
# removing them erased the apparent male predominance for Crohn's disease.
#
# Either that was an artefact specific to one litigated drug, or it is a general
# property of the database. This script settles which, by repeating the entire
# drug x event scan with lawyer-filed reports removed and asking how many of the
# 2,981 significant sex interactions survive.
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

MIN_A  <- 3L     # matches 03_analysis.R
MIN_PT <- 50L

demo <- fread(file.path(PROC, "demo.csv.gz"))
drug <- fread(file.path(PROC, "drug.csv.gz"))
reac <- fread(file.path(PROC, "reac.csv.gz"))
indi <- fread(file.path(PROC, "indi.csv.gz"))

demo <- demo[sex %in% c("M", "F"), .(primaryid, sex, occp_cod)]
keep <- intersect(demo$primaryid,
                  unique(indi[grepl(DERM, indi_pt, ignore.case = TRUE), primaryid]))
demo <- demo[primaryid %in% keep]
# uniqueness is taken on (primaryid, drug_group) exactly as in 03_analysis.R,
# with the class label carried as a separate lookup
drug_ps <- unique(drug[primaryid %in% keep & toupper(role_cod) == "PS",
                       .(primaryid, drug_group)])
CAT_DG <- unique(drug[, .(drug_group, category)])
reac <- unique(reac[primaryid %in% keep][
  !is.na(pt) & pt != "" & !grepl(NON_AE, pt, ignore.case = TRUE), .(primaryid, pt)])
cat("cohort:", nrow(demo), "reports\n")

demo[, lawyer := toupper(occp_cod) == "LW"]

# ============================================================ who reaches FDA
cat("\n== litigation channel in the whole dermatologic cohort ==\n")
lw_tab <- demo[, .(reports = .N, lawyer = sum(lawyer)), by = sex]
lw_tab[, pct_lawyer := round(100 * lawyer / reports, 1)]
print(lw_tab)
fwrite(lw_tab, file.path(SUP, "TableS13_lawyer_channel_by_sex.csv"))

# by therapeutic class, on the reports' primary suspect
ds <- drug_ps[demo, on = "primaryid", nomatch = NULL]
ds <- CAT_DG[ds, on = "drug_group", nomatch = NULL][
  , .(primaryid, sex, lawyer, category)]
ds <- unique(ds)
by_cat <- ds[, .(reports = .N, lawyer = sum(lawyer)), by = .(category, sex)]
by_cat[, pct_lawyer := round(100 * lawyer / reports, 1)]
by_cat <- by_cat[reports >= 500]
setorder(by_cat, category, sex)
fwrite(by_cat, file.path(SUP, "TableS14_lawyer_channel_by_category_sex.csv"))
cat("\nby therapeutic class (classes with >=500 per sex):\n")
print(by_cat, nrows = 40)

# ================================================ full scan without lawyers
# ------------------------------------------------- frequency floor (see below)
scan <- function(dv, label, pt_keep = NULL) {
  dp <- drug_ps[primaryid %in% dv$primaryid]
  rc <- reac[primaryid %in% dv$primaryid]
  de <- dp[, .(primaryid, drug_group)][rc, on = "primaryid",
                                       allow.cartesian = TRUE, nomatch = NULL][
    dv[, .(primaryid, sex)], on = "primaryid", nomatch = NULL]
  de <- unique(de[, .(primaryid, sex, drug_group, pt)])
  N_sex <- dv[, .(N = .N), by = sex]
  Nd  <- de[, .(Nd = uniqueN(primaryid)), by = .(drug_group, sex)]
  Ne  <- de[, .(Ne = uniqueN(primaryid)), by = .(pt, sex)]
  Nde <- de[, .(a = uniqueN(primaryid)), by = .(drug_group, pt, sex)]

  build <- function(sx) {
    dd <- Nde[sex == sx][Nd[sex == sx], on = "drug_group", nomatch = NULL][
      Ne[sex == sx], on = "pt", nomatch = NULL]
    dd[, `:=`(b = Nd - a, c = Ne - a, dcell = N_sex[sex == sx, N] - Nd - Ne + a)]
    dd <- dd[a >= MIN_A & b >= 0 & c >= 0 & dcell >= 0]
    dd[, `:=`(A = as.numeric(a), B = as.numeric(b),
              C = as.numeric(c), D = as.numeric(dcell))]
    dd[, `:=`(ror = (A * D) / (B * C), se = sqrt(1/A + 1/B + 1/C + 1/D))]
    dd[, .(drug_group, pt, a, ror, ln = log(ror), se)]
  }
  # frequency floor taken from the UNRESTRICTED cohort, exactly as 03_analysis.R
  # does, and reused for the restricted scan. Deriving it separately inside each
  # scan would quietly drop different events in each arm and break the pairing.
  if (is.null(pt_keep)) {
    pt_keep <- de[, .(tot = uniqueN(primaryid)), by = pt][tot >= MIN_PT, pt]
  }
  Fs <- build("F")[pt %in% pt_keep]
  Ms <- build("M")[pt %in% pt_keep]
  R  <- merge(Fs, Ms, by = c("drug_group", "pt"), suffixes = c("_F", "_M"))
  R[, `:=`(lorr = ln_F - ln_M, se_r = sqrt(se_F^2 + se_M^2))]
  R[, `:=`(rorr = exp(lorr), rorr_lci = exp(lorr - 1.96 * se_r),
           rorr_uci = exp(lorr + 1.96 * se_r))]
  R[, p_int := 2 * pnorm(-abs(lorr / se_r))]
  R[, q_int := p.adjust(p_int, method = "fdr")]
  sig <- R[p_int < 0.05 & (rorr_uci < 1 | rorr_lci > 1) &
             a_F >= MIN_A & a_M >= MIN_A & (a_F + a_M) >= 20]
  cat(sprintf("[%s] pairs=%d  significant=%d\n", label, nrow(R), nrow(sig)))
  list(R = R, sig = sig, pt_keep = pt_keep)
}

cat("\n== full drug x event scan, with and without litigation reports ==\n")
main <- scan(demo, "all reports")
clean <- scan(demo[lawyer == FALSE], "no lawyer", pt_keep = main$pt_keep)

# ================================================== how many signals survive
key <- function(dt) dt[, .(drug_group, pt)]
s1 <- main$sig[, .(drug_group, pt, rorr_main = rorr,
                   lci_main = rorr_lci, uci_main = rorr_uci)]
comp <- clean$R[key(s1), on = c("drug_group", "pt"),
                .(drug_group, pt, rorr_clean = rorr, lci_clean = rorr_lci,
                  uci_clean = rorr_uci, p_clean = p_int)]
s1 <- comp[s1, on = c("drug_group", "pt")]
# use the columns sitting on s1 itself: referring back to comp$ after the join
# pairs rows in the wrong order, which silently zeroed this flag
s1[, survived := p_clean < 0.05 & (uci_clean < 1 | lci_clean > 1) &
     !is.na(rorr_clean)]
s1[, same_direction := (rorr_main > 1) == (rorr_clean > 1)]
CATLBL <- fread(file.path(ROOT, "code", "drug_catalog.csv"))[
  , .(drug_group, drug_label, category)][!duplicated(drug_group)]
s1 <- CATLBL[s1, on = "drug_group", nomatch = NULL]
s1[, `:=`(rorr_main = round(rorr_main, 2), lci_main = round(lci_main, 2),
          uci_main = round(uci_main, 2), rorr_clean = round(rorr_clean, 2),
          lci_clean = round(lci_clean, 2), uci_clean = round(uci_clean, 2),
          p_clean = signif(p_clean, 3))]
s1[, rankd := -abs(log(rorr_main))]
setorder(s1, -survived, rankd)
s1[, rankd := NULL]
fwrite(s1, file.path(SUP, "TableS15_significant_pairs_with_without_lawyers.csv"))

n <- nrow(s1); ns <- s1[, sum(survived, na.rm = TRUE)]
cat(sprintf("\nsignificant sex interactions in the main scan : %d\n", n))
cat(sprintf("still significant without litigation reports   : %d (%.1f%%)\n",
            ns, 100 * ns / n))
cat(sprintf("retaining direction < 1 or > 1                 : %d (%.1f%%)\n",
            s1[, sum(same_direction, na.rm = TRUE)],
            100 * s1[, mean(same_direction, na.rm = TRUE)]))
cat(sprintf("crossed 1 (direction reversed)                 : %d\n",
            s1[, sum(!is.na(rorr_clean) & ((rorr_main > 1) != (rorr_clean > 1)))]))
cat(sprintf("no longer evaluable (too few reports)          : %d\n",
            s1[, sum(is.na(rorr_clean))]))

cat("\n-- largest signals that do NOT survive --\n")
print(head(s1[(!survived) & !is.na(rorr_clean)][
  order(-abs(log(rorr_main)))][
  , .(drug_label, pt, rorr_main, rorr_clean, survived)], 12))
cat("\n-- largest signals that DO survive --\n")
print(head(s1[(survived)][order(-abs(log(rorr_main)))][
  , .(drug_label, pt, rorr_main, rorr_clean)], 8))

# ------------------------------------------------- what kind of signal breaks
GI <- "bowel|colitis|crohn|diarrhoea|diarrhea|enteritis|gastrointestin| ileitis|"
GI <- paste0(GI, "rectal|gluten|colonic|abdominal|dyspepsia|nausea|diverticul")
s1[, is_gi := grepl(GI, pt, ignore.case = TRUE)]
fail <- s1[(!survived) & !is.na(rorr_clean)]
surv <- s1[(survived)]
cat(sprintf("\nsignals lost: %d, of which gastrointestinal terms: %d (%.0f%%)\n",
            nrow(fail), fail[, sum(is_gi)], 100 * fail[, mean(is_gi)]))
cat(sprintf("signals kept: %d, of which gastrointestinal terms: %d (%.0f%%)\n",
            nrow(surv), surv[, sum(is_gi)], 100 * surv[, mean(is_gi)]))
cat(sprintf("\nratio of GI odds: %.1f-fold enrichment among the lost signals\n",
            (fail[, mean(is_gi)] + 1e-9) / (surv[, mean(is_gi)] + 1e-9)))

cat("\n-- lost signals that are NOT gastrointestinal (for completeness) --\n")
print(head(fail[is_gi == FALSE][order(-abs(log(rorr_main)))][
  , .(drug_label, pt, rorr_main, rorr_clean)], 8))

cat("\n[done]\n")

# The README is written in full by 06_supplement.R; see the note in
# 07_reporting_evidence.R. The descriptions for TableS13-S15 live in 06.
if (!file.exists(file.path(SUP, "README.md"))) {
  cat("[warn] README.md missing: run 06_supplement.R first\n")
}
