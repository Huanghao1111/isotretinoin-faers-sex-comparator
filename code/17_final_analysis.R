#!/usr/bin/env Rscript
# ==============================================================================
# Final primary analysis, after two design faults were found and fixed.
#
# Fault 1 (Simpson). The comparator pool "all other dermatologic drugs" has a
# completely different age structure from the index arm (39.9% vs 2.9% aged 46+),
# so every age stratum gives a depression RORR above 1 while the pool gives 0.82.
#
# Fault 2 (sex-specific comparator drugs). The male comparator pool contains
# finasteride: 2.8% of its reports, but 56% of its depression reports (302.9 per
# 1,000 against 2.4 for dupilumab). That single drug inflates the male comparator
# depression rate and manufactures a female signal for the index drug.
#
# Both are design faults, not estimation faults, so both are fixed in the design:
#   - the comparator pool is restricted to drugs whose reports are not
#     overwhelmingly from one sex (20-80% female among drugs with >=200 reports);
#   - the pooled estimate is reported alongside the age-stratified one.
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

MIN_A <- 3L; MIN_PT <- 50L; INDEX <- "isotretinoin"
MIN_DRUG_N <- 200L; LO_SX <- 20; HI_SX <- 80

KEY <- c("Depression", "Anxiety", "Suicidal ideation", "Completed suicide",
         "Major depression", "Anhedonia", "Loss of libido",
         "Sexual dysfunction", "Libido decreased", "Lip dry",
         "Crohn's disease", "Inflammatory bowel disease", "Dry skin",
         "Epistaxis", "Acne")

demo <- fread(file.path(PROC, "demo.csv.gz"))[sex %in% c("M", "F"),
                                              .(primaryid, sex, age_yr, occp_cod)]
ps   <- unique(fread(file.path(PROC, "drug.csv.gz"))[, .(primaryid, role_cod, drug_group)])
rc   <- unique(fread(file.path(PROC, "reac.csv.gz"))[, .(primaryid, pt)])
indi <- unique(fread(file.path(PROC, "indi.csv.gz"))[, .(primaryid, indi_pt)])

demo <- demo[primaryid %in% unique(indi[grepl(DERM, indi_pt, ignore.case = TRUE),
                                        primaryid])]
cat("dermatologic-indication reports with known sex:", nrow(demo), "\n")
ps <- ps[primaryid %in% demo$primaryid & toupper(role_cod) == "PS"]
rc <- rc[primaryid %in% demo$primaryid & !is.na(pt) & pt != "" &
           !grepl(NON_AE, pt, ignore.case = TRUE)]
pop <- intersect(ps$primaryid, rc$primaryid)
demo <- demo[primaryid %in% pop]; ps <- ps[primaryid %in% pop]
rc <- rc[primaryid %in% pop]
cat("analysis population:", nrow(demo), "\n")

arm <- ps[, .(arm = if (any(drug_group == INDEX)) "index" else "comparator"),
          by = primaryid]
arm <- arm[!duplicated(primaryid)]
demo <- arm[demo, on = "primaryid", nomatch = NULL]
demo[, band := fifelse(is.na(age_yr), "unknown",
                fifelse(age_yr < 12, "<12", fifelse(age_yr <= 45, "12-45", "46+")))]

# ---------------------------------------------- which comparator drugs to keep
cmp <- ps[demo[arm == "comparator"], on = "primaryid", nomatch = NULL]
# every (report, drug) pair counts once: keeping one arbitrary drug per
# report made the mix depend on row order, and 17 and 18 then disagreed
cmp <- unique(cmp[!is.na(drug_group), .(primaryid, drug_group)])
cmp <- merge(cmp[, .(primaryid, drug_group)], demo[, .(primaryid, sex)],
             by = "primaryid", all.x = FALSE)
mix <- dcast(cmp[, .(n = .N), by = .(drug_group, sex)], drug_group ~ sex,
             value.var = "n", fill = 0)
mix[, tot := F + M]
mix[, pct_female := round(100 * F / tot, 1)]
keep_drug <- mix[tot >= MIN_DRUG_N & pct_female >= LO_SX & pct_female <= HI_SX,
                 drug_group]
drop_drug <- mix[tot >= MIN_DRUG_N & (pct_female < LO_SX | pct_female > HI_SX),
                 drug_group]
cat("\ncomparator drugs retained (sex-neutral use, >= ", MIN_DRUG_N, " reports): ",
    length(keep_drug), "\n", sep = "")
cat("comparator drugs dropped (used >", HI_SX, "% or <", LO_SX,
    "% by one sex):\n", sep = "")
print(mix[drug_group %in% drop_drug, .(drug_group, tot, pct_female)][order(-tot)],
      nrows = 30)

# a comparator report is usable only if every primary-suspect drug on it is one
# of the sex-neutral drugs; otherwise it still carries the contaminating drug
allps <- ps[, .(n = .N, bad = sum(!(drug_group %in% keep_drug))), by = primaryid]
demo <- merge(demo, allps[, .(primaryid, bad)], by = "primaryid", all.x = TRUE)
demo[, usable := arm == "index" | (is.na(bad) | bad == 0)]

# Size of the pool the restriction acts on. The Brief Report quotes the number
# removed; this is the file it comes from, so the fraction retained is quoted
# from the same object rather than from a count done separately.
pool <- data.table(
  step = c("comparator reports before restriction", "comparator reports removed",
           "comparator reports retained", "isotretinoin reports"),
  n = c(demo[arm == "comparator", .N],
        demo[arm == "comparator" & usable == FALSE, .N],
        demo[arm == "comparator" & usable == TRUE, .N],
        demo[arm == "index", .N])
)
pool[, pct_of_before := round(100 * n / n[1], 1)]
fwrite(pool, file.path(TAB, "comparator_pool_size.csv"))
cat("\n== comparator pool ==\n"); print(pool)

fin <- rc[pt == "Depression", unique(primaryid)]
fin_drug <- ps[drug_group == "finasteride", unique(primaryid)]
dup_drug <- ps[drug_group == "dupilumab", unique(primaryid)]
cat("\n== influence of finasteride on the male comparator depression rate ==\n")
infl <- rbindlist(lapply(c("F", "M"), function(sx) {
  ids <- demo[arm == "comparator" & sex == sx & band == "12-45", primaryid]
  nf <- demo[arm == "comparator" & sex == sx & band == "12-45" &
               primaryid %in% fin_drug, .N]
  nd <- sum(ids %in% fin)
  ndf <- sum(intersect(ids, fin) %in% fin_drug)
  dupn <- demo[arm == "comparator" & sex == sx & band == "12-45" &
                 primaryid %in% dup_drug, .N]
  dupd <- sum(intersect(ids, fin) %in% dup_drug)
  cat(sprintf("  %s: comparator n = %s, finasteride reports = %s, depression = %s, of which on finasteride = %s\n",
              sx, format(length(ids), big.mark = ","), format(nf, big.mark = ","),
              format(nd, big.mark = ","), format(ndf, big.mark = ",")))
  data.table(sex = sx, n_comparator = length(ids), n_finasteride = nf,
             n_depression = nd, n_depression_finasteride = ndf,
             depression_per_1000_finasteride = round(1000 * ndf / nf, 1),
             n_dupilumab = dupn, n_depression_dupilumab = dupd,
             depression_per_1000_dupilumab = round(1000 * dupd / dupn, 1))
}))
fwrite(infl, file.path(TAB, "finasteride_influence.csv"))

# ------------------------------------------------------------------ estimator
est <- function(d, rcc, label) {
  rcc <- rcc[primaryid %in% d$primaryid]
  freq <- rcc[, .(tot = .N), by = pt][tot >= MIN_PT, pt]
  rcc <- rcc[pt %in% freq]
  out <- rbindlist(lapply(KEY, function(ev) {
    pid <- rcc[pt == ev, unique(primaryid)]
    one <- function(sx) {
      ni <- d[arm == "index" & sex == sx, .N]
      ai <- d[arm == "index" & sex == sx & primaryid %in% pid, .N]
      nc <- d[arm == "comparator" & sex == sx, .N]
      ac <- d[arm == "comparator" & sex == sx & primaryid %in% pid, .N]
      if (min(ai, ni - ai, ac, nc - ac) <= 0) return(NULL)
      A <- as.numeric(ai); B <- as.numeric(ni - ai)
      C <- as.numeric(ac); D <- as.numeric(nc - ac)
      list(a = ai, n = ni, c = ac, ror = (A * D) / (B * C), ln = log((A * D) / (B * C)),
           se = sqrt(1/A + 1/B + 1/C + 1/D), oi = A / B, oc = C / D)
    }
    f <- one("F"); m <- one("M")
    if (is.null(f) || is.null(m) || min(f$a, m$a) < MIN_A) return(NULL)
    lo <- f$ln - m$ln; se <- sqrt(f$se^2 + m$se^2)
    data.table(pt = ev, a_F = f$a, a_M = m$a, n_iF = f$n, n_iM = m$n,
               c_F = f$c, c_M = m$c,
               ROR_F = round(f$ror, 2), ROR_M = round(m$ror, 2),
               RORR = round(exp(lo), 2), lci = round(exp(lo - 1.96 * se), 2),
               uci = round(exp(lo + 1.96 * se), 2),
               p = 2 * pnorm(-abs(lo / se)),
               index_OR_F_vs_M = round(f$oi / m$oi, 2),
               comp_OR_M_vs_F  = round(m$oc / f$oc, 2))
  }))
  out[, q := p.adjust(p, method = "fdr")]
  cat("\n== ", label, " ==\n")
  print(out[order(-RORR), .(pt, a_F, a_M, ROR_F, ROR_M, RORR, lci, uci,
                            p = signif(p, 3), q = signif(q, 3),
                            index_OR_F_vs_M, comp_OR_M_vs_F)], nrows = 20)
  out[]
}

d_all   <- demo[usable | arm == "index"]
d_young <- d_all[band == "12-45"]

A_all <- est(d_all,   rc, "PRIMARY comparator (sex-neutral drugs), all ages")
A_y   <- est(d_young, rc, "PRIMARY comparator (sex-neutral drugs), ages 12-45")

# the same thing with the untrimmed pool, for the record
d_raw <- demo
B_all <- est(d_raw, rc, "untrimmed comparator, all ages (retained for comparison)")
B_y   <- est(d_raw[band == "12-45"], rc, "untrimmed comparator, ages 12-45")

fwrite(A_all[, .(pt, a_F, a_M, n_iF, n_iM, c_F, c_M, ROR_F, ROR_M, RORR, lci, uci, p, q,
                 index_OR_F_vs_M, comp_OR_M_vs_F)],
       file.path(TAB, "final_rorr_all_ages.csv"))
fwrite(A_y[, .(pt, a_F, a_M, n_iF, n_iM, c_F, c_M, ROR_F, ROR_M, RORR, lci, uci, p, q,
               index_OR_F_vs_M, comp_OR_M_vs_F)],
       file.path(TAB, "final_rorr_age12_45.csv"))
fwrite(B_all[, .(pt, a_F, a_M, RORR, lci, uci, p)],
       file.path(TAB, "final_rorr_untrimmed_all_ages.csv"))
fwrite(B_y[, .(pt, a_F, a_M, RORR, lci, uci, p)],
       file.path(TAB, "final_rorr_untrimmed_age12_45.csv"))
fwrite(mix[tot >= MIN_DRUG_N,
           .(drug_group, n_female = F, n_male = M, tot, pct_female,
             retained = drug_group %in% keep_drug)][order(-tot)],
       file.path(TAB, "comparator_drug_sex_mix.csv"))
cat("\n[written] final_rorr_all_ages.csv, final_rorr_age12_45.csv,\n")
cat("          final_rorr_untrimmed_*.csv, comparator_drug_sex_mix.csv\n")
cat("[done]\n")
