#!/usr/bin/env Rscript
# ==============================================================================
# The four pre-existing restriction scenarios, re-run on the corrected design.
#
# The deposited sensitivity table (TableS4) was built twice wrong: it used the
# old three-denominator 2x2, and it used the comparator pool before sex-skewed
# drugs were removed. A reviewer will ask whether the male predominance survives
# the restrictions that the previous draft applied, so the same restrictions are
# applied here to the corrected population and the sex-neutral comparator. No new
# scenario is invented.
#
# Also produces:
#   - the formal age-band interaction under the primary comparator (the one in
#     TableS27 was computed on the uncorrected comparator and is now stale);
#   - the internal control: whether the correction moved mucocutaneous and
#     gastrointestinal terms too, or only the psychiatric ones. If it moved
#     everything, "correction" would be a global shift and the paper would have
#     no specific claim.
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

MIN_A <- 3L; INDEX <- "isotretinoin"; MIN_DRUG_N <- 200L
KEY <- c("Depression", "Anxiety", "Suicidal ideation", "Completed suicide",
         "Major depression", "Anhedonia",
         "Loss of libido", "Sexual dysfunction", "Libido decreased",
         "Lip dry", "Crohn's disease", "Inflammatory bowel disease",
         "Dry skin", "Epistaxis", "Acne")
NEURO <- c("Depression", "Anxiety", "Suicidal ideation", "Completed suicide",
           "Major depression", "Anhedonia")
SOMATIC <- c("Lip dry", "Crohn's disease", "Inflammatory bowel disease",
             "Dry skin", "Epistaxis", "Acne")

demo <- fread(file.path(PROC, "demo.csv.gz"))[sex %in% c("M", "F"),
                                              .(primaryid, sex, age_yr, occp_cod)]
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
demo[, band := fifelse(is.na(age_yr), "unknown",
                fifelse(age_yr < 12, "<12", fifelse(age_yr <= 45, "12-45", "46+")))]
ndrug <- ps[, .(ndrug = uniqueN(drug_group)), by = primaryid]
demo <- merge(demo, ndrug, by = "primaryid", all.x = TRUE)

# comparator trimming, as in 17_final_analysis.R
cmp <- ps[demo[arm == "comparator"], on = "primaryid", nomatch = NULL]
cmp <- cmp[!duplicated(primaryid)][!is.na(drug_group)]
cmp <- merge(cmp[, .(primaryid, drug_group)], demo[, .(primaryid, sex)],
             by = "primaryid", all.x = FALSE)
mix <- dcast(cmp[, .(n = .N), by = .(drug_group, sex)], drug_group ~ sex,
             value.var = "n", fill = 0)
mix[, `:=`(tot = F + M, pct_female = 100 * F / (F + M))]
keep <- mix[tot >= MIN_DRUG_N & pct_female >= 20 & pct_female <= 80, drug_group]
bad <- ps[, .(nbad = sum(!(drug_group %in% keep))), by = primaryid]
demo <- merge(demo, bad, by = "primaryid", all.x = TRUE)
demo[, usable := arm == "index" | (is.na(nbad) | nbad == 0)]

ev <- rc[pt %in% KEY, .(primaryid, pt)]

# one pass per scenario: denominators once, event numerators in one grouping
tabulate <- function(ids, label, trimmed = TRUE) {
  d <- if (trimmed) demo[primaryid %in% ids & usable == TRUE] else
                    demo[primaryid %in% ids]
  if (d[arm == "index", .N] < 50 || d[arm == "comparator", .N] < 50) return(NULL)
  base <- d[, .(n = .N), by = .(arm, sex)]
  num  <- merge(ev, d[, .(primaryid, arm, sex)], by = "primaryid",
                all.x = FALSE)[, .(a = .N), by = .(pt, arm, sex)]
  setkey(base, arm, sex); setkey(num, pt, arm, sex)
  # the parameter must not be called `a`: `num` has its own `a` column, and
  # data.table resolves i-expressions in the frame of x first
  cell <- function(p, armx, s) {
    n <- base[.(armx, s), n]; z <- num[.(p, armx, s), a]
    list(n = n, a = if (is.na(z)) 0L else z)
  }
  out <- rbindlist(lapply(KEY, function(x) {
    fi <- cell(x, "index", "F"); mi <- cell(x, "index", "M")
    fc <- cell(x, "comparator", "F"); mc <- cell(x, "comparator", "M")
    if (min(fi$a, mi$a, fc$a, mc$a) < MIN_A) return(NULL)
    oiF <- fi$a / (fi$n - fi$a); oiM <- mi$a / (mi$n - mi$a)
    ocF <- fc$a / (fc$n - fc$a); ocM <- mc$a / (mc$n - mc$a)
    if (min(oiF, oiM, ocF, ocM) <= 0) return(NULL)
    A <- as.numeric(fi$a); B <- as.numeric(fi$n - fi$a)
    C <- as.numeric(fc$a); D <- as.numeric(fc$n - fc$a)
    rf <- (A * D) / (B * C); sf <- sqrt(1/A + 1/B + 1/C + 1/D)
    A2 <- as.numeric(mi$a); B2 <- as.numeric(mi$n - mi$a)
    C2 <- as.numeric(mc$a); D2 <- as.numeric(mc$n - mc$a)
    rm <- (A2 * D2) / (B2 * C2); sm <- sqrt(1/A2 + 1/B2 + 1/C2 + 1/D2)
    lo <- log(rf) - log(rm); se <- sqrt(sf^2 + sm^2)
    data.table(pt = x, a_F = fi$a, a_M = mi$a, n_iF = fi$n, n_iM = mi$n,
               RORR = round(exp(lo), 2), lci = round(exp(lo - 1.96 * se), 2),
               uci = round(exp(lo + 1.96 * se), 2),
               p = 2 * pnorm(-abs(lo / se)),
               index_OR_F_vs_M = round(oiF / oiM, 2),
               comp_OR_M_vs_F = round(ocM / ocF, 2))
  }))
  cat(sprintf("[%-18s] index=%-7s comparator=%-9s\n", label,
              format(d[arm == "index", .N], big.mark = ","),
              format(d[arm == "comparator", .N], big.mark = ",")))
  out[, `:=`(scenario = label, n_reports = d[, .N])]
  out[]
}

single <- demo[!is.na(ndrug) & ndrug == 1, primaryid]
hp     <- demo[toupper(occp_cod) %in% c("MD", "PH", "HP", "RN"), primaryid]
young  <- demo[!is.na(age_yr) & age_yr >= 12 & age_yr <= 45, primaryid]
nolaw  <- demo[toupper(occp_cod) != "LW", primaryid]

cat("\n== restriction scenarios on the corrected design ==\n")
res <- rbindlist(list(
  tabulate(demo$primaryid,                  "A_main"),
  tabulate(single,                          "B_single_drug"),
  tabulate(hp,                              "C_healthprof"),
  tabulate(young,                           "D_age12_45"),
  tabulate(intersect(single, hp),           "E_single+hp"),
  tabulate(nolaw,                           "F_no_lawyer"),
  tabulate(intersect(nolaw, single),        "G_no_lawyer_single")
))

disp <- dcast(res[, .(scenario, pt,
                      v = sprintf("%.2f [%.2f-%.2f]", RORR, lci, uci))],
              pt ~ scenario, value.var = "v")
disp[, pt := factor(pt, levels = KEY)]
setorder(disp, pt)
cat("\n== ratios of ratios by scenario (corrected design, sex-neutral comparator) ==\n")
print(disp, nrows = 20)

cat("\n== direction stability by event ==\n")
for (x in KEY) {
  s <- res[pt == x]
  if (!nrow(s)) next
  rr <- s$RORR
  cat(sprintf("%-28s main=%6.2f  range %6.2f to %6.2f  %s\n", x,
              s[scenario == "A_main", RORR], min(rr), max(rr),
              ifelse(all(rr > 1) | all(rr < 1), "STABLE", "VARIES")))
}

# ------------------------------------- age-band interaction, primary comparator
cat("\n== age-band interaction under the primary comparator ==\n")
bands <- rbindlist(lapply(c("12-45", "46+", "unknown"), function(bd) {
  z <- tabulate(demo[band == bd, primaryid], paste0("band:", bd))
  if (is.null(z)) return(NULL)
  z[, band := bd]
  z[]
}))
if (!is.null(bands) && nrow(bands)) {
  y <- bands[band == "12-45", .(pt, rorr_y = RORR, ln_y = log(RORR),
                                se_y = (log(uci) - log(lci)) / (2 * 1.96))]
  o <- bands[band == "46+", .(pt, rorr_o = RORR, ln_o = log(RORR),
                              se_o = (log(uci) - log(lci)) / (2 * 1.96))]
  it <- merge(y, o, by = "pt")
  it[, `:=`(dln = ln_y - ln_o, se_d = sqrt(se_y^2 + se_o^2))]
  it[, `:=`(ratio = round(exp(dln), 2), r_lci = round(exp(dln - 1.96 * se_d), 2),
            r_uci = round(exp(dln + 1.96 * se_d), 2))]
  it[, p_int := round(2 * pnorm(-abs(dln / se_d)), 4)]
  it[, q_int := round(p.adjust(p_int, method = "fdr"), 4)]
  print(it[, .(pt, rorr_12_45 = rorr_y, rorr_46plus = rorr_o,
               ratio_12_45_vs_46plus = ratio, r_lci, r_uci, p_int, q_int)],
        nrows = 20)
  fwrite(it, file.path(TAB, "primary_age_band_interaction.csv"))
  fwrite(bands, file.path(TAB, "primary_rorr_by_age_band.csv"))
}

# --------------------------------------------------------------- internal control
cat("\n== internal control: did the correction move the non-psychiatric terms? ==\n")
allc <- tabulate(demo$primaryid, "all_drug_comparator", trimmed = FALSE)
pric <- tabulate(demo$primaryid, "sex_neutral_comparator", trimmed = TRUE)
ctl <- merge(allc[, .(pt, rorr_all = RORR, lci_all = lci, uci_all = uci)],
             pric[, .(pt, rorr_pri = RORR, lci_pri = lci, uci_pri = uci)],
             by = "pt")
ctl[, panel := fifelse(pt %in% NEURO, "neuropsychiatric",
                fifelse(pt %in% SOMATIC, "mucocutaneous / GI", "sexual"))]
ctl[, abs_log_shift := round(abs(log(rorr_pri) - log(rorr_all)), 2)]
print(ctl[order(panel, -abs_log_shift),
          .(pt, panel, rorr_all, rorr_pri, abs_log_shift)], nrows = 20)
cat("\nmedian |log shift| by panel:\n")
print(ctl[, .(median_shift = median(abs_log_shift), n = .N), by = panel])

fwrite(res, file.path(TAB, "primary_sensitivity_scenarios.csv"))
fwrite(ctl, file.path(TAB, "comparator_shift_by_panel.csv"))
cat("\n[written] primary_sensitivity_scenarios.csv, primary_age_band_interaction.csv,\n")
cat("          primary_rorr_by_age_band.csv, comparator_shift_by_panel.csv\n")
cat("[done]\n")
