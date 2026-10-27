#!/usr/bin/env Rscript
# ==============================================================================
# Single source of truth for every number quoted in the manuscript body.
#
# Primary estimand (fixed in response to review): ONE analysis population --
# reports carrying >=1 primary-suspect dermatologic drug AND >=1 retained
# adverse-event term -- partitioned by whether the primary suspect is
# isotretinoin (index) or another dermatologic drug (comparator). The earlier
# implementation drew cells from three differently restricted denominators; it
# is retained as a sensitivity analysis and differs in the third decimal.
#
# Also produces the formal age-band interaction that replaced the post-hoc
# "age is an effect modifier" assertion: 12-45 years versus 46+ years,
# compared as a ratio of ratios with a z test.
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

KEY <- c("Depression", "Anxiety", "Suicidal ideation", "Completed suicide",
         "Major depression", "Anhedonia", "Loss of libido",
         "Sexual dysfunction", "Libido decreased", "Lip dry",
         "Crohn's disease", "Inflammatory bowel disease", "Dry skin",
         "Epistaxis", "Acne")

dv0 <- fread(file.path(PROC, "demo.csv.gz"))[sex %in% c("M", "F"),
                                             .(primaryid, sex, age_yr)]
ps0 <- unique(fread(file.path(PROC, "drug.csv.gz"))[
  , .(primaryid, role_cod, drug_group)])
rc0 <- unique(fread(file.path(PROC, "reac.csv.gz"))[, .(primaryid, pt)])
in0 <- unique(fread(file.path(PROC, "indi.csv.gz"))[, .(primaryid, indi_pt)])

dv0 <- dv0[primaryid %in% unique(in0[grepl(DERM, indi_pt, ignore.case = TRUE),
                                     primaryid])]
ps0 <- ps0[primaryid %in% dv0$primaryid & toupper(role_cod) == "PS"]
rc0 <- rc0[primaryid %in% dv0$primaryid & !is.na(pt) & pt != "" &
             !grepl(NON_AE, pt, ignore.case = TRUE)]

# estimate the ROR ratio inside an arbitrary subset of report ids
rors <- function(ids) {
  dvv <- dv0[primaryid %in% ids]
  pss <- ps0[primaryid %in% ids]
  rcc <- rc0[primaryid %in% ids]
  pop <- intersect(pss$primaryid, rcc$primaryid)
  dvv <- dvv[primaryid %in% pop]
  pss <- pss[primaryid %in% pop]
  rcc <- rcc[primaryid %in% pop]
  arm <- pss[, .(arm = if (any(drug_group == INDEX)) "index" else "comparator"),
             by = primaryid]
  arm <- arm[!duplicated(primaryid)]
  dvv <- arm[dvv, on = "primaryid", nomatch = NULL]
  freq <- rcc[, .(tot = .N), by = pt][tot >= MIN_PT, pt]
  rcc <- rcc[pt %in% freq]
  out <- rbindlist(lapply(KEY, function(ev) {
    pid <- rcc[pt == ev, unique(primaryid)]
    one <- function(sx) {
      ni <- dvv[arm == "index" & sex == sx, .N]
      ai <- dvv[arm == "index" & sex == sx & primaryid %in% pid, .N]
      nc <- dvv[arm == "comparator" & sex == sx, .N]
      ac <- dvv[arm == "comparator" & sex == sx & primaryid %in% pid, .N]
      if (min(ai, ni - ai, ac, nc - ac) <= 0) return(NULL)
      A <- as.numeric(ai); B <- as.numeric(ni - ai)
      C <- as.numeric(ac); D <- as.numeric(nc - ac)
      list(a = ai, ror = (A * D) / (B * C), ln = log((A * D) / (B * C)),
           se = sqrt(1/A + 1/B + 1/C + 1/D))
    }
    f <- one("F"); m <- one("M")
    if (is.null(f) || is.null(m)) return(NULL)
    lo <- f$ln - m$ln; se <- sqrt(f$se^2 + m$se^2)
    data.table(pt = ev, a_F = f$a, a_M = m$a,
               ROR_F = round(f$ror, 2), ROR_M = round(m$ror, 2),
               RORR = round(exp(lo), 2), lci = round(exp(lo - 1.96 * se), 2),
               uci = round(exp(lo + 1.96 * se), 2),
               p = 2 * pnorm(-abs(lo / se)), ln = lo, se = se)
  }))
  if (nrow(out)) out[, q := p.adjust(p, method = "fdr")]
  out[]
}

all_ids <- dv0$primaryid
young   <- dv0[!is.na(age_yr) & age_yr >= 12 & age_yr <= 45, primaryid]
older   <- dv0[!is.na(age_yr) & age_yr > 45, primaryid]

Rall <- rors(all_ids); Ry <- rors(young); Ro <- rors(older)

# The all-ages estimate (0.97 for depression) sits below both the 12-45 and the
# 46+ estimate, so something in the remainder pulls it down. That remainder is
# made of reports with no recorded age. It is checked here rather than assumed.
child <- dv0[!is.na(age_yr) & age_yr < 12, primaryid]
unk   <- dv0[is.na(age_yr), primaryid]
cat(sprintf("\nreport counts by age band: <12 = %d, 12-45 = %d, 46+ = %d, unknown = %d\n",
            length(child), length(young), length(older), length(unk)))
Rc <- try(rors(child), silent = TRUE)
Ru <- try(rors(unk), silent = TRUE)
show_band <- function(x, ttl) {
  cat("\n==", ttl, "==\n")
  if (inherits(x, "try-error") || is.null(x) || !nrow(x)) {
    cat("(no event met the minimum-count rule)\n"); return(invisible(NULL))
  }
  print(x[order(-RORR), .(pt, a_F, a_M, RORR, lci, uci, p = signif(p, 3))],
        nrows = 20)
}
show_band(Rc, "ages <12")
show_band(Ru, "age not recorded")
cat("\n== PRIMARY: corrected ROR ratio, all ages ==\n")
print(Rall[order(-RORR), .(pt, a_F, a_M, ROR_F, ROR_M, RORR, lci, uci,
                           p = signif(p, 3), q = signif(q, 3))], nrows = 20)
cat("\n== ages 12-45 ==\n")
print(Ry[order(-RORR), .(pt, a_F, a_M, RORR, lci, uci, p = signif(p, 3))],
      nrows = 20)
cat("\n== ages 46+ ==\n")
print(Ro[order(-RORR), .(pt, a_F, a_M, RORR, lci, uci, p = signif(p, 3))],
      nrows = 20)

# ------------------------------------------------- formal age-band interaction
it <- merge(Ry[, .(pt, a_F_y = a_F, a_M_y = a_M, rorr_y = RORR,
                   lci_y = lci, uci_y = uci, ln_y = ln, se_y = se)],
            Ro[, .(pt, a_F_o = a_F, a_M_o = a_M, rorr_o = RORR,
                   lci_o = lci, uci_o = uci, ln_o = ln, se_o = se)],
            by = "pt")
it[, `:=`(dln = ln_y - ln_o, se_d = sqrt(se_y^2 + se_o^2))]
it[, `:=`(ratio = round(exp(dln), 2), r_lci = round(exp(dln - 1.96 * se_d), 2),
          r_uci = round(exp(dln + 1.96 * se_d), 2))]
it[, p_int := round(2 * pnorm(-abs(dln / se_d)), 4)]
it[, q_int := round(p.adjust(p_int, method = "fdr"), 4)]
cat("\n== age-band interaction: ratio of ratios, 12-45 versus 46+ ==\n")
print(it[order(p_int), .(pt, rorr_12_45 = rorr_y, rorr_46plus = rorr_o,
                         ratio, r_lci, r_uci, p_int, q_int)], nrows = 20)

fwrite(Rall[, .(pt, a_F, a_M, ROR_F, ROR_M, RORR, lci, uci, p, q)],
       file.path(TAB, "primary_rorr_all_ages.csv"))
fwrite(Ry[, .(pt, a_F, a_M, RORR, lci, uci, p)],
       file.path(TAB, "primary_rorr_age12_45.csv"))
fwrite(Ro[, .(pt, a_F, a_M, RORR, lci, uci, p)],
       file.path(TAB, "primary_rorr_age46plus.csv"))
fwrite(it[, .(pt, a_F_y, a_M_y, rorr_12_45 = rorr_y, lci_y, uci_y,
              a_F_o, a_M_o, rorr_46plus = rorr_o, lci_o, uci_o,
              ratio_12_45_vs_46plus = ratio, r_lci, r_uci, p_int, q_int)],
       file.path(TAB, "age_band_interaction.csv"))
cat("\n[written] primary_rorr_all_ages.csv, primary_rorr_age12_45.csv,\n")
cat("          primary_rorr_age46plus.csv, age_band_interaction.csv\n")
cat("[done]\n")
