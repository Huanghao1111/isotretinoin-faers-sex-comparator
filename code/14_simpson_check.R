#!/usr/bin/env Rscript
# ==============================================================================
# Why does the pooled estimate disagree with every stratum?
#
# All-ages depression RORR = 0.97, yet:
#   12-45  = 1.28,  46+ = 1.51,  age not recorded = 1.26.
# Every stratum points the same way and the pool says null. That is Simpson's
# paradox, and it happens here because isotretinoin reports and comparator
# reports sit in completely different age bands. If so, the headline null is an
# artefact of pooling and must not be reported as the primary result.
#
# This script (a) confirms it, (b) gives a stratum-standardised estimate
# (Mantel-Haenszel odds ratio of the 2x2, computed within each age band and
# pooled), and (c) reports how far apart the two arms are by age band.
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

MIN_A <- 3L; MIN_PT <- 50L; INDEX <- "isotretinoin"

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
cat("analysis population:", nrow(demo), "\n")

# (c) how far apart are the two arms by age band?
demo[, band := fifelse(is.na(age_yr), "unknown",
                fifelse(age_yr < 12, "<12",
                fifelse(age_yr <= 45, "12-45", "46+")))]
mix <- demo[, .(n = .N), by = .(band, arm)]
mix <- dcast(mix, band ~ arm, value.var = "n", fill = 0)
mix[, pct_index := round(100 * index / sum(mix$index), 1)]
mix[, pct_comp  := round(100 * comparator / sum(mix$comparator), 1)]
cat("\n== age-band composition of the two arms ==\n")
print(mix)

# (a) confirm: pooled over known-age reports only
cells <- function(d, ev) {
  pid <- rc[pt == ev, unique(primaryid)]
  rbindlist(lapply(c("F", "M"), function(sx) {
    ni <- d[arm == "index" & sex == sx, .N]
    ai <- d[arm == "index" & sex == sx & primaryid %in% pid, .N]
    nc <- d[arm == "comparator" & sex == sx, .N]
    ac <- d[arm == "comparator" & sex == sx & primaryid %in% pid, .N]
    data.table(sex = sx, a = ai, b = ni - ai, c = ac, d = nc - ac)
  }))
}
rorr <- function(cc) {
  if (min(cc$a, cc$b, cc$c, cc$d) <= 0) return(NULL)
  f <- cc[sex == "F"]; m <- cc[sex == "M"]
  rf <- (f$a * f$d) / (f$b * f$c); rm <- (m$a * m$d) / (m$b * m$c)
  lo <- log(rf) - log(rm)
  se <- sqrt(1/f$a + 1/f$b + 1/f$c + 1/f$d + 1/m$a + 1/m$b + 1/m$c + 1/m$d)
  data.table(RORR = round(exp(lo), 2), lci = round(exp(lo - 1.96 * se), 2),
             uci = round(exp(lo + 1.96 * se), 2),
             p = signif(2 * pnorm(-abs(lo / se)), 3))
}

KEY <- c("Depression", "Anxiety", "Suicidal ideation", "Completed suicide",
         "Major depression", "Anhedonia", "Loss of libido",
         "Sexual dysfunction", "Libido decreased", "Lip dry",
         "Crohn's disease", "Inflammatory bowel disease", "Dry skin",
         "Epistaxis", "Acne")

known <- demo[!is.na(age_yr), primaryid]
cat("\n== pooled over reports with a recorded age ==\n")
out <- rbindlist(lapply(KEY, function(ev) {
  cc <- cells(demo[primaryid %in% known], ev)
  if (is.null(cc) || min(cc$a) < MIN_A) return(NULL)
  rr <- rorr(cc); if (is.null(rr)) return(NULL)
  cbind(data.table(pt = ev, a_F = cc[sex == "F", a], a_M = cc[sex == "M", a]), rr)
}))
print(out[order(-RORR)], nrows = 20)

# (b) Mantel-Haenszel across the three recorded age bands, per sex, then ratio
cat("\n== age-standardised (Mantel-Haenszel within band, 12-45 / 46+ / <12) ==\n")
mh <- rbindlist(lapply(KEY, function(ev) {
  pid <- rc[pt == ev, unique(primaryid)]
  per <- rbindlist(lapply(c("12-45", "46+", "<12"), function(bd) {
    d <- demo[band == bd]
    rbindlist(lapply(c("F", "M"), function(sx) {
      ni <- d[arm == "index" & sex == sx, .N]
      ai <- d[arm == "index" & sex == sx & primaryid %in% pid, .N]
      nc <- d[arm == "comparator" & sex == sx, .N]
      ac <- d[arm == "comparator" & sex == sx & primaryid %in% pid, .N]
      if (min(ai, ni - ai, ac, nc - ac) <= 0) return(NULL)
      data.table(band = bd, sex = sx, a = ai, b = ni - ai, c = ac, d = nc - ac)
    }))
  }))
  if (is.null(per) || nrow(per) < 4) return(NULL)
  mh_one <- function(sx) {
    x <- per[sex == sx]
    num <- sum(x$a * x$d / (x$a + x$b + x$c + x$d))
    den <- sum(x$b * x$c / (x$a + x$b + x$c + x$d))
    if (den <= 0 || num <= 0) return(NULL)
    # Robins-Breslow-Greenland variance for the MH odds ratio
    ss <- per[sex == sx]
    P <- (ss$a + ss$d) / (ss$a + ss$b + ss$c + ss$d)
    Q <- (ss$b + ss$c) / (ss$a + ss$b + ss$c + ss$d)
    R <- (ss$a * ss$d) / (ss$a + ss$b + ss$c + ss$d)
    S <- (ss$b * ss$c) / (ss$a + ss$b + ss$c + ss$d)
    v <- sum(P * R) / (2 * sum(R)^2) +
         sum(P * S + Q * R) / (2 * sum(R) * sum(S)) +
         sum(Q * S) / (2 * sum(S)^2)
    list(or = num / den, var = v)
  }
  f <- mh_one("F"); m <- mh_one("M")
  if (is.null(f) || is.null(m)) return(NULL)
  lo <- log(f$or) - log(m$or); se <- sqrt(f$var / f$or^2 + m$var / m$or^2)
  data.table(pt = ev, a_F = per[sex == "F", sum(a)], a_M = per[sex == "M", sum(a)],
             MH_RORR = round(exp(lo), 2), lci = round(exp(lo - 1.96 * se), 2),
             uci = round(exp(lo + 1.96 * se), 2),
             p = signif(2 * pnorm(-abs(lo / se)), 3))
}))
print(mh[order(-MH_RORR)], nrows = 20)

fwrite(out, file.path(TAB, "pooled_known_age_rorr.csv"))
fwrite(mh,  file.path(TAB, "age_standardised_rorr.csv"))
fwrite(mix, file.path(TAB, "arm_age_composition.csv"))
cat("\n[written] pooled_known_age_rorr.csv, age_standardised_rorr.csv, arm_age_composition.csv\n")
cat("[done]\n")
