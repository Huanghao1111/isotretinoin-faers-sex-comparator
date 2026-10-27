#!/usr/bin/env Rscript
# ==============================================================================
# A comparator the index arm can actually be compared with.
#
# The pooled estimate cannot be read. Isotretinoin reports are 63.7% aged 12-45
# while the comparator pool ("all other dermatologic drugs") is 39.9% aged 46+,
# so every stratum gives a depression RORR above 1 (1.26 to 1.56) and the pool
# gives 0.82. That is Simpson's paradox, not a finding.
#
# Fix the design instead of the arithmetic: restrict the whole analysis to
# reports carrying an acne indication. Both arms then come from the same
# indication and, it turns out, a near-identical age distribution. If the RORR
# survives there, it is not an artefact of case mix.
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "C:/Users/lenovo/WorkBuddy/2026-10-17-00-49-45"
PROC <- file.path(ROOT, "data", "processed")
TAB  <- file.path(ROOT, "output", "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

ACNE <- "acne"
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

MIN_A <- 3L; INDEX <- "isotretinoin"

demo <- fread(file.path(PROC, "demo.csv.gz"))[sex %in% c("M", "F"),
                                              .(primaryid, sex, age_yr)]
ps   <- unique(fread(file.path(PROC, "drug.csv.gz"))[, .(primaryid, role_cod, drug_group)])
rc   <- unique(fread(file.path(PROC, "reac.csv.gz"))[, .(primaryid, pt)])
indi <- unique(fread(file.path(PROC, "indi.csv.gz"))[, .(primaryid, indi_pt)])

acne_pid <- unique(indi[grepl(ACNE, indi_pt, ignore.case = TRUE), primaryid])
demo <- demo[primaryid %in% acne_pid]
cat("acne-indication reports with known sex:", nrow(demo), "\n")

ps <- ps[primaryid %in% demo$primaryid & toupper(role_cod) == "PS"]
rc <- rc[primaryid %in% demo$primaryid & !is.na(pt) & pt != "" &
           !grepl(NON_AE, pt, ignore.case = TRUE)]
pop <- intersect(ps$primaryid, rc$primaryid)
demo <- demo[primaryid %in% pop]
ps <- ps[primaryid %in% pop]
rc <- rc[primaryid %in% pop]
arm <- ps[, .(arm = if (any(drug_group == INDEX)) "index" else "comparator"),
          by = primaryid]
arm <- arm[!duplicated(primaryid)]
demo <- arm[demo, on = "primaryid", nomatch = NULL]
cat("analysis population:", nrow(demo), "\n")
print(demo[, .(n = .N), by = .(arm, sex)])

demo[, band := fifelse(is.na(age_yr), "unknown",
                fifelse(age_yr < 12, "<12",
                fifelse(age_yr <= 45, "12-45", "46+")))]
mix <- dcast(demo[, .(n = .N), by = .(band, arm)], band ~ arm,
             value.var = "n", fill = 0)
mix[, pct_index := round(100 * index / sum(index), 1)]
mix[, pct_comp  := round(100 * comparator / sum(comparator), 1)]
cat("\n== age-band composition, acne-indication population ==\n")
print(mix)
cat("\nmedian age, index:", median(demo[arm == "index" & !is.na(age_yr), age_yr]),
    "| comparator:", median(demo[arm == "comparator" & !is.na(age_yr), age_yr]), "\n")

# what is in the comparator arm, now that the indication is fixed
top <- merge(ps[demo[arm == "comparator"], on = "primaryid"][
  !duplicated(primaryid), .(primaryid, drug_group)],
  demo[, .(primaryid, arm)], by = "primaryid")
cat("\n== comparator arm drug mix (top 15) ==\n")
print(top[!is.na(drug_group), .(n = .N), by = drug_group][order(-n)][1:15])

KEY <- c("Depression", "Anxiety", "Suicidal ideation", "Completed suicide",
         "Major depression", "Anhedonia", "Loss of libido",
         "Sexual dysfunction", "Libido decreased", "Lip dry",
         "Crohn's disease", "Inflammatory bowel disease", "Dry skin",
         "Epistaxis", "Acne")

est <- function(d, rcc, pss, label) {
  pss <- pss[primaryid %in% d$primaryid]
  rcc <- rcc[primaryid %in% d$primaryid]
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
      list(a = ai, n = ni, ror = (A * D) / (B * C), ln = log((A * D) / (B * C)),
           se = sqrt(1/A + 1/B + 1/C + 1/D),
           nc = nc, ac = ac, oddsr = A / B, oddsc = C / D)
    }
    f <- one("F"); m <- one("M")
    if (is.null(f) || is.null(m)) return(NULL)
    if (min(f$a, m$a) < MIN_A) return(NULL)
    lo <- f$ln - m$ln; se <- sqrt(f$se^2 + m$se^2)
    data.table(pt = ev, a_F = f$a, a_M = m$a, n_index_F = f$n, n_index_M = m$n,
               ROR_F = round(f$ror, 2), ROR_M = round(m$ror, 2),
               RORR = round(exp(lo), 2), lci = round(exp(lo - 1.96 * se), 2),
               uci = round(exp(lo + 1.96 * se), 2),
               p = 2 * pnorm(-abs(lo / se)),
               index_OR_F_vs_M = round((f$oddsr) / (m$oddsr), 2),
               comparator_OR_M_vs_F = round((m$oddsc) / (f$oddsc), 2))
  }))
  out[, q := p.adjust(p, method = "fdr")]
  cat("\n== ", label, " ==\n")
  print(out[order(-RORR), .(pt, a_F, a_M, ROR_F, ROR_M, RORR, lci, uci,
                            p = signif(p, 3), q = signif(q, 3),
                            index_OR_F_vs_M, comparator_OR_M_vs_F)], nrows = 20)
  out[]
}

Rall <- est(demo, rc, ps, "acne-indication population, all ages")
young <- demo[!is.na(age_yr) & age_yr >= 12 & age_yr <= 45, primaryid]
Ry <- est(demo[primaryid %in% young], rc, ps, "acne-indication population, 12-45")

fwrite(Rall[, .(pt, a_F, a_M, ROR_F, ROR_M, RORR, lci, uci, p, q,
                n_index_F, n_index_M, index_OR_F_vs_M, comparator_OR_M_vs_F)],
       file.path(TAB, "acne_comparator_rorr.csv"))
fwrite(Ry[, .(pt, a_F, a_M, ROR_F, ROR_M, RORR, lci, uci, p, q)],
       file.path(TAB, "acne_comparator_rorr_age12_45.csv"))
fwrite(mix, file.path(TAB, "acne_arm_age_composition.csv"))
cat("\n[written] acne_comparator_rorr.csv, acne_comparator_rorr_age12_45.csv, acne_arm_age_composition.csv\n")
cat("[done]\n")
