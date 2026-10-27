#!/usr/bin/env Rscript
# ==============================================================================
# Robustness of the primary finding.
#
# Trimming sex-skewed drugs out of the comparator pool moves depression from
# 0.97 (null) to 0.59 (male predominance). A jump that large has to survive
# alternative trimming thresholds, dropping finasteride alone rather than all
# five flagged drugs, and the age-stratified analysis, since the untrimmed pool
# was also Simpson-contaminated.
#
# Counts are tabulated in a single pass per rule. The previous implementation
# re-grouped the reaction table once per rule-by-event cell and never finished.
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

MIN_A <- 3L; INDEX <- "isotretinoin"; MIN_DRUG_N <- 200L
EVENTS <- c("Depression", "Anxiety", "Suicidal ideation", "Completed suicide",
            "Major depression", "Loss of libido")

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
cat("analysis population:", nrow(demo), "\n")

arm <- ps[, .(arm = if (any(drug_group == INDEX)) "index" else "comparator"),
          by = primaryid]
arm <- arm[!duplicated(primaryid)]
demo <- arm[demo, on = "primaryid", nomatch = NULL]
demo[, band := fifelse(is.na(age_yr), "unknown",
                fifelse(age_yr < 12, "<12", fifelse(age_yr <= 45, "12-45", "46+")))]
cat("index arm:", demo[arm == "index", .N], "\n")

# comparator drug sex mix
cmp <- ps[demo[arm == "comparator"], on = "primaryid", nomatch = NULL]
# every (report, drug) pair counts once: keeping one arbitrary drug per
# report made the mix depend on row order, and 17 and 18 then disagreed
cmp <- unique(cmp[!is.na(drug_group), .(primaryid, drug_group)])
cmp <- merge(cmp[, .(primaryid, drug_group)], demo[, .(primaryid, sex)],
             by = "primaryid", all.x = FALSE)
mix <- dcast(cmp[, .(n = .N), by = .(drug_group, sex)], drug_group ~ sex,
             value.var = "n", fill = 0)
mix[, tot := F + M]
mix[, pct_female := 100 * F / tot]

# one long event table, reused by every rule
ev <- rc[pt %in% EVENTS, .(primaryid, pt)]
ev <- ev[primaryid %in% demo$primaryid]

bad_of <- function(dropped) {
  # "kept" means kept on purpose: a drug with fewer than MIN_DRUG_N reports is
  # not retained, it is simply too small to judge. Leaving those in `ok` kept
  # their reports in the pool and made every rule remove 467 fewer reports than
  # the primary analysis.
  ok <- mix[tot >= MIN_DRUG_N & !(drug_group %in% dropped), drug_group]
  ps[, .(bad = sum(!(drug_group %in% ok))), by = primaryid]
}

sel <- function(dropped) {
  if (length(dropped) == 0) return(demo[, .(primaryid, arm, sex, band)])
  b <- bad_of(dropped)
  d <- merge(demo[, .(primaryid, arm, sex, band)], b, by = "primaryid", all.x = TRUE)
  d[arm == "index" | (is.na(bad) | bad == 0)][, .(primaryid, arm, sex, band)]
}

tabulate <- function(d) {
  base <- d[, .(n = .N), by = .(arm, sex)]
  num  <- merge(ev, d[, .(primaryid, arm, sex)], by = "primaryid",
                all.x = FALSE)[, .(a = .N), by = .(pt, arm, sex)]
  setkey(base, arm, sex); setkey(num, pt, arm, sex)
  cell <- function(ptx, ar, sx) {
    n <- base[.(ar, sx), n]
    a <- num[.(ptx, ar, sx), a]
    if (is.na(a)) a <- 0L
    list(n = n, a = a)
  }
  rbindlist(lapply(EVENTS, function(x) {
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
    data.table(pt = x, a_F = fi$a, a_M = mi$a,
               n_iF = fi$n, n_iM = mi$n, n_cF = fc$n, n_cM = mc$n,
               c_F = fc$a, c_M = mc$a,
               ROR_F = round(rf, 2), ROR_M = round(rm, 2),
               RORR = round(exp(lo), 2), lci = round(exp(lo - 1.96 * se), 2),
               uci = round(exp(lo + 1.96 * se), 2),
               index_OR_F_vs_M = round(oiF / oiM, 2),
               comp_OR_M_vs_F  = round(ocM / ocF, 2))
  }))
}

rules <- list(
  "no trimming"             = character(0),
  "finasteride only"        = "finasteride",
  "20-80% female (primary)" = "20-80",
  "30-70% female"           = "30-70",
  "40-60% female"           = "40-60"
)

rows <- rbindlist(lapply(names(rules), function(nm) {
  r <- rules[[nm]]
  dropped <- if (length(r) == 0) character(0) else
    if (grepl("-", r, fixed = TRUE)) {
      b <- as.numeric(strsplit(r, "-", fixed = TRUE)[[1]])
      mix[tot >= MIN_DRUG_N & (pct_female < b[1] | pct_female > b[2]), drug_group]
    } else r
  d <- sel(dropped)
  z <- tabulate(d)
  cat(sprintf("[%s] n_index=%s n_comparator=%s dropped=%d\n", nm,
              format(d[arm == "index", .N], big.mark = ","),
              format(d[arm == "comparator", .N], big.mark = ","),
              length(dropped)))
  z[, rule := nm]
  z[]
}))

for (x in EVENTS) {
  cat("\n== ", x, " ==\n")
  print(rows[pt == x, .(rule, n_iF, n_iM, n_cF, n_cM, a_F, a_M, RORR, lci, uci,
                        index_OR_F_vs_M, comp_OR_M_vs_F)], nrows = 10)
}

d <- sel(mix[tot >= MIN_DRUG_N & (pct_female < 20 | pct_female > 80), drug_group])
cat("\n== age-band composition after trimming ==\n")
m2 <- dcast(d[, .(n = .N), by = .(band, arm)], band ~ arm, value.var = "n", fill = 0)
m2[, pct_index := round(100 * index / sum(index), 1)]
m2[, pct_comp := round(100 * comparator / sum(comparator), 1)]
print(m2)

cat("\n== by age band, trimmed comparator ==\n")
for (bd in c("12-45", "46+", "unknown")) {
  z <- tabulate(d[band == bd])
  if (!is.null(z) && nrow(z)) {
    cat("-- ", bd, " --\n")
    print(z[, .(pt, a_F, a_M, RORR, lci, uci, index_OR_F_vs_M, comp_OR_M_vs_F)],
          nrows = 10)
  }
}

fwrite(rows, file.path(TAB, "robustness_trimming_rules.csv"))
fwrite(m2, file.path(TAB, "trimmed_arm_age_composition.csv"))
cat("\n[written] robustness_trimming_rules.csv, trimmed_arm_age_composition.csv\n")
cat("[done]\n")
