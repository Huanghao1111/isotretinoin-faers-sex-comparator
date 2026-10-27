#!/usr/bin/env Rscript
# ==============================================================================
# Leave-one-drug-out across the comparator pool.
#
# The paper's claim is that one comparator drug, finasteride, manufactured the
# female-predominant signal. The symmetric question is whether one comparator
# drug manufactures the corrected male-predominant signal: if the estimate
# depends on a single drug in the pool, the correction has simply moved the
# problem rather than removed it.
#
# Each retained comparator drug is removed in turn and the primary 2x2 is
# recomputed. Removing finasteride alone is included as the reference.
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
MIN_DRUG_N <- 200L; LO_SX <- 20; HI_SX <- 80
TERMS <- c("Depression", "Anxiety", "Suicidal ideation")

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

cmp <- ps[demo[arm == "comparator"], on = "primaryid", nomatch = NULL]
cmp <- cmp[!duplicated(primaryid)][!is.na(drug_group)]
cmp <- merge(cmp[, .(primaryid, drug_group)], demo[, .(primaryid, sex)],
             by = "primaryid", all.x = FALSE)
mix <- dcast(cmp[, .(n = .N), by = .(drug_group, sex)], drug_group ~ sex,
             value.var = "n", fill = 0)
mix[, tot := F + M]
mix[, pct_female := round(100 * F / tot, 1)]
keep_drug <- mix[tot >= MIN_DRUG_N & pct_female >= LO_SX & pct_female <= HI_SX,
                 drug_group]
allps <- ps[, .(bad = sum(!(drug_group %in% keep_drug))), by = primaryid]
demo <- merge(demo, allps, by = "primaryid", all.x = TRUE)
demo[, usable := arm == "index" | (is.na(bad) | bad == 0)]
d <- demo[usable == TRUE]
cat("comparator drugs retained:", length(keep_drug), "\n")

# cells for the full pool
cell <- function(ev) {
  pid <- rc[pt == ev, unique(primaryid)]
  one <- function(sx, dd) {
    ni <- dd[arm == "index" & sex == sx, .N]
    ai <- dd[arm == "index" & sex == sx & primaryid %in% pid, .N]
    nc <- dd[arm == "comparator" & sex == sx, .N]
    ac <- dd[arm == "comparator" & sex == sx & primaryid %in% pid, .N]
    list(ni = ni, ai = ai, nc = nc, ac = ac)
  }
  list(F = one("F", d), M = one("M", d), pid = pid)
}

rorr_of <- function(f, m) {
  if (min(f$ai, f$ni - f$ai, f$ac, f$nc - f$ac,
          m$ai, m$ni - m$ai, m$ac, m$nc - m$ac) <= 0) return(NULL)
  ror_f <- (f$ai / (f$ni - f$ai)) / (f$ac / (f$nc - f$ac))
  ror_m <- (m$ai / (m$ni - m$ai)) / (m$ac / (m$nc - m$ac))
  rorr <- ror_f / ror_m
  se <- sqrt(sum(1 / c(f$ai, f$ni - f$ai, f$ac, f$nc - f$ac,
                       m$ai, m$ni - m$ai, m$ac, m$nc - m$ac)))
  list(rorr = rorr, lci = exp(log(rorr) - 1.96 * se), uci = exp(log(rorr) + 1.96 * se))
}

# report counts and event counts per comparator drug, computed once
cmpset <- unique(ps[demo[arm == "comparator" & usable == TRUE],
                    on = "primaryid", nomatch = NULL,
                    .(primaryid, drug_group)])
cmpset <- merge(cmpset, d[, .(primaryid, sex)], by = "primaryid")
tot <- dcast(cmpset[, .N, by = .(drug_group, sex)], drug_group ~ sex,
             value.var = "N", fill = 0)
setnames(tot, c("F", "M"), c("n_F", "n_M"))

res <- rbindlist(lapply(TERMS, function(ev) {
  c0 <- cell(ev)
  base <- rorr_of(c0$F, c0$M)
  out <- data.table(pt = ev, removed = "none",
                    RORR = base$rorr, lci = base$lci, uci = base$uci)
  # cases attributable to each drug, among reports carrying the event
  ev_cnt <- cmpset[primaryid %in% c0$pid, .N, by = .(drug_group, sex)]
  ev_cnt <- dcast(ev_cnt, drug_group ~ sex, value.var = "N", fill = 0)
  setnames(ev_cnt, c("F", "M"), c("a_F", "a_M"))
  tab <- merge(tot, ev_cnt, by = "drug_group", all.x = TRUE)
  tab[is.na(a_F), a_F := 0L]; tab[is.na(a_M), a_M := 0L]
  for (i in seq_len(nrow(tab))) {
    f <- c0$F; m <- c0$M
    f$nc <- f$nc - tab$n_F[i]; f$ac <- f$ac - tab$a_F[i]
    m$nc <- m$nc - tab$n_M[i]; m$ac <- m$ac - tab$a_M[i]
    z <- rorr_of(f, m)
    if (is.null(z)) next
    out <- rbind(out, data.table(pt = ev, removed = tab$drug_group[i],
                                 RORR = z$rorr, lci = z$lci, uci = z$uci))
  }
  out
}))

res[, `:=`(RORR = round(RORR, 3), lci = round(lci, 3), uci = round(uci, 3))]
fwrite(res, file.path(TAB, "leave_one_drug_out.csv"))

cat("\n== spread of the estimate across leave-one-out refits ==\n")
print(res[removed != "none", .(refits = .N, min = min(RORR), max = max(RORR),
                               median = median(RORR)), by = pt])
cat("\nmost influential removals for depression:\n")
print(res[pt == "Depression"][order(-abs(log(RORR) - log(RORR[removed == "none"])))][1:6])
cat("\n[written] leave_one_drug_out.csv\n[done]\n")
