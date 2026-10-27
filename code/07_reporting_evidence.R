#!/usr/bin/env Rscript
# ==============================================================================
# Produces the quantitative evidence needed to satisfy RECORD / RECORD-PE items
# that cannot be answered with prose alone:
#
#   RECORD 13.1   detail of cohort selection (data-quality / availability filters)
#   RECORD-PE 8.a who generates the exposure records, and reporting route by sex
#   RECORD-PE 12a did the "reporting propensity cancels" assumption hold?
#   RECORD-PE 20.a confounding by indication: does the same drug carry the same
#                  indication in men and in women? This is the one threat that a
#                  within-sex comparator does NOT remove, so it is measured.
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "C:/Users/lenovo/WorkBuddy/2026-10-17-00-49-45"
PROC <- file.path(ROOT, "data", "processed")
TAB  <- file.path(ROOT, "output", "tables")
SUP  <- file.path(ROOT, "supplement")
dir.create(SUP, recursive = TRUE, showWarnings = FALSE)

DERM <- "psoriasis|psoriatic|dermatitis|eczema|acne|urticaria|angioedema|hidradenitis|alopecia|hair loss|vitiligo|melanoma|basal cell carcinoma|squamous cell carcinoma|actinic keratosis|skin cancer|skin neoplasm|rosacea|pemphigus|pemphigoid|lupus|pruritus|itch|onychomycosis|tinea|lichen planus|morphea|scleroderma|pyoderma|photo|seborrhoe|hyperhidrosis|cutane|skin|wart|cellulitis|impetigo|scabies|herpes|erythema nodosum|granuloma annulare|sarcoid|keratosis"
NON_AE <- "^(product use issue|off label use|unapproved drug use|product used? for unapproved indication|wrong (patient|drug|dose)|drug (dispensing|administration) error|medication error|incorrect dose administered|intercepted (drug|medication)|intentional product (misuse|use issue)|product (misuse|quality issue)|no adverse (event|drug effect)|treatment noncompliance|therapy (interrupted|cessation)|drug (withdrawn|dose omission)|circumstance or information capable of leading to|product (substitution|label|packaging) issue|expired product)"

demo <- fread(file.path(PROC, "demo.csv.gz"))
drug <- fread(file.path(PROC, "drug.csv.gz"))
indi <- fread(file.path(PROC, "indi.csv.gz"))
reac <- fread(file.path(PROC, "reac.csv.gz"))

pct <- function(x, d = 1) round(100 * x, d)

# ============================================================ RECORD 13.1 flow
# Counts before shipment come from the extraction log; everything from the
# "sex recorded" step onward is recomputed here so the table cannot go stale.
flow_pre <- data.table(
  step = 1:3,
  filter_applied = c(
    "Quarterly FAERS ASCII extracts available at the FDA public endpoint, 2012Q4 to 2026Q2",
    "Report lists at least one dermatologic agent as a suspect drug",
    "Deduplicated: one record per caseid, highest caseversion retained"),
  reports = c(NA_integer_, NA_integer_, 3897025L),
  dropped = NA_integer_,
  note = c("FDA retains no earlier quarters at this endpoint; 55 quarters used",
           "Identified by word-bounded name match on DRUG.drugname and DRUG.prod_ai",
           "Earlier quarterly files carry multiple versions of the same case")
)

n_sex <- nrow(demo)
indi_u <- unique(indi[!is.na(indi_pt) & indi_pt != "", .(primaryid, indi_pt)])
derm_pid <- indi_u[grepl(DERM, indi_pt, ignore.case = TRUE), unique(primaryid)]
n_derm <- length(intersect(demo$primaryid, derm_pid))

keep <- intersect(demo$primaryid, derm_pid)
d2 <- demo[primaryid %in% keep]
g2 <- drug[primaryid %in% keep]
r2t <- reac[primaryid %in% keep]

n_ps <- uniqueN(g2[toupper(role_cod) == "PS", primaryid])
ps_pid <- unique(g2[toupper(role_cod) == "PS", primaryid])
# retained terms are counted only among reports that still have a primary
# suspect after step 6, otherwise this row adds reports back in
r2t <- r2t[primaryid %in% ps_pid]
r2t <- r2t[!is.na(pt) & pt != "" & !grepl(NON_AE, pt, ignore.case = TRUE)]
n_pt <- uniqueN(r2t$primaryid)
drug_ps <- unique(g2[toupper(role_cod) == "PS", .(primaryid, drug_group)])
pairs <- drug_ps[r2t[, .(primaryid, pt)], on = "primaryid",
                 allow.cartesian = TRUE, nomatch = NULL][, .(primaryid, drug_group, pt)]
n_pair_reports <- uniqueN(pairs$primaryid)
n_pairs <- uniqueN(pairs[, .(drug_group, pt)])

flow_post <- data.table(
  step = 4:9,
  filter_applied = c(
    "Sex recorded as male or female",
    "At least one dermatologic indication recorded (INDI.indi_pt)",
    "Index agent recorded as primary suspect (DRUG.role_cod = PS)",
    "At least one adverse event preferred term retained after removing non-AE terms",
    "Report contributes at least one drug-event pair",
    "Distinct drug-event pairs before applying the frequency filter"),
  reports = c(n_sex, n_derm, n_ps, n_pt, n_pair_reports, n_pairs),
  dropped = c(3897025L - n_sex, n_sex - n_derm, n_derm - n_ps, n_ps - n_pt,
              n_pt - n_pair_reports, NA_integer_),
  note = c("Sex missing or uncodable in 11.1% of deduplicated reports",
           "Removes rheumatologic and gastroenterologic use of the same agents",
           "Concomitant-suspect reports excluded from the primary analysis",
           "Removed terms describe handling events, not toxicity",
           "No further loss: every report reaching this step contributes",
           "Unit of analysis hereafter; the 50-report floor leaves 21,565 pairs")
)

flow <- rbindlist(list(flow_pre, flow_post), fill = TRUE)
fwrite(flow, file.path(SUP, "TableS8_cohort_selection_flow.csv"))
cat("[S8] cohort selection flow\n"); print(flow, nrows = 20)

# ============================================================ RECORD-PE 12a
# Assumption under test: within a sex, index and comparator reports share the
# same reporting route mix, so differential reporting propensity cancels in the
# ratio of ratios. If men and women reach FAERS by different routes, that is
# visible here and the health-professional restriction becomes the arbiter.
iso <- g2[toupper(role_cod) == "PS" & drug_group == "isotretinoin", unique(primaryid)]
route_sex <- d2[primaryid %in% iso, .(n = .N), by = .(sex, occp_cod)]
route_sex[, pct_within_sex := pct(n / sum(n)), by = sex]
setorder(route_sex, sex, -n)
fwrite(route_sex, file.path(SUP, "TableS9_reporter_by_sex_isotretinoin.csv"))
cat("\n[S9] reporter occupation by sex, isotretinoin reports\n")
print(route_sex, nrows = 20)

# ============================================================ RECORD-PE 20.a
# Confounding by indication survives a within-sex comparator: it only fails if
# the SAME drug is used for different reasons in men and women. Measured, not
# assumed, for the headline drug.
iso_indi <- indi_u[primaryid %in% iso]
iso_indi <- d2[, .(primaryid, sex)][iso_indi, on = "primaryid", nomatch = NULL]
iso_indi[, indi_pt := tolower(trimws(indi_pt))]
# NOTE: pregnancy prevention is checked before the administrative catch-all,
# because "product used for unknown indication" and contraception mean entirely
# different things and conflating them hides the teratogenicity signal
iso_indi[, indi_grp := fifelse(grepl("contracept|pregnan", indi_pt),
                               "Pregnancy prevention",
                       fifelse(grepl("acne", indi_pt), "Acne",
                       fifelse(grepl("hidradenitis|rosacea|seborrhoe", indi_pt),
                               "Other inflammatory dermatosis",
                       fifelse(grepl("product used|prophylaxis|off label", indi_pt),
                               "Product-use / administrative", "Other"))))]
# a report can carry several indication terms, so both units are given
nx <- iso_indi[, .(reports_total = uniqueN(primaryid)), by = sex]
ind_tab <- iso_indi[, .(indication_rows = .N, reports = uniqueN(primaryid)),
                    by = .(sex, indi_grp)]
ind_tab <- nx[ind_tab, on = "sex"]
ind_tab[, pct_of_sex := pct(reports / reports_total)][, reports_total := NULL]
setorder(ind_tab, sex, -reports)
fwrite(ind_tab, file.path(SUP, "TableS10_indication_by_sex_isotretinoin.csv"))
cat("\n[S10] indication profile by sex, isotretinoin reports\n")
print(ind_tab, nrows = 12)

# top literal preferred terms, so no post-hoc grouping is being hidden
top_raw <- iso_indi[, .(n = .N), by = .(sex, indi_pt)]
top_raw <- top_raw[, .SD[order(-n)][1:8], by = sex]
fwrite(top_raw, file.path(SUP, "TableS11_top_indication_terms_isotretinoin.csv"))
cat("\n[S11] most frequent literal indication terms, by sex\n")
print(top_raw, nrows = 20)

# ======================================================= cohort availability
yr <- d2[!is.na(event_year), .(reports = .N, female = sum(sex == "F"),
                               male = sum(sex == "M")), by = event_year]
yr[, pct_female := pct(female / reports)]
setorder(yr, event_year)
fwrite(yr, file.path(SUP, "TableS12_reports_by_year.csv"))
cat("\n[S12] reports by calendar year\n"); print(yr, nrows = 20)

# The README is written in full by 06_supplement.R and is not extended here.
# Appending to it produced duplicated and out-of-order entries, and the marker
# file that guarded the append survived a README rewrite, so the entries could
# also be lost. The descriptions for TableS8-S12 live in 06.
if (!file.exists(file.path(SUP, "README.md"))) {
  cat("[warn] README.md missing: run 06_supplement.R first\n")
}

cat("\n[done] reporting-evidence tables written to", SUP, "\n")
