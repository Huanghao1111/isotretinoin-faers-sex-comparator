suppressPackageStartupMessages({
  library(data.table)
})
ROOT <- "C:/Users/lenovo/WorkBuddy/2026-10-17-00-49-45"
PROC <- file.path(ROOT, "data", "processed")
DERM <- "psoriasis|psoriatic|dermatitis|eczema|acne|urticaria|angioedema|hidradenitis|alopecia|hair loss|vitiligo|melanoma|basal cell carcinoma|squamous cell carcinoma|actinic keratosis|skin cancer|skin neoplasm|rosacea|pemphigus|pemphigoid|lupus|pruritus|itch|onychomycosis|tinea|lichen planus|morphea|scleroderma|pyoderma|photo|seborrhoe|hyperhidrosis|cutane|skin|wart|cellulitis|impetigo|scabies|herpes|erythema nodosum|granuloma annulare|sarcoid|keratosis"
NON_AE <- "^(product use issue|off label use|unapproved drug use|product used? for unapproved indication|wrong (patient|drug|dose)|drug (dispensing|administration) error|medication error|incorrect dose administered|intercepted (drug|medication)|intentional product (misuse|use issue)|product (misuse|quality issue)|no adverse (event|drug effect)|treatment noncompliance|therapy (interrupted|cessation)|drug (withdrawn|dose omission)|circumstance or information capable of leading to|product (substitution|label|packaging) issue|expired product)"

demo <- fread(file.path(PROC, "demo.csv.gz"))[sex %in% c("M", "F"),
          .(primaryid, sex, age_yr)]
drug <- fread(file.path(PROC, "drug.csv.gz"))
reac <- fread(file.path(PROC, "reac.csv.gz"))
indi <- fread(file.path(PROC, "indi.csv.gz"))
keep <- intersect(demo$primaryid,
                  unique(indi[grepl(DERM, indi_pt, ignore.case = TRUE), primaryid]))
reac <- unique(reac[primaryid %in% keep & !is.na(pt) & pt != "" &
                      !grepl(NON_AE, pt, ignore.case = TRUE), .(primaryid, pt)])
iso <- drug[toupper(role_cod) == "PS" & drug_group == "isotretinoin" &
              primaryid %in% keep, unique(primaryid)]

for (restrict in c("all ages", "12-45")) {
  dv <- if (restrict == "all ages") demo[primaryid %in% keep] else
          demo[primaryid %in% keep & age_yr >= 12 & age_yr <= 45]
  ev <- "Depression"
  pid <- reac[pt == ev, unique(primaryid)]
  cat("\n===", restrict, "| event:", ev, "===\n")
  cell <- function(arm, sx) {
    n <- dv[primaryid %in% arm & sex == sx, .N]
    a <- dv[primaryid %in% arm & sex == sx & primaryid %in% pid, .N]
    c(n = n, a = a, odds = ifelse(n - a > 0, a / (n - a), NA))
  }
  f_i <- cell(iso, "F"); m_i <- cell(iso, "M")
  comp <- setdiff(dv$primaryid, iso)
  f_c <- cell(comp, "F"); m_c <- cell(comp, "M")
  cat(sprintf("index  F: n=%d a=%d odds=%.4f\n", f_i["n"], f_i["a"], f_i["odds"]))
  cat(sprintf("index  M: n=%d a=%d odds=%.4f\n", m_i["n"], m_i["a"], m_i["odds"]))
  cat(sprintf("compar F: n=%d a=%d odds=%.5f\n", f_c["n"], f_c["a"], f_c["odds"]))
  cat(sprintf("compar M: n=%d a=%d odds=%.5f\n", m_c["n"], m_c["a"], m_c["odds"]))
  ror_f <- f_i["odds"] / f_c["odds"]; ror_m <- m_i["odds"] / m_c["odds"]
  cat(sprintf("ROR_F=%.3f  ROR_M=%.3f  RORR=%.3f\n", ror_f, ror_m, ror_f / ror_m))
  cat(sprintf("check: index odds ratio F/M = %.3f\n", f_i["odds"] / m_i["odds"]))
  cat(sprintf("check: comparator odds ratio M/F = %.3f\n", m_c["odds"] / f_c["odds"]))
  cat(sprintf("product = %.3f\n",
              (f_i["odds"] / m_i["odds"]) * (m_c["odds"] / f_c["odds"])))
}
