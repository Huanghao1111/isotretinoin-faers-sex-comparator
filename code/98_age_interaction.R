suppressPackageStartupMessages({library(data.table)})
ROOT <- "C:/Users/lenovo/WorkBuddy/2026-10-17-00-49-45"
TAB  <- file.path(ROOT, "output", "tables")
DERM <- "psoriasis|psoriatic|dermatitis|eczema|acne|urticaria|angioedema|hidradenitis|alopecia|hair loss|vitiligo|melanoma|basal cell carcinoma|squamous cell carcinoma|actinic keratosis|skin cancer|skin neoplasm|rosacea|pemphigus|pemphigoid|lupus|pruritus|itch|onychomycosis|tinea|lichen planus|morphea|scleroderma|pyoderma|photo|seborrhoe|hyperhidrosis|cutane|skin|wart|cellulitis|impetigo|scabies|herpes|erythema nodosum|granuloma annulare|sarcoid|keratosis"
NON_AE <- "^(product use issue|off label use|unapproved drug use|product used? for unapproved indication|wrong (patient|drug|dose)|drug (dispensing|administration) error|medication error|incorrect dose administered|intercepted (drug|medication)|intentional product (misuse|use issue)|product (misuse|quality issue)|no adverse (event|drug effect)|treatment noncompliance|therapy (interrupted|cessation)|drug (withdrawn|dose omission)|circumstance or information capable of leading to|product (substitution|label|packaging) issue|expired product)"
PROC <- file.path(ROOT, "data", "processed")
demo <- fread(file.path(PROC, "demo.csv.gz"))[sex %in% c("M","F"), .(primaryid, sex, age_yr)]
drug <- fread(file.path(PROC, "drug.csv.gz")); reac <- fread(file.path(PROC, "reac.csv.gz"))
indi <- fread(file.path(PROC, "indi.csv.gz"))
keep <- intersect(demo$primaryid, unique(indi[grepl(DERM, indi_pt, ignore.case=TRUE), primaryid]))
demo <- demo[primaryid %in% keep]
ps <- unique(drug[primaryid %in% keep & toupper(role_cod)=="PS", .(primaryid, drug_group)])
rc <- unique(reac[primaryid %in% keep & !is.na(pt) & pt!="" & !grepl(NON_AE, pt, ignore.case=TRUE), .(primaryid, pt)])
pop <- intersect(ps$primaryid, rc$primaryid)
demo <- demo[primaryid %in% pop]; ps <- ps[primaryid %in% pop]; rc <- rc[primaryid %in% pop]
arm <- ps[, .(arm = if(any(drug_group=="isotretinoin")) "index" else "comparator"), by=primaryid]
arm <- arm[!duplicated(primaryid)]
demo <- arm[demo, on="primaryid", nomatch=NULL]

one <- function(dv, rcc, ev) {
  pid <- rcc[pt==ev, unique(primaryid)]
  orr <- function(sx) {
    ni <- dv[arm=="index" & sex==sx, .N]; ai <- dv[arm=="index" & sex==sx & primaryid %in% pid, .N]
    nc <- dv[arm=="comparator" & sex==sx, .N]; ac <- dv[arm=="comparator" & sex==sx & primaryid %in% pid, .N]
    if (min(ni-ai, nc-ac, ai, ac) <= 0) return(c(a=ai, ror=NA_real_, se=NA_real_))
    A<-ai; B<-ni-ai; C<-ac; D<-nc-ac
    r <- (A*D)/(B*C); s <- sqrt(1/A+1/B+1/C+1/D)
    c(a=ai, ror=r, se=s)
  }
  f <- orr("F"); m <- orr("M")
  if (any(is.na(c(f["ror"], m["ror"])))) return(NULL)
  lo <- log(f["ror"]) - log(m["ror"]); se <- sqrt(f["se"]^2 + m["se"]^2)
  data.table(pt=ev, a_F=f["a"], a_M=m["a"], rorr=round(exp(lo),2),
             lci=round(exp(lo-1.96*se),2), uci=round(exp(lo+1.96*se),2),
             ln=lo, se=se)
}
ev <- "Depression"
bands <- list("12-45"=demo[age_yr>=12 & age_yr<=45, primaryid],
              "46+"=demo[age_yr>45, primaryid])
res <- rbindlist(lapply(names(bands), function(b) one(demo[primaryid %in% bands[[b]]],
                                                       rc[primaryid %in% bands[[b]]], ev)[, band:=b]))
print(res[, .(band, a_F, a_M, rorr, lci, uci)])
if (nrow(res)==2) {
  dln <- res[band=="12-45", ln] - res[band=="46+", ln]
  se  <- sqrt(res[band=="12-45", se]^2 + res[band=="46+", se]^2)
  z <- dln/se
  cat(sprintf("\nratio of ratios, 12-45 vs 46+: %.2f (95%% CI %.2f-%.2f), p = %.3g\n",
              exp(dln), exp(dln-1.96*se), exp(dln+1.96*se), 2*pnorm(-abs(z))))
}
