suppressPackageStartupMessages({library(data.table)})
ROOT <- "C:/Users/lenovo/WorkBuddy/2026-10-17-00-49-45"
PROC <- file.path(ROOT, "data", "processed")
DERM <- paste0("psoriasis|psoriatic|dermatitis|eczema|acne|urticaria|angioedema|",
  "hidradenitis|alopecia|hair loss|vitiligo|melanoma|basal cell carcinoma|",
  "squamous cell carcinoma|actinic keratosis|skin cancer|skin neoplasm|",
  "rosacea|pemphigus|pemphigoid|lupus|pruritus|itch|onychomycosis|tinea|",
  "lichen planus|morphea|scleroderma|pyoderma|photo|seborrhoe|",
  "hyperhidrosis|cutane|skin|wart|cellulitis|impetigo|scabies|herpes|",
  "erythema nodosum|granuloma annulare|sarcoid|keratosis")
NON_AE <- paste0("^(product use issue|off label use|unapproved drug use|",
  "product used? for unapproved indication|wrong (patient|drug|dose)|",
  "drug (dispensing|administration) error|medication error|",
  "incorrect dose administered|intercepted (drug|medication)|",
  "intentional product (misuse|use issue)|product (misuse|quality issue)|",
  "no adverse (event|drug effect)|treatment noncompliance|",
  "therapy (interrupted|cessation)|drug (withdrawn|dose omission)|",
  "circumstance or information capable of leading to|",
  "product (substitution|label|packaging) issue|expired product)")
demo <- fread(file.path(PROC,"demo.csv.gz"))[sex %in% c("M","F"), .(primaryid,sex,age_yr)]
ps <- unique(fread(file.path(PROC,"drug.csv.gz"))[,.(primaryid,role_cod,drug_group)])
rc <- unique(fread(file.path(PROC,"reac.csv.gz"))[,.(primaryid,pt)])
indi <- unique(fread(file.path(PROC,"indi.csv.gz"))[,.(primaryid,indi_pt)])
demo <- demo[primaryid %in% unique(indi[grepl(DERM,indi_pt,ignore.case=TRUE),primaryid])]
ps <- ps[primaryid %in% demo$primaryid & toupper(role_cod)=="PS"]
rc <- rc[primaryid %in% demo$primaryid & !is.na(pt) & pt!="" & !grepl(NON_AE,pt,ignore.case=TRUE)]
pop <- intersect(ps$primaryid, rc$primaryid)
demo <- demo[primaryid %in% pop]; ps <- ps[primaryid %in% pop]; rc <- rc[primaryid %in% pop]
arm <- ps[,.(arm=if(any(drug_group=="isotretinoin")) "index" else "comparator"), by=primaryid]
arm <- arm[!duplicated(primaryid)]
demo <- arm[demo, on="primaryid", nomatch=NULL]
demo[, band := fifelse(is.na(age_yr),"unknown", fifelse(age_yr<12,"<12",
      fifelse(age_yr<=45,"12-45","46+")))]
ev <- "Depression"
pid <- rc[pt==ev, unique(primaryid)]
tab <- demo[, .(n=.N, a=sum(primaryid %in% pid)), by=.(band,arm,sex)]
setorder(tab, band, arm, sex)
print(tab, nrows=40)
cat("\n--- ROR by band and sex ---\n")
w <- dcast(tab, band+sex ~ arm, value.var=c("n","a"))
w[, ROR := round((a_index/(n_index-a_index))/(a_comparator/(n_comparator-a_comparator)),2)]
print(w[,.(band,sex,n_index,a_index,n_comparator,a_comparator,ROR)], nrows=40)
cat("\n--- RORR by band ---\n")
rr <- dcast(w[,.(band,sex,ROR)], band ~ sex, value.var="ROR")
rr[, RORR := round(F/M,2)]
print(rr, nrows=20)
