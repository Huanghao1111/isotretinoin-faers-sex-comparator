#!/usr/bin/env Rscript
# ==============================================================================
# The comparator pool is not sex-neutral, and TableS16 shows how far that
# matters: within the dermatologic pool, reports of loss of libido are 102 times
# more frequent among men than women, and reports of dry lip are 8 times more
# frequent among women. Because ROR conditions on the comparator pool, that
# difference propagates straight into the ROR ratio and can manufacture a sex
# difference that has nothing to do with the index drug.
#
# The defence is an active comparator: drugs given to both sexes for the same
# indication, without the sex-skewed toxicity profile of finasteride and related
# agents. For acne, that is the oral tetracycline class. Isotretinoin is compared
# against doxycycline, minocycline and sarecycline directly, within each sex.
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "C:/Users/lenovo/WorkBuddy/2026-10-17-00-49-45"
PROC <- file.path(ROOT, "data", "processed")
TAB  <- file.path(ROOT, "output", "tables")
SUP  <- file.path(ROOT, "supplement")

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

ACNE <- "acne|hidradenitis|rosacea"
INDEX <- "isotretinoin"
COMPARE <- c("doxycycline", "minocycline", "sarecycline")

EVENTS <- c("Loss of libido", "Sexual dysfunction", "Libido decreased",
            "Depression", "Anxiety", "Completed suicide", "Suicidal ideation",
            "Major depression", "Anhedonia", "Lip dry", "Dry skin", "Epistaxis",
            "Acne", "Crohn's disease", "Inflammatory bowel disease")

demo <- fread(file.path(PROC, "demo.csv.gz"))
drug <- fread(file.path(PROC, "drug.csv.gz"))
reac <- fread(file.path(PROC, "reac.csv.gz"))
indi <- fread(file.path(PROC, "indi.csv.gz"))

demo <- demo[sex %in% c("M", "F"), .(primaryid, sex)]
acne_pid <- unique(indi[grepl(ACNE, indi_pt, ignore.case = TRUE), primaryid])
demo <- demo[primaryid %in% acne_pid]
cat("reports with an acne-spectrum indication:", nrow(demo), "\n")

dp <- unique(drug[primaryid %in% demo$primaryid & toupper(role_cod) == "PS",
                  .(primaryid, drug_group)])
grp <- dp[, .(primaryid,
              arm = fifelse(drug_group == INDEX, "isotretinoin",
                    fifelse(drug_group %in% COMPARE, "tetracycline",
                            NA_character_)))]
grp <- grp[!is.na(arm)]
grp <- grp[, .(arm = arm[1]), by = primaryid]     # one arm per report
reac <- unique(reac[primaryid %in% demo$primaryid &
                      !is.na(pt) & pt != "" &
                      !grepl(NON_AE, pt, ignore.case = TRUE), .(primaryid, pt)])

d <- demo[grp, on = "primaryid", nomatch = NULL]
cat("isotretinoin arm:", d[arm == "isotretinoin", .N],
    "| tetracycline arm:", d[arm == "tetracycline", .N], "\n")
print(d[, .(n = .N), by = .(arm, sex)])

res <- rbindlist(lapply(EVENTS, function(ev) {
  pid <- reac[pt == ev, unique(primaryid)]
  out <- data.table(pt = ev)
  for (sx in c("F", "M")) {
    a <- d[arm == "isotretinoin" & sex == sx & primaryid %in% pid, .N]
    b <- d[arm == "isotretinoin" & sex == sx, .N] - a
    c <- d[arm == "tetracycline" & sex == sx & primaryid %in% pid, .N]
    dd <- d[arm == "tetracycline" & sex == sx, .N] - c
    set(out, j = paste0("a_", sx), value = a)
    set(out, j = paste0("c_", sx), value = c)
    set(out, j = paste0("n_index_", sx), value = a + b)
    if (a > 0 && c > 0 && b > 0 && dd > 0) {
      r <- (a * dd) / (b * c)
      s <- sqrt(1/a + 1/b + 1/c + 1/dd)
      set(out, j = paste0("ror_", sx), value = round(r, 2))
      set(out, j = paste0("lci_", sx), value = round(exp(log(r) - 1.96 * s), 2))
      set(out, j = paste0("uci_", sx), value = round(exp(log(r) + 1.96 * s), 2))
    } else {
      set(out, j = paste0("ror_", sx), value = NA_real_)
      set(out, j = paste0("lci_", sx), value = NA_real_)
      set(out, j = paste0("uci_", sx), value = NA_real_)
    }
  }
  if (!is.na(out$ror_F) && !is.na(out$ror_M)) {
    lo <- log(out$ror_F) - log(out$ror_M)
    se <- sqrt((log(out$uci_F) - log(out$ror_F))^2 +
               (log(out$uci_M) - log(out$ror_M))^2) / 1.96
    out[, `:=`(rorr = round(exp(lo), 2),
               rorr_lci = round(exp(lo - 1.96 * se), 2),
               rorr_uci = round(exp(lo + 1.96 * se), 2),
               p_int = signif(2 * pnorm(-abs(lo / se)), 3))]
  } else {
    out[, `:=`(rorr = NA_real_, rorr_lci = NA_real_, rorr_uci = NA_real_,
               p_int = NA_real_)]
  }
  out
}))

res[, ci := ifelse(is.na(rorr), "",
                   sprintf("%.2f [%.2f-%.2f]", rorr, rorr_lci, rorr_uci))]
res[, srt := -abs(log(ifelse(is.na(rorr) | rorr <= 0, 1, rorr)))]
setorder(res, srt)
res[, srt := NULL]
fwrite(res, file.path(SUP, "TableS17_active_comparator_tetracycline.csv"))

cat("\n== isotretinoin versus oral tetracyclines, within each sex ==\n")
cat("reference = tetracycline arm; ROR ratio > 1 = stronger signal in women\n\n")
print(res[, .(pt, a_F, a_M, c_F, c_M, ror_F, ror_M, ci, p_int)], nrows = 20)

cat("\n[done]\n")
