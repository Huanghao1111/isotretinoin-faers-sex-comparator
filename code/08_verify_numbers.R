#!/usr/bin/env Rscript
# ==============================================================================
# Verifies every number quoted in the manuscript against the file that produced
# it. A transcription slip between a result table and the text is invisible to
# proof-reading and fatal under review, so it is checked here.
#
# Exit status is non-zero if any check fails, so this can run before deposit.
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "C:/Users/lenovo/WorkBuddy/2026-10-17-00-49-45"
TAB  <- file.path(ROOT, "output", "tables")
PROC <- file.path(ROOT, "data", "processed")
MS   <- file.path(ROOT, "manuscript", "JAAD_BriefReport_isotretinoin_sex.md")

ok <- 0L; bad <- 0L
check <- function(label, truth, claimed, tol = 0.005) {
  pass <- isTRUE(all.equal(as.numeric(truth), as.numeric(claimed), tolerance = tol))
  if (pass) ok <<- ok + 1L else bad <<- bad + 1L
  cat(sprintf("[%s] %-54s file=%-12s manuscript=%s\n",
              ifelse(pass, "PASS", "FAIL"), label, truth, claimed))
}
check_txt <- function(label, pattern) {
  pass <- any(grepl(pattern, readLines(MS, warn = FALSE), fixed = TRUE))
  if (pass) ok <<- ok + 1L else bad <<- bad + 1L
  cat(sprintf("[%s] %-54s \"%s\" present\n",
              ifelse(pass, "PASS", "FAIL"), label, pattern))
}
nums <- function(s) as.numeric(regmatches(s, gregexpr("[0-9]+\\.?[0-9]*", s))[[1]])

t1   <- fread(file.path(TAB, "table1_isotretinoin.csv"))
pri  <- fread(file.path(TAB, "final_rorr_all_ages.csv"))
yng  <- fread(file.path(TAB, "final_rorr_age12_45.csv"))
rob  <- fread(file.path(TAB, "robustness_trimming_rules.csv"))
fin  <- fread(file.path(TAB, "finasteride_influence.csv"))
A    <- readRDS(file.path(PROC, "analysis_objects.rds"))
dv   <- as.data.table(A$demo)

cat("=== cohort ===\n")
check("reports after deduplication", nrow(dv), 1148083, tol = 0)
check("female reports", dv[sex == "F", .N], 709690, tol = 0)
check("female percent", round(100 * dv[sex == "F", .N] / nrow(dv), 1), 61.8)

cat("\n=== analysis population and index arm ===\n")
nz <- rob[rule == "no trimming"][1]
pop <- nz$n_iF + nz$n_iM + nz$n_cF + nz$n_cM
check("analysis population", pop, 1091621, tol = 0)
check("isotretinoin reports, total", nz$n_iF + nz$n_iM, 13829, tol = 0)
check("isotretinoin reports, female", nz$n_iF, 7496, tol = 0)
check("isotretinoin reports, male", nz$n_iM, 6333, tol = 0)

cat("\n=== Table I, rows quoted in the body ===\n")
# preferred term, cases F, cases M, uncorrected RORR, lci, uci,
# primary RORR, lci, uci, comparator factor, index OR
quoted <- list(
  list("Depression",        1366, 1440, 0.97, 0.89, 1.07, 0.59, 0.54, 0.65, 0.78, 0.76),
  list("Anxiety",            436,  394, 1.09, 0.94, 1.27, 0.65, 0.56, 0.76,   NA,   NA),
  list("Suicidal ideation",  350,  463, 1.30, 1.10, 1.53, 0.64, 0.54, 0.76,   NA,   NA),
  list("Completed suicide",   22,   68, 6.15, 3.22, 11.76, 3.26, 1.68, 6.33, 12.03, NA),
  list("Loss of libido",      22,   40, 64.26, 30.62, 134.85, 6.71, 2.83, 15.91, 14.48, NA),
  list("Major depression",    29,   21, 2.49, 1.35, 4.58, 0.89, 0.47, 1.70,   NA,   NA),
  list("Lip dry",            475,  484, 0.21, 0.17, 0.26, 0.21, 0.16, 0.26,   NA,   NA),
  list("Crohn's disease",    553,  840, 0.39, 0.33, 0.45, 0.38, 0.33, 0.44,   NA,   NA)
)
for (q in quoted) {
  ptv <- q[[1]]
  r <- t1[Event == ptv]
  if (!nrow(r)) { bad <- bad + 1L
    cat(sprintf("[FAIL] %-54s row missing\n", ptv)); next }
  check(paste0(ptv, ": cases female"), r$`Cases, female`, q[[2]], tol = 0)
  check(paste0(ptv, ": cases male"),   r$`Cases, male`,   q[[3]], tol = 0)
  u <- nums(r$`ROR ratio, all-drug comparator (95% CI)`)
  check(paste0(ptv, ": all-drug RORR"), u[1], q[[4]])
  check(paste0(ptv, ": all-drug lower"), u[2], q[[5]])
  check(paste0(ptv, ": all-drug upper"), u[3], q[[6]])
  v <- nums(r$`ROR ratio, sex-neutral comparator (95% CI)`)
  check(paste0(ptv, ": sex-neutral RORR"), v[1], q[[7]])
  check(paste0(ptv, ": sex-neutral lower"), v[2], q[[8]])
  check(paste0(ptv, ": sex-neutral upper"), v[3], q[[9]])
  if (!is.na(q[[10]])) check(paste0(ptv, ": comparator factor"),
                             r$`Comparator OR, men vs women`, q[[10]])
  if (!is.na(q[[11]])) check(paste0(ptv, ": isotretinoin OR, women vs men"),
                             r$`Isotretinoin OR, women vs men`, q[[11]])
}

cat("\n=== ages 12-45, quoted in the body ===\n")
check("depression 12-45, all-drug",
      fread(file.path(TAB, "final_rorr_untrimmed_age12_45.csv"))[pt == "Depression", RORR],
      1.28)
check("depression 12-45, all-drug lower",
      fread(file.path(TAB, "final_rorr_untrimmed_age12_45.csv"))[pt == "Depression", lci],
      1.12)
check("depression 12-45, all-drug upper",
      fread(file.path(TAB, "final_rorr_untrimmed_age12_45.csv"))[pt == "Depression", uci],
      1.45)
check("depression 12-45, sex-neutral", yng[pt == "Depression", RORR], 0.58)
check("depression 12-45, sex-neutral lower", yng[pt == "Depression", lci], 0.50)
check("depression 12-45, sex-neutral upper", yng[pt == "Depression", uci], 0.67)

cat("\n=== comparator trimming rules quoted in the body ===\n")
for (z in list(c("finasteride only", 0.59), c("30-70% female", 0.56),
               c("40-60% female", 0.46))) {
  check(paste0("depression under '", z[1], "'"),
        rob[rule == z[1] & pt == "Depression", RORR], as.numeric(z[2]))
}
check("depression index OR is invariant to trimming",
      rob[pt == "Depression", index_OR_F_vs_M][1], 0.76)

cat("\n=== finasteride mechanism ===\n")
m <- fin[sex == "M"]
check("male comparator reports, 12-45", m$n_comparator, 92976, tol = 0)
check("of which finasteride", m$n_finasteride, 2572, tol = 0)
check("male comparator depression reports", m$n_depression, 1385, tol = 0)
check("of which on finasteride", m$n_depression_finasteride, 779, tol = 0)
check("depression per 1000 on finasteride", m$depression_per_1000_finasteride, 302.9)
check("depression per 1000 on dupilumab", m$depression_per_1000_dupilumab, 2.4)

cat("\n=== restriction scenarios on the corrected design ===\n")
sen2 <- fread(file.path(TAB, "primary_sensitivity_scenarios.csv"))
dep <- sen2[pt == "Depression"]
check("depression, main", dep[scenario == "A_main", RORR], 0.59)
check("depression, single drug", dep[scenario == "B_single_drug", RORR], 0.59)
check("depression, health professional", dep[scenario == "C_healthprof", RORR], 0.55)
check("depression, ages 12-45", dep[scenario == "D_age12_45", RORR], 0.58)
check("depression, no lawyer", dep[scenario == "F_no_lawyer", RORR], 0.57)
check("depression, lowest across scenarios", min(dep$RORR), 0.55)
check("depression, highest across scenarios", max(dep$RORR), 0.59)
pass <- all(dep$RORR < 1) && all(dep$uci < 1)
if (pass) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] %-54s every scenario below 1, upper bound too\n",
            ifelse(pass, "PASS", "FAIL"), "depression direction stability"))

cat("\n=== age-band interaction on the primary comparator ===\n")
age <- fread(file.path(TAB, "primary_age_band_interaction.csv"))
dr <- age[pt == "Depression"]
# column names come from the object written by 20_primary_sensitivity.R, not
# from the subset printed to the console; a wrong name yields NULL, sprintf then
# emits nothing, and the check fails silently. Assert the columns exist first.
need <- c("rorr_y", "rorr_o", "ratio", "r_lci", "r_uci", "p_int", "q_int")
miss <- setdiff(need, names(dr))
if (length(miss)) {
  bad <- bad + 1L
  cat(sprintf("[FAIL] %-54s missing columns: %s\n", "age interaction table",
              paste(miss, collapse = ", ")))
}
check("depression 12-45 RORR", dr$rorr_y, 0.58)
check("depression 46+ RORR", dr$rorr_o, 1.40)
check("interaction ratio", dr$ratio, 0.41)
check("interaction lower", dr$r_lci, 0.23)
check("interaction upper", dr$r_uci, 0.76)
check("interaction p", dr$p_int, 0.0041, tol = 0.02)
pass <- dr$q_int < 0.05
if (pass) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] %-54s q = %s (survives FDR)\n",
            ifelse(pass, "PASS", "FAIL"), "interaction survives BH", dr$q_int))

cat("\n=== internal control: the correction is targeted ===\n")
sh <- fread(file.path(TAB, "comparator_shift_by_panel.csv"))
med <- sh[, .(m = median(abs_log_shift)), by = panel]
check("median shift, mucocutaneous / GI",
      med[panel == "mucocutaneous / GI", m], 0.02)
check("median shift, neuropsychiatric",
      med[panel == "neuropsychiatric", m], 0.67)
check("median shift, sexual", med[panel == "sexual", m], 2.26)
pass <- med[panel == "mucocutaneous / GI", m] <
        med[panel == "neuropsychiatric", m]
if (pass) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] %-54s %s < %s\n", ifelse(pass, "PASS", "FAIL"),
            "somatic terms moved less than psychiatric",
            round(med[panel == "mucocutaneous / GI", m], 2),
            round(med[panel == "neuropsychiatric", m], 2)))

cat("\n=== calendar-time stability (quoted in the technical note) ===\n")
NOTE <- file.path(ROOT, "supplement", "NOTE_comparator_design.md")
ntxt <- paste(readLines(NOTE, warn = FALSE, encoding = "UTF-8"), collapse = " ")
txt  <- paste(readLines(MS, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
ts <- fread(file.path(TAB, "time_stability_rorr.csv"))
th <- fread(file.path(TAB, "time_stability_heterogeneity.csv"))
dep <- ts[pt == "Depression" & band == "all"][order(era)]
for (i in seq_len(nrow(dep))) {
  s <- sprintf("%.2f (%.2f-%.2f)", dep$RORR[i], dep$lci[i], dep$uci[i])
  pass <- grepl(s, ntxt, fixed = TRUE)
  if (pass) ok <- ok + 1L else bad <- bad + 1L
  cat(sprintf("[%s] %-54s note=%s  file=%s\n", ifelse(pass, "PASS", "FAIL"),
              paste("depression", dep$era[i]), ifelse(pass, "quoted", "MISSING"), s))
}
dep45 <- ts[pt == "Depression" & band == "12-45"][order(era)]
for (i in seq_len(nrow(dep45))) {
  s <- sprintf("%.2f (%.2f-%.2f)", dep45$RORR[i], dep45$lci[i], dep45$uci[i])
  pass <- grepl(s, ntxt, fixed = TRUE)
  if (pass) ok <- ok + 1L else bad <- bad + 1L
  cat(sprintf("[%s] %-54s note=%s  file=%s\n", ifelse(pass, "PASS", "FAIL"),
              paste("depression 12-45,", dep45$era[i]),
              ifelse(pass, "quoted", "MISSING"), s))
}
era_s <- sprintf("%.2f, against %.2f and %.2f earlier",
                 dep45[era == "2022-2026", RORR], dep45[era == "2012-2016", RORR],
                 dep45[era == "2017-2021", RORR])
pass <- grepl(era_s, txt, fixed = TRUE)
if (pass) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] %-54s body=%s  file=%s\n", ifelse(pass, "PASS", "FAIL"),
            "era split quoted in the body", ifelse(pass, "quoted", "MISSING"), era_s))

for (bd in c("all", "12-45")) {
  row <- th[pt == "Depression" & band == bd]
  s <- sprintf("%.2f (%.2f-%.2f)", row$pooled_RORR, row$pooled_lci, row$pooled_uci)
  pass <- grepl(s, ntxt, fixed = TRUE)
  if (pass) ok <- ok + 1L else bad <- bad + 1L
  cat(sprintf("[%s] %-54s note=%s  file=%s\n", ifelse(pass, "PASS", "FAIL"),
              paste("era-pooled estimate,", bd), ifelse(pass, "quoted", "MISSING"), s))
}
for (tm in c("Anxiety", "Suicidal ideation")) {
  v <- ts[pt == tm & band == "12-45"][order(era), sprintf("%.2f", RORR)]
  s <- sprintf("%s, %s and %s", v[1], v[2], v[3])
  pass <- grepl(s, ntxt, fixed = TRUE)
  if (pass) ok <- ok + 1L else bad <- bad + 1L
  cat(sprintf("[%s] %-54s note=%s  file=%s\n", ifelse(pass, "PASS", "FAIL"),
              paste(tm, "by era at 12-45"), ifelse(pass, "quoted", "MISSING"), s))
}
pass <- grepl("0.79 when the eras are pooled", txt, fixed = TRUE) &&
  abs(th[pt == "Depression" & band == "12-45", pooled_RORR] - 0.788) < 0.005
if (pass) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] %-54s file=%.3f\n", ifelse(pass, "PASS", "FAIL"),
            "era-pooled 12-45 value quoted in the body",
            th[pt == "Depression" & band == "12-45", pooled_RORR]))

p45 <- th[pt == "Depression" & band == "12-45", p_heterogeneity]
pass <- grepl(sprintf("p = %.3f", p45), ntxt, fixed = TRUE)
if (pass) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] %-54s \n", ifelse(pass, "PASS", "FAIL"),
            "12-45 heterogeneity p = 0.002 quoted"))
shares <- dep45[, round(100 * n_iM / (n_iF + n_iM), 1)]
pass <- all(vapply(sprintf("%.1f%%", shares), function(x) grepl(x, ntxt, fixed = TRUE),
                   logical(1)))
if (pass) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] %-54s file=%s\n", ifelse(pass, "PASS", "FAIL"),
            "male share of the 12-45 index arm quoted",
            paste(sprintf("%.1f%%", shares), collapse = ", ")))

p_het <- th[pt == "Depression" & band == "all", p_heterogeneity]
pass <- grepl(sprintf("p = %.3f", p_het), ntxt, fixed = TRUE)
if (pass) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] %-54s note=%s  file=%.3f\n", ifelse(pass, "PASS", "FAIL"),
            "era heterogeneity p quoted", ifelse(pass, "quoted", "MISSING"), p_het))

fs <- fread(file.path(TAB, "finasteride_share_by_era.csv"))[band == "12-45"][order(era)]
shr <- sprintf("%.1f%%", fs$pct_depression_on_finasteride)
pass <- all(vapply(shr, function(x) grepl(x, ntxt, fixed = TRUE), logical(1)))
if (pass) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] %-54s note=%s  file=%s\n", ifelse(pass, "PASS", "FAIL"),
            "finasteride share by era (12-45)", ifelse(pass, "quoted", "MISSING"),
            paste(shr, collapse = ", ")))
pass <- length(unique(dep$era)) == 3L && nrow(fs) == 3L
if (pass) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] %-54s %d eras x 3 finasteride rows\n", ifelse(pass, "PASS", "FAIL"),
            "three eras present", length(unique(dep$era))))

cat("\n=== pregnancy-prevention channel (quoted in the technical note) ===\n")
ms <- fread(file.path(TAB, "monitoring_sensitivity.csv"))
mc <- fread(file.path(TAB, "monitoring_report_counts.csv"))
idx <- mc[arm == "index"]
pass <- all(vapply(sprintf("%.1f%%", idx$pct), function(x) grepl(x, ntxt, fixed = TRUE),
                   logical(1)))
if (pass) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] %-54s file=%s\n", ifelse(pass, "PASS", "FAIL"),
            "monitoring-only share, both sexes",
            paste(sprintf("%s %.1f%%", idx$sex, idx$pct), collapse = ", ")))
for (bd in c("all", "12-45")) {
  for (tm in c("Depression", "Anxiety", "Suicidal ideation")) {
    # the note quotes depression and suicidal ideation at 12-45, not anxiety
    if (bd == "12-45" && tm == "Anxiety") next
    row <- ms[pt == tm & band == bd & grepl("^monitoring", scenario)]
    s <- sprintf("%.2f (%.2f-%.2f)", row$RORR, row$lci, row$uci)
    hit <- grepl(s, ntxt, fixed = TRUE)
    if (hit) ok <- ok + 1L else bad <- bad + 1L
    cat(sprintf("[%s] %-54s file=%s\n", ifelse(hit, "PASS", "FAIL"),
                paste(tm, bd, "without monitoring"), s))
  }
}

cat("\n=== age strata, including the unrecorded band ===\n")
bands <- fread(file.path(TAB, "primary_rorr_by_age_band.csv"))
for (bd in c("12-45", "46+", "unknown")) {
  row <- bands[pt == "Depression" & band == bd]
  s <- sprintf("%.2f (%.2f-%.2f)", row$RORR, row$lci, row$uci)
  hit <- grepl(s, ntxt, fixed = TRUE)
  if (hit) ok <- ok + 1L else bad <- bad + 1L
  cat(sprintf("[%s] %-54s note=%s  file=%s\n", ifelse(hit, "PASS", "FAIL"),
              paste("depression, band", bd), ifelse(hit, "quoted", "MISSING"), s))
}

cat("\n=== leave-one-drug-out ===\n")
loo <- fread(file.path(TAB, "leave_one_drug_out.csv"))
for (tm in c("Depression", "Anxiety", "Suicidal ideation")) {
  z <- loo[pt == tm & removed != "none"]
  # the note writes some ranges with "to" and some with a hyphen
  cand <- c(sprintf("%.3f-%.3f", min(z$RORR), max(z$RORR)),
            sprintf("%.3f to %.3f", min(z$RORR), max(z$RORR)))
  s <- cand[1]
  hit <- any(vapply(cand, function(x) grepl(x, ntxt, fixed = TRUE), logical(1)))
  if (hit) ok <- ok + 1L else bad <- bad + 1L
  cat(sprintf("[%s] %-54s file=%s (%d refits)\n", ifelse(hit, "PASS", "FAIL"),
              paste("leave-one-out range,", tm), s, nrow(z)))
}
base <- loo[pt == "Depression" & removed == "none", RORR]
pri_dep <- fread(file.path(TAB, "final_rorr_all_ages.csv"))[pt == "Depression", RORR]
pass <- abs(base - pri_dep) < 0.01
if (pass) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] %-54s file=%.3f  primary=0.588\n", ifelse(pass, "PASS", "FAIL"),
            "leave-one-out baseline equals the primary estimate", base))

cat("\n=== indication-matched active comparator ===\n")
ac <- fread(file.path(TAB, "active_comparator_acne.csv"))
for (tm in c("Depression", "Anxiety")) {
  row <- ac[pt == tm]
  # the body writes "0.88, 0.52-1.49"; the note writes "0.88 (0.52-1.49)"
  cand <- c(sprintf("%.2f, %.2f-%.2f", row$RORR, row$lci, row$uci),
            sprintf("%.2f (%.2f-%.2f)", row$RORR, row$lci, row$uci))
  s <- cand[1]
  hit <- any(vapply(cand, function(x) grepl(x, txt, fixed = TRUE) ||
                      grepl(x, ntxt, fixed = TRUE), logical(1)))
  if (hit) ok <- ok + 1L else bad <- bad + 1L
  cat(sprintf("[%s] %-54s file=%s\n", ifelse(hit, "PASS", "FAIL"),
              paste(tm, "versus tetracyclines in acne"), s))
}

cat("\n=== figure panel A (quoted in the figure legend) ===\n")
txt <- paste(readLines(MS, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
pa  <- fread(file.path(TAB, "panel_a_drug_rates.csv"))
pool <- fread(file.path(TAB, "panel_a_pooled.csv"))
poolM <- pool[sex == "M", rate_per_1000]
fin_r <- pa[drug_group == "finasteride", rate_M]
pass <- grepl(sprintf("%.1f per 1,000", poolM), txt, fixed = TRUE)
if (pass) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] %-54s file=%.1f  legend=%s\n", ifelse(pass, "PASS", "FAIL"),
            "pooled male comparator rate quoted", poolM,
            ifelse(pass, "quoted", "MISSING")))
pass <- grepl(sprintf("%.1f per 1,000", fin_r), txt, fixed = TRUE)
if (pass) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] %-54s file=%.1f  legend=%s\n", ifelse(pass, "PASS", "FAIL"),
            "finasteride male rate quoted", fin_r,
            ifelse(pass, "quoted", "MISSING")))
pass <- grepl(sprintf("%d drugs carrying at least 200 male reports", nrow(pa)), txt,
              fixed = TRUE)
if (pass) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] %-54s file=%d  legend=%s\n", ifelse(pass, "PASS", "FAIL"),
            "drug count in panel A quoted", nrow(pa),
            ifelse(pass, "quoted", "MISSING")))
second <- pa[drug_group != "finasteride"][order(-rate_M)][1]
pass <- grepl(sprintf("%.1f per 1,000", second$rate_M), ntxt, fixed = TRUE)
if (pass) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] %-54s note=%s  file=%.1f (%s)\n", ifelse(pass, "PASS", "FAIL"),
            "second-highest drug rate quoted in the note",
            ifelse(pass, "quoted", "MISSING"), second$rate_M, second$drug_group))
omitted <- pa[drug_group != "finasteride"][order(-rate_M)][-(1:14)]
pass <- nrow(omitted) == nrow(pa) - 15L && max(omitted$rate_M) < 10 &&
  grepl(sprintf("the other %d all below 10 per 1,000", nrow(omitted)), txt, fixed = TRUE)
if (pass) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] %-54s %d drugs below 10 per 1,000\n", ifelse(pass, "PASS", "FAIL"),
            "omitted drugs and their ceiling", nrow(omitted)))
fin <- pa[drug_group == "finasteride"]
rate_ex <- 1000 * (pool[sex == "M", d] - fin$dep_M) / (pool[sex == "M", n] - fin$n_M)
pass <- grepl(sprintf("without finasteride it is %.1f", rate_ex), txt, fixed = TRUE) &&
  grepl(sprintf("%.1f to %.1f depression", poolM, rate_ex), ntxt, fixed = TRUE)
if (pass) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] %-54s file=%.1f per 1,000 without finasteride\n",
            ifelse(pass, "PASS", "FAIL"), "pool rate without finasteride quoted", rate_ex))
pass <- grepl("15 comparator drugs", txt, fixed = TRUE)
if (pass) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] %-54s \n", ifelse(pass, "PASS", "FAIL"),
            "panel A shows the top 15"))

cat("\n=== numbers in the body that no other block covers ===\n")
si <- fread(file.path(TAB, "primary_rorr_by_age_band.csv"))[pt == "Suicidal ideation" & band == "12-45"]
si_s <- sprintf("%.2f (%.2f-%.2f)", si$RORR, si$lci, si$uci)
pass <- grepl(si_s, txt, fixed = TRUE)
if (pass) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] %-54s body=%s  file=%s\n", ifelse(pass, "PASS", "FAIL"),
            "suicidal ideation at 12-45 quoted", ifelse(pass, "quoted", "MISSING"), si_s))

cat("\n=== in-text citations ===\n")
for (k in 1:5) {
  sup <- c("\u00b9", "\u00b2", "\u00b3", "\u2074", "\u2075")[k]
  hit <- grepl(sup, txt, fixed = TRUE)
  if (hit) ok <- ok + 1L else bad <- bad + 1L
  cat(sprintf("[%s] %-54s citation %d present in the body\n",
              ifelse(hit, "PASS", "FAIL"), "superscript marker", k))
}

cat("\n=== word limit and display budget ===\n")
txt <- paste(readLines(MS, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
body <- sub(".*\\*To the Editor:\\*", "", txt)
body <- sub("Word count.*", "", body)
nw <- length(strsplit(trimws(body), "\\s+")[[1]])
pass <- nw <= 500
if (pass) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] %-54s body words = %d (limit 500)\n",
            ifelse(pass, "PASS", "FAIL"), "word count", nw))
declared <- as.integer(sub(".*Word count \\(body\\): ([0-9]+).*", "\\1", txt))
pass <- !is.na(declared) && abs(declared - nw) <= 1
if (pass) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] %-54s declared %s, computed %d\n", ifelse(pass, "PASS", "FAIL"),
            "declared word count matches", as.character(declared), nw))

nref <- length(grep("^\\d+\\. ", readLines(MS, warn = FALSE)))
pass <- nref <= 5
if (pass) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] %-54s references = %d (limit 5)\n",
            ifelse(pass, "PASS", "FAIL"), "reference count", nref))

cat("\n=== comparator pool size, body number against the file ===\n")
pool <- fread(file.path(TAB, "comparator_pool_size.csv"))
rm_txt <- pool[grepl("^comparator reports removed$", step), n]
bs <- regmatches(txt, regexpr("[0-9,]+ comparator reports removed", txt))
body_rm <- if (length(bs)) as.numeric(gsub("[^0-9]", "", bs)) else NA_real_
pass <- length(body_rm) == 1 && !is.na(body_rm) && body_rm == rm_txt
if (pass) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] %-54s file=%s  manuscript=%s\n", ifelse(pass, "PASS", "FAIL"),
            "reports removed from comparator pool",
            format(rm_txt, big.mark = ","),
            ifelse(is.na(body_rm), "NA", format(body_rm, big.mark = ","))))

cat("\n=== strings that must survive editing ===\n")

# the 46+ interaction was removed from the body: the age contrast belongs to the
# comparator, and the within-drug contrast is what is now reported
s5 <- readLines(MS, warn = FALSE, encoding = "UTF-8")
body_txt <- paste(s5, collapse = " ")
pass <- !grepl("0.41, 0.23-0.76; p = 0.004", body_txt, fixed = TRUE)
if (pass) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] %-54s\n", ifelse(pass, "PASS", "FAIL"),
            "the comparator-driven 46+ interaction is not in the body"))

for (s in c("1,148,083", "61.8%",
            "2,572 of 92,976", "779 of 1,385", "302.9 per 1,000",
            "0.97 (95% CI 0.89-1.07)", "1.28 (1.12-1.45)", "0.59 (0.54-0.65)",
            "0.58 (0.50-0.67)", "0.65 (0.56-0.76)", "0.64 (0.54-0.76)",
            "6.71 (2.83-15.91)", "4.7 to 14.5", "finasteride", "spironolactone",
            "comparator reports removed",
            "The comparator exclusion was added after the first complete analysis",
            "cohort data",
            "The female predominance in isotretinoin depression reporting is a property of the comparison, not of the drug.",
            "The age contrast belongs to the comparator",
            "0.55 to 0.64"))

  check_txt("body or legend contains", s)

# Table I is deposited, so its denominator is asserted against the supplement README
rtxt <- paste(readLines(file.path(ROOT, "supplement", "README.md"), warn = FALSE, encoding = "UTF-8"),
              collapse = " ")
for (s in c("13,829", "7,496", "6,333")) {
  pass <- grepl(s, rtxt, fixed = TRUE)
  if (pass) ok <- ok + 1L else bad <- bad + 1L
  cat(sprintf("[%s] %-54s supplement README contains \"%s\"\n",
              ifelse(pass, "PASS", "FAIL"), "deposited table denominator", s))
}

cat(sprintf("\n===== %d passed, %d failed =====\n", ok, bad))
if (bad > 0) quit(status = 1)
