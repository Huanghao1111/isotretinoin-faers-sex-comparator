#!/usr/bin/env Rscript
# ==============================================================================
# The six supplemental tables a reader will actually open.
#
# Forty-three tables were built, and all forty-three stay in the repository, but a
# deposit is not an archive. These six each carry a claim the manuscript makes:
# which drugs left the pool and why, what the corrected estimate is, how it moves
# with the exclusion threshold, how it moves across eras, what the monitoring
# channel does, and what an indication-matched comparator gives.
#
# Everything else - the flow counts, the exposure vocabulary, the per-drug rates,
# the leave-one-out refits - is process documentation, and process documentation
# belongs in the repository, where the code that produced it lives.
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})

ROOT <- "C:/Users/lenovo/WorkBuddy/2026-10-17-00-49-45"
SUP  <- file.path(ROOT, "supplement")
OUT  <- file.path(SUP, "JAAD_supplemental_tables.xlsx")

PICK <- data.table(
  sheet = c("S1", "S2", "S3", "S4", "S5", "S6"),
  file  = c("TableS24_comparator_drug_sex_mix.csv",
            "TableS41_primary_estimates_all_terms.csv",
            "TableS23_comparator_trimming_rules.csv",
            "TableS35_time_stability_rorr.csv",
            "TableS40_monitoring_sensitivity.csv",
            "TableS42_active_comparator_acne.csv"),
  title = c("Comparator drugs by sex composition, before and after trimming",
            "Ratios of ratios for 15 events, both comparators, with the decomposition",
            "The estimate under five comparator exclusion rules",
            "The estimate by calendar era, all ages and ages 12-45",
            "The estimate with pregnancy-programme monitoring reports removed",
            "Isotretinoin versus tetracyclines, both arms restricted to acne"),
  shows = c("Which drugs left the pool and why; finasteride supplies 2.1% of its reports from women",
            "The primary result: depression 0.59 (0.54-0.65); negative controls barely move",
            "0.59 at the 20-80% rule, 0.56 at 30-70%, 0.46 at 40-60%: direction holds, magnitude does not",
            "0.77 / 0.65 / 0.48 by era all ages; 0.87 / 0.90 / 0.37 at 12-45",
            "Depression 0.64 (0.58-0.70) after removing monitoring-only reports",
            "Depression 0.88 (0.52-1.49), anxiety 0.52 (0.30-0.92) on an acne-matched comparison")
)

for (i in seq_len(nrow(PICK))) {
  f <- file.path(SUP, PICK$file[i])
  if (!file.exists(f)) stop("missing source table: ", PICK$file[i])
}

wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", border = "bottom")
hdr <- createStyle(textDecoration = "bold", fontSize = 12)

## ------------------------------------------------------------------ index
addWorksheet(wb, "Contents")
writeData(wb, "Contents", "Supplemental tables", colNames = FALSE, startRow = 1)
addStyle(wb, "Contents", hdr, rows = 1, cols = 1)
writeData(wb, "Contents", data.table(
  Sheet = PICK$sheet,
  Table = PICK$title,
  `What it shows` = PICK$shows,
  `Source file in the code repository` = PICK$file
), startRow = 3, headerStyle = hs)
setColWidths(wb, "Contents", 1:4, c(8, 52, 78, 48))
writeData(wb, "Contents", paste(
  "Counts are at the report level (distinct primaryid), never at the report-by-drug level.",
  "A blank cell means the quantity is not estimable for that row, usually because one arm has no case.",
  "The full set of 43 tables, including cohort flow, exposure vocabulary, per-drug rates and the",
  "leave-one-out refits, is in the code repository that accompanies the manuscript."),
  startRow = 12, colNames = FALSE)

## ------------------------------------------------------------------ tables
for (i in seq_len(nrow(PICK))) {
  d <- fread(file.path(SUP, PICK$file[i]))
  addWorksheet(wb, PICK$sheet[i])
  writeData(wb, PICK$sheet[i], PICK$title[i], colNames = FALSE, startRow = 1)
  addStyle(wb, PICK$sheet[i], hdr, rows = 1, cols = 1)
  writeData(wb, PICK$sheet[i], d, startRow = 3, headerStyle = hs)
  setColWidths(wb, PICK$sheet[i], seq_len(ncol(d)), pmax(10, pmin(30, nchar(names(d)) + 4)))
  freezePane(wb, PICK$sheet[i], firstActiveRow = 4)
}

saveWorkbook(wb, OUT, overwrite = TRUE)
cat(sprintf("[written] %s (%.2f MB, %d sheets)\n", basename(OUT),
            file.size(OUT) / 1e6, nrow(PICK) + 1))
for (i in seq_len(nrow(PICK))) cat(sprintf("  %s <- %s\n", PICK$sheet[i], PICK$file[i]))
