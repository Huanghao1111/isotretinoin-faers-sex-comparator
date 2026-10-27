suppressPackageStartupMessages(library(pdftools))
f <- file.path("C:/Users/lenovo/WorkBuddy/2026-10-17-00-49-45",
               "output/figures/fig1_isotretinoin_rorr.pdf")
i <- pdf_info(f)
cat("pages:", i$pages, "| pdf version:", i$version, "\n")
for (nm in names(i$fonts)[1:min(4, length(names(i$fonts)))]) {
  cat("  font:", nm, "->", i$fonts[[nm]]$name,
      "| embedded:", isTRUE(i$fonts[[nm]]$embedded), "\n")
}
txt <- pdf_text(f)
cat("\n=== extracted text layer (first 700 chars) ===\n")
cat(substr(txt[1], 1, 700), "\n")
cat("\n=== total text chars:", nchar(txt[1]), "===\n")

# The figure must be editable vector art: every label has to come back as text.
must <- c("Depression", "Anxiety", "Suicidal ideation", "Completed suicide",
          "Major depression", "Anhedonia",
          "Comparator: all other dermatologic drugs", "Comparator: excluding sex-skewed drugs",
          "Ratio of ratios", "32",   # "32" proves the axis reaches the widest interval
          # panel A: the mechanism has to be in the figure, not only in the text
          "Finasteride", "Doxycycline", "302.9",
          "Depression reports per 1,000 male reports",
          "A", "B")                   # the two panel tags
ok <- 0L; bad <- 0L
for (s in must) {
  pass <- grepl(s, txt[1], fixed = TRUE)
  if (pass) ok <- ok + 1L else bad <- bad + 1L
  cat(sprintf("[%s] text layer contains \"%s\"\n",
              ifelse(pass, "PASS", "FAIL"), s))
}
# The figure is drawn from four tables. If any of them is newer than the PDF,
# the figure is showing numbers nobody can reproduce from the repository.
ROOT <- dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))))
inputs <- c("panel_a_drug_rates.csv", "panel_a_pooled.csv",
            "final_rorr_all_ages.csv", "final_rorr_untrimmed_all_ages.csv")
mt <- file.mtime(file.path(ROOT, "output/tables", inputs))
mt <- mt[!is.na(mt)]
pdf_mt <- file.mtime(f)
fresh <- length(mt) == 0 || pdf_mt > max(mt)
if (fresh) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] figure is newer than its four source tables (%s)\n",
            ifelse(fresh, "PASS", "FAIL - rerun 04_figures.R"),
            format(max(mt), "%H:%M:%S")))

# and the TIFF that is actually uploaded must be the same vintage
tif <- file.path(ROOT, "output/figures/fig1_isotretinoin_rorr.tif")
same <- file.exists(tif) && file.mtime(tif) >= pdf_mt - 5
if (same) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] uploaded TIFF matches the current PDF\n",
            ifelse(same, "PASS", "FAIL - rerun 04_figures.R")))

# the plotted values must exist as a file, or the figure cannot be audited
src <- file.path(ROOT, "output/tables/fig1_source_data.csv")
has_src <- file.exists(src) && nrow(read.csv(src)) >= 45
if (has_src) ok <- ok + 1L else bad <- bad + 1L
cat(sprintf("[%s] figure source data present (%s rows)\n",
            ifelse(has_src, "PASS", "FAIL - rerun 32"), 
            ifelse(has_src, nrow(read.csv(src)), 0)))

cat(sprintf("\n===== %d passed, %d failed =====\n", ok, bad))
if (bad > 0) quit(status = 1)
