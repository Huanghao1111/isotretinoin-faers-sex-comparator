#!/usr/bin/env Rscript
# ==============================================================================
# Submission-format Word file: title page, double-spaced body with continuous
# line numbers, references, figure legend, three-line Table I.
#
# The text is read from the markdown rather than retyped, so the docx cannot
# drift from the file that 08_verify_numbers.R checks. Line numbering is not
# exposed by officer, so it is injected into the section properties afterwards by
# 23b_line_numbers.py: journals ask for it and a manuscript without it looks
# unfinished.
# ==============================================================================
suppressPackageStartupMessages({
  library(officer)
  library(flextable)
  library(data.table)
})

ROOT <- "C:/Users/lenovo/WorkBuddy/2026-10-17-00-49-45"
MS   <- file.path(ROOT, "manuscript", "JAAD_BriefReport_isotretinoin_sex.md")
TAB  <- file.path(ROOT, "output", "tables")
OUT  <- file.path(ROOT, "manuscript", "JAAD_BriefReport_isotretinoin_sex.docx")

# encoding has to be declared: on a Chinese Windows locale readLines otherwise
# returns native-encoded bytes, and an em dash then truncates its paragraph
lines <- readLines(MS, warn = FALSE, encoding = "UTF-8")

body_start <- grep("^\\*To the Editor:\\*$", lines)
body_end   <- grep("^\\*\\*Word count", lines)
if (length(body_start) != 1L || length(body_end) != 1L) stop("cannot locate the body")
body <- lines[(body_start + 1L):(body_end - 1L)]
body <- paste(body, collapse = "\n")
paras <- trimws(strsplit(body, "\n\\s*\n")[[1]])
paras <- paras[nzchar(paras)]
# markdown emphasis is a typing convention, not part of the text
paras <- gsub("\\*\\*", "", paras)
cat("body paragraphs:", length(paras), "\n")

refs <- grep("^\\d+\\. ", lines, value = TRUE)
cat("references:", length(refs), "\n")

legend_line <- grep("^\\*\\*Figure 1\\.\\*\\*", lines, value = TRUE)
legend <- gsub("^\\*\\*|\\*\\*$", "", legend_line)
legend <- gsub("\\*\\*", "", legend)

# the 15-term table is deposited as TableS41; this file carries no table

FONT <- "Times New Roman"
base_txt <- fp_text(font.family = FONT, font.size = 12)
bold_txt <- fp_text(font.family = FONT, font.size = 12, bold = TRUE)
ital_txt <- fp_text(font.family = FONT, font.size = 12, italic = TRUE)
small_txt <- fp_text(font.family = FONT, font.size = 10)
small_ital <- fp_text(font.family = FONT, font.size = 10, italic = TRUE)

dbl <- fp_par(line_spacing = 2, text.align = "left", padding.bottom = 0)
sgl <- fp_par(line_spacing = 1, text.align = "left", padding.bottom = 0)

add_par <- function(x, text, props = dbl, fmt = base_txt) {
  x <- body_add_fpar(x, fpar(ftext(text, fmt), fp_p = props))
  x
}

doc <- read_docx()

# ------------------------------------------------------------------ title page
doc <- add_par(doc, "JAAD Brief Reports", sgl, small_ital)
doc <- add_par(doc, "", sgl)
doc <- add_par(doc, sub("^# ", "", gsub("\\*\\*", "", lines[grep("^# ", lines)[1]])), sgl, bold_txt)
doc <- add_par(doc, "", sgl)
doc <- add_par(doc, "Short title: Rebuilding the comparator for isotretinoin depression reporting", sgl, ital_txt)
doc <- add_par(doc, "", sgl)
doc <- add_par(doc, "Haozhang Huang\u00b9, Jiajia Li\u00b2, Hao Xiong\u00b2", dbl)
doc <- add_par(doc, "\u00b9 Guangdong Provincial Geriatrics Institute, Guangdong Provincial People's Hospital (Guangdong Academy of Medical Sciences), Southern Medical University, Guangzhou 510080, China", dbl)
doc <- add_par(doc, "\u00b2 Department of Dermatology, Guangdong Provincial People's Hospital (Guangdong Academy of Medical Sciences), Southern Medical University, Guangzhou 510080, China", dbl)
doc <- add_par(doc, "", sgl)
doc <- add_par(doc, "Corresponding author: Haozhang Huang, Guangdong Provincial Geriatrics Institute, Guangdong Provincial People's Hospital (Guangdong Academy of Medical Sciences), Southern Medical University, Guangzhou 510080, China. Email: 1017988724@qq.com", dbl)
doc <- add_par(doc, "", sgl)

meta <- c(
  "Word count, body: 496 (limit 500)",
  "References: 5 (limit 5)",
  "Displays: one figure (two panels); the 15-event table is sheet S2 of the supplemental workbook",
  "Abstract: none (Brief Reports do not carry one)",
  "Keywords: isotretinoin; pharmacovigilance; disproportionality analysis; sex differences; FAERS",
  "",
  "Funding: none.",
  "Competing interests: none declared.",
  "IRB approval: not required; de-identified public data (FAERS).",
  "Patient consent: not applicable; no patients were involved.",
  "Data availability: the data are publicly available from the FDA FAERS quarterly ASCII extracts (https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html), 2012 Q4 to 2026 Q2; supplemental material (six tables) is submitted as one workbook and will be deposited at Mendeley Data; the full set of 43 tables and the analysis code are in the accompanying repository.",
  "Author contributions: all authors meet the ICMJE criteria for authorship - each contributed to the design or to the analysis and interpretation of data, drafted or critically revised the manuscript, approved the version submitted, and agrees to be accountable for the work.",
  "Acknowledgements: none.",
  "Use of artificial intelligence: no AI tool was used to generate data or conclusions. An AI-assisted writing tool was used for language editing of the text; the authors reviewed and take responsibility for all content."
)
for (m in meta) doc <- add_par(doc, m, sgl, if (nzchar(m)) base_txt else base_txt)

doc <- body_add_break(doc)

# ------------------------------------------------------------------------ body
# read the key message from the manuscript so the two cannot diverge again
km <- grep("Key Message:", lines, value = TRUE)[1]
km <- gsub("**", "", km, fixed = TRUE)
km <- sub("Key Message: ", "", km, fixed = TRUE)
doc <- add_par(doc, paste0("Key Message: ", km), dbl, bold_txt)
doc <- add_par(doc, "", sgl)
doc <- add_par(doc, "To the Editor:", dbl, ital_txt)
doc <- add_par(doc, "", sgl)
for (p in paras) doc <- add_par(doc, p, dbl)
doc <- add_par(doc, "", sgl)

# ------------------------------------------------------------------ references
doc <- add_par(doc, "References", dbl, bold_txt)
for (r in refs) doc <- add_par(doc, gsub("\\*", "", r), dbl)

doc <- add_par(doc, "", sgl)
doc <- add_par(doc, "Figure legend", dbl, bold_txt)
doc <- add_par(doc, legend, dbl)

print(doc, target = OUT)
cat("[written]", OUT, "\n[done]\n")
