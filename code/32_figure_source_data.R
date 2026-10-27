#!/usr/bin/env Rscript
# ==============================================================================
# The numbers behind Figure 1, in a file.
#
# A reader can see 302.9 on the page but cannot check it without the tables the
# figure was drawn from. This writes exactly the values that go into the two
# panels, one row per plotted mark, so the figure can be reconstructed and
# audited without running ggplot.
# ==============================================================================
suppressPackageStartupMessages(library(data.table))

ROOT <- "C:/Users/lenovo/WorkBuddy/2026-10-17-00-49-45"
TAB  <- file.path(ROOT, "output", "tables")

pa  <- fread(file.path(TAB, "panel_a_drug_rates.csv"))
pool <- fread(file.path(TAB, "panel_a_pooled.csv"))

# panel A: every drug that carries at least 200 male reports at ages 12-45
a <- pa[, .(panel = "A",
            label = drug_label,
            value = round(rate_M, 1),
            n_male_reports = n_M,
            n_depression = dep_M,
            drug_group)]
a <- a[order(-value)]
a[, plotted := seq_len(.N) <= 15]          # the figure draws the 15 highest
a[, note := fifelse(drug_group == "finasteride", "highlighted",
             fifelse(plotted, "plotted", "not plotted"))]
a[, drug_group := NULL]
# the legend promises the unplotted drugs are all below 10 per 1,000
if (any(a[plotted == FALSE, value] >= 10)) stop("legend claim fails: an unplotted drug reaches 10")

# panel B: both comparators for the six terms drawn
pri <- fread(file.path(TAB, "final_rorr_all_ages.csv"))
unc <- fread(file.path(TAB, "final_rorr_untrimmed_all_ages.csv"))
terms <- c("Depression", "Anxiety", "Suicidal ideation",
           "Completed suicide", "Major depression", "Anhedonia")

b1 <- pri[pt %in% terms, .(panel = "B", label = pt,
                           comparator = "sex-neutral (primary)",
                           value = RORR, lo = lci, hi = uci)]
b2 <- unc[pt %in% terms, .(panel = "B", label = pt,
                           comparator = "all other dermatologic drugs",
                           value = RORR, lo = lci, hi = uci)]
b <- rbind(b1, b2, fill = TRUE)

out <- rbind(a, b, fill = TRUE)
setcolorder(out, c("panel", "label", "comparator", "value", "lo", "hi",
                   "n_male_reports", "n_depression", "plotted", "note"))
fwrite(out, file.path(TAB, "fig1_source_data.csv"))

cat(sprintf("[written] fig1_source_data.csv: %d rows (%d panel A, %d panel B)\n",
            nrow(out), nrow(a), nrow(b)))
cat("  pooled table columns:", paste(names(pool), collapse = ", "), "\n")
