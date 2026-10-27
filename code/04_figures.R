#!/usr/bin/env Rscript
# ==============================================================================
# Figure 1 : neuropsychiatric ratios of ratios under the two comparator
#            definitions, side by side. The point of the figure is that the
#            answer flips when a comparator drug used almost only by men is
#            removed, so both are plotted rather than the primary alone.
# Table I  : all fifteen pre-specified terms, both comparators, plus the two
#            factors the ratio of ratios decomposes into.
# Graphics are editable vector PDF (text stays text, never rasterised).
#
# Inputs are produced by 17_final_analysis.R, which must run first.
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

ROOT <- "C:/Users/lenovo/WorkBuddy/2026-10-17-00-49-45"
PROC <- file.path(ROOT, "data", "processed")
TAB  <- file.path(ROOT, "output", "tables")
FIG  <- file.path(ROOT, "output", "figures")
for (d in c(TAB, FIG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

NEURO <- c("Depression", "Anxiety", "Suicidal ideation",
           "Completed suicide", "Major depression", "Anhedonia")
SEXUAL <- c("Loss of libido", "Sexual dysfunction", "Libido decreased")
SOMATIC <- c("Lip dry", "Crohn's disease", "Inflammatory bowel disease",
             "Dry skin", "Epistaxis", "Acne")

pri <- fread(file.path(TAB, "final_rorr_all_ages.csv"))
unc <- fread(file.path(TAB, "final_rorr_untrimmed_all_ages.csv"))
stopifnot(nrow(pri) > 0, nrow(unc) > 0)

T <- merge(pri[, .(pt, a_F, a_M, ROR_F, ROR_M,
                   rorr_pri = RORR, lci_pri = lci, uci_pri = uci, p_pri = p,
                   index_OR_F_vs_M, comp_OR_M_vs_F)],
           unc[, .(pt, rorr_all = RORR, lci_all = lci, uci_all = uci)],
           by = "pt")
T <- T[pt %in% c(NEURO, SEXUAL, SOMATIC)]
T[, grp := fifelse(pt %in% NEURO, "Neuropsychiatric",
            fifelse(pt %in% SEXUAL, "Sexual function",
                    "Mucocutaneous / gastrointestinal"))]
T[, grp := factor(grp, levels = c("Neuropsychiatric", "Sexual function",
                                  "Mucocutaneous / gastrointestinal"))]

# Index-drug denominator for the table title. Taken from the primary output so
# it is the same population the ratios were computed in: reports with >=1
# primary-suspect dermatologic drug AND >=1 retained adverse-event term.
N_F <- pri[, max(n_iF)]; N_M <- pri[, max(n_iM)]; N_ISO <- N_F + N_M
n_iso <- data.table(sex = c("F", "M"), N = c(N_F, N_M))
cat(sprintf("isotretinoin PS reports: %s total (%s F, %s M)\n",
            format(N_ISO, big.mark = ","), format(N_F, big.mark = ","),
            format(N_M, big.mark = ",")))
fwrite(n_iso, file.path(TAB, "table1_denominator.csv"))

# ------------------------------------------------------------ Figure 1 data
# Panel A answers "why is the ratio of ratios wrong?" and Panel B answers "what
# is the ratio of ratios?". The paper's claim is the first, so it gets a panel
# rather than a sentence: the male comparator pool's depression rate is set by
# one drug, and every other comparator is an order of magnitude below it.
TOP <- 15L
pa <- fread(file.path(TAB, "panel_a_drug_rates.csv"))
poolM <- fread(file.path(TAB, "panel_a_pooled.csv"))[sex == "M", rate_per_1000]
stopifnot(nrow(pa) >= TOP)
pa <- pa[order(-rate_M)][seq_len(TOP)]
pa[, lab := factor(drug_label, levels = rev(drug_label))]
pa[, fin := drug_label == "Finasteride"]
cat("panel A drugs:", nrow(pa), " pooled male rate:", poolM, "\n")

pA <- ggplot(pa, aes(y = lab)) +
  geom_vline(xintercept = poolM, linetype = "dashed", colour = "grey35",
             linewidth = 0.4) +
  geom_point(aes(x = rate_M, colour = fin, size = n_M)) +
  geom_text(data = pa[fin == TRUE], aes(x = rate_M, label = sprintf("%.1f", rate_M)),
            hjust = -0.45, size = 2.4, colour = "#1A1A1A") +
  scale_x_log10(breaks = c(1, 3, 10, 30, 100, 300),
                labels = c("1", "3", "10", "30", "100", "300"),
                limits = c(0.55, 900)) +
  scale_colour_manual(values = c(`FALSE` = "grey35", `TRUE` = "#1A1A1A"),
                      guide = "none") +
  scale_size_area(max_size = 2.8, guide = "none") +
  annotation_logticks(sides = "b", colour = "grey70", linewidth = 0.25) +
  labs(x = "Depression reports per 1,000 male reports", y = NULL) +
  theme_classic(base_size = 9) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        axis.text        = element_text(colour = "black", size = 7.5),
        axis.title.x     = element_text(colour = "black", size = 8.5),
        plot.margin      = margin(10, 14, 6, 6))

fig <- T[pt %in% NEURO]
setorder(fig, -rorr_pri)
fig[, lab := factor(pt, levels = rev(pt))]
cat("events in figure:", nrow(fig), "\n")

fd <- rbind(
  fig[, .(lab, est = "Comparator: all other dermatologic drugs",
          x = rorr_all, lo = lci_all, hi = uci_all)],
  fig[, .(lab, est = "Comparator: excluding sex-skewed drugs",
          x = rorr_pri, lo = lci_pri, hi = uci_pri)]
)
fd[, est := factor(est, levels = c("Comparator: all other dermatologic drugs",
                                   "Comparator: excluding sex-skewed drugs"))]
fd[, lab := factor(lab, levels = levels(fig$lab))]

pd <- position_dodge(width = 0.55)
pB <- ggplot(fd, aes(x = x, y = lab, colour = est, group = est)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45", linewidth = 0.5) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.12, linewidth = 0.45,
                alpha = 0.85, position = pd, orientation = "y") +
  geom_point(size = 2.1, position = pd) +
  # range has to cover the widest interval in the panel: anhedonia reaches 23.9
  # on the all-drug comparator, which a 0.1-8 axis would silently cut off
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 8, 32),
                labels = c("0.1", "0.25", "0.5", "1", "2", "8", "32"),
                limits = c(0.08, 40)) +
  scale_colour_manual(values = c(
    "Comparator: all other dermatologic drugs" = "#9AA0A6",
    "Comparator: excluding sex-skewed drugs" = "#B3222B")) +
  annotation_logticks(sides = "b", colour = "grey70", linewidth = 0.25) +
  labs(x = "Ratio of ratios (female / male), logarithmic scale",
       y = NULL, colour = NULL) +
  theme_classic(base_size = 9) +
  theme(panel.background  = element_rect(fill = "white", colour = NA),
        plot.background   = element_rect(fill = "white", colour = NA),
        axis.text         = element_text(colour = "black", size = 8),
        axis.title        = element_text(colour = "black", size = 9),
        legend.position   = "top",
        legend.direction  = "vertical",
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key        = element_rect(fill = NA, colour = NA),
        legend.text       = element_text(size = 8),
        plot.margin       = margin(10, 12, 6, 6)) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))

suppressPackageStartupMessages(library(patchwork))
p <- pA + pB +
  plot_layout(widths = c(1, 1.25)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

ggsave(file.path(FIG, "fig1_isotretinoin_rorr.pdf"), p,
       device = cairo_pdf, width = 7.4, height = 4.7, units = "in")

# JAAD takes figures as .tif/.eps/.jpg, never .pdf. 600 dpi is enough for this
# line-art composite at the submitted width.
if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(file.path(FIG, "fig1_isotretinoin_rorr.tif"), p,
         device = ragg::agg_tiff, width = 7.4, height = 4.7, units = "in",
         dpi = 600, compression = "lzw", bg = "white")
}
cat("[fig] fig1_isotretinoin_rorr.pdf (two panels)\n")

# -------------------------------------------------- Table I (three-line, SCI)
fmt_ci <- function(x, l, u) sprintf("%.2f (%.2f–%.2f)", x, l, u)
pfmt <- function(p) ifelse(p < 1e-4, "<0.0001", sprintf("%.4f", p))

t1 <- copy(T)[order(grp, -rorr_pri)]
t1 <- t1[, .(Event = pt,
             Group = grp,
             `Cases, female` = a_F,
             `Cases, male`   = a_M,
             `ROR ratio, all-drug comparator (95% CI)` = fmt_ci(rorr_all, lci_all, uci_all),
             `ROR ratio, sex-neutral comparator (95% CI)`     = fmt_ci(rorr_pri, lci_pri, uci_pri),
             `P`                               = pfmt(p_pri),
             `Comparator OR, men vs women`     = comp_OR_M_vs_F,
             `Isotretinoin OR, women vs men`   = index_OR_F_vs_M)]
fwrite(t1, file.path(TAB, "table1_isotretinoin.csv"))
cat("[tab] table1_isotretinoin.csv  (", nrow(t1), "rows )\n")

if (requireNamespace("flextable", quietly = TRUE) &&
    requireNamespace("officer", quietly = TRUE)) {
  ft <- flextable::flextable(t1)
  ft <- flextable::border_remove(ft)
  bd <- officer::fp_border(color = "black", width = 1.2)
  bd2 <- officer::fp_border(color = "black", width = 0.6)
  ft <- flextable::hline_top(ft, border = bd)
  ft <- flextable::hline(ft, i = 1, part = "header", border = bd2)
  ft <- flextable::hline_bottom(ft, border = bd)
  ft <- flextable::autofit(ft)
  flextable::save_as_docx(ft, path = file.path(TAB, "table1_isotretinoin.docx"))
  cat("[tab] table1_isotretinoin.docx (three-line)\n")
}

print(t1[, .(Event, `Cases, female`, `Cases, male`,
             `ROR ratio, all-drug comparator (95% CI)`, `ROR ratio, sex-neutral comparator (95% CI)`, P,
             `Comparator OR, men vs women`)], nrows = 30)

cat("\n[done]\n")
