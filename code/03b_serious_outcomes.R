#!/usr/bin/env Rscript
# ==============================================================================
# Serious outcomes: sex contrast with drug as stratification variable.
# The naive sex comparison is confounded because men receive intrinsically
# more toxic drugs (checkpoint inhibitors, targeted oncology agents). Each drug
# contributes one stratum and the stratum-specific ORs are pooled with a
# random-effects model, which yields a valid CI.
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(metafor)
})

ROOT <- "C:/Users/lenovo/WorkBuddy/2026-10-17-00-49-45"
PROC <- file.path(ROOT, "data", "processed")
TAB  <- file.path(ROOT, "output", "tables")

A <- readRDS(file.path(PROC, "analysis_objects.rds"))
ser <- as.data.table(A$ser)

MIN_STRATUM_N <- 40   # need enough cases per drug for a usable stratum

# NOTE: iterate over the column codes actually present in `ser`
out <- rbindlist(lapply(c("DE", "HO", "LT", "DS"),
  function(ev) {
    d <- ser[, .(m_n = sum(get(ev) > 0 & sex == "M"),
                 m_t = sum(sex == "M"),
                 f_n = sum(get(ev) > 0 & sex == "F"),
                 f_t = sum(sex == "F")), by = drug_group]
    d <- d[m_t >= MIN_STRATUM_N & f_t >= MIN_STRATUM_N]
    # continuity correction only where a cell is empty
    cc <- 0.5 * (d$m_n == 0 | d$f_n == 0 |
                 d$m_n == d$m_t | d$f_n == d$f_t)
    d[, `:=`(a = m_n + cc, b = m_t - m_n + cc,
             c = f_n + cc, dd = f_t - f_n + cc)]
    d <- d[a > 0 & b > 0 & c > 0 & dd > 0]
    if (nrow(d) < 3) return(NULL)
    d[, `:=`(yi = log((a * dd) / (b * c)),
             vi = 1/a + 1/b + 1/c + 1/dd)]
    m <- rma(yi, vi, data = d, method = "REML", slab = d$drug_group)
    data.table(outcome = ev, strata = nrow(d),
               OR_M_vs_F = round(exp(m$b[1]), 3),
               lci = round(exp(m$ci.lb), 3),
               uci = round(exp(m$ci.ub), 3),
               p = signif(m$pval, 3),
               I2 = round(m$I2, 1), tau2 = round(m$tau2, 4))
  }))

cat("\n== serious outcomes, drug-stratified random-effects pooling ==\n")
cat("reference = female; OR > 1 means more frequent in men\n\n")
print(out)

nms <- c(DE = "Death", HO = "Hospitalisation",
         LT = "Life-threatening", DS = "Disability")
out[, outcome := nms[outcome]]
fwrite(out, file.path(TAB, "serious_outcome_MH_drugstratified.csv"))

# crude, unadjusted contrast for comparison
crude <- ser[, .(n = .N,
                 Death = round(100 * sum(DE > 0) / .N, 1),
                 Hospitalisation = round(100 * sum(HO > 0) / .N, 1),
                 `Life-threatening` = round(100 * sum(LT > 0) / .N, 1),
                 Disability = round(100 * sum(DS > 0) / .N, 1)), by = sex]
cat("\n== crude rates (%) ==\n"); print(crude)
fwrite(crude, file.path(TAB, "serious_outcome_crude.csv"))
cat("\n[done]\n")
