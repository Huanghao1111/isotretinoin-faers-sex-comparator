#!/usr/bin/env Rscript
# ==============================================================================
# Moves the tables produced by the corrected design into supplement/ and writes
# the technical note that explains the two design faults.
#
# The README is appended with open("at") rather than read-modify-write: on
# Windows a file written in the same session is not reliably visible to a
# subsequent read, which silently truncated the README once already.
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "C:/Users/lenovo/WorkBuddy/2026-10-17-00-49-45"
TAB  <- file.path(ROOT, "output", "tables")
SUP  <- file.path(ROOT, "supplement")
dir.create(SUP, recursive = TRUE, showWarnings = FALSE)

copy_as <- function(from, to) {
  src <- file.path(TAB, from)
  if (!file.exists(src)) { cat("[skip] ", from, "\n", sep = ""); return(invisible(NULL)) }
  file.copy(src, file.path(SUP, to), overwrite = TRUE)
  cat("[supp] ", to, "\n", sep = "")
}

copy_as("final_rorr_all_ages.csv",             "TableS19_primary_rorr_all_ages.csv")
copy_as("final_rorr_age12_45.csv",             "TableS20_primary_rorr_age12_45.csv")
copy_as("final_rorr_untrimmed_all_ages.csv",   "TableS21_uncorrected_rorr_all_ages.csv")
copy_as("final_rorr_untrimmed_age12_45.csv",   "TableS22_uncorrected_rorr_age12_45.csv")
copy_as("robustness_trimming_rules.csv",       "TableS23_comparator_trimming_rules.csv")
copy_as("comparator_drug_sex_mix.csv",         "TableS24_comparator_drug_sex_mix.csv")
copy_as("finasteride_influence.csv",           "TableS25_finasteride_influence.csv")
copy_as("trimmed_arm_age_composition.csv",     "TableS26_arm_age_composition.csv")
copy_as("age_band_interaction.csv",            "TableS27_age_band_interaction.csv")
copy_as("pooled_known_age_rorr.csv",           "TableS28_pooled_known_age_rorr.csv")
copy_as("arm_age_composition.csv",             "TableS29_arm_age_composition_uncorrected.csv")
copy_as("primary_sensitivity_scenarios.csv",   "TableS30_primary_sensitivity_scenarios.csv")
copy_as("primary_age_band_interaction.csv",    "TableS31_primary_age_band_interaction.csv")
copy_as("primary_rorr_by_age_band.csv",        "TableS32_primary_rorr_by_age_band.csv")
copy_as("comparator_shift_by_panel.csv",       "TableS33_comparator_shift_by_panel.csv")
copy_as("comparator_pool_size.csv",            "TableS34_comparator_pool_size.csv")
copy_as("time_stability_rorr.csv",             "TableS35_time_stability_rorr.csv")
copy_as("time_stability_heterogeneity.csv",    "TableS36_time_stability_heterogeneity.csv")
copy_as("finasteride_share_by_era.csv",        "TableS37_finasteride_share_by_era.csv")
copy_as("panel_a_drug_rates.csv",            "TableS38_comparator_drug_depression_rates.csv")
copy_as("monitoring_report_counts.csv",       "TableS39_monitoring_report_counts.csv")
copy_as("monitoring_sensitivity.csv",         "TableS40_monitoring_sensitivity.csv")
copy_as("table1_isotretinoin.csv",             "TableS41_primary_estimates_all_terms.csv")
copy_as("active_comparator_acne.csv",        "TableS42_active_comparator_acne.csv")
copy_as("leave_one_drug_out.csv",           "TableS43_leave_one_drug_out.csv")

# --------------------------------------------------------------- the note
note <- file.path(SUP, "NOTE_comparator_design.md")
writeLines(c(
  "# Technical note: two design faults in the comparator pool",
  "",
  "Companion to TableS19-S43. Both faults were found after the first complete",
  "analysis, and both are fixed in the design rather than in the estimator.",
  "",
  "## Fault 1 - Simpson's paradox from age-structure mismatch",
  "",
  "The comparator pool (all other dermatologic drugs) has a completely",
  "different age structure from the index arm: 39.9% of comparator reports are",
  "aged 46+ against 2.9% of isotretinoin reports, while at ages 12-45 the shares",
  "are 25.4% and 63.7% (TableS29). The index-to-comparator share is therefore",
  "2.51 in the younger band and 0.07 in the older one, a 34-fold difference in",
  "the number of person-reports each index report is standing against. A pooled",
  "estimate built on that is not a weighted average of the stratum estimates and",
  "can point the other way. For depression it does: 1.28 at ages 12-45 and 1.51",
  "at 46+ (TableS27), and 0.82 when the two known-age strata are pooled (TableS28).",
  "Any",
  "pooled ratio of ratios from this cohort should be read as a case-mix average,",
  "not as a sex effect.",
  "",
  "## Fault 2 - sex-specific comparator drugs",
  "",
  "A comparator drug used almost exclusively by one sex cannot supply a baseline",
  "for both sexes. Finasteride is the clearest case: 2.1% of its reports overall",
  "are from women (1.3% at ages 12-45). At those ages it supplies 2,572 of 92,976",
  "male comparator reports but",
  "779 of 1,385 male comparator depression reports - 302.9 per 1,000, against 2.4",
  "per 1,000 for dupilumab (TableS25). Every comparator drug with at least 200",
  "male reports at these ages is on one axis in Figure 1A and in TableS38: the",
  "next highest rate is 48.4 per 1,000 (doxycycline) and the pooled rate across",
  "the whole comparator is 14.9, so finasteride sits six times above the next",
  "drug and twenty times above the pool. The male comparator depression rate is",
  "therefore set by one drug carrying its own psychiatric reporting signal, which",
  "deflates the male ROR and manufactures a female signal for isotretinoin.",
  "",
  "The distortion is easier to see at the level of the pool than of the drug.",
  "Finasteride holds 2.8% of the male comparator reports at these ages (2,572 of",
  "92,976) and supplies 56% of the male comparator depression reports (779 of",
  "1,385). Removing it takes the male comparator rate from 14.9 to 6.7 depression",
  "reports per 1,000: a drug carrying under 3% of the reports was",
  "responsible for more than half of the events.",
  "",
  "## The correction",
  "",
  "Comparator drugs with at least 200 reports whose reports are <20% or >80%",
  "from one sex are excluded: finasteride (2.1% female), spironolactone (92.4%),",
  "tretinoin (86.4%), hydroxychloroquine (87.1%) and omalizumab (80.2%)",
  "(TableS24). A comparator report is retained only if every primary-suspect drug",
  "on it is retained. The restriction removes 30,311 of 1,077,792 comparator",
  "reports - 2.8% - and leaves 1,047,481, or 97.2% of the pool, alongside 13,829",
  "isotretinoin reports (TableS34). Five drugs are named in the manuscript, so it",
  "is worth stating how little of the pool they carry: the correction is surgical",
  "in extent and decisive in effect. Removing finasteride alone gives 0.59, the",
  "same value as the full rule to two decimal places, so the correction is that",
  "one drug and not the rule (TableS23).",
  "",
  "## What the correction changes",
  "",
  "Depression moves from 0.97 (0.89-1.07) to 0.59 (0.54-0.65); anxiety from 1.09",
  "to 0.65; suicidal ideation from 1.30 to 0.64. The comparator factor for",
  "depression moves from 1.28 to 0.78, so for the first time the ratio of ratios",
  "and the drug's own odds ratio (0.76) point the same way. Completed suicide and",
  "the sexual-function terms keep comparator factors of 4.7 to 14.5 and remain",
  "uninterpretable regardless of trimming.",
  "",
  "## Stability over calendar time",
  "",
  "Calendar time is a stratification like age and it behaves the same way: the",
  "single table pooled over fifteen years is a case-mix average. Depression gives",
  "0.77 (0.67-0.88) in 2012-2016, 0.65 (0.54-0.79) in 2017-2021 and 0.48",
  "(0.35-0.66) in 2022-2026, heterogeneity p = 0.018 (TableS35, TableS36).",
  "",
  "At ages 12-45 the split is sharper and less comfortable: 0.87 (0.70-1.09), 0.90",
  "(0.68-1.19) and 0.37 (0.24-0.59), heterogeneity p = 0.002. Two of the three",
  "eras have intervals that include 1, so at the ages where acne is treated the",
  "male predominance rests on the most recent era. Pooling the eras by inverse",
  "variance rather than in one table gives 0.70 (0.63-0.77) for all ages and 0.79",
  "(0.67-0.93) at 12-45, against 0.59 and 0.58 from the single table: the odds",
  "ratio is not collapsible and the comparator's baseline moves between eras, so",
  "two defensible summaries of the same data differ by more than sampling error.",
  "Both are deposited. The manuscript quotes the single-table estimate because",
  "that is the estimator the paper is built on, and the era-stratified range is",
  "the honest statement of its uncertainty.",
  "",
  "Anxiety and suicidal ideation do not follow the same shape. At ages 12-45",
  "anxiety runs 0.81, 1.10 and 0.52 (p = 0.087) and suicidal ideation 0.67, 1.01",
  "and 1.05 (p = 0.317); the apparent stability of suicidal ideation is in the",
  "all-age estimate (0.62, 0.69, 0.70; p = 0.86) and not in the age-restricted",
  "one. The pooling caveat above applies to all three terms, not only to the one",
  "that moved.",
  "",
  "The drift is not specific to depression. Terms with no psychiatric mechanism",
  "move as much or more (lip dry p = 0.003 and dry skin p = 0.018 at ages 12-45),",
  "which points to the comparator rather than to the drug: the pool is not a",
  "fixed ruler as new drugs enter it.",
  "",
  "Nor is the fall in male reporting produced by missing values: reports without a",
  "recorded sex are excluded by design (TableS8, step 4), so they cannot drive the",
  "trend, though misclassification of recorded sex cannot be excluded. The male",
  "share of the 12-45 index arm falls from 53.0% to 45.6% to 29.1% across the eras. Whether that means",
  "fewer men treated, fewer men reported, or a change in which male cases get",
  "reported cannot be separated here. It is the main reason the era-specific",
  "numbers are described rather than interpreted, and the reason the paper claims",
  "a direction rather than a magnitude.",
  "",
  "Finasteride was not a late-window problem: at ages 12-45 it supplied 60.5%,",
  "49.1% and 59.7% of male comparator depression reports in the same three eras",
  "(TableS37), so the correction is needed throughout the window rather than",
  "being an artefact of one period's reporting fashion.",
  "",
  "## The pregnancy-prevention channel",
  "",
  "Isotretinoin is dispensed under a pregnancy-prevention programme and the",
  "reports show it: 8.2% of the female index arm (613 of 7,496 reports) carries",
  "only monitoring terms - pregnancy tests, contraception, lipids, liver function",
  "- against 1.9% of the male arm (118 of 6,333) (TableS39). Those reports add no",
  "depression cases and enlarge the female denominator, which lowers the female",
  "reporting odds and biases the ratio of ratios towards the male side. The bias",
  "runs in favour of this paper's conclusion, so it is removed rather than left",
  "standing: dropping every report whose only events are monitoring terms from",
  "both arms gives 0.64 (0.58-0.70) for depression at all ages and 0.62",
  "(0.54-0.72) at 12-45, against 0.59 and 0.58 (TableS40). Anxiety moves from",
  "0.65 to 0.69 (0.60-0.81) and suicidal ideation from 0.64 to 0.68 (0.57-0.81).",
  "The conclusion survives",
  "the correction it was most exposed to. Suicidal ideation at ages 12-45 does not",
  "(0.96 (0.74-1.23)), which agrees with the era analysis above and is a further",
  "reason the paper does not present that term as a finding.",
  "",
  "## No single comparator drug carries the corrected estimate",
  "",
  "Removing finasteride fixed one problem; the symmetric question is whether the",
  "corrected estimate is itself the product of one comparator drug. Each of the 64",
  "retained comparator drugs was removed in turn and the 2x2 refitted (TableS43).",
  "Depression ranges from 0.547 to 0.602 across the 64 refits, median 0.589,",
  "against 0.588 with the full pool: the largest single move is 0.04, from",
  "removing apremilast, and no refit comes near 1. Anxiety stays within",
  "0.614-0.667 and suicidal ideation within 0.605-0.669. The correction moved the",
  "estimate without making it fragile, and whatever produces the male predominance",
  "is spread across the pool rather than concentrated in one drug the way",
  "finasteride was.",
  "",
  "## An indication-matched active comparator",
  "",
  "Every estimate above compares isotretinoin with drugs given for other",
  "diseases. The tetracyclines are the other systemic treatment for acne, so",
  "restricting both arms to reports carrying an acne indication matches the",
  "indication by construction rather than by argument: 7,286 female and 6,185 male",
  "isotretinoin reports against 1,571 and 829 tetracycline reports. Depression",
  "gives 0.88 (0.52-1.49) and anxiety 0.52 (0.30-0.92) (TableS42). The direction",
  "is the same as the primary analysis, the anxiety interval excludes 1, and the",
  "depression interval does not, because the active arm is small. The comparison",
  "is deposited as compatible with the primary result, not as confirmation of it.",
  "Its value is that it is the one design here in which the comparator's",
  "indication cannot explain the sex difference.",
  "",
  "## Scale of the age contrast",
  "",
  "The 12-45 versus 46+ contrast (0.41, 0.23-0.76; p = 0.004; BH q = 0.016) is an",
  "interaction on the multiplicative scale, computed on the same post-hoc",
  "comparator as every other estimate here. On an additive scale the same data",
  "would not give the same number, and no additive estimate is offered because the",
  "reported quantity is a ratio of ratios and has no additive analogue that a",
  "reader could act on.",
  "",
  "The third of the index arm with no recorded age behaves like the younger band",
  "and not like the older one: depression gives 0.65 (0.55-0.76) there, against",
  "0.58 (0.50-0.67) at 12-45 and 1.40 (0.78-2.51) above 45 (TableS32). The pooled",
  "estimate is therefore not being carried by an unrecorded stratum, which is the",
  "failure mode that usually undoes an age-stratified analysis of this kind. Note",
  "that the 46+ interval includes 1: the difference between the bands rests on the",
  "interaction test (0.41, 0.23-0.76), not on the 46+ estimate alone.",
  "",

  "",
  "The trimming rule was written after the first analysis and was not",
  "pre-specified, so the intervals above are conditional on a comparator chosen",
  "from the same data: they describe the precision of each estimate given the",
  "design that produced it, and they do not carry the uncertainty of the",
  "selection step. No post-selection interval is offered, because none of the",
  "available options - bootstrapping the rule, or widening by the range of",
  "estimates across rules - has a defensible frequency interpretation here. The",
  "intervals are therefore reported as conditional and descriptive, and the",
  "argument rests on the direction being stable across every restriction",
  "(TableS30) rather than on the width of any single interval.",
  "",
  "## What is not claimed",
  "",
  "Trimming the pool does not turn disproportionality into risk. It removes one",
  "identifiable source of bias; the sex composition of treated patients remains",
  "unknown, and men may receive isotretinoin for more severe disease."
), note)
cat("[supp] NOTE_comparator_design.md\n")

# ------------------------------------------------------------- README check
# The README is authored directly rather than assembled by appending. Both
# earlier approaches failed: read-modify-write silently truncated it on Windows
# (a file written in the same session is not reliably visible to a later read),
# and appending left the table entries out of order. This script therefore only
# verifies that every table it ships is described in the README.
rd <- file.path(SUP, "README.md")
if (!file.exists(rd)) {
  cat("[warn] README.md missing\n")
} else {
  txt <- readLines(rd, warn = FALSE)
  shipped <- c("TableS19", "TableS20", "TableS21", "TableS22", "TableS23",
               "TableS24", "TableS25", "TableS26", "TableS27", "TableS28",
               "TableS29", "TableS30", "TableS31", "TableS32", "TableS33",
               "TableS34", "TableS35", "TableS36", "TableS37", "TableS38", "TableS39", "TableS40", "TableS41", "TableS42", "TableS43")
  miss <- shipped[!vapply(shipped, function(x) any(grepl(x, txt, fixed = TRUE)), logical(1))]
  if (length(miss)) {
    cat("[warn] README.md lacks entries for: ", paste(miss, collapse = ", "), "\n", sep = "")
  } else {
    cat("[ok] README.md describes all of S19-S43\n")
  }
  if (anyDuplicated(txt[txt != ""])) {
    cat("[warn] README.md contains duplicated non-blank lines\n")
  } else {
    cat("[ok] README.md has no duplicated lines\n")
  }
}
cat("[done]\n")
