# Isotretinoin and sex differences in depression reporting

Analysis code for a FAERS case/non-case study of sex differences in isotretinoin
depression reporting, 2012 Q4 to 2026 Q2.

The pipeline builds a dermatology-indication cohort from FAERS quarterly ASCII
files, estimates within-sex reporting odds ratios (ROR) for isotretinoin against
other dermatologic drugs, compares them as a ratio (RORR), and rebuilds the
comparator pool to exclude drugs used almost exclusively by one sex.

## Requirements

- R: data.table, ggplot2, openxlsx, officer, pdftools
- Python 3: pandas, python-docx

## Data

FAERS quarterly ASCII extracts, 55 quarters (2012 Q4 to 2026 Q2):

```
https://fis.fda.gov/content/Exports/faers_ascii_{YYYY}q{Q}.zip
```

Data and generated tables are not included in this repository.

## Main pipeline

```
python  code/01b_download_parallel.py
python  code/02_build_cohort.py
Rscript code/12_correct_2x2.R
Rscript code/17_final_analysis.R
Rscript code/18_robustness.R
Rscript code/22_time_stability.R
Rscript code/24_panel_a_data.R
Rscript code/25_monitoring_reports.R
Rscript code/26_active_comparator.R
Rscript code/29_leave_one_out.R
Rscript code/20_primary_sensitivity.R
Rscript code/19_supplement_update.R
Rscript code/04_figures.R
Rscript code/06_supplement.R
Rscript code/07_reporting_evidence.R
Rscript code/08_verify_numbers.R
Rscript code/99_check_pdf.R
```

The remaining scripts are exploratory scans and document builders. Output
tables, figures and manuscripts are not tracked here.
