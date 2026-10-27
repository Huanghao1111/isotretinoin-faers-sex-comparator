suppressPackageStartupMessages({library(data.table)})
ROOT <- "C:/Users/lenovo/WorkBuddy/2026-10-17-00-49-45"
A <- readRDS(file.path(ROOT, "data", "processed", "analysis_objects.rds"))
R <- as.data.table(A$R)
r <- R[drug_group == "isotretinoin" & pt == "Depression"]
print(r[, .(pt, a_F, b_F, c_F, d_F, ror_F, a_M, b_M, c_M, d_M, ror_M, rorr)])
cat("\n-- reconstruct the 2x2 from the saved cells --\n")
cat(sprintf("F: a=%d b=%d c=%d d=%d  N=%d  Nd=a+b=%d  Ne=a+c=%d\n",
            r$a_F, r$b_F, r$c_F, r$d_F, sum(unlist(r[,.(a_F,b_F,c_F,d_F)])),
            r$a_F + r$b_F, r$a_F + r$c_F))
cat(sprintf("M: a=%d b=%d c=%d d=%d  N=%d  Nd=a+b=%d  Ne=a+c=%d\n",
            r$a_M, r$b_M, r$c_M, r$d_M, sum(unlist(r[,.(a_M,b_M,c_M,d_M)])),
            r$a_M + r$b_M, r$a_M + r$c_M))
cat("\n-- cohort totals for comparison --\n")
dv <- as.data.table(A$demo)
print(dv[, .(N = .N), by = sex])
