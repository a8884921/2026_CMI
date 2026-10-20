.release_script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
.release_root <- normalizePath(file.path(.release_script_dir, "..", ".."))
source(file.path(.release_root, "privacy.R"), encoding = "UTF-8")
.d <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f)) dirname(normalizePath(f)) else getwd()})
d <- validate_anonymous_data(read.csv(file.path(.d, "inputs", "patient_five_score_medians_robust_minmax.csv")))
d$subtype <- factor(d$subtype, c("BL", "DLBCL", "FL", "MZL"))
fmt <- function(p) ifelse(p < 1e-300, "< 1.00e-300", sprintf("%.2e", p))
out <- do.call(rbind, lapply(c("CNV", "MTP", "Prol", "Diff", "Ig"), function(v) {
  med <- tapply(d[[v]], d$subtype, median)
  data.frame(dimension = v, BL_median = round(med[["BL"]], 3), DLBCL_median = round(med[["DLBCL"]], 3),
             FL_median = round(med[["FL"]], 3), MZL_median = round(med[["MZL"]], 3),
             KruskalWallis_P = fmt(kruskal.test(d[[v]] ~ d$subtype)$p.value),
             one_way_ANOVA_P = fmt(summary(aov(d[[v]] ~ d$subtype))[[1]][1, 5]))
}))
write.csv(out, file.path(.d, "results", "subtype_Kruskal_Wallis_sensitivity.csv"), row.names = FALSE, fileEncoding = "UTF-8")
print(out)
