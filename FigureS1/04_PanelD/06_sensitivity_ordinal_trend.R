.release_script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
.release_root <- normalizePath(file.path(.release_script_dir, "..", ".."))
source(file.path(.release_root, "privacy.R"), encoding = "UTF-8")
.d <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f)) dirname(normalizePath(f)) else getwd()})
suppressMessages(library(lme4))
P <- validate_anonymous_data(read.csv(file.path(.d, "inputs", "patient_five_score_medians_robust_minmax.csv")))
b <- validate_anonymous_data(read.csv(file.path(.d, "inputs", "data", "five_score_matrix.csv")))
stopifnot(all(c("sample_id", "patient_id", "malignancy_label") %in% names(b)))
b <- b[b$malignancy_label == "tumor", ]; b$pt <- b$patient_id
jt <- function(x, g) { s <- 0; lv <- sort(unique(g))
  for (i in lv) for (j in lv) if (i < j) { a <- x[g == i]; c <- x[g == j]; s <- s + sum(outer(a, c, "<")) + 0.5 * sum(outer(a, c, "==")) }
  s }
jt_p <- function(x, g, B = 100000) { set.seed(20261012); o <- jt(x, g); r <- replicate(B, jt(x, sample(g)))
  (sum(abs(r - mean(r)) >= abs(o - mean(r))) + 1) / (B + 1) }
fmt <- function(p, B = 100000) ifelse(p <= 1 / (B + 1) + 1e-12, sprintf("< %.2e", 1 / B), sprintf("%.2e", p))
ord <- list("four_levels_MZL_FL_DLBCL_BL" = c(MZL = 1, FL = 2, DLBCL = 3, BL = 4),
            "three_levels_MZL_FL_DLBCL_BL_tied" = c(MZL = 1, FL = 2, DLBCL = 3, BL = 3))
v <- c(CNV = "cnv_score", MTP = "mtp_score", Prol = "prol_score", Diff = "diff_score", Ig = "Ig_score")
out <- do.call(rbind, lapply(names(ord), function(nm) do.call(rbind, lapply(names(v), function(k) {
  g <- ord[[nm]][P$subtype]; x <- P[[k]]; kt <- cor.test(x, g, method = "kendall", exact = FALSE)
  b$y <- rank(b[[v[k]]]) / nrow(b); b$g <- ord[[nm]][b$subtype]
  m1 <- lmer(y ~ g + (1 | pt), b, REML = FALSE); m0 <- lmer(y ~ 1 + (1 | pt), b, REML = FALSE)
  data.frame(ordering = nm, dimension = k, patient_JT_P = fmt(jt_p(x, g)), Kendall_tau = round(unname(kt$estimate), 2),
             Kendall_P = sprintf("%.2e", kt$p.value), cell_mixed_model_trend_P = sprintf("%.2e", anova(m0, m1)[2, 8]))
}))))
write.csv(out, file.path(.d, "results", "ordinal_trend_sensitivity.csv"), row.names = FALSE, fileEncoding = "UTF-8")
print(out)
