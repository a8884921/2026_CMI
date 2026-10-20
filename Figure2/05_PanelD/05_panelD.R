.release_script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
.release_root <- normalizePath(file.path(.release_script_dir, "..", ".."))
source(file.path(.release_root, "privacy.R"), encoding = "UTF-8")
.script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE); f <- gsub("~+~", " ", sub("^--file=", "", a[grep("^--file=", a)]), fixed = TRUE)
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
PANEL_DIR <- normalizePath(.script_dir)
source_release_dependency(file.path(.release_root, "shared", "figure2_config.R"))
DATASETS_TXT <- "Analysis cohorts"
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })
W_MM <- 64; H_MM <- 73
roc  <- validate_anonymous_data(utils::read.csv(P("panelD_ROC_curves.csv"), check.names = FALSE))
perf <- validate_anonymous_data(utils::read.csv(P("panelD_discrimination.csv"), check.names = FALSE))
g <- function(item, col = "value") perf[[col]][perf$item == item]
auc_s <- g("Ig-score"); auc_b <- g("Ig_burden"); p_pair <- g("paired DeLong (Ig-score vs Ig_burden)")
n_cells <- g("Ig-score", "n")
stopifnot(setequal(unique(roc$curve), c("Ig-score", "Ig_burden")), auc_s > auc_b, p_pair < 0.05, n_cells == 48745)
rob <- validate_anonymous_data(utils::read.csv(P("panelD_robustness.csv"), check.names = FALSE))
rv  <- function(item) rob$value[rob$item == item]
stopifnot(rv("ds_wins") == rv("ds_evaluable"), rv("ds_max_p") < 0.05, rv("oof_vs_burden_p") < 0.05,
          (is.na(rv("pt_boot_p")) && rv("pt_boot_n_le0") == 0) || isTRUE(rv("pt_boot_p") < 0.05))
stopifnot(auc_b > 0.5)
CURVE_COLS <- c("Ig-score" = IG_COL, "Ig_burden" = DIR_COLS[["higher in tumor"]])
CURVE_LW   <- c("Ig-score" = LW_DATA, "Ig_burden" = LW_DATA * 0.8)
roc$curve <- factor(roc$curve, levels = c("Ig_burden", "Ig-score"))
ci_txt <- function(k) sprintf("%s  AUC %.3f (%.3f, %.3f)", k, g(k), g(k, "lo"), g(k, "hi"))
stat_lines <- c(ci_txt("Ig-score"), ci_txt("Ig_burden"), P_txt(p_pair),
                sprintf("%d/%d datasets: P \u2264 %s", as.integer(rv("ds_wins")), as.integer(rv("ds_evaluable")), fmtP(rv("ds_max_p"))))
pD <- ggplot(roc, aes(fpr, tpr)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = GREY_MID, linewidth = LW_AXIS * 0.8) +
  geom_path(aes(colour = curve, linewidth = curve), lineend = "round") +
  scale_colour_manual(NULL, values = CURVE_COLS, breaks = c("Ig-score", "Ig_burden")) +
  scale_linewidth_manual(NULL, values = CURVE_LW, breaks = c("Ig-score", "Ig_burden")) +
  scale_x_continuous("1 \u2212 Specificity", limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"),
                     expand = c(0.01, 0)) +
  scale_y_continuous("Sensitivity", limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"),
                     expand = c(0.02, 0)) +
  coord_cartesian(clip = "off") + theme_fig2() + row3_legend_theme()
build_D <- function(pad_pt) {
  q <- (stat_strip(stat_lines, cols = c(CURVE_COLS[["Ig-score"]], CURVE_COLS[["Ig_burden"]], "black", "black")) /
          (pD + theme(plot.margin = ggplot2::margin(0, SZ_TICK * 0.9, 0, SZ_TICK * 0.6),
                      legend.box.margin = ggplot2::margin(-SZ_TICK * 0.5 + pad_pt, 0, 0, 0)))) +
    patchwork::plot_layout(heights = grid::unit(c(STAT_STRIP_MM, 1), c("mm", "null")))
  dataset_line(q, sprintf("%s \u00b7 %s B cells", DATASETS_TXT, format(n_cells, big.mark = ",")),
               pw_theme = theme(plot.margin = ggplot2::margin(0, 0, 0, 0)))
}
pad_mm <- ROW3_BOTTOM_MM - panel_bottom_mm(build_D(0), W_MM, H_MM)
stopifnot(pad_mm >= 0)
pD <- build_D(pad_mm * 72.27 / 25.4)
stopifnot(abs(panel_bottom_mm(pD, W_MM, H_MM) - ROW3_BOTTOM_MM) < 0.01)
pD <- PANEL_TAG(pD, "D")
save_panel(pD, W_MM, H_MM)
out <- data.frame(
  metric = c("Igscore_AUC", "Ig_burden_AUC", "Igscore_vs_Ig_burden_paired_DeLong_P",
           "Igscore_cell_level_five_fold_out_of_fold_AUC",
           "out_of_fold_vs_Ig_burden_paired_DeLong_P",
           "patient_block_bootstrap_delta_AUC", "patient_bootstrap_percentile_P",
           sprintf("Igscore_higher_datasets_normal_ge%d", as.integer(rv("ds_min_normal"))), "maximum_dataset_paired_P",
           "training_five_fold_CV_AUC"),
  value = c(sprintf("%.4f", auc_s), sprintf("%.4f", auc_b), fmtP(p_pair),
         sprintf("%.4f", rv("oof_auc")), fmtP(rv("oof_vs_burden_p")),
         sprintf("%.4f", rv("pt_boot_delta")), if (is.na(rv("pt_boot_p"))) "< 2.00e-05" else fmtP(rv("pt_boot_p")),
         sprintf("%d / %d", as.integer(rv("ds_wins")), as.integer(rv("ds_evaluable"))), fmtP(rv("ds_max_p")),
         sprintf("%.4f", g("model 5-fold CV AUC (training set)"))),
  `95%CI` = c(sprintf("%.3f\u2013%.3f", g("Ig-score", "lo"), g("Ig-score", "hi")),
              sprintf("%.3f\u2013%.3f", g("Ig_burden", "lo"), g("Ig_burden", "hi")), "",
              sprintf("%.3f\u2013%.3f", rob$lo[rob$item == "oof_auc"], rob$hi[rob$item == "oof_auc"]), "",
              sprintf("%.3f-%.3f (percentile)", rv("pt_boot_delta_lo"), rv("pt_boot_delta_hi")), "", "", "",
              sprintf("%.3f\u2013%.3f", g("model 5-fold CV AUC (training set)", "lo"), g("model 5-fold CV AUC (training set)", "hi"))),
  n = c(n_cells, n_cells, n_cells, n_cells, n_cells, n_cells, n_cells, NA, NA, g("model 5-fold CV AUC (training set)", "n")),
  check.names = FALSE)
utils::write.csv(out, O("Igscore_discrimination.csv"), row.names = FALSE)
