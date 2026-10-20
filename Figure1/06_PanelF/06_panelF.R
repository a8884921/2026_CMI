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
source_release_dependency(file.path(.release_root, "shared", "figure1_config.R"))
DATASETS_TXT <- "Analysis cohorts"
suppressMessages({ library(readr); library(dplyr); library(ggplot2) })
W_MM <- 68; H_MM <- 62
curves <- validate_anonymous_data(readr::read_csv(P("panelF_ROC_curves.csv"), show_col_types = FALSE))
aucs   <- validate_anonymous_data(readr::read_csv(P("panelF_AUC_summary.csv"), show_col_types = FALSE))
stopifnot(all(c("method", "fpr", "tpr", "tpr_lo", "tpr_hi") %in% names(curves)),
          all(c("method", "auc", "auc_lo", "auc_hi") %in% names(aucs)),
          setequal(aucs$method, SCORE_COLS), setequal(unique(curves$method), SCORE_COLS))
n_cells <- unique(aucs$n_cells); stopifnot(length(n_cells) == 1)
aucs <- aucs |>
  dplyr::mutate(name = lab_of(method),
                label = sprintf("%s  %.3f (%.3f\u2013%.3f)", name, auc, auc_lo, auc_hi)) |>
  dplyr::arrange(dplyr::desc(auc))
readr::write_csv(aucs[, c("method", "name", "auc", "auc_lo", "auc_hi", "n_cells",
                          "n_tumor", "n_normal", "label")], O("AUC_and_DeLong_confidence_intervals.csv"))
curves$method <- factor(curves$method, levels = aucs$method)
LINE_GAP <- 0.085
txt <- aucs |>
  dplyr::mutate(x = 1.00, y = 0.030 + (rev(seq_len(dplyr::n())) - 1) * LINE_GAP)
hdr <- data.frame(x = 1.00, y = max(txt$y) + LINE_GAP, label = "AUC (95% CI)")
pF <- ggplot(curves, aes(fpr, tpr, colour = method, fill = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = SEM[["neutral"]], linewidth = LW_AXIS) +
  geom_ribbon(aes(ymin = tpr_lo, ymax = tpr_hi), alpha = 0.16, colour = NA) +
  geom_line(aes(linewidth = method == "Ig_score")) +
  scale_linewidth_manual(values = c(`TRUE` = LW_DATA * 1.25, `FALSE` = LW_DATA * 0.85), guide = "none") +
  geom_text(data = txt, aes(x = x, y = y, label = label, colour = method), inherit.aes = FALSE,
            hjust = 1, size = SZ_TICK / PT, family = FIG_FONT) +
  geom_text(data = hdr, aes(x = x, y = y, label = label), inherit.aes = FALSE,
            hjust = 1, size = SZ_TICK / PT, family = FIG_FONT, colour = "black") +
  scale_colour_manual(values = ROC_COLS, guide = "none") +
  scale_fill_manual(values = ROC_COLS, guide = "none") +
  scale_x_continuous("1 \u2212 Specificity", limits = c(0, 1),
                     breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), expand = c(0.01, 0.01)) +
  scale_y_continuous("Sensitivity", limits = c(0, 1),
                     breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), expand = c(0.01, 0.01)) +
  theme_fig1() +
  theme(plot.margin = ggplot2::margin(SZ_TICK * 0.9, SZ_TICK * 0.6, SZ_TICK * 0.4, SZ_TICK * 0.4),
        axis.title.x = element_text(margin = ggplot2::margin(t = SZ_TICK * 0.45 + 2.07)))
pF <- PANEL_TAG(pF, "F")
pF <- dataset_line(pF, sprintf("Datasets: 6 studies \u00b7 %s B cells", format(n_cells, big.mark = ",")))
save_panel(pF, W_MM, H_MM)
print(as.data.frame(aucs[, c("method", "name", "auc", "auc_lo", "auc_hi")]), row.names = FALSE)
stopifnot(aucs$name[aucs$method == "prol_score"] == "Prol",
          aucs$name[aucs$method == "diff_score"] == "Diff",
          aucs$name[which.min(aucs$auc)] == "Prol",
          aucs$name[which.max(aucs$auc)] == "Ig")
