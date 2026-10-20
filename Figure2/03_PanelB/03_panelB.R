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
W_MM <- 180; H_MM <- 120
COLW <- c(lab = 34, forest = 44, dir = 9, p = 15, n = 12, gap = 3, imp = 52)
source_release_dependency(file.path(.release_root, "shared", "panelBC_layout.R"))
d   <- validate_anonymous_data(utils::read.csv(P("panelB_metric_discrimination.csv"), check.names = FALSE))
imp <- validate_anonymous_data(utils::read.csv(P("panelC_feature_importance.csv"), check.names = FALSE))
ord <- validate_anonymous_data(utils::read.csv(P("panelBC_row_order.csv"), check.names = FALSE))
stopifnot(nrow(d) == 24, setequal(d$feature, FEATURE24), all(d$igscore_auc_same_cells > d$auc_strength),
          setequal(imp$feature, FEATURE24), abs(sum(imp$pct) - 100) < 1e-6)
rows <- bc_rows(ord)
stopifnot(identical(rows$metrics$direction, d$direction[match(rows$metrics$feature, d$feature)]))
m <- dplyr::left_join(rows$metrics, dplyr::select(d, -category, -top4, -direction), by = "feature") |>
  dplyr::mutate(sig = ifelse(q_auc_BH < 0.05, "Single metric (q < 0.05)", "Single metric (n.s.)"),
                dir = ifelse(direction == "higher in tumor", "Higher in tumor", "Lower in tumor"),
                p_lab = fmtP(p_igscore_vs_metric), n_lab = format(n, big.mark = ",", trim = TRUE))
th <- bc_theme()
pLab <- ggplot() + bc_bands(rows) +
  geom_text(data = rows$headers, aes(x = 0.02, y = y, label = label), hjust = 0, fontface = "bold",
            size = SZ_TICK / PT, family = FIG_FONT) +
  geom_text(data = m, aes(x = 0.06, y = y, label = disp), hjust = 0, size = SZ_TICK / PT, family = FIG_FONT) +
  annotate("text", x = 0.02, y = rows$top, label = "Ig metric", hjust = 0, size = SZ_TICK / PT,
           family = FIG_FONT, colour = GREY_DARK) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + bc_y(rows) + th + bc_blank_x() +
  theme(plot.margin = ggplot2::margin(SZ_TICK * 1.6, 1.5 * 72.27 / 25.4, 0, 0))
pFor <- ggplot(m) + bc_bands(rows) +
  annotate("segment", x = 0.5, xend = 0.5, y = 1.5, yend = rows$n_rows + 0.6, linetype = "dashed", colour = GREY_MID,
           linewidth = LW_AXIS * 0.8) +
  geom_segment(aes(x = auc_strength, xend = igscore_auc_same_cells, y = y, yend = y),
               colour = GREY_LIGHT, linewidth = 1.1) +
  geom_errorbar(aes(xmin = strength_lo, xmax = strength_hi, y = y, colour = direction), orientation = "y", width = 0,
                linewidth = LW_ERR * 0.8) +
  scale_colour_manual(values = DIR_COLS, guide = "none") +
  geom_point(aes(x = auc_strength, y = y, shape = sig, fill = sig), size = 1.5, stroke = 0.5, colour = GREY_DARK,
             alpha = 0) +
  lapply(list(c("higher in tumor", TRUE), c("higher in tumor", FALSE), c("lower in tumor", TRUE), c("lower in tumor", FALSE)),
         function(k) {
           sel <- m$direction == k[1] & (m$q_auc_BH < 0.05) == as.logical(k[2])
           geom_point(data = m[sel, ], aes(x = auc_strength, y = y), shape = 21, size = 1.5, stroke = 0.5,
                      colour = DIR_COLS[[k[1]]], fill = if (as.logical(k[2])) DIR_COLS[[k[1]]] else "white",
                      show.legend = FALSE)
         }) +
  geom_point(aes(x = igscore_auc_same_cells, y = y, shape = "Ig-score, same cells",
                 fill = "Ig-score, same cells"), size = 1.9, stroke = 0.4, colour = IG_COL) +
  scale_shape_manual(NULL, values = c("Single metric (q < 0.05)" = 21, "Single metric (n.s.)" = 21,
                                      "Ig-score, same cells" = 23),
                     breaks = c("Single metric (q < 0.05)", "Single metric (n.s.)", "Ig-score, same cells")) +
  scale_fill_manual(NULL, values = c("Single metric (q < 0.05)" = GREY_DARK, "Single metric (n.s.)" = "white",
                                     "Ig-score, same cells" = IG_COL),
                    breaks = c("Single metric (q < 0.05)", "Single metric (n.s.)", "Ig-score, same cells")) +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  scale_x_continuous("Discrimination (AUC)", breaks = c(0.5, 0.6, 0.7, 0.8, 0.9),
                     labels = c("0.5", "0.6", "0.7", "0.8", "0.9"), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.44, 0.93)) +
  bc_y(rows) + th + bc_x_axis() + theme(legend.position = "bottom")
pDir <- ggplot(m) + bc_bands(rows) +
  geom_point(data = m[m$q_auc_BH < 0.05, ], aes(x = 0.5, y = y, shape = dir, fill = dir), size = 1.5, stroke = 0.3,
             colour = "white") +
  geom_point(data = m[m$q_auc_BH >= 0.05, ], aes(x = 0.5, y = y, shape = dir, colour = dir), size = 1.25, stroke = 0.55,
             fill = "white", show.legend = FALSE) +
  scale_shape_manual(NULL, values = c("Higher in tumor" = 24, "Lower in tumor" = 25)) +
  scale_fill_manual(NULL, values = c("Higher in tumor" = BOX_COLS[["tumor"]], "Lower in tumor" = BOX_COLS[["normal"]])) +
  scale_colour_manual(values = c("Higher in tumor" = BOX_COLS[["tumor"]], "Lower in tumor" = BOX_COLS[["normal"]]),
                      guide = "none") +
  annotate("text", x = 0.5, y = rows$top, label = "Tumor", size = SZ_TICK / PT, family = FIG_FONT, colour = GREY_DARK) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + bc_y(rows) + th + bc_blank_x() +
  theme(legend.position = "bottom")
pP <- ggplot(m) + bc_bands(rows) +
  geom_text(aes(x = 1, y = y, label = p_lab), hjust = 1, size = SZ_TICK / PT, family = FIG_FONT) +
  annotate("text", x = 1, y = rows$top, label = "Paired P", hjust = 1, size = SZ_TICK / PT, family = FIG_FONT,
           colour = GREY_DARK) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + bc_y(rows) + th + bc_blank_x()
pN <- ggplot(m) + bc_bands(rows) +
  geom_text(aes(x = 1, y = y, label = n_lab), hjust = 1, size = SZ_TICK / PT, family = FIG_FONT) +
  annotate("text", x = 1, y = rows$top, label = "n", hjust = 1, size = SZ_TICK / PT, family = FIG_FONT,
           colour = GREY_DARK) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + bc_y(rows) + th + bc_blank_x()
pImp <- bc_importance_col(rows, imp)
n_train <- unique(imp$n_train); stopifnot(length(n_train) == 1)
pB <- (pLab | pFor | pDir | pP | pN | patchwork::plot_spacer() | pImp) +
  patchwork::plot_layout(widths = grid::unit(COLW, "mm"), guides = "collect")
pB <- PANEL_TAG(pB, "B")
n_cells <- d$n[d$feature == "Ig_burden"]; stopifnot(n_cells == max(d$n))
pB <- dataset_line(pB, sprintf("%s \u00b7 %s B cells | C: GBM training set, %s B cells (lymphoma, tonsil, Raji, pan-cancer)",
                               DATASETS_TXT, format(n_cells, big.mark = ","), format(n_train, big.mark = ",")),
                   pw_theme = theme(legend.position = "bottom", legend.box = "horizontal",
                                    legend.box.margin = ggplot2::margin(-SZ_TICK * 1.1, 0, -SZ_TICK * 0.4, 0),
                                    legend.spacing.y = grid::unit(0, "pt"),
                                    legend.margin = ggplot2::margin(0, 0, 0, 0),
                                    legend.text = element_text(size = SZ_TICK, family = FIG_FONT),
                                    legend.key.size = grid::unit(8, "pt"),
                                    plot.margin = ggplot2::margin(0, 0, 0, 0)))
C_X <- 0.6 + panel_box_mm(pB, W_MM, H_MM, "panel-7")[["x0"]] - panel_box_mm(pB, W_MM, H_MM, "panel-1")[["x0"]]
STARS <- bc_star_pos(pB, rows, "panel-1", COLW[["lab"]], W_MM, H_MM)
save_panel(pB, W_MM, H_MM, draw_fun = function() {
  print(pB)
  grid::grid.text("C", x = grid::unit(C_X, "mm"), y = grid::unit(H_MM - 0.6, "mm"), just = c("left", "top"),
                  gp = grid::gpar(fontface = "bold", fontsize = SZ_TAG, fontfamily = FIG_FONT, col = "black"))
  draw_stars(STARS)
})
out <- m |> dplyr::arrange(y) |>
  dplyr::transmute(category = category, feature = feature, top4_driver = top4, n, n_normal, n_tumor,
                   direction = direction, AUC_tumor_positive = round(auc_tumor, 4),
                   discrimination_strength = round(auc_strength, 4), discrimination_strength_95CI = sprintf("%.3f\u2013%.3f", strength_lo, strength_hi),
                   P_AUC_vs_0.5 = fmtP(p_auc_vs_0.5), q_BH = fmtP(q_auc_BH),
                   Igscore_AUC_same_cells = round(igscore_auc_same_cells, 4),
                   Igscore_95CI = sprintf("%.3f\u2013%.3f", igscore_lo, igscore_hi),
                   P_Igscore_vs_metric_paired_DeLong = fmtP(p_igscore_vs_metric),
                   P_Welch = fmtP(p_welch), q_Welch_BH = fmtP(q_welch_BH))
utils::write.csv(out, O("metric_discrimination_and_Igscore_comparison.csv"), row.names = FALSE)
