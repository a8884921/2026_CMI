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
source_release_dependency(file.path(.release_root, "shared", "panelBC_layout.R"))
W_MM <- 180; H_MM <- 120
COLW <- c(lab = 36, box = 100, gap = 3, p = 16, n = 14)
bx  <- validate_anonymous_data(utils::read.csv(P("panelA2_box_statistics.csv"), check.names = FALSE))
sg  <- validate_anonymous_data(utils::read.csv(P("panelA2_significance.csv"), check.names = FALSE))
ord <- validate_anonymous_data(utils::read.csv(P("panelBC_row_order.csv"), check.names = FALSE))
stopifnot(nrow(bx) == 48, setequal(bx$group, c("normal", "tumor")), nrow(sg) == 24, setequal(sg$feature, FEATURE24),
          all(bx$p05 <= bx$p25 & bx$p25 <= bx$p50 & bx$p50 <= bx$p75 & bx$p75 <= bx$p95),
          all(bx$p25 >= -2.0 & bx$p75 <= 2.9))
rows <- bc_rows(ord)
OFF <- c(normal = -0.215, tumor = 0.215); HH <- 0.195
b <- dplyr::left_join(bx, rows$metrics[, c("feature", "y")], by = "feature") |>
  dplyr::mutate(yy = y + OFF[group], grp = GRP_LABELS[group])
stopifnot(identical(rows$metrics$direction, sg$direction[match(rows$metrics$feature, sg$feature)]))
s <- dplyr::left_join(dplyr::select(rows$metrics, -direction), dplyr::select(sg, -category, -top4), by = "feature") |>
  dplyr::mutate(p_lab = fmtP(p_welch), n_lab = format(n, big.mark = ",", trim = TRUE))
n_grp <- c(normal = unique(sg$n_normal[sg$feature == "Ig_burden"]), tumor = unique(sg$n_tumor[sg$feature == "Ig_burden"]))
LEG <- setNames(sprintf("%s (n = %s)", GRP_LABELS, format(n_grp[names(GRP_LABELS)], big.mark = ",", trim = TRUE)),
                GRP_LABELS)
th <- bc_theme()
pLab <- ggplot() + bc_bands(rows) +
  geom_text(data = rows$headers, aes(x = 0.02, y = y, label = label), hjust = 0, fontface = "bold",
            size = SZ_TICK / PT, family = FIG_FONT) +
  geom_text(data = rows$metrics, aes(x = 0.06, y = y, label = disp), hjust = 0, size = SZ_TICK / PT, family = FIG_FONT) +
  annotate("text", x = 0.02, y = rows$top, label = "Ig metric", hjust = 0, size = SZ_TICK / PT,
           family = FIG_FONT, colour = GREY_DARK) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + bc_y(rows) + th + bc_blank_x() +
  theme(plot.margin = ggplot2::margin(SZ_TICK * 1.6, 1.5 * 72.27 / 25.4, 0, 0))
XL <- c(-2.0, 2.9)
cut_hi <- b[b$p95 > XL[2], ]; cut_lo <- b[b$p05 < XL[1], ]
AR <- grid::arrow(length = grid::unit(1.1, "mm"), type = "closed")
pBox <- ggplot(b) + bc_bands(rows) +
  annotate("segment", x = 0, xend = 0, y = 1.5, yend = rows$n_rows + 0.6, linetype = "dashed", colour = GREY_MID,
           linewidth = LW_AXIS * 0.8) +
  geom_segment(aes(x = p05, xend = p95, y = yy, yend = yy, colour = grp), linewidth = LW_ERR * 0.8) +
  geom_rect(aes(xmin = p25, xmax = p75, ymin = yy - HH, ymax = yy + HH, fill = grp), colour = NA) +
  geom_segment(aes(x = p50, xend = p50, y = yy - HH, yend = yy + HH), colour = "white", linewidth = 0.6) +
  geom_segment(data = b[b$p05 == b$p95, ], aes(x = p50, xend = p50, y = yy - HH, yend = yy + HH, colour = grp),
               linewidth = 0.9) +
  geom_segment(data = cut_hi, aes(x = XL[2] - 0.15, xend = XL[2], y = yy, yend = yy, colour = grp), arrow = AR,
               linewidth = LW_ERR * 0.8, show.legend = FALSE) +
  geom_segment(data = cut_lo, aes(x = XL[1] + 0.15, xend = XL[1], y = yy, yend = yy, colour = grp), arrow = AR,
               linewidth = LW_ERR * 0.8, show.legend = FALSE) +
  scale_colour_manual(NULL, values = setNames(unname(BOX_COLS[names(GRP_LABELS)]), GRP_LABELS), labels = LEG) +
  scale_fill_manual(NULL, values = setNames(unname(BOX_COLS[names(GRP_LABELS)]), GRP_LABELS), labels = LEG) +
  scale_x_continuous("Value (z-score)", breaks = c(-1, 0, 1, 2), expand = c(0, 0)) +
  coord_cartesian(xlim = XL) +
  bc_y(rows) + th + bc_x_axis() +
  theme(legend.position = "bottom", legend.key.size = grid::unit(SZ_TICK, "pt"),
        legend.text = element_text(size = SZ_TICK, family = FIG_FONT))
pP <- ggplot(s) + bc_bands(rows) +
  geom_text(aes(x = 1, y = y, label = p_lab), hjust = 1, size = SZ_TICK / PT, family = FIG_FONT) +
  annotate("text", x = 1, y = rows$top, label = "P", hjust = 1, size = SZ_TICK / PT, family = FIG_FONT,
           colour = GREY_DARK) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + bc_y(rows) + th + bc_blank_x()
pN <- ggplot(s) + bc_bands(rows) +
  geom_text(aes(x = 1, y = y, label = n_lab), hjust = 1, size = SZ_TICK / PT, family = FIG_FONT) +
  annotate("text", x = 1, y = rows$top, label = "n", hjust = 1, size = SZ_TICK / PT, family = FIG_FONT,
           colour = GREY_DARK) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + bc_y(rows) + th + bc_blank_x()
pA2 <- (pLab | pBox | patchwork::plot_spacer() | pP | pN) +
  patchwork::plot_layout(widths = grid::unit(COLW, "mm"), guides = "collect")
pA2 <- dataset_line(pA2, sprintf("%s \u00b7 %s B cells", DATASETS_TXT, format(sum(n_grp), big.mark = ",")),
                    pw_theme = theme(legend.position = "bottom", legend.box = "horizontal",
                                     legend.box.margin = ggplot2::margin(-SZ_TICK * 1.1, 0, -SZ_TICK * 0.4, 0),
                                     legend.margin = ggplot2::margin(0, 0, 0, 0),
                                     legend.text = element_text(size = SZ_TICK, family = FIG_FONT),
                                     legend.key.size = grid::unit(8, "pt"),
                                     plot.margin = ggplot2::margin(0, 0, 0, 0)))
STARS <- bc_star_pos(pA2, rows, "panel-1", COLW[["lab"]], W_MM, H_MM)
save_panel(pA2, W_MM, H_MM, draw_fun = function() { print(pA2); draw_stars(STARS) })
wide <- bx |> dplyr::select(feature, group, p25, p50, p75) |>
  tidyr::pivot_wider(names_from = group, values_from = c(p25, p50, p75))
out <- s |> dplyr::arrange(y) |> dplyr::left_join(wide, by = "feature") |>
  dplyr::transmute(category = category, feature = feature, top4_driver = top4, n, n_normal, n_tumor, direction = direction,
                   normal_median_z = round(p50_normal, 3), normal_IQR_z = sprintf("%.3f \u2013 %.3f", p25_normal, p75_normal),
                   tumor_median_z = round(p50_tumor, 3), tumor_IQR_z = sprintf("%.3f \u2013 %.3f", p25_tumor, p75_tumor),
                   P_Welch = fmtP(p_welch), q_BH = fmtP(q_welch_BH), q_below_0.05 = q_welch_BH < 0.05)
utils::write.csv(out, O("metric_quantiles_and_Welch_tests.csv"), row.names = FALSE)
