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
W_MM <- 180; H_MM <- 62
COLW <- c(lab = 42, forest = 78, gap = 3, p = 22, n = 26)
ds   <- validate_anonymous_data(utils::read.csv(P("dataset_discrimination.csv"), check.names = FALSE))
ds$label <- anonymous_cohorts(ds$label)
perf <- validate_anonymous_data(utils::read.csv(P("panelD_discrimination.csv"), check.names = FALSE))
g <- function(item, col = "value") perf[[col]][perf$item == item]
ev <- ds[ds$evaluable, ]
stopifnot(nrow(ds) == 9, nrow(ev) >= 1, all(ev$auc_ig > ev$auc_burden), all(ev$p_paired < 0.05))
ev <- ev[order(-ev$auc_ig), ]
all_row <- data.frame(label = "All cells (9 datasets)", subtype = "", n_normal = sum(ds$n_normal), n_tumor = sum(ds$n_tumor),
                      auc_ig = g("Ig-score"), ig_lo = g("Ig-score", "lo"), ig_hi = g("Ig-score", "hi"),
                      auc_burden = g("Ig_burden"), burden_lo = g("Ig_burden", "lo"), burden_hi = g("Ig_burden", "hi"),
                      p_paired = g("paired DeLong (Ig-score vs Ig_burden)"))
stopifnot(all_row$n_normal == 777, all_row$n_tumor == 47968)
r <- rbind(transform(ev[, names(all_row)], label = sprintf("%s (%s)", ev$label, ev$subtype)), all_row)
r$y <- c(seq_len(nrow(ev)) + 1, nrow(ev) + 2.4)
n_rows <- max(r$y)
r$p_lab <- fmtP(r$p_paired)
r$n_lab <- sprintf("%s / %s", format(r$n_normal, big.mark = ",", trim = TRUE), format(r$n_tumor, big.mark = ",", trim = TRUE))
SEP_Y <- nrow(ev) + 1.7
yrev <- scale_y_reverse(limits = c(n_rows + 0.6, 0.4), expand = c(0, 0))
th <- theme_fig2() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
                           axis.title.y = element_blank(), plot.margin = ggplot2::margin(SZ_TICK * 1.6, 0, 0, 0))
blank_x <- theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(),
                 axis.title.x = element_blank())
hdr <- function(x, label, hj) annotate("text", x = x, y = 1, label = label, hjust = hj, size = SZ_TICK / PT,
                                        family = FIG_FONT, colour = GREY_DARK)
stopifnot(all(r$subtype[r$subtype != ""] %in% names(SUBTYPE_COLS)))
TAG_W <- 2.6
pLab <- ggplot(r) + hdr(0, "Dataset (subtype)", 0) +
  geom_point(data = r[r$subtype != "", ], aes(x = 1.0 / COLW[["lab"]], y = y, colour = subtype), shape = 15, size = 2.3) +
  scale_colour_manual(values = SUBTYPE_COLS, guide = "none") +
  geom_text(aes(x = ifelse(subtype != "", TAG_W / COLW[["lab"]], 0), y = y, label = label), hjust = 0,
            size = SZ_TICK / PT, family = FIG_FONT) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + yrev + th + blank_x +
  theme(plot.margin = ggplot2::margin(SZ_TICK * 1.6, 1.5 * 72.27 / 25.4, 0, 0))
KEY_B <- "Ig_burden"; KEY_I <- "Ig-score"
pFor <- ggplot(r) +
  annotate("segment", x = 0.5, xend = 0.5, y = 1.5, yend = n_rows + 0.6, linetype = "dashed", colour = GREY_MID,
           linewidth = LW_AXIS * 0.8) +
  geom_segment(aes(x = auc_burden, xend = auc_ig, y = y, yend = y), colour = GREY_LIGHT, linewidth = 1.1) +
  geom_errorbar(aes(xmin = burden_lo, xmax = burden_hi, y = y), orientation = "y", width = 0,
                colour = DIR_COLS[["higher in tumor"]], linewidth = LW_ERR * 0.8) +
  geom_errorbar(aes(xmin = ig_lo, xmax = ig_hi, y = y), orientation = "y", width = 0, colour = IG_COL, linewidth = LW_ERR * 0.8) +
  geom_point(aes(x = auc_burden, y = y, shape = KEY_B, fill = KEY_B), size = 1.6, stroke = 0.4,
             colour = DIR_COLS[["higher in tumor"]]) +
  geom_point(aes(x = auc_ig, y = y, shape = KEY_I, fill = KEY_I), size = 1.9, stroke = 0.4, colour = IG_COL) +
  scale_shape_manual(NULL, values = setNames(c(21, 23), c(KEY_B, KEY_I)), breaks = c(KEY_I, KEY_B)) +
  scale_fill_manual(NULL, values = setNames(c(DIR_COLS[["higher in tumor"]], IG_COL), c(KEY_B, KEY_I)), breaks = c(KEY_I, KEY_B)) +
  scale_x_continuous("AUC (95% CI)", breaks = c(0.4, 0.6, 0.8, 1), labels = c("0.4", "0.6", "0.8", "1"), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.4, 1)) + yrev + th +
  theme(axis.line.x = element_line(linewidth = LW_AXIS, colour = "black"), legend.position = "bottom")
pP <- ggplot(r) + hdr(1, "Paired P", 1) +
  geom_text(aes(x = 1, y = y, label = p_lab), hjust = 1, size = SZ_TICK / PT, family = FIG_FONT) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + yrev + th + blank_x
pN <- ggplot(r) + hdr(1, "n (normal / tumor)", 1) +
  geom_text(aes(x = 1, y = y, label = n_lab), hjust = 1, size = SZ_TICK / PT, family = FIG_FONT) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + yrev + th + blank_x
pS <- (pLab | pFor | patchwork::plot_spacer() | pP | pN) +
  patchwork::plot_layout(widths = grid::unit(COLW, "mm"), guides = "collect")
pS <- dataset_line(pS, sprintf("%s \u00b7 %s B cells", DATASETS_TXT, format(all_row$n_normal + all_row$n_tumor, big.mark = ",")),
                   pw_theme = row3_legend_theme() + theme(legend.box = "horizontal",
                                                          legend.box.margin = ggplot2::margin(-SZ_TICK * 1.0, 0, -SZ_TICK * 0.4, 0),
                                                          plot.margin = ggplot2::margin(0, 0, 0, 0)))
SEP <- local({
  tf <- tempfile(fileext = ".svg"); svglite::svglite(tf, width = W_MM / 25.4, height = H_MM / 25.4)
  nm <- grep("^panel-[0-9]+$", patchwork::patchworkGrob(pS)$layout$name, value = TRUE)
  grDevices::dev.off(); unlink(tf)
  bx <- vapply(nm, function(k) panel_box_mm(pS, W_MM, H_MM, k), numeric(4))
  stopifnot(length(nm) == 4, diff(range(bx["y0", ])) < 1e-6, diff(range(bx["y1", ])) < 1e-6)
  y0 <- bx["y0", 1]; y1 <- bx["y1", 1]
  c(x0 = min(bx["x0", ]), x1 = max(bx["x1", ]), y = y0 + (y1 - y0) * ((n_rows + 0.6) - SEP_Y) / ((n_rows + 0.6) - 0.4))
})
stopifnot(SEP[["x0"]] > 0.3, SEP[["x1"]] < W_MM - 0.3)
draw_08 <- function() {
  print(pS)
  grid::grid.segments(grid::unit(SEP[["x0"]], "mm"), grid::unit(SEP[["y"]], "mm"),
                      grid::unit(SEP[["x1"]], "mm"), grid::unit(SEP[["y"]], "mm"),
                      gp = grid::gpar(col = GREY_LIGHT, lwd = 0.4 * ggplot2::.pt, lineend = "butt"))
}
save_panel(pS, W_MM, H_MM, draw_fun = draw_08)
out <- ds |> dplyr::transmute(dataset = label, subtype = subtype, normal_cells = n_normal, tumor_cells = n_tumor,
                              evaluable_normal_ge20 = evaluable,
                              Igscore_AUC = round(auc_ig, 4), Igscore_95CI = ifelse(evaluable, sprintf("%.3f\u2013%.3f", ig_lo, ig_hi), ""),
                              Ig_burden_AUC = round(auc_burden, 4),
                              Ig_burden_95CI = ifelse(evaluable, sprintf("%.3f\u2013%.3f", burden_lo, burden_hi), ""),
                              paired_DeLong_P = ifelse(evaluable, fmtP(p_paired), ""))
utils::write.csv(out, O("dataset_Igscore_vs_Ig_burden.csv"), row.names = FALSE)
