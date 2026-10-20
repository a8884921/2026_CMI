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
W_MM <- 116; H_MM <- 73
bins  <- validate_anonymous_data(utils::read.csv(P("panelE_binned_trends.csv"), check.names = FALSE))
tests <- validate_anonymous_data(utils::read.csv(P("panelE_trend_tests.csv"), check.names = FALSE))
TOP4 <- tests$feature[order(tests$rank)]
stopifnot(length(TOP4) == 4, setequal(unique(bins$feature), TOP4), all(tests$p < 0.05),
          all(tests$direction %in% names(DIR_COLS)), all((tests$rho > 0) == (tests$direction == "higher in tumor")),
          all(tests$rho_lo < tests$rho & tests$rho < tests$rho_hi), all(sign(tests$rho_lo) == sign(tests$rho_hi)),
          all(bins$y_lo <= bins$y_mean & bins$y_mean <= bins$y_hi),
          all(table(bins$feature[!bins$tie_bin]) == 10))
X_SCALE <- list(
  IgD_IgM_ratio     = scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.002, base = 10),
                                         limits = c(0, 1), breaks = c(0, 0.01, 0.1, 1), labels = c("0", "0.01", "0.1", "1"),
                                         expand = c(0.08, 0)),
  mucosal_secretory = scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"),
                                         expand = c(0.08, 0)),
  Ig_burden         = scale_x_log10(limits = c(5e-4, 0.3), breaks = c(0.001, 0.01, 0.1), labels = c("0.001", "0.01", "0.1"),
                                    expand = c(0.02, 0)),
  J_chain_secretory = scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"),
                                         expand = c(0.08, 0)))
stopifnot(setequal(names(X_SCALE), TOP4),
          all(bins$x_mean[bins$feature == "Ig_burden"] > 5e-4))
cells <- validate_anonymous_data(utils::read.csv(R("labeled_cells_24_indices_and_Igscore.csv"), check.names = FALSE))
TR <- list(IgD_IgM_ratio = list(f = function(x) asinh(x / (2 * 0.002)) / log(10), b = function(t) 2 * 0.002 * sinh(t * log(10))),
           Ig_burden = list(f = log10, b = function(t) 10^t),
           mucosal_secretory = list(f = identity, b = identity), J_chain_secretory = list(f = identity, b = identity))
fits <- do.call(rbind, lapply(TOP4, function(f) {
  d <- data.frame(x = cells[[f]], y = cells$ig_score); d <- d[is.finite(d$x) & is.finite(d$y), ]
  stopifnot(nrow(d) == tests$n[tests$feature == f])
  d$tx <- TR[[f]]$f(d$x); d <- d[is.finite(d$tx), ]
  m <- stats::lm(y ~ tx, data = d)
  stopifnot(sign(coef(m)[["tx"]]) == sign(tests$rho[tests$feature == f]))
  xr <- TR[[f]]$f(range(bins$x_mean[bins$feature == f]))
  g <- data.frame(tx = seq(xr[1], xr[2], length.out = 100))
  pr <- stats::predict(m, g, se.fit = TRUE)
  data.frame(feature = f, x = TR[[f]]$b(g$tx), fit = pr$fit, lo = pr$fit - 1.96 * pr$se.fit, hi = pr$fit + 1.96 * pr$se.fit,
             n_fit = nrow(d), slope = coef(m)[["tx"]], r2 = summary(m)$r.squared)
}))
KEY_TIE <- "Tied value (\u2265 5% of cells)"
DIR_LABS <- c("higher in tumor" = "Higher in tumor", "lower in tumor" = "Lower in tumor")
mk <- function(f, first) {
  b <- bins[bins$feature == f, ]; nb <- b[!b$tie_bin, ]; tb <- b[b$tie_bin, ]
  st <- tests[tests$feature == f, ]
  dcol <- DIR_COLS[[st$direction]]
  fi <- fits[fits$feature == f, ]
  lo <- min(c(nb$y_lo, tb$y_lo, fi$lo)); hi <- max(c(nb$y_hi, tb$y_hi, fi$hi)); pad <- 0.06 * (hi - lo)
  ybrk_f <- pretty(c(lo, hi), n = 3); ybrk_f <- ybrk_f[ybrk_f >= lo - pad & ybrk_f <= hi + pad]
  ylim_f <- range(c(lo - pad, hi + pad, ybrk_f))
  leg <- data.frame(x_mean = rep(nb$x_mean[1], 4), y_mean = rep(nb$y_mean[1], 4), dir = rep(names(DIR_COLS), each = 2))
  mn <- function(v) gsub("-", "\u2212", sprintf("%.3f", v), fixed = TRUE)
  stat_lines <- c(sprintf("\u03c1 = %s", mn(st$rho)), sprintf("(%s, %s)", mn(st$rho_lo), mn(st$rho_hi)), P_txt(st$p),
                  sprintf("n = %s", format(st$n, big.mark = ",")))
  p <- ggplot() +
    geom_ribbon(data = fi, aes(x = x, ymin = lo, ymax = hi), fill = dcol, alpha = 0.20) +
    geom_line(data = fi, aes(x = x, y = fit), colour = dcol, linewidth = LW_DATA * 2.0, alpha = 0.35, lineend = "round") +
    geom_errorbar(data = nb, aes(x = x_mean, ymin = y_lo, ymax = y_hi), width = 0, colour = dcol, linewidth = LW_ERR) +
    geom_line(data = nb, aes(x = x_mean, y = y_mean), colour = dcol, linewidth = LW_DATA * 0.75) +
    geom_point(data = nb, aes(x = x_mean, y = y_mean), colour = dcol, size = 1.0) +
    geom_line(data = leg, aes(x = x_mean, y = y_mean, colour = dir, group = dir), linewidth = LW_DATA * 0.75, alpha = 0) +
    geom_point(data = leg, aes(x = x_mean, y = y_mean, colour = dir), size = 1.0, alpha = 0) +
    scale_colour_manual(NULL, values = DIR_COLS, limits = names(DIR_COLS), labels = DIR_LABS[names(DIR_COLS)],
                        guide = guide_legend(order = 1, override.aes = list(alpha = 1))) +
    X_SCALE[[f]] +
    scale_y_continuous(if (first) "Ig-score" else NULL, limits = ylim_f, breaks = ybrk_f,
                       labels = function(v) sub("\\.?0+$", "", sprintf("%.2f", v)), expand = c(0.03, 0)) +
    labs(x = f) + coord_cartesian(clip = "off") + theme_fig2() + row3_legend_theme() +
    theme(plot.margin = ggplot2::margin(0, SZ_TICK * 0.30, 0, 0),
          axis.title.x = element_text(size = SZ_TICK))
  if (nrow(tb))
    p <- p + geom_errorbar(data = tb, aes(x = x_mean, ymin = y_lo, ymax = y_hi), width = 0,
                           colour = dcol, linewidth = LW_ERR) +
      geom_point(data = tb, aes(x = x_mean, y = y_mean), shape = 21, fill = "white", colour = dcol,
                 size = 1.9, stroke = 0.6) +
      geom_point(data = tb, aes(x = x_mean, y = y_mean, shape = KEY_TIE), fill = "white", colour = GREY_DARK,
                 size = 1.9, stroke = 0.6, alpha = 0) +
      scale_shape_manual(NULL, values = setNames(21, KEY_TIE), guide = guide_legend(order = 2, override.aes = list(alpha = 1)))
  if (first) p <- p + theme(plot.margin = ggplot2::margin(0, SZ_TICK * 0.30, 0, SZ_TICK * 0.25))
  list(plot = p, strip = stat_strip(stat_lines, cols = c(dcol, dcol, "black", "black")))
}
pl <- lapply(seq_along(TOP4), function(i) mk(TOP4[i], i == 1))
n_cells <- tests$n[tests$feature == "Ig_burden"]; stopifnot(length(n_cells) == 1, n_cells == max(tests$n))
build_E <- function(pad_pt) {
  q <- patchwork::wrap_plots(c(lapply(pl, `[[`, "strip"), lapply(pl, `[[`, "plot")), nrow = 2,
                             heights = grid::unit(c(STAT_STRIP_MM, 1), c("mm", "null"))) +
    patchwork::plot_layout(guides = "collect")
  dataset_line(q, sprintf("%s \u00b7 %s B cells", DATASETS_TXT, format(n_cells, big.mark = ",")),
               pw_theme = row3_legend_theme() + theme(legend.box = "horizontal", plot.margin = ggplot2::margin(0, 0, 0, 0),
                                                      legend.box.margin = ggplot2::margin(-SZ_TICK * 0.5 + pad_pt, 0, 0, 0)))
}
pad_mm <- ROW3_BOTTOM_MM - panel_bottom_mm(build_E(0), W_MM, H_MM)
stopifnot(pad_mm >= 0)
pE <- build_E(pad_mm * 72.27 / 25.4)
stopifnot(abs(panel_bottom_mm(pE, W_MM, H_MM) - ROW3_BOTTOM_MM) < 0.01)
pE <- PANEL_TAG(pE, "E")
save_panel(pE, W_MM, H_MM)
out <- bins |> dplyr::left_join(tests[, c("feature", "n", "rho", "rho_lo", "rho_hi", "p", "rank", "importance_pct")],
                                by = "feature", suffix = c("_bin", "_feature")) |>
  dplyr::arrange(rank, x_mean) |>
  dplyr::transmute(importance_rank = rank, feature = feature, importance_pct = round(importance_pct, 2),
                   bin_type = ifelse(tie_bin, "tied_value_bin", "decile_bin"), n_cells_in_bin = n_bin,
                   x_mean = signif(x_mean, 4), x_min = signif(x_min, 4), x_max = signif(x_max, 4),
                   Igscore_mean = round(y_mean, 4), Igscore_95CI_lower = round(y_lo, 4), Igscore_95CI_upper = round(y_hi, 4),
                   Igscore_median = round(y_med, 4), Igscore_Q25 = round(y_q25, 4), Igscore_Q75 = round(y_q75, 4),
                   n_available_cells = n_feature, Spearman_rho = round(rho, 4),
                   Spearman_rho_95CI = sprintf("%.3f \u2013 %.3f", rho_lo, rho_hi), Spearman_P = fmtP(p))
utils::write.csv(out, O("top4_drivers_and_binned_trends.csv"), row.names = FALSE)
utils::write.csv(transform(fits, x = signif(x, 5), fit = round(fit, 4), lo = round(lo, 4), hi = round(hi, 4), slope = signif(slope, 4), r2 = round(r2, 4)),
                 O("correlation_regression_lines.csv"), row.names = FALSE)
for (f in TOP4) invisible(NULL)
for (i in seq_len(nrow(tests)))
  invisible(NULL)
