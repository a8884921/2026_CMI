# =============================================================================
# Figure library: one set of ggplot builders for the live interface and the exported files.
# Color roles: CMI grades use an ordinal single-hue ramp (Mal-1 #F2A26A, Mal-2 #D45624,
# Mal-3 #861D21); not malignant is gray (#D3D8E2); CMI-low / CMI-high are blue / dark red.
# Typeface: Liberation Sans (metric-compatible with Arial). Text is 7-10.5 pt at the export width
# (7.2 in); the page shows a 200 dpi PNG of the same figure.
# =============================================================================

CMI_COL <- list(
  ink = "#1C2340", ink2 = "#4B536E", muted = "#6E7590", grid = "#E7E9EF", axis = "#C5CAD5", surface = "#FFFFFF",
  below = "#D3D8E2", mal1 = "#F2A26A", mal2 = "#D45624", mal3 = "#861D21", cd8 = "#2160C4", other = "#D3D8E2",
  low = "#2160C4", high = "#861D21", accent = "#1C2340", outline = "#8A90A2",
  div_low = "#1C5CAB", div_mid = "#F1F0ED", div_high = "#A8322F")
STATE_COLS <- c("Below threshold" = CMI_COL$below, "Mal-1" = CMI_COL$mal1, "Mal-2" = CMI_COL$mal2, "Mal-3" = CMI_COL$mal3)
GROUP_COLS <- c("CMI-low" = CMI_COL$low, "CMI-high" = CMI_COL$high)
FIG_FAMILY <- "Liberation Sans"
FIG_BASE <- 9
# Export sizes in inches (width, height); NA heights are set from the data (see *_height helpers).
FIG_SIZES <- list(fig1_cmi = c(7.2, 4.2), fig2_grades = c(7.2, 2.8), fig1_riskscore = c(7.2, 4.0), fig1_grade_map = c(7.2, 6.8))
fig_width <- function(name) FIG_SIZES[[name]][1]

# Text rendering needs a UTF-8 locale (Delta, minus and dash glyphs). Hosted R usually runs in
# en_US.UTF-8; minimal containers may only have C.UTF-8.
ensure_utf8 <- function() {
  if (isTRUE(l10n_info()$`UTF-8`)) return(invisible(TRUE))
  for (loc in c("en_US.UTF-8", "C.UTF-8", "C.utf8")) if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) return(invisible(TRUE))
  invisible(FALSE)
}

# Make the bundled Liberation Sans visible to every device: a private fontconfig file adds the app's
# font folder for cairo (PDF) and systemfonts (ragg PNG, svglite). Call before any device opens.
register_figure_fonts <- function(dir) {
  ensure_utf8()
  dir <- normalizePath(dir, mustWork = FALSE)
  if (dir.exists(dir) && .Platform$OS.type == "unix" && !nzchar(Sys.getenv("CMI_FONTCONFIG_DONE"))) {
    conf <- file.path(tempdir(), "cmi-fonts.conf")
    writeLines(c('<?xml version="1.0"?>', '<!DOCTYPE fontconfig SYSTEM "fonts.dtd">', '<fontconfig>',
                 '  <include ignore_missing="yes">/etc/fonts/fonts.conf</include>', sprintf('  <dir>%s</dir>', dir),
                 sprintf('  <cachedir>%s</cachedir>', file.path(tempdir(), "fontconfig-cache")), '</fontconfig>'), conf)
    Sys.setenv(FONTCONFIG_FILE = conf, CMI_FONTCONFIG_DONE = "1")
  }
  ok <- FALSE
  if (requireNamespace("systemfonts", quietly = TRUE)) {
    have <- tryCatch(any(tolower(systemfonts::system_fonts()$family) == tolower(FIG_FAMILY)), error = function(e) FALSE)
    if (!have && dir.exists(dir)) {
      have <- tryCatch({ systemfonts::register_font(name = FIG_FAMILY, plain = file.path(dir, "LiberationSans-Regular.ttf"),
                                                    bold = file.path(dir, "LiberationSans-Bold.ttf"), italic = file.path(dir, "LiberationSans-Italic.ttf"),
                                                    bolditalic = file.path(dir, "LiberationSans-Bold.ttf")); TRUE },
                       error = function(e) grepl("already exists", conditionMessage(e)))
    }
    ok <- have
  }
  if (!ok) { FIG_FAMILY <<- "sans"; message("Figure font Liberation Sans could not be registered; using the default sans-serif font.") }
  invisible(ok)
}

# ---- text helpers -------------------------------------------------------------------------------
lab_size <- function(points) points / ggplot2::.pt                    # geom_text size (mm) from points
wrap_chars <- function(width_in, points) max(30L, floor((width_in - 0.45) * 72 / (0.5 * points)))
wrap_txt <- function(x, width = 110) paste(vapply(strsplit(x, "\n")[[1]], function(l) paste(strwrap(l, width = width), collapse = "\n"), character(1)), collapse = "\n")
labs_wrapped <- function(..., subtitle = NULL, caption = NULL, width_in = 7.2, base_size = FIG_BASE) {
  ggplot2::labs(..., subtitle = if (!is.null(subtitle) && nzchar(subtitle)) wrap_txt(subtitle, wrap_chars(width_in, base_size - 0.5)),
                caption = if (!is.null(caption) && nzchar(caption)) wrap_txt(caption, wrap_chars(width_in, base_size - 2)))
}
annotate_set <- function(title, subtitle = NULL, caption = NULL, width_in = 7.2, base_size = FIG_BASE) {
  patchwork::plot_annotation(title = title, subtitle = if (!is.null(subtitle) && nzchar(subtitle)) wrap_txt(subtitle, wrap_chars(width_in, base_size - 0.5)),
                             caption = if (!is.null(caption) && nzchar(caption)) wrap_txt(caption, wrap_chars(width_in, base_size - 2)),
                             theme = theme_cmi(base_size))
}
has_patchwork <- function() requireNamespace("patchwork", quietly = TRUE)

theme_cmi <- function(base_size = FIG_BASE, grid = "y") {
  g <- ggplot2::element_line(colour = CMI_COL$grid, linewidth = 0.3)
  ggplot2::theme_minimal(base_size = base_size, base_family = FIG_FAMILY) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = CMI_COL$ink2),
      plot.title = ggplot2::element_text(colour = CMI_COL$ink, face = "bold", size = base_size + 1.5, hjust = 0, margin = ggplot2::margin(b = 3)),
      plot.subtitle = ggplot2::element_text(colour = CMI_COL$ink2, size = base_size - 0.5, hjust = 0, margin = ggplot2::margin(b = 8), lineheight = 1.15),
      plot.caption = ggplot2::element_text(colour = CMI_COL$muted, size = base_size - 2, hjust = 0, margin = ggplot2::margin(t = 8), lineheight = 1.15),
      plot.title.position = "plot", plot.caption.position = "plot",
      axis.title = ggplot2::element_text(colour = CMI_COL$ink2, size = base_size - 0.5),
      axis.text = ggplot2::element_text(colour = CMI_COL$ink2, size = base_size - 1),
      axis.line.x = ggplot2::element_line(colour = CMI_COL$axis, linewidth = 0.35),
      axis.ticks = ggplot2::element_line(colour = CMI_COL$axis, linewidth = 0.35), axis.ticks.length = grid::unit(2.5, "pt"),
      panel.grid.major.x = if (grid %in% c("x", "xy")) g else ggplot2::element_blank(),
      panel.grid.major.y = if (grid %in% c("y", "xy")) g else ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(colour = CMI_COL$ink, face = "bold", size = base_size - 0.5, hjust = 0, margin = ggplot2::margin(b = 4)),
      legend.position = "top", legend.justification = "left", legend.title = ggplot2::element_text(size = base_size - 1, colour = CMI_COL$ink2),
      legend.text = ggplot2::element_text(size = base_size - 1, colour = CMI_COL$ink2), legend.key.size = grid::unit(9, "pt"),
      legend.margin = ggplot2::margin(0, 0, 0, 0), legend.box.spacing = grid::unit(4, "pt"),
      plot.margin = ggplot2::margin(10, 14, 8, 8), plot.background = ggplot2::element_rect(fill = CMI_COL$surface, colour = NA))
}
panel_title <- function(base_size = FIG_BASE) ggplot2::theme(plot.title = ggplot2::element_text(colour = CMI_COL$ink, face = "bold", size = base_size, margin = ggplot2::margin(b = 4)))

minus <- function(x) gsub("-", "\u2212", x, fixed = TRUE)
num_lab <- function(x) minus(format(x, trim = TRUE, drop0trailing = TRUE))
fmt1 <- function(v) { r <- round(v, 1); r[!is.na(r) & r == 0] <- 0; minus(sprintf("%.1f", r)) }
# Tick labels with as many decimals as the tick spacing needs.
tick_lab <- function(x) {
  ok <- is.finite(x); if (sum(ok) < 2) return(num_lab(x))
  d <- max(0, -floor(log10(min(diff(sort(unique(x[ok]))))) + 1e-9)); out <- rep("", length(x)); out[ok] <- minus(formatC(x[ok], format = "f", digits = d)); out
}
# Two or three ticks inside the limits (narrow facets): try finer pretty() grids, then thin them.
pretty_in <- function(n = 3) function(l) {
  inside <- function(b) b[b >= l[1] - 1e-12 & b <= l[2] + 1e-12]
  b <- inside(pretty(l, n)); k <- n
  while (length(b) < 2 && k < n + 6) { k <- k + 1; b <- inside(pretty(l, k)) }
  if (length(b) > 3) b <- b[seq(1, length(b), by = ceiling((length(b) - 1) / 2))]
  b
}
integer_breaks <- function(lims) {
  b <- pretty(lims, n = 5); b <- b[abs(b - round(b)) < 1e-8 & b >= 0]
  if (length(b) < 2) b <- unique(round(seq(max(0, lims[1]), lims[2], length.out = 3)))
  b
}
# P values in figures avoid superscript glyphs (not in every figure font); exact values are in the tables.
fmt_p_fig <- function(p) {
  if (is.null(p) || !length(p) || is.na(p)) return("P not available")
  if (p < 0.001) return("P < 0.001")
  sprintf("P = %s", formatC(p, digits = if (p < 0.01) 1 else 2, format = "fg"))
}
sup_digits <- function(x) chartr("-0123456789", "\u207B\u2070\u00B9\u00B2\u00B3\u2074\u2075\u2076\u2077\u2078\u2079", x)
fmt_p_txt <- function(p) {
  if (is.null(p) || !length(p) || is.na(p)) return("P not available")
  if (p < 1e-3) { e <- floor(log10(p)); m <- p / 10^e; return(sprintf("P = %.1f \u00D7 10%s", m, sup_digits(as.character(e)))) }
  sprintf("P = %s", formatC(p, digits = 2, format = "fg"))
}
synthetic_tag <- function(synthetic) if (isTRUE(synthetic)) "Synthetic demonstration data \u2014 illustrates the software, not biology. " else ""
# Readable names for deconvolution classes: CD8_T -> CD8+ T cells, Normal_B -> Normal B cells.
pretty_class <- function(x) {
  y <- gsub("[_.]+", " ", x)
  y <- sub("^(CD[0-9]+) ?T( cells?)?$", "\\1+ T cells", y, ignore.case = TRUE)
  y <- sub("^(Normal|Naive|Memory|Plasma) B$", "\\1 B cells", y)
  y <- sub("^(NK|Mast|Dendritic|Stromal|Endothelial|Epithelial|Myeloid)$", "\\1 cells", y)
  y <- sub("^(Macrophage|Monocyte|Neutrophil|Fibroblast)$", "\\1s", y)
  ifelse(nchar(y) > 24, paste0(substr(y, 1, 23), "\u2026"), y)
}
empty_fig <- function(title, message) {
  ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = wrap_txt(message, 80), size = lab_size(9), colour = CMI_COL$ink2, family = FIG_FAMILY, lineheight = 1.1) +
    ggplot2::labs(title = title) + ggplot2::theme_void(base_family = FIG_FAMILY) +
    ggplot2::theme(plot.title = ggplot2::element_text(colour = CMI_COL$ink, face = "bold", size = FIG_BASE + 1.5, margin = ggplot2::margin(b = 6)),
                   plot.title.position = "plot", plot.margin = ggplot2::margin(10, 14, 10, 8), plot.background = ggplot2::element_rect(fill = CMI_COL$surface, colour = NA))
}
# Equal-width bins aligned to each state segment, so no bin straddles a state boundary.
seg_breaks <- function(lo, hi, cuts, w = 0.1) {
  e <- c(lo, cuts[cuts > lo & cuts < hi], hi)
  unique(unlist(lapply(seq_along(e)[-1], function(i) seq(e[i - 1], e[i], length.out = max(1, round((e[i] - e[i - 1]) / w)) + 1))))
}

# ---- single-cell -------------------------------------------------------------------------------
# CMI values crowd just below 1, so the histogram axis is logit-spaced and labeled in CMI values.
CMI_CUTS <- c(CMI_THRESHOLD, CMI_INDEX2, CMI_INDEX3)
CMI_CUT_LABELS <- c("0.9755", "0.98926", "0.99620")
cmi_axis <- function(lo, hi) {
  cand <- c(0.5, 0.8, 0.9, 0.95, 0.97, 0.98, 0.99, 0.995, 0.998, 0.999, 0.9995, 0.9999)
  at <- stats::qlogis(cand); ok <- at >= lo & at <= hi; at <- at[ok]; cand <- cand[ok]
  keep <- integer(); last <- -Inf
  for (i in seq_along(at)) if (at[i] - last >= 0.6) { keep <- c(keep, i); last <- at[i] }
  list(at = at[keep], lab = format(cand[keep], drop0trailing = TRUE, scientific = FALSE, trim = TRUE))
}
fig_sc_cmi <- function(cells, synthetic = FALSE, width_in = fig_width("fig1_cmi")) {
  d <- cells; d$state <- factor(d$state, levels = names(STATE_COLS))
  d$lg <- stats::qlogis(pmin(pmax(d$CMI, 1e-6), 1 - 1e-6))
  cuts <- stats::qlogis(CMI_CUTS)
  q <- stats::quantile(d$lg, c(0.01, 0.99), names = FALSE)
  lo <- min(stats::qlogis(0.9), floor(q[1] * 2) / 2); hi <- max(stats::qlogis(0.999), ceiling(q[2] * 2) / 2)
  bw <- if (nrow(d) >= 5000) 0.05 else 0.1
  brks <- seg_breaks(lo, hi, cuts, bw)
  n_lo <- sum(d$lg < lo); n_hi <- sum(d$lg >= hi)
  dh <- d[d$lg >= lo & d$lg < hi, , drop = FALSE]
  ax <- cmi_axis(lo, hi)
  bands <- data.frame(mid = (c(lo, cuts) + c(cuts, hi)) / 2, lab = c("Below threshold", "Mal-1", "Mal-2", "Mal-3"))
  nm <- sum(d$malignant)
  ggplot2::ggplot(dh, ggplot2::aes(x = .data$lg, fill = .data$state)) +
    ggplot2::geom_histogram(breaks = brks, colour = CMI_COL$surface, linewidth = 0.15, closed = "left") +
    ggplot2::geom_vline(xintercept = cuts, colour = CMI_COL$ink2, linewidth = 0.3, linetype = "33") +
    ggplot2::geom_text(data = bands, ggplot2::aes(x = .data$mid, y = Inf, label = .data$lab), inherit.aes = FALSE, vjust = 1.4,
                       size = lab_size(8), fontface = "bold", colour = CMI_COL$ink2, family = FIG_FAMILY) +
    ggplot2::annotate("text", x = cuts, y = Inf, label = CMI_CUT_LABELS, hjust = -0.1, vjust = 3.6, size = lab_size(7), colour = CMI_COL$muted, family = FIG_FAMILY) +
    (if (n_hi > 0) ggplot2::annotate("text", x = hi, y = Inf, label = sprintf("%s %s above \u2192", format(n_hi, big.mark = ","), if (n_hi == 1) "cell" else "cells"), hjust = 1.02, vjust = 5.6, size = lab_size(7), colour = CMI_COL$ink2, family = FIG_FAMILY)) +
    (if (n_lo > 0) ggplot2::annotate("text", x = lo, y = Inf, label = sprintf("\u2190 %s %s below", format(n_lo, big.mark = ","), if (n_lo == 1) "cell" else "cells"), hjust = -0.02, vjust = 5.6, size = lab_size(7), colour = CMI_COL$ink2, family = FIG_FAMILY)) +
    ggplot2::scale_fill_manual(values = STATE_COLS, drop = FALSE, guide = "none") +
    ggplot2::scale_x_continuous(limits = c(lo, hi), expand = c(0, 0), breaks = ax$at, labels = ax$lab, name = "CMI") +
    ggplot2::scale_y_continuous(name = "Cells", breaks = integer_breaks, expand = ggplot2::expansion(mult = c(0, 0.3))) +
    ggplot2::coord_cartesian(clip = "off") +
    labs_wrapped(title = "CMI of each cell", width_in = width_in,
                 subtitle = sprintf("%s cells; %s (%.1f%%) at or above the malignancy threshold", format(nrow(d), big.mark = ","), format(nm, big.mark = ","), 100 * nm / nrow(d)),
                 caption = paste0(synthetic_tag(synthetic), "Dashed lines: malignancy threshold (CMI 0.9755) and grade boundaries (0.98926, 0.99620). The axis is logit-spaced because CMI values crowd just below 1.")) +
    theme_cmi() + ggplot2::theme(plot.margin = ggplot2::margin(6, 14, 6, 8))
}

fig_sc_grades <- function(cells, synthetic = FALSE, width_in = fig_width("fig2_grades")) {
  lv <- names(STATE_COLS); n <- as.integer(table(factor(cells$state, levels = lv))); tot <- max(1L, sum(n))
  d <- data.frame(state = factor(lv, levels = rev(lv)), n = n, lab = sprintf("%s (%.1f%%)", format(n, big.mark = ",", trim = TRUE), 100 * n / tot))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$n, y = .data$state, fill = .data$state)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = .data$lab), hjust = -0.1, size = lab_size(8), colour = CMI_COL$ink, family = FIG_FAMILY) +
    ggplot2::scale_fill_manual(values = STATE_COLS, guide = "none") +
    ggplot2::scale_x_continuous(name = "Cells", breaks = integer_breaks, expand = ggplot2::expansion(mult = c(0, 0.25))) +
    labs_wrapped(y = NULL, title = "Cells by CMI grade", width_in = width_in, subtitle = sprintf("%s cells", format(sum(n), big.mark = ",")),
                 caption = paste0(synthetic_tag(synthetic), "Below threshold: CMI < 0.9755; Mal-1: up to 0.98926; Mal-2: up to 0.99620; Mal-3: above 0.99620.")) +
    theme_cmi(grid = "x") + ggplot2::theme(axis.text.y = ggplot2::element_text(colour = CMI_COL$ink, face = "bold"))
}

# ---- bulk ---------------------------------------------------------------------------------------
fig_bulk_ranked <- function(scores, synthetic = FALSE, width_in = fig_width("fig1_riskscore")) {
  d <- scores[order(scores$riskscore), , drop = FALSE]; d$rank <- seq_len(nrow(d))
  d$cmi_group <- factor(d$cmi_group, levels = names(GROUP_COLS))
  med <- stats::median(d$riskscore); big <- nrow(d) > 300; ids <- nrow(d) <= 40
  ng <- as.integer(table(d$cmi_group))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$rank, y = .data$riskscore, colour = .data$cmi_group)) +
    ggplot2::geom_hline(yintercept = med, colour = CMI_COL$ink2, linewidth = 0.4) +
    ggplot2::annotate("text", x = 0.5, y = med, label = sprintf("Cohort median %.2f", med), hjust = 0, vjust = -0.7, size = lab_size(7.5), colour = CMI_COL$ink2, family = FIG_FAMILY) +
    (if (!big) ggplot2::geom_segment(ggplot2::aes(xend = .data$rank, y = med, yend = .data$riskscore), linewidth = 0.5, alpha = 0.55)) +
    ggplot2::geom_point(size = if (big) 1.1 else 2.2) +
    ggplot2::scale_colour_manual(values = GROUP_COLS, name = NULL, drop = FALSE, labels = sprintf("%s (n = %s)", names(GROUP_COLS), format(ng, big.mark = ","))) +
    labs_wrapped(x = if (ids) NULL else "Samples, ranked by CMI-riskscore", y = "CMI-riskscore", title = "CMI-riskscore of each sample", width_in = width_in,
                 subtitle = sprintf("%s samples; CMI-high: above the cohort median", format(nrow(d), big.mark = ",")),
                 caption = paste0(synthetic_tag(synthetic), "The CMI-riskscore is computed for each sample from its own gene ranks. ",
                                  "Higher values mean higher predicted risk; it is not a probability. The CMI-high/low split uses the cohort median.")) +
    theme_cmi()
  if (ids) p + ggplot2::scale_x_continuous(breaks = d$rank, labels = d$sample_id, expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = FIG_BASE - 2.5))
  else p + ggplot2::scale_x_continuous(breaks = integer_breaks, expand = ggplot2::expansion(mult = 0.02))
}

# ---- spatial ------------------------------------------------------------------------------------
GRADE_COLS <- c(STATE_COLS[c("Mal-1", "Mal-2", "Mal-3")], "Non-malignant" = CMI_COL$other)
# Point size (ggplot units) so that spots of the given pitch nearly touch on a panel of the given size.
spot_size_for <- function(xy, pitch, panel_w_in, panel_h_in) {
  rx <- diff(range(xy[, 1])); ry <- diff(range(xy[, 2]))
  if (!is.finite(pitch) || pitch <= 0 || rx <= 0 || ry <= 0) return(1.2)
  mm_u <- min(panel_w_in * 25.4 / (rx + pitch), panel_h_in * 25.4 / (ry + pitch))
  max(0.25, min(4, 0.92 * pitch * mm_u / 0.78))
}
fig_sp_map <- function(res, reverse_y = TRUE, synthetic = FALSE, width_in = fig_width("fig1_grade_map"), height_in = FIG_SIZES$fig1_grade_map[2]) {
  s <- res$spots; pitch <- res$pitch %||% NA_real_
  s$grade <- factor(s$grade, levels = names(GRADE_COLS)); s <- s[order(s$grade != "Non-malignant"), ]
  n <- as.integer(table(s$grade))
  pt <- spot_size_for(cbind(s$x, s$y), pitch, width_in - 0.5, height_in - 1.6)
  p <- ggplot2::ggplot(s, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$grade), size = pt, shape = 16) +
    ggplot2::scale_colour_manual(values = GRADE_COLS, labels = sprintf("%s (%s)", names(GRADE_COLS), format(n, big.mark = ",", trim = TRUE)), name = NULL, drop = FALSE) +
    ggplot2::coord_fixed(clip = "off") + ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3), nrow = 1)) +
    labs_wrapped(x = NULL, y = NULL, title = "CMI grade of each spot", width_in = width_in, subtitle = sprintf("%s spots", format(nrow(s), big.mark = ",")),
                 caption = paste0(synthetic_tag(synthetic), "Each spot takes the class of highest deconvolution weight; spots won by a non-malignant class are gray.")) +
    theme_cmi(grid = "none") + ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.line.x = ggplot2::element_blank())
  if (reverse_y) p <- p + ggplot2::scale_y_reverse()
  p
}

# ---- export --------------------------------------------------------------------------------------
# Each figure: SVG (editable text), PDF (vector, embedded fonts), 600 dpi PNG and a 200 dpi PNG for the
# page; plus a "panel" version without title, subtitle and caption (PDF and SVG) for assembling
# your own figures. The legend text is returned as attribute "legend" (written to legends.txt).
fig_text <- function(p) {
  a <- if (inherits(p, "patchwork")) p$patches$annotation else tryCatch(p$labels, error = function(e) list())
  clean <- function(x) if (is.null(x) || inherits(x, "waiver") || !length(x)) NULL else gsub("\\s*\n\\s*", " ", as.character(x))
  list(title = clean(a$title), subtitle = clean(a$subtitle), caption = clean(a$caption))
}
strip_text <- function(p) {
  if (inherits(p, "patchwork")) {            # plot_annotation(title = NULL) keeps the old text; clear the fields
    for (k in c("title", "subtitle", "caption")) p$patches$annotation[k] <- list(NULL)
    p
  } else p + ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL)
}
write_svg <- function(p, f, width, height) {
  args <- list(filename = f, width = width, height = height, bg = "white")
  if ("fix_text_size" %in% names(formals(svglite::svglite))) args$fix_text_size <- FALSE
  do.call(svglite::svglite, args); print(p); grDevices::dev.off()
  svg <- readLines(f, warn = FALSE, encoding = "UTF-8")
  writeLines(gsub("font-family: *['\"]?Liberation Sans['\"]?", "font-family: Arial, Liberation Sans, Helvetica, sans-serif", svg), f, useBytes = TRUE)
}
write_pdf <- function(p, f, width, height) {
  if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf(f, width = width, height = height, family = if (FIG_FAMILY == "sans") "sans" else FIG_FAMILY)
  else grDevices::pdf(f, width = width, height = height, family = "Helvetica", useDingbats = FALSE)
  print(p); grDevices::dev.off()
}
save_figure_set <- function(p, path_base, width = 7.2, height = 4.6) {
  out <- character()
  if (requireNamespace("svglite", quietly = TRUE)) { f <- paste0(path_base, ".svg"); write_svg(p, f, width, height); out <- c(out, f) }
  f <- paste0(path_base, ".pdf"); write_pdf(p, f, width, height); out <- c(out, f)
  for (spec in list(c("_600dpi.png", 600), c("_web.png", 200))) {
    f <- paste0(path_base, spec[1]); r <- as.numeric(spec[2])
    if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png(f, width = width, height = height, units = "in", res = r, background = "white")
    else grDevices::png(f, width = width, height = height, units = "in", res = r, bg = "white")
    print(p); grDevices::dev.off(); out <- c(out, f)
  }
  txt <- fig_text(p)
  if (!is.null(txt$title)) {
    pdir <- file.path(dirname(path_base), "panels"); dir.create(pdir, showWarnings = FALSE)
    pb <- file.path(pdir, basename(path_base)); q <- strip_text(p)
    h <- max(2, height - 0.22 * (1 + !is.null(txt$subtitle) + 2 * !is.null(txt$caption)))
    if (requireNamespace("svglite", quietly = TRUE)) write_svg(q, paste0(pb, ".svg"), width, h)
    write_pdf(q, paste0(pb, ".pdf"), width, h)
  }
  attr(out, "legend") <- txt
  out
}
