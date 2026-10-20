.release_script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
.release_root <- normalizePath(file.path(.release_script_dir, "..", ".."))
source(file.path(.release_root, "privacy.R"), encoding = "UTF-8")
if (!exists("FIG4", inherits = TRUE)) {
  .d <- local({ a <- commandArgs(FALSE); f <- gsub("~+~", " ", sub("^--file=", "", a[grep("^--file=", a)]), fixed = TRUE)
                if (length(f)) dirname(normalizePath(f)) else getwd() })
  source_release_dependency(file.path(.release_root, "shared", "figure4_config.R"))
  source_release_dependency(file.path(.release_root, "shared", "figure4_plot_utils.R"))
}
if (!exists("build_only", inherits = TRUE)) build_only <- FALSE
PANEL <- "PanelF"
Fd <- validate_anonymous_data(read_panel_csv(PANEL, "Fig4F_within_subtype_effect_sizes.csv"))
if (!all(c("MHC_I", "MHC_II") %in% Fd$prog)) stop("Panel F requires prepared MHC-I and MHC-II plot data.")
ORD  <- intersect(c("MHC_I", "MHC_II"), Fd$prog)
LAB  <- setNames(Fd$label[match(ORD, Fd$prog)], ORD)
ROWS <- c(intersect(ENT_ORD, unique(Fd$row)), "All")
stopifnot(setequal(unique(Fd$row), ROWS), Fd[row == "All", uniqueN(n) == 1L])
NR <- length(ROWS)
Fd[, `:=`(blk = match(prog, ORD), lab = factor(LAB[prog], LAB), all = row == "All")]
Fd[, y := (length(ORD) - blk) * (NR + 1) + (NR + 1 - match(row, ROWS))]
Fd[, txt := sprintf("%d/%d %s", n_consistent, n, stars)]
Fd[, dirn := fifelse(hl_pct < 0, "down", "up")]
COL <- c(down = DIV_LOW, up = DIV_HIGH)
lo <- min(c(Fd$lo_pct, 0)); hi <- max(c(Fd$hi_pct, 0))
XR <- c(floor((lo - 3) / 10) * 10, ceiling((hi + 3) / 10) * 10)
BRK <- pretty(XR, n = 3); BRK <- BRK[BRK >= XR[1] & BRK <= XR[2]]
YB <- Fd[, .(y, row, txt)][order(-y)]
p <- ggplot(Fd, aes(y = y)) +
  geom_vline(xintercept = 0, linetype = "22", linewidth = LW_AXIS, colour = "grey45") +
  geom_segment(aes(x = lo_pct, xend = hi_pct, yend = y, colour = dirn), linewidth = LW_DATA) +
  geom_point(data = Fd[all == FALSE], aes(x = hl_pct, fill = dirn), shape = 21, size = 2.1, stroke = 0.35, colour = "white") +
  geom_point(data = Fd[all == TRUE], aes(x = hl_pct, fill = dirn), shape = 23, size = 3.0, stroke = 0.35, colour = "white") +
  scale_colour_manual(values = COL, guide = "none") +
  scale_fill_manual(values = COL, guide = "none") +
  facet_wrap(~ lab, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = BRK, limits = XR, expand = c(0, 0)) +
  scale_y_continuous(breaks = YB$y, labels = YB$row, expand = expansion(add = 0.6),
                     sec.axis = sec_axis(~ ., breaks = YB$y, labels = YB$txt)) +
  labs(x = "Malig-3 vs Malig-1 (% change)", y = NULL) +
  theme_fig4() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y.right = element_text(size = SZ_SMALL, colour = "black", hjust = 0, margin = margin(l = SZ * 0.3)),
        panel.spacing.y = unit(SZ * 0.5, "pt"),
        strip.text = element_text(size = SZ, colour = "black", hjust = 0.5, margin = margin(b = SZ * 0.3, t = 0)))
p <- preserve_header(p)
if (!build_only) {
  save_figure(draw_single_panel(p, "F", DS_TAG$F, PANEL_MM$F[2]), RES(PANEL), "Figure4_PanelF_within_subtype_effects",
         PANEL_MM$F[1], PANEL_MM$F[2],
         highlight_color = DIV_LOW)
}
