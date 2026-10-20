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
PANEL <- "PanelS4D"
R <- validate_anonymous_data(read_panel_csv(PANEL, "FigS4D_gradient_robustness.csv"))
ORD <- c("main", "no_cut", "excl_multi", "no_cns", "no_liquid", "cnv_high", "upper_tiers", "lmm_dataset", "lodo", "include_excluded_sample")
OPT <- c("no_liquid", "cnv_high")
if (!all(setdiff(ORD, OPT) %in% R$id)) stop("S4D missing analysis: ", paste(setdiff(setdiff(ORD, OPT), R$id), collapse = ", "))
if (length(setdiff(ORD, R$id))) warning("S4D missing optional analysis: ", paste(setdiff(ORD, R$id), collapse = ", "))
stopifnot(all(R$id %in% ORD))
ORD <- ORD[ORD %in% R$id]
if (any(R$p >= 0.05)) warning("S4D contains analyses with P >= 0.05.")
ROW_LABELS <- c(main = "Main analysis", no_cut = "No threshold", excl_multi = "Exclude repeated specimens", no_cns = "Exclude CNS samples", no_liquid = "Exclude liquid samples", cnv_high = "High CNV", upper_tiers = "Upper tiers", lmm_dataset = "Cohort-adjusted mixed model", lodo = "Leave one cohort out", include_excluded_sample = "Include excluded sample")
R[, row := unname(ROW_LABELS[id])]
R[, `:=`(y = length(ORD) + 1 - match(id, ORD), lp = -log10(p), main = id == "main", lab = fmtP(p))]
XR <- c(0, 9.0)
stopifnot(max(R$lp) < XR[2] - 0.3)
R[, inside := lp >= 3.8]
p <- ggplot(R, aes(y = y)) +
  geom_vline(xintercept = -log10(0.05), linetype = "22", linewidth = LW_AXIS, colour = "grey30") +
  geom_rect(aes(xmin = 0, xmax = lp, ymin = y - 0.34, ymax = y + 0.34, fill = main), colour = NA) +
  scale_fill_manual(values = c(`TRUE` = HI_COL, `FALSE` = LYM_FILL), guide = "none") +
  geom_text(data = R[inside == TRUE], aes(x = lp - 0.15, label = lab), hjust = 1, size = SZ_SMALL / PT, family = "Arial",
            colour = ifelse(R[inside == TRUE]$main, "white", "black")) +
  geom_text(data = R[inside == FALSE], aes(x = lp + 0.15, label = lab), hjust = 0, size = SZ_SMALL / PT, family = "Arial",
            colour = "black") +
  geom_text(aes(x = XR[2] + 0.3, label = n_lab), hjust = 0, size = SZ_SMALL / PT, family = "Arial", colour = NEU_COL) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8), expand = c(0, 0)) +
  scale_y_continuous(breaks = R$y, labels = R$row, expand = c(0, 0)) +
  coord_cartesian(xlim = XR, ylim = c(0.45, length(ORD) + 0.55), clip = "off") +
  labs(x = "\u2212log10 P", y = NULL) +
  theme_fig4() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = SZ_SMALL + 0.5, colour = "black"),
        plot.margin = margin(HEADER_MM / 25.4 * 72, 19.5 / 25.4 * 72, SZ * 0.35, SZ * 0.35))
if (!build_only) {
  save_figure(draw_single_panel(p, "D", DS_TAG$S4D, PANEL_MM$S4D[2]), RES(PANEL), "FigureS2_PanelD_gradient_robustness",
         PANEL_MM$S4D[1], PANEL_MM$S4D[2])
}
