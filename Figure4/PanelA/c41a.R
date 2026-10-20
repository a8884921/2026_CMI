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
PANEL <- "PanelA"
dens <- validate_anonymous_data(read_panel_csv(PANEL, "Fig4A_population_density.csv"))
med  <- validate_anonymous_data(read_panel_csv(PANEL, "Fig4A_population_medians_and_cell_counts.csv"))
cuts <- validate_anonymous_data(read_panel_csv(PANEL, "Fig4A_tier_thresholds.csv"))
stopifnot(setequal(unique(dens$pop), POP_ORD), nrow(cuts) == 3L)
XL  <- c(-2.75, 4.75); H <- 0.86
c12 <- cuts$dCMI[2]; c23 <- cuts$dCMI[3]
dens <- dens[x >= XL[1] & x <= XL[2]]
dens[, y0 := RIDGE_POS[pop]][, y1 := y0 + H * density_scaled]
dens[, cls := fifelse(pop %in% ENT_ORD, "lym", pop)]
TOP <- max(RIDGE_POS) + 1.02
mk <- med[pop %in% ENT_ORD]
mk[, y0 := RIDGE_POS[pop]]
mk[, ytop := mapply(function(g, m) { d <- dens[pop == g]; RIDGE_POS[[g]] + H * approx(d$x, d$density_scaled, m)$y }, pop, median_dCMI)]
mk[, lab := sprintf("%.2f", median_dCMI)]
band <- data.table(tier = c("Baseline", TIER_LV), xmin = c(XL[1], 0, c12, c23), xmax = c(0, c12, c23, XL[2]))
band_lab <- data.table(x = c(-1.38, c12 / 2, (c12 + c23) / 2, (c23 + XL[2]) / 2),
                       lab = c("Baseline", "M1", "M2", "M3"),
                       col = c("#6B7280", "#C9790F", TIER_COL[["Malig-2"]], TIER_COL[["Malig-3"]]))
ncell <- med[, .(pop, y = RIDGE_POS[pop] + 0.50, lab = paste0("n = ", format_thousands(n_cells)))]
p <- ggplot() +
  lapply(seq_len(nrow(band)), function(i) annotate("rect", xmin = band$xmin[i], xmax = band$xmax[i],
         ymin = RIDGE_YLIM[1], ymax = TOP - 0.10, fill = TIER_COL[[band$tier[i]]], alpha = 0.20, colour = NA)) +
  annotate("segment", x = c(c12, c23), xend = c(c12, c23), y = RIDGE_YLIM[1], yend = TOP - 0.10,
           linetype = "22", linewidth = LW_AXIS, colour = "grey35") +
  annotate("segment", x = 0, xend = 0, y = RIDGE_YLIM[1], yend = TOP - 0.10, linewidth = LW_ERR, colour = "grey15") +
  geom_ribbon(data = dens[cls == "Reactive B"], aes(x, ymin = y0, ymax = y1), fill = BEN_COL[["Reactive B"]], alpha = 0.55) +
  geom_ribbon(data = dens[cls == "Pan-cancer B"], aes(x, ymin = y0, ymax = y1), fill = BEN_COL[["Pan-cancer B"]], alpha = 0.60) +
  geom_line(data = dens[cls != "lym"], aes(x, y1, group = pop), colour = BEN_COL[["Reactive B"]], linewidth = LW_ERR) +
  geom_ribbon(data = dens[cls == "lym"], aes(x, ymin = y0, ymax = y1, group = pop), fill = LYM_FILL, alpha = 0.85) +
  geom_line(data = dens[cls == "lym"], aes(x, y1, group = pop), colour = LYM_LINE, linewidth = LW_ERR) +
  geom_segment(data = mk, aes(x = median_dCMI, xend = median_dCMI, y = y0, yend = ytop), colour = HI_COL, linewidth = LW_DATA) +
  white_label(mk, aes(x = median_dCMI + 0.07, y = y0 + 0.21, label = lab), colour = HI_COL, hjust = 0, vjust = 0.5,
           size = SZ_SMALL / PT) +
  geom_text(data = ncell, aes(x = XL[2] - 0.04, y = y, label = lab), hjust = 1, colour = NEU_COL,
            size = SZ_SMALL / PT, family = "Arial") +
  geom_text(data = band_lab, aes(x = x, y = TOP + 0.02, label = lab), colour = band_lab$col, vjust = 0,
            size = SZ_SMALL / PT, family = "Arial") +
  scale_x_continuous(limits = XL, breaks = c(-2, 0, 2, 4), expand = c(0, 0)) +
  scale_y_continuous(limits = c(RIDGE_YLIM[1], TOP + 0.55), breaks = RIDGE_POS + 0.30, labels = POP_LAB[names(RIDGE_POS)],
                     expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = "\u0394CMI (logit)", y = NULL) +
  theme_fig4() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = SZ, colour = "black", margin = margin(r = SZ * 0.2)))
p <- preserve_header(p)
if (!build_only) {
  save_figure(draw_single_panel(p, "A", DS_TAG$A, PANEL_MM$A[2]), RES(PANEL), "Figure4_PanelA_malignancy_distribution",
         PANEL_MM$A[1], PANEL_MM$A[2])
}
