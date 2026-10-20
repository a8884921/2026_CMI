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
PANEL <- "PanelB"
comp <- validate_anonymous_data(read_panel_csv(PANEL, "Fig4B_mean_sample_tier_composition.csv"))
st   <- validate_anonymous_data(read_panel_csv(PANEL, "Fig4B_statistical_tests.csv"))
stopifnot(setequal(unique(comp$entity), ENT_ORD))
comp[, tier := factor(tier, TIER_LV)]
comp[, entity := factor(entity, ENT_ORD)]; setorder(comp, entity, tier)
comp[, `:=`(xmax = cumsum(frac) * 100, xmin = (cumsum(frac) - frac) * 100), by = entity]
comp[, `:=`(y0 = RIDGE_POS[as.character(entity)] + 0.05, y1 = RIDGE_POS[as.character(entity)] + 0.55)]
m3 <- comp[tier == "Malig-3"]
m3[, `:=`(inside = frac >= 0.20, lab = sprintf("%.1f%%", frac * 100))]
nn <- unique(comp[, .(entity, n_samples, y = RIDGE_POS[as.character(entity)] + 0.30)])
p_trend <- st[id == "trend3", p]
TOP <- max(RIDGE_POS) + 1.02
XMAX <- 123
right_axis <- sec_axis(~ ., breaks = nn$y, labels = sprintf("n = %d", nn$n_samples))
key <- data.table(tier = TIER_LV, x = c(0, 37, 74), y = RIDGE_POS[["Reactive B"]] + 0.30)
p <- ggplot(comp) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = y0, ymax = y1, fill = tier), colour = "white", linewidth = 0.25 / PT) +
  scale_fill_manual(values = TIER_COL, guide = "none") +
  geom_text(data = m3[inside == TRUE], aes(x = 98.5, y = (y0 + y1) / 2, label = lab), hjust = 1, colour = "white",
            size = SZ_SMALL / PT, family = "Arial") +
  geom_text(data = m3[inside == FALSE], aes(x = 101.5, y = (y0 + y1) / 2, label = lab), hjust = 0, colour = HI_COL,
            size = SZ_SMALL / PT, family = "Arial") +
  geom_rect(data = key, aes(xmin = x, xmax = x + 6, ymin = y - 0.18, ymax = y + 0.18, fill = tier), colour = NA) +
  geom_text(data = key, aes(x = x + 8, y = y, label = tier), hjust = 0, size = SZ_SMALL / PT, family = "Arial") +
  annotate("text", x = 0, y = RIDGE_POS[["Pan-cancer B"]] + 0.30, hjust = 0, vjust = 0.5, colour = HI_COL,
           size = SZ_SMALL / PT, family = "Arial", label = sprintf("Malig-3 trend P = %s", fmtP(p_trend))) +
  annotate("segment", x = 0, xend = 100, y = RIDGE_YLIM[1], yend = RIDGE_YLIM[1], linewidth = LW_AXIS, colour = "black",
           lineend = "square") +
  scale_x_continuous(limits = c(0, XMAX), breaks = c(0, 50, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(RIDGE_YLIM[1], TOP + 0.55), breaks = RIDGE_POS + 0.30, labels = NULL, expand = c(0, 0),
                     sec.axis = right_axis) +
  coord_cartesian(clip = "off") +
  labs(x = "Cells per sample (%)", y = NULL) +
  theme_fig4() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y.right = element_text(size = SZ_SMALL, colour = NEU_COL, margin = margin(l = SZ * 0.5)),
        axis.title.x = element_text(hjust = 100 / XMAX / 2, margin = margin(t = SZ * 0.45)))
p <- preserve_header(p)
if (!build_only) {
  p_single <- suppressMessages(p + scale_y_continuous(limits = c(RIDGE_YLIM[1], TOP + 0.55), breaks = RIDGE_POS[ENT_ORD] + 0.30,
                                                     labels = ENT_ORD, expand = c(0, 0), sec.axis = right_axis)) +
    theme(axis.text.y = element_text(size = SZ, colour = "black", margin = margin(r = SZ * 0.2)))
  save_figure(draw_single_panel(p_single, "B", DS_TAG$B, PANEL_MM$B[2]), RES(PANEL), "Figure4_PanelB_tier_composition",
         PANEL_MM$B[1], PANEL_MM$B[2])
}
