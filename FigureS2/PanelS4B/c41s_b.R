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
PANEL <- "PanelS4B"
S  <- validate_anonymous_data(read_panel_csv(PANEL, "FigS4B_sample_tier_proportions.csv"))
st <- validate_anonymous_data(read_panel_csv(PANEL, "FigS4B_statistical_tests.csv"))
gm <- validate_anonymous_data(read_panel_csv(PANEL, "FigS4B_group_medians.csv"))
stopifnot(nrow(S) > 0, nrow(gm) == 3L)
S[, `:=`(x = match(entity, ENT_ORD), y = f_m3 * 100)]
ns <- S[, .N, by = entity]
xlab <- sprintf("%s\n(%d)", ENT_ORD, ns$N[match(ENT_ORD, ns$entity)])
XR <- S[, .(xlo = min(x), xhi = max(x)), by = .(grp = g3)]
gm <- merge(gm, XR, by = "grp")[, m := m * 100]
stair <- step_data(gm)
ymax <- max(S$y); DROP <- 3
br <- data.table(a1 = c(3, 1), a2 = c(3, 2), b1 = c(4, 3), b2 = c(5, 5), y = ymax + 7 + c(0, 16),
                 lab = sprintf("P = %s", fmtP(c(st[id == "aggr_vs_fl", p], st[id == "gc_vs_nongc", p]))))
YTOP <- max(br$y) + 13
p <- ggplot(S, aes(x, y)) +
  geom_boxplot(aes(group = x), width = 0.62, outlier.shape = NA, fill = "#E7EAEE", colour = "grey35",
               linewidth = LW_ERR, fatten = 1.6) +
  geom_point(position = position_jitter(width = 0.17, height = 0, seed = 20261011), size = 0.55, alpha = 0.55,
             colour = "grey30", stroke = 0) +
  geom_segment(data = stair, aes(x = x, xend = xend, y = y, yend = yend), colour = HI_COL, linewidth = LW_DATA,
               lineend = "square", inherit.aes = FALSE) +
  group_brackets(br, drop = DROP) +
  scale_x_continuous(breaks = 1:5, labels = xlab, limits = c(0.45, 5.55), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 25, 50, 75), limits = c(-1.5, YTOP), expand = c(0, 0)) +
  labs(x = "Lymphoma (n samples)", y = "Malig-3 cells per sample (%)",
       subtitle = sprintf("Trend P = %s", fmtP(st[id == "trend3", p]))) +
  theme_fig4() +
  theme(plot.subtitle = element_text(size = SZ_SMALL, colour = HI_COL, hjust = 1, family = "Arial",
                                     margin = margin(b = SZ * 0.25)))
p <- preserve_header(p)
if (!build_only) {
  save_figure(draw_single_panel(p, "B", DS_TAG$S4B, PANEL_MM$S4B[2]), RES(PANEL), "FigureS2_PanelB_sample_Malig3_proportion",
         PANEL_MM$S4B[1], PANEL_MM$S4B[2])
}
