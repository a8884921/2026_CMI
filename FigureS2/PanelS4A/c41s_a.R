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
PANEL <- "PanelS4A"
cp <- validate_anonymous_data(read_panel_csv(PANEL, "FigS4A_subtype_stage_composition.csv"))
kw <- validate_anonymous_data(read_panel_csv(PANEL, "FigS4A_statistical_tests.csv"))
STG <- c("Naive_B", "GC_DZ", "GC_LZ", "Memory_B", "Unassigned")
key <- c(f_naive = "Naive_B", f_dz = "GC_DZ", f_lz = "GC_LZ", f_mem = "Memory_B", f_un = "Unassigned")
kw[, stage := key[stage_fraction]]
COLLAB <- c(Naive_B = "Naive\nB", GC_DZ = "GC\nDZ", GC_LZ = "GC\nLZ", Memory_B = "Memory\nB", Unassigned = "Un-\nassigned")
cp <- cp[CJ(entity = ENT_ORD, stage = STG), on = .(entity, stage)][is.na(cells), `:=`(cells = 0L, frac = 0)]
cp[, `:=`(xi = match(stage, STG), yi = 6 - match(entity, ENT_ORD))]
cp[, lab := ifelse(frac >= 0.005, sprintf("%.0f", frac * 100), "0")]
cp[, txt := ifelse(frac > 0.42, "white", "black")]
nent <- cp[, .(n = sum(cells)), by = entity]
xlab <- paste0(COLLAB[STG], "\n", stars(kw$KW_p[match(STG, kw$stage)]))
p <- ggplot(cp, aes(xi, yi)) +
  geom_tile(aes(fill = frac * 100), colour = "white", linewidth = 1 / PT, width = 0.96, height = 0.92) +
  geom_text(aes(label = lab), colour = cp$txt, size = SZ_SMALL / PT, family = "Arial") +
  scale_fill_gradient(low = "#F4F1F8", high = "#54278F", limits = c(0, 100), guide = "none") +
  annotate("rect", xmin = 1.52, xmax = 2.48, ymin = 0.54, ymax = 5.46, fill = NA, colour = HI_COL, linewidth = LW_DATA) +
  scale_x_continuous(breaks = 1:5, labels = xlab, position = "top", expand = c(0, 0), limits = c(0.5, 5.5)) +
  scale_y_continuous(breaks = 5:1, labels = sprintf("%s (%s)", ENT_ORD, format_thousands(nent$n[match(ENT_ORD, nent$entity)])),
                     expand = c(0, 0), limits = c(0.5, 5.5)) +
  labs(x = NULL, y = NULL) +
  theme_fig4() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x.top = element_text(size = SZ_SMALL, lineheight = 0.9, margin = margin(b = SZ * 0.3)))
p <- preserve_header(p)
if (!build_only) {
  save_figure(draw_single_panel(p, "A", DS_TAG$S4A, PANEL_MM$S4A[2]), RES(PANEL), "FigureS2_PanelA_stage_composition",
         PANEL_MM$S4A[1], PANEL_MM$S4A[2],
         highlight_color = HI_COL)
}
