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
PANEL <- "PanelD"
D  <- validate_anonymous_data(read_panel_csv(PANEL, "Fig4D_sample_stage_median_dCMI.csv"))[include == TRUE]
tt <- validate_anonymous_data(read_panel_csv(PANEL, "Fig4D_within_stage_tests.csv"))[plotted == TRUE]
mm <- validate_anonymous_data(read_panel_csv(PANEL, "Fig4D_stage_adjusted_model.csv"))
gm <- validate_anonymous_data(read_panel_csv(PANEL, "Fig4D_group_medians.csv"))
STG <- c("Naive_B", "GC_DZ", "Memory_B")
stopifnot(setequal(unique(D$stage), STG), setequal(tt$stage, STG), max(D$median_dCMI) < CD_YLIM[2] - 0.45)
D[, ent_i := match(entity, ENT_ORD)]
D[, slot := frank(ent_i, ties.method = "dense"), by = stage]
D[, n_se := .N, by = .(stage, entity)]
keys <- unique(D[, .(stage, slot, entity, n_se)])[order(match(stage, STG), slot)]
keys[, key := sprintf("%s|%d", stage, slot)]
D[, key := factor(sprintf("%s|%d", stage, slot), keys$key)]
xl <- setNames(sprintf("%s\n(%d)", keys$entity, keys$n_se), keys$key)
XR <- D[, .(xlo = min(slot), xhi = max(slot)), by = .(stage, grp = g3)]
gm <- merge(gm, XR, by = c("stage", "grp"))
stair <- rbindlist(lapply(STG, function(s) step_data(gm[stage == s])[, stage := s]))
lab_stage <- function(v) factor(v, STG, labels = STAGE_LAB[STG])
D[, stage := lab_stage(stage)]; stair[, stage := lab_stage(stage)]; tt[, stage := lab_stage(stage)]
p <- ggplot(D, aes(key, median_dCMI)) +
  geom_hline(yintercept = 0, linetype = "22", linewidth = LW_AXIS, colour = "grey45") +
  geom_boxplot(width = 0.62, outlier.shape = NA, fill = "#E7EAEE", colour = "grey35", linewidth = LW_ERR, fatten = 1.6) +
  geom_point(position = position_jitter(width = 0.17, height = 0, seed = 20261011), size = 0.55, alpha = 0.55,
             colour = "grey30", stroke = 0) +
  geom_segment(data = stair, aes(x = x, xend = xend, y = y, yend = yend), colour = HI_COL, linewidth = LW_DATA,
               lineend = "square", inherit.aes = FALSE) +
  geom_text(data = tt, aes(x = -Inf, y = CD_YLIM[2], label = sprintf("P = %s", fmtP(trend_p))), inherit.aes = FALSE,
            hjust = -0.08, vjust = 1.25, colour = HI_COL, size = SZ_SMALL / PT, family = "Arial") +
  facet_grid(. ~ stage, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = xl) +
  scale_y_continuous(breaks = 0:3, limits = CD_YLIM, expand = c(0, 0)) +
  labs(x = "Lymphoma (n samples)", y = "Median \u0394CMI per sample",
       subtitle = sprintf("All stages, stage-adjusted: P = %s", fmtP(mm$p))) +
  theme_fig4() +
  theme(plot.subtitle = element_text(size = SZ_SMALL, colour = HI_COL, hjust = 0.5, family = "Arial",
                                     margin = margin(b = SZ * 0.25)),
        panel.spacing.x = unit(SZ * 0.9, "pt"),
        strip.text = element_text(size = SZ, colour = "black", margin = margin(b = SZ * 0.3, t = 0)))
p <- preserve_header(p)
if (!build_only) {
  save_figure(draw_single_panel(p, "D", DS_TAG$D, PANEL_MM$D[2]), RES(PANEL), "Figure4_PanelD_stage_gradient",
         PANEL_MM$D[1], PANEL_MM$D[2])
}
