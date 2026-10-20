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
PANEL <- "PanelE"
E  <- validate_anonymous_data(read_panel_csv(PANEL, "Fig4E_sample_tier_program_means.csv"))
tt <- validate_anonymous_data(read_panel_csv(PANEL, "Fig4E_paired_tests.csv"))[shown_in_E == TRUE]
if (!all(c("label", "expect", "p_31_fmt") %in% names(tt)) || !all(c("MHC_I", "MHC_II") %in% tt$prog))
  stop("Panel E requires prepared MHC-I and MHC-II plot data.")
ORD <- intersect(c("MHC_I", "MHC_II"), tt$prog)
LAB <- setNames(tt$label[match(ORD, tt$prog)], ORD)
stopifnot(setequal(unique(E$prog), ORD), nrow(tt) == length(ORD), all(tt$n == uniqueN(E$sample)))
E[, `:=`(lab = factor(LAB[prog], LAB), xi = match(tier, TIER_LV))]
tt[, lab := factor(LAB[prog], LAB)]
tt[, cnt := fifelse(expect == "down", sprintf("M3 < M1: %d/%d", down_31, n), sprintf("M3 > M1: %d/%d", up_31, n))]
rng <- E[, .(lo = min(score), hi = max(score)), by = lab][, r := hi - lo]
rng[, `:=`(yb = hi + 0.06 * r, drop = 0.035 * r, yp = hi + 0.075 * r, yc = hi + 0.165 * r, top = hi + 0.27 * r, bot = lo - 0.03 * r)]
tt <- merge(tt, rng, by = "lab")
p <- ggplot(E, aes(xi, score)) +
  geom_line(aes(group = sample), colour = "grey55", alpha = 0.20, linewidth = 0.30 / PT) +
  geom_boxplot(aes(group = xi, fill = tier), width = 0.56, outlier.shape = NA, colour = "grey20",
               linewidth = LW_ERR, fatten = 1.6, alpha = 0.95) +
  scale_fill_manual(values = TIER_COL, guide = "none") +
  geom_blank(data = rng, aes(x = 2, y = top), inherit.aes = FALSE) +
  geom_blank(data = rng, aes(x = 2, y = bot), inherit.aes = FALSE) +
  geom_segment(data = tt, aes(x = 1, xend = 3, y = yb, yend = yb), linewidth = LW_AXIS, inherit.aes = FALSE) +
  geom_segment(data = tt, aes(x = 1, xend = 1, y = yb - drop, yend = yb), linewidth = LW_AXIS, inherit.aes = FALSE) +
  geom_segment(data = tt, aes(x = 3, xend = 3, y = yb - drop, yend = yb), linewidth = LW_AXIS, inherit.aes = FALSE) +
  geom_text(data = tt, aes(x = 2, y = yp, label = sprintf("P = %s", p_31_fmt)), inherit.aes = FALSE,
            vjust = 0, size = SZ_SMALL / PT, family = "Arial") +
  geom_text(data = tt, aes(x = 2, y = yc, label = cnt), inherit.aes = FALSE, vjust = 0, colour = HI_COL,
            size = SZ / PT, family = "Arial") +
  facet_wrap(~ lab, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = 1:3, labels = c("M1", "M2", "M3"), limits = c(0.45, 3.55), expand = c(0, 0)) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 4), labels = function(v) sprintf("%.2f", v), expand = c(0, 0)) +
  labs(x = sprintf("Malignancy tier (%d paired samples)", uniqueN(E$sample)), y = "UCell score") +
  theme_fig4() +
  theme(panel.spacing.x = unit(SZ * 1.6, "pt"),
        strip.text = element_text(size = SZ, colour = "black", margin = margin(b = SZ * 0.3, t = 0)))
p <- preserve_header(p)
if (!build_only) {
  save_figure(draw_single_panel(p, "E", DS_TAG$E, PANEL_MM$E[2]), RES(PANEL), "Figure4_PanelE_immune_escape_programs",
         PANEL_MM$E[1], PANEL_MM$E[2],
         highlight_color = HI_COL)
}
