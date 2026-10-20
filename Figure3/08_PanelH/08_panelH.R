.release_script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
.release_root <- normalizePath(file.path(.release_script_dir, "..", ".."))
source(file.path(.release_root, "privacy.R"), encoding = "UTF-8")
PANEL_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
PANEL_DIR <- normalizePath(PANEL_DIR)
source_release_dependency(file.path(.release_root, "shared", "figure3_config.R"))
suppressMessages({ library(readr); library(dplyr); library(ggplot2); library(patchwork) })
auc <- validate_anonymous_data(readr::read_csv(plot_csv("Figure3H_cohort_AUC.csv"), show_col_types = FALSE))
DS_LEVELS <- anonymous_cohorts(DS_LEVELS)
DS_LABELS <- DS_LEVELS
stopifnot(all(DS_LEVELS %in% unique(auc$dataset)))
ROW_LEVELS <- c("prol_score", "cnv_score", "mtp_score", "diff_score", "Ig_score", "CMI")
fst <- auc |>
  dplyr::mutate(dataset = factor(dataset, levels = DS_LEVELS),
                method  = factor(method,  levels = ROW_LEVELS),
                is_cmi  = ifelse(method == "CMI", "CMI", "single"))
n_lab <- auc |> dplyr::distinct(dataset, n_cells, n_tumor) |>
  dplyr::mutate(dataset = factor(dataset, levels = DS_LEVELS),
                lab = sprintf("n = %s\n(%s tumor)",
                              format(n_cells, big.mark = ","),
                              format(n_tumor, big.mark = ",")))
W_MM <- 205; H_MM <- 60
bg <- data.frame(dataset = factor(DS_LEVELS, levels = DS_LEVELS))
pH <- ggplot(fst, aes(AUC, method)) +
  geom_rect(data = bg, aes(fill = dataset), xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha = 0.55) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 5.5, ymax = 6.5,
           fill = SEM[["case"]], alpha = 0.10) +
  geom_vline(xintercept = seq(0.5, 1.0, 0.1), colour = "white", linewidth = LW_ERR) +
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = is_cmi),
                 height = 0.26, linewidth = LW_ERR) +
  geom_point(aes(colour = is_cmi, size = is_cmi)) +
  geom_text(aes(label = sprintf("%.2f", AUC), colour = is_cmi),
            vjust = -1.05, size = PT(SZ_TICK), family = FIG_FONT,
            fontface = "bold", show.legend = FALSE) +
  geom_text(data = n_lab, aes(x = 0.55, y = 0.55, label = lab),
            inherit.aes = FALSE, hjust = 0, vjust = 0,
            size = PT(SZ_SMALL), family = FIG_FONT, colour = SEM[["neutral"]],
            lineheight = 0.95) +
  facet_wrap(~dataset, nrow = 1,
             labeller = as_labeller(setNames(DS_LABELS, DS_LEVELS))) +
  scale_fill_manual(values = PANEL_FILL, guide = "none") +
  scale_colour_manual(values = c(CMI = SEM[["case"]], single = "grey35"), guide = "none") +
  scale_size_manual(values = c(CMI = 2.4, single = 1.5), guide = "none") +
  scale_y_discrete(labels = c(prol_score = "Prol", cnv_score = "CNV", mtp_score = "MTP",
                              diff_score = "Diff", Ig_score = "Ig", CMI = "CMI")) +
  scale_x_continuous("AUC (95% CI)", limits = c(0.5, 1.04),
                     breaks = c(0.5, 0.9), expand = c(0, 0)) +
  labs(y = NULL) +
  theme_fig3() +
  theme(panel.grid       = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text       = element_text(size = SZ_TICK, face = "bold", family = FIG_FONT),
        panel.spacing   = grid::unit(2.2, "pt"),
        axis.text.x     = element_text(size = SZ_SMALL),
        plot.margin     = ggplot2::margin(SZ_TICK * 0.5, SZ_TICK * 0.5,
                                          SZ_TICK * 0.4, SZ_TICK * 0.4))
recolor_strips <- function(p, fills) {
  g  <- ggplot2::ggplotGrob(p)
  st <- which(grepl("strip-t", g$layout$name))
  st <- st[order(g$layout$l[st])]
  for (i in seq_along(st)) {
    rect <- g$grobs[[st[i]]]$grobs[[1]]$children[[1]]
    rect$gp$fill <- fills[i]
    g$grobs[[st[i]]]$grobs[[1]]$children[[1]] <- rect
  }
  g
}
pH <- tag_plot(pH, "H")
gH <- recolor_strips(pH, unname(STRIP_FILL[DS_LEVELS]))
gH_plot <- patchwork::wrap_elements(gH)
DS_TXT <- "4 cohorts \u00b7 internal train 38,995 / test 9,750 cells \u00b7 Cohort_S 121,214 \u00b7 Cohort_T 446,179 cells (both external cohorts include the 47,968 internal tumour cells)"
gH_plot <- dataset_line(gH_plot, DS_TXT)
save_panel(gH_plot, "Figure3H_cohort_AUC_forest", W_MM, H_MM,
  dataset = DS_TXT)
for (dd in DS_LEVELS) {
  sub <- auc[auc$dataset == dd & auc$method != "CMI", ]
  b   <- sub$method[which.max(sub$AUC)]
  cmi <- auc$AUC[auc$dataset == dd & auc$method == "CMI"]
  n   <- auc$n_cells[auc$dataset == dd][1]
  invisible(NULL)
}
