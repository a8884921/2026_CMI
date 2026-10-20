.release_script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
.release_root <- normalizePath(file.path(.release_script_dir, "..", ".."))
source(file.path(.release_root, "privacy.R"), encoding = "UTF-8")
PANEL_DIR <- local({ a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd() })
PANEL_DIR <- normalizePath(PANEL_DIR)
source_release_dependency(file.path(.release_root, "shared", "figure3_config.R"))
suppressMessages({ library(readr); library(ggplot2) })
imp <- validate_anonymous_data(readr::read_csv(plot_csv("Figure3D_feature_importance.csv"), show_col_types = FALSE))
stopifnot(all(c("feature", "importance") %in% names(imp)))
stopifnot(identical(as.character(imp$feature), SCORE_LEVELS))
imp$feature <- factor(imp$feature, levels = rev(SCORE_LEVELS))
pD <- ggplot(imp, aes(x = importance, y = feature)) +
  geom_col(aes(fill = importance), width = 0.62) +
  geom_text(aes(label = sprintf("%.2f", importance), colour = importance),
            hjust = -0.15, size = PT(SZ_TICK), family = FIG_FONT) +
  scale_fill_gradient(low = IMP_LOW, high = IMP_HIGH, guide = "none") +
  scale_colour_gradient(low = "#4A90C4", high = "#12456F", guide = "none") +
  scale_x_continuous("Variable Importance", limits = c(0, 1.18),
                     breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  scale_y_discrete(labels = SCORE_LABELS[rev(SCORE_LEVELS)]) +
  labs(y = NULL) +
  theme_fig3() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = ggplot2::margin(SZ_TICK * 0.5, SZ_TICK * 0.5,
                                      SZ_TICK * 0.4, SZ_TICK * 0.9))
pD <- tag_plot(pD, "D")
pD <- pD + ggplot2::theme(plot.tag.position = c(0, 0.96))
DS_D <- "GBM g1-19 \u00b7 SHAP importance \u00b7 n = 9,750 test cells"
pD <- dataset_line(pD, DS_D)
if (requireNamespace("systemfonts", quietly = TRUE)) {
  w_cap  <- systemfonts::string_width(DS_D, family = FIG_FONT, size = SZ_DATA) * 0.3528
  w_lim  <- 61
  invisible(NULL)
  stopifnot(w_cap <= w_lim)
}
save_panel(pD, "Figure3D_CMI_feature_importance", 84, 46,
           dataset = DS_D)
