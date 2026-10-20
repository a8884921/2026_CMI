.release_script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
.release_root <- normalizePath(file.path(.release_script_dir, "..", ".."))
source(file.path(.release_root, "privacy.R"), encoding = "UTF-8")
.script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE); f <- gsub("~+~", " ", sub("^--file=", "", a[grep("^--file=", a)]), fixed = TRUE)
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
PANEL_DIR <- normalizePath(.script_dir)
source_release_dependency(file.path(.release_root, "shared", "figure2_config.R"))
DATASETS_TXT <- "Analysis cohorts"
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })
source_release_dependency(file.path(.release_root, "shared", "panelBC_layout.R"))
W_MM <- 92; H_MM <- 120
imp <- validate_anonymous_data(utils::read.csv(P("panelC_feature_importance.csv"), check.names = FALSE))
ord <- validate_anonymous_data(utils::read.csv(P("panelBC_row_order.csv"), check.names = FALSE))
stopifnot(setequal(imp$feature, FEATURE24), abs(sum(imp$pct) - 100) < 1e-6, sum(imp$top4) == 4)
rows <- bc_rows(ord)
pLab <- ggplot() + bc_bands(rows) +
  geom_text(data = rows$headers, aes(x = 0.02, y = y, label = label), hjust = 0, fontface = "bold",
            size = SZ_TICK / PT, family = FIG_FONT) +
  geom_text(data = rows$metrics, aes(x = 0.06, y = y, label = disp), hjust = 0, size = SZ_TICK / PT, family = FIG_FONT) +
  annotate("text", x = 0.02, y = rows$top, label = "Ig metric", hjust = 0, size = SZ_TICK / PT,
           family = FIG_FONT, colour = GREY_DARK) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + bc_y(rows) + bc_theme() + bc_blank_x() +
  theme(plot.margin = ggplot2::margin(SZ_TICK * 1.6, 1.5 * 72.27 / 25.4, 0, 0))
pC <- (pLab | bc_importance_col(rows, imp)) + patchwork::plot_layout(widths = grid::unit(c(36, 50), "mm"))
pC <- PANEL_TAG(pC, "C")
pC <- dataset_line(pC, sprintf("Datasets: GBM training set \u00b7 %s B cells", format(unique(imp$n_train), big.mark = ",")),
                   pw_theme = theme(plot.margin = ggplot2::margin(0, 0, 0, 0)))
STARS <- bc_star_pos(pC, rows, "panel-1", 36, W_MM, H_MM)
save_panel(pC, W_MM, H_MM, draw_fun = function() { print(pC); draw_stars(STARS) })
out <- imp |> dplyr::transmute(rank = rank, feature = feature, category = category,
                               importance_unnormalized = round(importance, 2), percentage = round(pct, 2), top4_driver = top4,
                               model_training_cells = n_train)
utils::write.csv(out, O("GBM_feature_importance.csv"), row.names = FALSE)
