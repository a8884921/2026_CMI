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
source_release_dependency(file.path(.release_root, "shared", "figure1_config.R"))
DATASETS_TXT <- "Analysis cohorts"
suppressMessages({ library(readr); library(dplyr); library(ggplot2) })
W_MM <- 60; H_MM <- 56
UMAP_AR <- 0.72
um <- validate_anonymous_data(readr::read_csv(P("panelB_UMAP.csv"), show_col_types = FALSE))
stopifnot(all(c("UMAP_1", "UMAP_2", "malignancy_label") %in% names(um)), nrow(um) == 48745)
n_tum <- sum(um$malignancy_label == "tumor"); n_nor <- sum(um$malignancy_label == "normal")
stopifnot(n_tum == 47968, n_nor == 777)
umap_lims <- function(x, y, ar = UMAP_AR) {
  yr <- range(y) + c(-0.03, 0.03) * diff(range(y)); xr <- range(x)
  pad <- max(0, (diff(yr) / ar - diff(xr)) / 2)
  list(x = xr + c(-pad, pad), y = yr)
}
L <- umap_lims(um$UMAP_1, um$UMAP_2)
pur <- um |> dplyr::group_by(seurat_clusters) |>
  dplyr::summarise(n_cells = dplyr::n(), n_normal = sum(malignancy_label == "normal"),
                   tumor_pct = round(100 * mean(malignancy_label == "tumor"), 2),
                   samples = paste(sort(unique(sample_id)), collapse = "+"), .groups = "drop") |>
  dplyr::arrange(tumor_pct)
readr::write_csv(pur, O("normal_cell_cluster_purity.csv"))
in_pure <- sum(pur$n_normal[pur$tumor_pct >= 95])
stopifnot(in_pure == 676, nrow(pur) == 23, sum(pur$tumor_pct >= 95) == 21)
pB <- ggplot() +
  geom_point(data = dplyr::filter(um, malignancy_label == "tumor"),
             aes(UMAP_1, UMAP_2, colour = "Tumor"), size = 0.12, stroke = 0, shape = 16) +
  geom_point(data = dplyr::filter(um, malignancy_label == "normal"),
             aes(UMAP_1, UMAP_2, colour = "Normal"), size = 0.40, stroke = 0, shape = 16) +
  scale_colour_manual(values = UMAP_GROUP_COLS, name = NULL,
                      breaks = UMAP_GROUP_LEVELS,
                      labels = sprintf("%s (%s)", UMAP_GROUP_LEVELS, format(c(n_nor, n_tum), big.mark = ",", trim = TRUE))) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 2.2))) +
  coord_equal(xlim = L$x, ylim = L$y, expand = FALSE) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_fig1() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank(),
        legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center",
        legend.key.size = grid::unit(7, "pt"), legend.key.height = grid::unit(SZ_TICK * 1.35, "pt"),
        legend.text = element_text(size = SZ_TICK, margin = ggplot2::margin(l = 2, r = SZ_TICK * 0.7, unit = "pt")),
        legend.margin = ggplot2::margin(0, 0, 0, 0), legend.box.margin = ggplot2::margin(-SZ_TICK * 0.4, 0, 0, 0),
        legend.spacing.x = grid::unit(1, "pt"), legend.spacing.y = grid::unit(0, "pt"),
        plot.margin = ggplot2::margin(SZ_TICK * 1.6, SZ_TICK * 0.3, SZ_TICK * 0.3, SZ_TICK * 0.3))
pB <- PANEL_TAG(pB, "B")
pB <- dataset_line(pB, sprintf("Datasets: 6 studies \u00b7 %s B cells", format(nrow(um), big.mark = ",")))
save_panel(pB, W_MM, H_MM)
