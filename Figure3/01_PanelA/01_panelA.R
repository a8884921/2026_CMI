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
suppressMessages({
  library(readr); library(dplyr); library(ComplexHeatmap); library(circlize); library(grid)
})
W_MM  <- as.numeric(Sys.getenv("FIG3_PANEL_A_W", "107"))
H_MM  <- as.numeric(Sys.getenv("FIG3_PANEL_A_H", "245"))
N_ROW <- suppressWarnings(as.integer(Sys.getenv("FIG3_PANEL_A_ROWS", "0")))
OUT_SUFFIX <- Sys.getenv("FIG3_PANEL_A_SUFFIX", "")
CELL_TEXT <- Sys.getenv("FIG3_PANEL_A_CELLTEXT", "1") == "1"
mat_df <- validate_anonymous_data(readr::read_csv(plot_csv("Figure3A_model_AUC_matrix.csv"), show_col_types = FALSE))
COHORT_LEVELS <- anonymous_cohorts(COHORT_LEVELS)
stopifnot(nrow(mat_df) == 101, all(COHORT_LEVELS %in% names(mat_df)))
stopifnot(all(mat_df$algorithm_family %in% FAMILY_LEVELS))
mat_df <- mat_df |> arrange(desc(mean_auc))
N_ALL <- nrow(mat_df)
if (!is.na(N_ROW) && N_ROW > 0 && N_ROW < N_ALL) {
  invisible(NULL)
  mat_df <- mat_df[seq_len(N_ROW), ]
}
mat      <- as.matrix(mat_df[, COHORT_LEVELS])
rownames(mat) <- mat_df$model_label
family   <- as.character(mat_df$algorithm_family)
mean_auc <- mat_df$mean_auc
col_fun <- circlize::colorRamp2(AUC_BREAKS, AUC_COLS)
top_anno <- HeatmapAnnotation(
  `Cohort type` = COHORT_TYPE[COHORT_LEVELS],
  col = list(`Cohort type` = COHORT_TYPE_COLS),
  annotation_name_gp = gpar(fontsize = SZ_SMALL, fontfamily = FIG_FONT),
  annotation_legend_param = list(
    `Cohort type` = list(
      title_gp  = gpar(fontsize = SZ_TICK, fontface = "bold", fontfamily = FIG_FONT),
      labels_gp = gpar(fontsize = SZ_TICK, fontfamily = FIG_FONT)))
)
right_anno <- rowAnnotation(
  `Mean AUC` = anno_barplot(mean_auc, baseline = 0.85, ylim = c(0.85, 0.97),
                            gp = gpar(fill = FAMILY_COLS[family], col = NA),
                            border = FALSE, bar_width = 0.85,
                            width = unit(13, "mm"),
                            axis_param = list(at = c(0.86, 0.90, 0.94),
                                              gp = gpar(fontsize = SZ_SMALL * 0.85,
                                                        fontfamily = FIG_FONT))),
  annotation_name_gp = gpar(fontsize = SZ_SMALL, fontfamily = FIG_FONT)
)
CELL_PT <- SZ_TICK * 0.56
layer_fun <- function(j, i, x, y, w, h, fill) {
  if (!CELL_TEXT) return(invisible(NULL))
  v <- pindex(mat, i, j)
  grid.text(sprintf("%.2f", v), x, y,
            gp = gpar(fontsize = CELL_PT, fontfamily = FIG_FONT,
                      col = ifelse(v > 0.965, "white", "grey15")))
}
ht <- Heatmap(
  mat, name = "AUC", col = col_fun,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = SZ_SMALL, fontfamily = FIG_FONT),
  column_names_side = "top", column_names_rot = 45,
  column_labels = c("train", "test", "Cohort_S", "Cohort_T"),
  column_names_gp = gpar(fontsize = SZ_TICK, fontface = "bold", fontfamily = FIG_FONT),
  top_annotation = top_anno, right_annotation = right_anno,
  layer_fun = layer_fun,
  rect_gp = gpar(col = "white", lwd = 0.35),
  width  = unit(if (N_ALL > 0 && !is.na(N_ROW) && N_ROW > 0 && N_ROW < N_ALL) 30 else 32, "mm"),
  height = unit(H_MM - if (!is.na(N_ROW) && N_ROW > 0 && N_ROW < N_ALL) 44 else 52, "mm"),
  heatmap_legend_param = list(
    title = "AUC", at = c(0.75, 0.85, 0.95, 1.0),
    legend_height = unit(16, "mm"),
    title_gp  = gpar(fontsize = SZ_TICK, fontface = "bold", fontfamily = FIG_FONT),
    labels_gp = gpar(fontsize = SZ_TICK, fontfamily = FIG_FONT))
)
fam_show <- c(GBM = "GBM", XGBoost = "XGBoost", StackedEnsemble = "StackedEns.",
              XRT = "XRT", DRF = "DRF", GLM = "GLM", DeepLearning = "DeepLearning")
fam_lgd <- Legend(title = "Algorithm family", labels = unname(fam_show[FAMILY_LEVELS]),
                  legend_gp = gpar(fill = FAMILY_COLS[FAMILY_LEVELS]),
                  title_gp  = gpar(fontsize = SZ_TICK, fontface = "bold", fontfamily = FIG_FONT),
                  labels_gp = gpar(fontsize = SZ_TICK, fontfamily = FIG_FONT))
DATASET_TXT <- "101 H2O AutoML models \u00b7 4 cohorts \u00b7 n = 38,995 / 9,750 internal \u00b7 n = 121,214 / 446,179 external cells"
RES_PDF <- file.path(RES_DIR, paste0("Figure3A_model_AUC_heatmap", OUT_SUFFIX, ".pdf"))
RES_SVG <- file.path(RES_DIR, paste0("Figure3A_model_AUC_heatmap", OUT_SUFFIX, ".svg"))
RES_PNG <- file.path(RES_DIR, paste0("Figure3A_model_AUC_heatmap", OUT_SUFFIX, ".png"))
draw_A <- function() {
  draw(ht, annotation_legend_list = list(fam_lgd), merge_legend = TRUE,
       heatmap_legend_side = "right", annotation_legend_side = "right",
       padding = unit(c(3, 2, 17, 2), "mm"))
  grid.text("A", x = unit(1.2, "mm"), y = unit(1, "npc") - unit(1.2, "mm"),
            just = c("left", "top"),
            gp = gpar(fontsize = SZ_TAG, fontface = "bold", fontfamily = FIG_FONT))
  grid.text(DATASET_TXT, x = unit(1.2, "mm"), y = unit(1.2, "mm"),
            just = c("left", "bottom"),
            gp = gpar(fontsize = SZ_DATA, col = SEM[["neutral"]], fontfamily = FIG_FONT))
}
open_pdf(RES_PDF, W_MM, H_MM); draw_A(); close_pdf()
svglite::svglite(RES_SVG, width = W_MM / 25.4, height = H_MM / 25.4); draw_A(); grDevices::dev.off()
grDevices::png(RES_PNG, width = W_MM / 25.4, height = H_MM / 25.4, units = "in",
               res = 600, family = FIG_FONT, type = "cairo")
draw_A(); grDevices::dev.off()
