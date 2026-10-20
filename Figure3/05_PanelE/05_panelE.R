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
suppressMessages({ library(readr); library(dplyr); library(ggplot2) })
d <- validate_anonymous_data(readr::read_csv(plot_csv("Figure3E_SHAP_beeswarm.csv"), show_col_types = FALSE))
FEATS <- SCORE_LEVELS
SHAP_C <- paste0("shap_", FEATS)
VAL_C  <- paste0("value_", FEATS)
stopifnot(all(c(SHAP_C, VAL_C) %in% names(d)))
n_cell <- nrow(d)
n_tumor <- sum(d$GROUP == "tumor")
set.seed(20260614)
long <- do.call(rbind, lapply(seq_along(FEATS), function(k) {
  data.frame(feature = FEATS[k], y = k,
             shap    = d[[SHAP_C[k]]],
             value   = d[[VAL_C[k]]])
}))
long$feature <- factor(long$feature, levels = rev(FEATS))
long$value_rank <- ave(long$value, long$feature,
                       FUN = function(v) (rank(v, ties.method = "average") - 1) / (length(v) - 1))
beeswarm_offset <- function(x, nbins = 220) {
  rng <- range(x, na.rm = TRUE)
  bw  <- diff(rng) / nbins
  bin <- round((x - rng[1]) / bw)
  off <- ave(seq_along(bin), bin, FUN = function(i) {
    n <- length(i); if (n <= 1) return(0)
    seq(-0.5, 0.5, length.out = n)
  })
  list(bin = bin, off = as.numeric(off))
}
bs <- beeswarm_offset(long$shap)
long$ypos <- long$y + bs$off * 0.62
W_MM <- 80; H_MM <- 62
pE <- ggplot(long, aes(shap, ypos)) +
  geom_vline(xintercept = 0, colour = "#D9D9D9", linewidth = LW_AXIS) +
  geom_hline(yintercept = seq_along(FEATS), colour = "#EFEFEF", linewidth = LW_AXIS) +
  geom_point(aes(colour = value_rank), size = 0.30, alpha = 0.75, shape = 16) +
  scale_colour_gradientn(colours = SHAP_COLS, limits = c(0, 1), guide = "none") +
  scale_y_continuous(breaks = seq_along(FEATS), labels = rev(SCORE_LABELS[FEATS]),
                     limits = c(0.4, length(FEATS) + 0.6), expand = c(0, 0)) +
  scale_x_continuous("SHAP value", breaks = c(-1, 0, 1)) +
  labs(y = NULL) +
  theme_fig3() +
  theme(axis.line.y  = element_blank(), axis.ticks.y = element_blank(),
        plot.margin  = ggplot2::margin(SZ_TICK * 0.5, SZ_TICK * 0.5,
                                       SZ_TICK * 0.4, SZ_TICK * 0.4))
pE <- tag_plot(pE, "E")
DS_TXT <- "In-house mature B-cell lymphoma scRNA-seq \u00b7 internal test set\nn = 9,750 cells (9,594 tumour)"
pE <- dataset_line(pE, DS_TXT)
save_panel(pE, "Figure3E_SHAP_beeswarm", W_MM, H_MM,
  dataset = DS_TXT)
for (f in FEATS) {
  s <- d[[paste0("shap_", f)]]; v <- d[[paste0("value_", f)]]
  invisible(NULL)
}
