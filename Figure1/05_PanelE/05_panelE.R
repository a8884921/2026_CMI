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
suppressMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
W_MM <- 112; H_MM <- 62
df <- validate_anonymous_data(readr::read_csv(P("panelE_five_score_matrix.csv"), show_col_types = FALSE))
stopifnot(all(SCORE_COLS %in% names(df)), "malignancy_label" %in% names(df))
dfz <- df |>
  dplyr::filter(malignancy_label %in% c("normal", "tumor")) |>
  dplyr::mutate(group = factor(ifelse(malignancy_label == "tumor", "Tumor B cell", "Normal B cell"),
                               levels = VIOLIN_LEVELS),
                dplyr::across(dplyr::all_of(SCORE_COLS), ~ as.numeric(scale(.x))))
n_tum <- sum(dfz$group == "Tumor B cell"); n_nor <- sum(dfz$group == "Normal B cell")
stopifnot(nrow(dfz) == 48745, n_tum == 47968, n_nor == 777)
long <- dfz |>
  dplyr::select(group, dplyr::all_of(SCORE_COLS)) |>
  tidyr::pivot_longer(dplyr::all_of(SCORE_COLS), names_to = "score", values_to = "z") |>
  dplyr::filter(is.finite(z)) |>
  dplyr::mutate(score = factor(score, levels = SCORE_COLS))
q <- long |>
  dplyr::group_by(score, group) |>
  dplyr::summarise(p25 = quantile(z, .25), p50 = median(z), p75 = quantile(z, .75), .groups = "drop") |>
  dplyr::mutate(x = as.numeric(score) + ifelse(group == "Normal B cell", -0.20, 0.20))
arr <- q |>
  dplyr::select(score, group, p50, x) |>
  tidyr::pivot_wider(names_from = group, values_from = c(p50, x), names_sep = "|") |>
  dplyr::rename(y0 = `p50|Normal B cell`, y1 = `p50|Tumor B cell`,
                x0 = `x|Normal B cell`,  x1 = `x|Tumor B cell`)
tum_v <- function(s) dfz[[s]][dfz$group == "Tumor B cell"]
nor_v <- function(s) dfz[[s]][dfz$group == "Normal B cell"]
p_raw <- sapply(SCORE_COLS, function(s) stats::wilcox.test(tum_v(s), nor_v(s))$p.value)
p_bh  <- stats::p.adjust(p_raw, method = "BH")
cliff <- sapply(SCORE_COLS, function(s) {
  U <- as.numeric(stats::wilcox.test(tum_v(s), nor_v(s), exact = FALSE)$statistic)
  2 * U / (n_tum * n_nor) - 1
})
med_diff <- arr$y1[match(SCORE_COLS, arr$score)] - arr$y0[match(SCORE_COLS, arr$score)]
stopifnot(all(med_diff > 0), all(p_bh < 0.05))
stat_tab <- data.frame(score = lab_of(SCORE_COLS), score_name = SCORE_LONG[SCORE_COLS],
                       n_tumor_cells = n_tum, n_normal_cells = n_nor,
                       z_median_tumor = round(arr$y1[match(SCORE_COLS, arr$score)], 3),
                       z_median_normal = round(arr$y0[match(SCORE_COLS, arr$score)], 3),
                       median_difference = round(med_diff, 3), Cliff_delta = round(cliff, 3),
                       P_raw = fmtP(p_raw), P_BH_adjusted = fmtP(p_bh),
                       test = "Wilcoxon rank-sum (two-sided, unpaired)", multiple_comparisons = "Benjamini-Hochberg (5 tests)",
                       check.names = FALSE)
pub <- validate_anonymous_data(readr::read_csv(P("panelE_p_value_annotations.csv"), show_col_types = FALSE))
stat_tab$previous_annotation <- pub$p_label[match(SCORE_COLS, pub$score)]
readr::write_csv(stat_tab, O("five_score_group_comparisons.csv"))
lab <- data.frame(score = factor(SCORE_COLS, levels = SCORE_COLS), y = 3.62, label = P_txt(p_bh))
ev <- df |> dplyr::group_by(sample_id) |>
  dplyr::summarise(k = dplyr::n_distinct(malignancy_label), .groups = "drop") |>
  dplyr::filter(k == 2) |> dplyr::pull(sample_id)
per_s <- do.call(rbind, lapply(SCORE_COLS, function(s) {
  m <- df |> dplyr::filter(sample_id %in% ev) |> dplyr::group_by(sample_id) |>
    dplyr::summarise(t = median(.data[[s]][malignancy_label == "tumor"]),
                     n = median(.data[[s]][malignancy_label == "normal"]), .groups = "drop")
  d <- m$t - m$n
  w <- suppressWarnings(stats::wilcox.test(m$t, m$n, paired = TRUE))
  data.frame(dimension = lab_of(s), n_samples = nrow(m), n_tumor_higher = sum(d > 0), n_normal_higher = sum(d < 0), n_ties = sum(d == 0),
             P_sample_paired = fmtP(w$p.value), check.names = FALSE)
}))
readr::write_csv(per_s, O("sample_paired_tests.csv"))
print(per_s, row.names = FALSE)
stopifnot(length(ev) == 17)
Z_SHOW <- 3.2
long_show <- dplyr::filter(long, abs(z) <= Z_SHOW)
pE <- ggplot(long_show, aes(x = score, y = z, fill = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey75", linewidth = LW_AXIS * 0.8) +
  geom_violin(position = position_dodge(0.80), width = 0.78, linewidth = 0.25,
              colour = "grey30", scale = "width", trim = TRUE) +
  geom_segment(data = q, aes(x = x, xend = x, y = p25, yend = p75),
               linewidth = 0.4, colour = "grey10", inherit.aes = FALSE) +
  geom_segment(data = arr, aes(x = x0, xend = x1, y = y0, yend = y1), inherit.aes = FALSE,
               colour = "grey15", linewidth = 0.55,
               arrow = grid::arrow(length = grid::unit(1.3, "mm"), type = "closed")) +
  geom_point(data = q, aes(x = x, y = p50), size = 1.1, shape = 21, fill = "white",
             colour = "grey10", stroke = 0.35, inherit.aes = FALSE) +
  geom_text(data = lab, aes(x = score, y = y, label = label), inherit.aes = FALSE,
            size = SZ_TICK / PT, family = FIG_FONT) +
  scale_fill_manual(values = VIOLIN_COLS, name = NULL, breaks = VIOLIN_LEVELS,
                    labels = sprintf("%s (n = %s)", VIOLIN_LEVELS, format(c(n_nor, n_tum), big.mark = ",", trim = TRUE))) +
  scale_x_discrete(labels = SCORE_AXIS_LABELS) +
  scale_y_continuous("Standardized score (z)", limits = c(-Z_SHOW, 3.9), breaks = c(-2, 0, 2)) +
  labs(x = NULL) +
  theme_fig1() +
  theme(legend.position = "bottom", legend.justification = "center",
        legend.direction = "horizontal", legend.key.size = grid::unit(8, "pt"),
        legend.text = element_text(size = SZ_TICK), legend.margin = ggplot2::margin(0, 0, 0, 0),
        legend.box.margin = ggplot2::margin(-SZ_TICK * 0.3, 0, 0, 0),
        plot.margin = ggplot2::margin(SZ_TICK * 0.9, SZ_TICK * 0.6, SZ_TICK * 0.4, SZ_TICK * 0.4))
pE <- PANEL_TAG(pE, "E")
pE <- dataset_line(pE, sprintf("Datasets: 6 studies \u00b7 %s B cells", format(nrow(dfz), big.mark = ",")))
save_panel(pE, W_MM, H_MM)
print(stat_tab[, c("score", "median_difference", "Cliff_delta", "P_raw", "P_BH_adjusted")], row.names = FALSE)
