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
suppressMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2); library(pROC); library(patchwork) })
W_MM <- 90; H_MM <- 68
SCORE2METHOD <- c(Ig_score = "Ig", cnv_score = "CNV", mtp_score = "MTP",
                  prol_score = "Prol", diff_score = "Diff")
METHODS <- unname(SCORE2METHOD[SCORE_COLS])
df  <- validate_anonymous_data(readr::read_csv(P("panelG_five_score_matrix.csv"), show_col_types = FALSE))
thr <- validate_anonymous_data(readr::read_csv(P("panelG_patient_Youden_thresholds.csv"), show_col_types = FALSE))
stopifnot(all(c(SCORE_COLS, "sample_id", "patient_id") %in% names(df)), all(c("sample_id", "method", "threshold") %in% names(thr)),
          all(unique(thr$method) %in% SCORE_COLS))
evaluable <- df |> dplyr::group_by(sample_id) |>
  dplyr::summarise(n_cls = dplyr::n_distinct(malignancy_label), .groups = "drop") |>
  dplyr::filter(n_cls >= 2) |> dplyr::pull(sample_id)
pool <- dplyr::filter(df, sample_id %in% evaluable)
youden_thr <- function(d, m)
  as.numeric(pROC::coords(pROC::roc(d$malignancy_label, d[[m]], levels = c("normal", "tumor"),
                                    direction = "<", quiet = TRUE),
                          "best", best.method = "youden", ret = "threshold")[1, 1])
global_thr <- sapply(SCORE_COLS, function(m) youden_thr(pool, m))
thr_wide <- thr |> dplyr::select(sample_id, method, threshold) |>
  tidyr::pivot_wider(names_from = method, values_from = threshold, names_glue = "{method}_thr")
stopifnot(all(paste0(SCORE_COLS, "_thr") %in% names(thr_wide)))
d2 <- df |>
  dplyr::left_join(thr_wide, by = "sample_id") |>
  dplyr::mutate(dplyr::across(dplyr::ends_with("_thr"),
                              ~ dplyr::coalesce(.x, global_thr[sub("_thr$", "", dplyr::cur_column())])))
calls <- d2 |>
  dplyr::transmute(sample_id, malignancy_label,
                   Ig = Ig_score > Ig_score_thr, CNV = cnv_score > cnv_score_thr,
                   MTP = mtp_score > mtp_score_thr, Prol = prol_score > prol_score_thr,
                   Diff = diff_score > diff_score_thr) |>
  dplyr::mutate(k = Ig + CNV + MTP + Prol + Diff,
                group = factor(ifelse(malignancy_label == "tumor", "Tumor B cell", "Normal B cell"),
                               levels = VIOLIN_LEVELS))
stopifnot(!any(is.na(calls[, METHODS])), nrow(calls) == 48745)
n_tum <- sum(calls$group == "Tumor B cell"); n_nor <- sum(calls$group == "Normal B cell")
dist <- calls |> dplyr::count(group, k) |>
  tidyr::complete(group, k = 0:5, fill = list(n = 0)) |>
  dplyr::group_by(group) |> dplyr::mutate(pct = 100 * n / sum(n)) |> dplyr::ungroup()
readr::write_csv(dist, O("abnormal_dimension_count_distribution.csv"))
ge3_t <- sum(calls$k[calls$group == "Tumor B cell"] >= 3)
ge3_n <- sum(calls$k[calls$group == "Normal B cell"] >= 3)
fis   <- stats::fisher.test(matrix(c(ge3_t, n_tum - ge3_t, ge3_n, n_nor - ge3_n), 2))
wil   <- stats::wilcox.test(calls$k[calls$group == "Tumor B cell"], calls$k[calls$group == "Normal B cell"])
pct_t <- 100 * ge3_t / n_tum; pct_n <- 100 * ge3_n / n_nor
per_dim <- do.call(rbind, lapply(METHODS, function(m) {
  a <- sum(calls[[m]][calls$group == "Tumor B cell"]); c <- sum(calls[[m]][calls$group == "Normal B cell"])
  ft <- stats::fisher.test(matrix(c(a, n_tum - a, c, n_nor - c), 2))
  data.frame(dimension = m, tumor_positive_pct = round(100 * a / n_tum, 2), normal_positive_pct = round(100 * c / n_nor, 2),
             OR = signif(unname(ft$estimate), 3), P_Fisher = fmtP(ft$p.value))
}))
readr::write_csv(per_dim, O("dimension_positive_rates.csv"))
tests <- data.frame(
  comparison = c("tumor_vs_normal_ge3_abnormal_dimensions", "tumor_vs_normal_abnormal_dimension_distribution", "tumor_vs_normal_all_dimensions_negative"),
  tumor = c(sprintf("%.1f%% (%s / %s)", pct_t, format(ge3_t, big.mark = ","), format(n_tum, big.mark = ",")),
           sprintf("mean %.2f", mean(calls$k[calls$group == "Tumor B cell"])),
           sprintf("%.2f%%", dist$pct[dist$group == "Tumor B cell" & dist$k == 0])),
  normal = c(sprintf("%.1f%% (%d / %d)", pct_n, ge3_n, n_nor),
           sprintf("mean %.2f", mean(calls$k[calls$group == "Normal B cell"])),
           sprintf("%.2f%%", dist$pct[dist$group == "Normal B cell" & dist$k == 0])),
  test = c(sprintf("Fisher exact test, OR = %.2f", unname(fis$estimate)), "Wilcoxon rank-sum test", "\u2014"),
  P = c(fmtP(fis$p.value), fmtP(wil$p.value), ""), check.names = FALSE)
readr::write_csv(tests, O("statistical_tests.csv"))
k_by <- function(d, thr_vec) Reduce(`+`, lapply(SCORE_COLS, function(m) d[[m]] > thr_vec[[m]]))
paired_test <- function(per) {
  w <- stats::wilcox.test(per$tumor, per$normal, paired = TRUE, exact = TRUE)
  list(n = nrow(per), up = sum(per$tumor > per$normal), p = w$p.value,
       med_t = median(per$tumor), med_n = median(per$normal))
}
per_sample <- calls |> dplyr::filter(sample_id %in% evaluable) |>
  dplyr::group_by(sample_id) |>
  dplyr::summarise(n_normal = sum(group == "Normal B cell"), n_tumor = sum(group == "Tumor B cell"),
                   normal = mean(k[group == "Normal B cell"] >= 3), tumor = mean(k[group == "Tumor B cell"] >= 3),
                   .groups = "drop") |>
  dplyr::mutate(patient_id = df$patient_id[match(sample_id, df$sample_id)],
                subtype = df$subtype[match(sample_id, df$sample_id)])
stopifnot(nrow(per_sample) == 17, !anyDuplicated(per_sample$patient_id))
pt_main <- paired_test(per_sample)
pt_ge10 <- paired_test(dplyr::filter(per_sample, n_normal >= 10))
readr::write_csv(dplyr::mutate(per_sample, normal_pct = round(100 * normal, 1), tumor_pct = round(100 * tumor, 1)),
                 O("sample_ge3_dimension_proportions.csv"))
df$patient_key <- df$patient_id
pool$patient_key <- pool$patient_id
k_lopo <- integer(nrow(df))
for (pk in unique(df$patient_key)) {
  idx <- which(df$patient_key == pk)
  thr_o <- sapply(SCORE_COLS, function(m) youden_thr(dplyr::filter(pool, patient_key != pk), m))
  k_lopo[idx] <- k_by(df[idx, ], thr_o)
}
k_glob <- k_by(df, global_thr)
rob_one <- function(kv, label) {
  is_t <- df$malignancy_label == "tumor"
  a <- sum(kv[is_t] >= 3); b <- sum(kv[!is_t] >= 3)
  ft <- stats::fisher.test(matrix(c(a, sum(is_t) - a, b, sum(!is_t) - b), 2))
  per <- data.frame(sample_id = df$sample_id, t = is_t, k = kv) |>
    dplyr::filter(sample_id %in% evaluable) |> dplyr::group_by(sample_id) |>
    dplyr::summarise(normal = mean(k[!t] >= 3), tumor = mean(k[t] >= 3), .groups = "drop")
  pp <- paired_test(per)
  data.frame(analysis = label, tumor_ge3 = sprintf("%.1f%%", 100 * a / sum(is_t)),
             normal_ge3 = sprintf("%.1f%%", 100 * b / sum(!is_t)), OR = round(unname(ft$estimate), 2),
             P_cell_Fisher = fmtP(ft$p.value),
             samples_tumor_higher = sprintf("%d / %d", pp$up, pp$n), P_sample_paired_Wilcoxon = fmtP(pp$p),
             check.names = FALSE)
}
robust <- rbind(
  data.frame(analysis = "patient_Youden_thresholds_with_pooled_fallback",
             tumor_ge3 = sprintf("%.1f%%", pct_t), normal_ge3 = sprintf("%.1f%%", pct_n),
             OR = round(unname(fis$estimate), 2), P_cell_Fisher = fmtP(fis$p.value),
             samples_tumor_higher = sprintf("%d / %d", pt_main$up, pt_main$n), P_sample_paired_Wilcoxon = fmtP(pt_main$p),
             check.names = FALSE),
  data.frame(analysis = "paired_samples_normal_ge10", tumor_ge3 = "\u2014", normal_ge3 = "\u2014", OR = NA,
             P_cell_Fisher = "\u2014", samples_tumor_higher = sprintf("%d / %d", pt_ge10$up, pt_ge10$n),
             P_sample_paired_Wilcoxon = fmtP(pt_ge10$p), check.names = FALSE),
  rob_one(k_lopo, "leave_one_patient_out_thresholds"),
  rob_one(k_glob, "pooled_Youden_thresholds"))
readr::write_csv(robust, O("robustness_tests.csv"))
print(robust, row.names = FALSE)
dist$kf <- factor(dist$k, levels = 0:5)
ymax  <- max(dist$pct)
y_top <- ymax * 1.62
lh    <- y_top * 0.088
ann <- data.frame(x = 5, y = y_top - lh * (0:3) - lh * 0.15,
                  label = c("\u22653 dimensions", sprintf("Tumor %.1f%%", pct_t),
                            sprintf("Normal %.1f%%", pct_n), P_txt(fis$p.value)),
                  col = c("black", VIOLIN_COLS[["Tumor B cell"]], VIOLIN_COLS[["Normal B cell"]], "black"))
pBar <- ggplot(dist, aes(kf, pct, fill = group)) +
  geom_blank() +
  annotate("rect", xmin = 3.5, xmax = 6.5, ymin = -Inf, ymax = Inf, fill = "#FBE3DF", colour = NA) +
  geom_col(position = position_dodge(0.78), width = 0.74) +
  geom_text(data = ann, aes(x = x, y = y, label = label), inherit.aes = FALSE, colour = ann$col,
            hjust = 0.5, vjust = 1, family = FIG_FONT, size = SZ_TICK / PT) +
  scale_fill_manual(values = VIOLIN_COLS, name = NULL, breaks = VIOLIN_LEVELS,
                    labels = sprintf("%s (n = %s)", VIOLIN_LEVELS,
                                     format(c(n_nor, n_tum), big.mark = ",", trim = TRUE))) +
  scale_y_continuous("Cells (%)", limits = c(0, y_top), breaks = c(0, 10, 20, 30),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Abnormal dimensions") +
  theme_fig1() +
  theme(legend.position = "bottom", legend.justification = "center",
        legend.direction = "horizontal", legend.key.size = grid::unit(8, "pt"),
        legend.text = element_text(size = SZ_TICK), legend.margin = ggplot2::margin(0, 0, 0, 0),
        plot.margin = ggplot2::margin(SZ_TICK * 0.9, SZ_TICK * 0.3, SZ_TICK * 0.2, SZ_TICK * 0.4))
ps_long <- per_sample |>
  tidyr::pivot_longer(c(normal, tumor), names_to = "grp", values_to = "v") |>
  dplyr::mutate(grp = factor(ifelse(grp == "normal", "Normal", "Tumor"), levels = c("Normal", "Tumor")),
                v = 100 * v)
ps_long <- ps_long |> dplyr::group_by(grp, v) |>
  dplyr::mutate(dx = if (dplyr::n() > 1) (seq_len(dplyr::n()) - (dplyr::n() + 1) / 2) * 0.035 else 0) |>
  dplyr::ungroup() |>
  dplyr::mutate(xpos = as.numeric(grp) + dx)
med_df <- data.frame(grp = factor(c("Normal", "Tumor"), levels = c("Normal", "Tumor")),
                     v = 100 * c(pt_main$med_n, pt_main$med_t))
med_df$xpos <- as.numeric(med_df$grp)
pSlope <- ggplot(ps_long, aes(xpos, v)) +
  geom_line(aes(group = sample_id), colour = "#B8B8B8", linewidth = LW_AXIS * 0.8) +
  geom_point(aes(colour = grp), size = 1.0, shape = 16) +
  geom_line(data = med_df, aes(group = 1), colour = "black", linewidth = LW_DATA) +
  geom_point(data = med_df, colour = "black", size = 1.6, shape = 18) +
  annotate("text", x = 1.5, y = 112, label = P_txt(pt_main$p), family = FIG_FONT,
           size = SZ_TICK / PT, vjust = 0) +
  scale_colour_manual(values = c(Normal = VIOLIN_COLS[["Normal B cell"]], Tumor = VIOLIN_COLS[["Tumor B cell"]]),
                      guide = "none") +
  scale_x_continuous(breaks = 1:2, labels = c("Normal", "Tumor"), limits = c(0.62, 2.38), expand = c(0, 0)) +
  scale_y_continuous("Cells \u22653 abnormal (%)", limits = c(0, 122), breaks = c(0, 50, 100),
                     expand = expansion(mult = c(0.02, 0))) +
  labs(x = sprintf("%d samples", pt_main$n)) +
  coord_cartesian(clip = "off") +
  theme_fig1() +
  theme(plot.margin = ggplot2::margin(SZ_TICK * 0.9, SZ_TICK * 0.3, SZ_TICK * 0.2, SZ_TICK * 0.6))
pG <- (pBar | pSlope) + patchwork::plot_layout(widths = c(2.15, 1), guides = "collect")
pG <- PANEL_TAG(pG, "G")
pG <- dataset_line(pG, sprintf("Datasets: 6 studies \u00b7 %s B cells \u00b7 %d paired samples",
                               format(nrow(calls), big.mark = ","), pt_main$n),
                   pw_theme = theme(legend.position = "bottom", plot.margin = ggplot2::margin(0, 0, 0, 0)))
save_panel(pG, W_MM, H_MM)
print(tests, row.names = FALSE); print(per_dim, row.names = FALSE)
pub <- validate_anonymous_data(readr::read_csv(P("panelG_reference_coverage_statistics.csv"), show_col_types = FALSE))
pub_v <- setNames(pub$value, pub$metric)
stopifnot(sum(calls$k[calls$group == "Tumor B cell"] == 0) == pub_v[["Blind spot (caught by none)"]],
          sum(calls$k[calls$group == "Tumor B cell"] >= 1) == pub_v[["Caught by union of 5 methods"]],
          sum(calls$k[calls$group == "Tumor B cell"] == 5) == pub_v[["Caught by all 5 methods"]],
          n_tum == pub_v[["Total evaluable tumor cells"]], ge3_t > 0.6 * n_tum, pct_t > pct_n,
          pt_main$up == pt_main$n, pt_main$p < 0.001, pt_ge10$up == pt_ge10$n)
