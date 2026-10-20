.release_script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
.release_root <- normalizePath(file.path(.release_script_dir, "..", ".."))
source(file.path(.release_root, "privacy.R"), encoding = "UTF-8")
.script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE); f <- gsub("~+~", " ", sub("^--file=", "", a[grep("^--file=", a)]), fixed = TRUE)
  if (length(f)) dirname(normalizePath(f[1]))
  else getwd()
})
source(file.path(.script_dir, "00_config.R"), encoding = "UTF-8")
suppressMessages({ library(readr); library(dplyr); library(ggplot2) })
W_MM <- 180; H_MM <- 45
MAR <- c(SZ_TICK * 0.35, SZ_TICK * 1.75, SZ_TICK * 0.3, SZ_TICK * 0.3)
N_SAMPLES <- 21L
A <- validate_anonymous_data(readr::read_csv(IN_("panelA_label_distribution_by_subtype.csv"), show_col_types = FALSE))
stopifnot(all(c("subtype", "n_malignant", "n_ambiguous", "n_non_malignant", "n_total", "pct_malignant") %in% names(A)),
          setequal(A$subtype, DISEASE_LEVELS))
A <- A[match(DISEASE_LEVELS, A$subtype), ]
cnt <- as.matrix(A[, c("n_malignant", "n_ambiguous", "n_non_malignant")])
dimnames(cnt) <- list(subtype = A$subtype, label = LABEL_LEVELS)
N_TOTAL <- sum(cnt)
pct_tum <- cnt[, "Tumor"] / rowSums(cnt) * 100
stopifnot(all(rowSums(cnt) == A$n_total),
          all(A$n_total == c(2736, 26396, 20716, 4913)),
          N_TOTAL == 54761,
          all(abs(round(pct_tum, 1) - A$pct_malignant) < 1e-9),
          all(abs(round(pct_tum, 1) - c(88.2, 82.2, 92.0, 97.6)) < 1e-9),
          sum(cnt[, "Tumor"]) == 47968)
stopifnot(all(pct_tum >= 82))
chi <- stats::chisq.test(cnt, correct = FALSE)
chi_stat <- unname(chi$statistic); chi_df <- unname(chi$parameter); chi_p <- chi$p.value
chi_log10p <- stats::pchisq(chi_stat, chi_df, lower.tail = FALSE, log.p = TRUE) / log(10)
stopifnot(chi_df == 6, min(chi$expected) >= 5)
readr::write_csv(data.frame(
  test = "Pearson chi-square (subtype by classification)", n_cells_total = N_TOTAL,
  table_dimensions = "4 x 3 (subtype by Tumor/Ambiguous/Normal cell counts)",
  chi2_statistic = round(chi_stat, 3), degrees_of_freedom = chi_df,
  P = chi_p, P_text = fmtP(chi_p),
  log10_P = round(chi_log10p, 1), min_expected_count = round(min(chi$expected), 1),
  check.names = FALSE), OUT_("subtype_label_composition_chi_square.csv"))
BAR_H <- 0.78
ypos  <- setNames(rev(seq_along(DISEASE_LEVELS)), DISEASE_LEVELS)
rect <- do.call(rbind, lapply(DISEASE_LEVELS, function(s) {
  p <- cnt[s, ] / sum(cnt[s, ]) * 100
  data.frame(subtype = s, y = ypos[[s]], label = LABEL_LEVELS,
             xmin = c(0, cumsum(p)[-3]), xmax = pmin(cumsum(p), 100), row.names = NULL)
}))
rect$label <- factor(rect$label, levels = LABEL_LEVELS)
lab <- data.frame(subtype = DISEASE_LEVELS, y = unname(ypos[DISEASE_LEVELS]),
                  x = unname(pct_tum[DISEASE_LEVELS]) / 2,
                  txt = sprintf("%.1f%% (n = %s)", pct_tum[DISEASE_LEVELS],
                                format(rowSums(cnt)[DISEASE_LEVELS], big.mark = ",", trim = TRUE)))
pA <- ggplot(rect) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = y - BAR_H / 2, ymax = y + BAR_H / 2, fill = label),
            colour = "white", linewidth = 0.3) +
  geom_text(data = lab, aes(x = x, y = y, label = txt), colour = "white",
            size = SZ_TICK / PT, family = FIG_FONT) +
  scale_fill_manual(values = LABEL_COLS, breaks = LABEL_LEVELS, name = NULL) +
  scale_x_continuous("% of cells", limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
                     labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  scale_y_continuous(breaks = unname(ypos), labels = names(ypos), limits = c(0.5, 4.5), expand = c(0, 0)) +
  labs(y = NULL) +
  theme_figS1() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.key.size = grid::unit(8, "pt"), legend.text = element_text(size = SZ_TICK),
        legend.margin = ggplot2::margin(0, 0, 0, 0),
        legend.box.margin = ggplot2::margin(0, 0, 0, 0),
        legend.box.spacing = grid::unit(SZ_TICK * 0.55, "pt"),
        plot.margin = mar_pt(MAR))
pA <- tag_plot(pA, "A", W_MM, H_MM, MAR)
pA <- dataset_line(pA, sprintf("%s \u00b7 %d samples \u00b7 n = %s cells", DATASETS_TXT, N_SAMPLES,
                               format(N_TOTAL, big.mark = ",", trim = TRUE)))
save_panel(pA, "panelA_label_distribution", W_MM, H_MM)
print(cbind(as.data.frame(cnt), n_total = rowSums(cnt), pct_tumor = round(pct_tum, 1)))
