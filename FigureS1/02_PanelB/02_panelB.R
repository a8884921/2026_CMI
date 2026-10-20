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
suppressMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
W_MM <- 180; H_MM <- 62
MAR <- c(SZ_TICK * 0.35, SZ_TICK * 0.5, SZ_TICK * 0.3, SZ_TICK * 0.5)
B <- validate_anonymous_data(readr::read_csv(IN_("panelB_sample_label_distribution.csv"), show_col_types = FALSE))
stopifnot(all(c("sample_id", "patient_id", "include_in_analysis") %in% names(B)))
B$include_in_analysis <- validate_inclusion(B$include_in_analysis)
EXCLUDED <- as.character(B$sample_id[!B$include_in_analysis])
REPEATED <- duplicated(B$patient_id) | duplicated(B$patient_id, fromLast = TRUE)
SAME_PT <- as.character(B$sample_id[REPEATED])
stopifnot(all(c("plot_order", "sample_id", "subtype", "n_malignant", "n_ambiguous", "n_non_malignant",
                "n_cells", "pct_malignant", "pct_ambiguous", "pct_non_malignant") %in% names(B)),
          nrow(B) == 21L, !anyDuplicated(B$sample_id), all(B$subtype %in% DISEASE_LEVELS),
          all(B$n_malignant + B$n_ambiguous + B$n_non_malignant == B$n_cells),
          sum(B$n_cells) == 54761,
          all(EXCLUDED %in% B$sample_id), all(SAME_PT %in% B$sample_id), !any(SAME_PT %in% EXCLUDED))
B <- B |>
  dplyr::mutate(pct_tum = n_malignant / n_cells * 100, pct_amb = n_ambiguous / n_cells * 100,
                pct_nor = n_non_malignant / n_cells * 100)
stopifnot(all(abs(B$pct_malignant - B$pct_tum) < 0.05), all(abs(B$pct_ambiguous - B$pct_amb) < 0.05),
          all(abs(B$pct_non_malignant - B$pct_nor) < 0.05))
stopifnot(setequal(B$sample_id[B$n_ambiguous == B$n_cells], EXCLUDED))
B <- B |>
  dplyr::mutate(subtype = factor(subtype, levels = DISEASE_LEVELS)) |>
  dplyr::arrange(subtype, dplyr::desc(pct_tum), plot_order) |>
  dplyr::mutate(x = dplyr::row_number(), sample_id = factor(sample_id, levels = sample_id))
sub_n <- table(B$subtype)
long <- B |>
  dplyr::select(sample_id, x, subtype, Tumor = pct_tum, Ambiguous = pct_amb, Normal = pct_nor) |>
  tidyr::pivot_longer(dplyr::all_of(LABEL_LEVELS), names_to = "grp", values_to = "pct") |>
  dplyr::mutate(grp = factor(grp, levels = LABEL_LEVELS))
Y_DAG <- 111; Y_B0 <- 112.5; Y_B1 <- 130; Y_TOP <- 131
band <- B |>
  dplyr::group_by(subtype) |>
  dplyr::summarise(lo = min(x) - 0.5, hi = max(x) + 0.5, n = dplyr::n(), .groups = "drop") |>
  dplyr::mutate(fill = unname(DISEASE_COLS[as.character(subtype)]), txt_col = text_on(fill),
                lab = sprintf("%s (n = %d)", subtype, n))
dag <- B |> dplyr::filter(as.character(sample_id) %in% EXCLUDED)
dbl <- B |> dplyr::filter(as.character(sample_id) %in% SAME_PT)
stopifnot(nrow(dag) == length(EXCLUDED), nrow(dbl) == length(SAME_PT))
pB <- ggplot() +
  geom_col(data = long, aes(x = x, y = pct, fill = grp), width = 0.86, colour = "white", linewidth = 0.2,
           position = position_stack(reverse = TRUE)) +
  geom_rect(data = band, aes(xmin = lo + 0.04, xmax = hi - 0.04, ymin = Y_B0, ymax = Y_B1), fill = band$fill) +
  geom_text(data = band, aes(x = (lo + hi) / 2, y = (Y_B0 + Y_B1) / 2, label = lab), colour = band$txt_col,
            size = SZ_TICK / PT, family = FIG_FONT) +
  geom_text(data = dag, aes(x = x, y = 101.5, label = "\u2020"), vjust = 0,
            size = SZ_TICK / PT, family = FIG_FONT, colour = "black") +
  geom_text(data = dbl, aes(x = x, y = 101.5, label = "\u2021"), vjust = 0,
            size = SZ_TICK / PT, family = FIG_FONT, colour = "black") +
  annotate("segment", x = -Inf, xend = -Inf, y = 0, yend = 100, linewidth = LW_AXIS, colour = "black") +
  scale_fill_manual(values = LABEL_COLS, breaks = LABEL_LEVELS, name = NULL) +
  scale_x_continuous(breaks = B$x, labels = as.character(B$sample_id), expand = expansion(add = c(0.55, 0.55))) +
  scale_y_continuous("% of cells", limits = c(0, Y_TOP), breaks = c(0, 50, 100), expand = c(0, 0)) +
  labs(x = NULL) +
  theme_figS1() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = SZ_TICK,
                                   margin = ggplot2::margin(t = SZ_TICK * 0.2)),
        axis.line.x = element_blank(), axis.line.y = element_blank(),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.key.size = grid::unit(8, "pt"), legend.text = element_text(size = SZ_TICK),
        legend.margin = ggplot2::margin(0, 0, 0, 0),
        legend.box.margin = ggplot2::margin(SZ_TICK * 0.1, 0, 0, 0),
        plot.margin = mar_pt(MAR))
pB <- tag_plot(pB, "B", W_MM, H_MM, MAR)
pB <- dataset_line(pB, sprintf("%s \u00b7 %d samples \u00b7 n = %s cells", DATASETS_TXT, nrow(B),
                               format(sum(B$n_cells), big.mark = ",", trim = TRUE)))
save_panel(pB, "panelB_sample_distribution", W_MM, H_MM)
