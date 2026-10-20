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
W_MM <- 112; H_MM <- 84
MAR <- c(SZ_TICK * 0.35, SZ_TICK * 0.5, SZ_TICK * 0.3, SZ_TICK * 0.3)
JITTER_SEED <- 20261010L
JITTER_HALF <- 0.30
MED_HALF    <- 0.36
X_LIM       <- c(0.45, 4.55)
PT_SIZE     <- 1.5
P <- validate_anonymous_data(readr::read_csv(IN_("patient_five_score_medians_robust_minmax.csv"), show_col_types = FALSE))
OLD <- validate_anonymous_data(readr::read_csv(IN_("panelD_radar_medians.csv"), show_col_types = FALSE))
stopifnot(all(c("patient_id", "subtype", "n_tumor_cells", "samples", AX, "group") %in% names(P)),
          all(c("subtype", "n_cells", AX) %in% names(OLD)),
          nrow(P) > 0, !anyDuplicated(P$patient_id), all(P$subtype %in% DISEASE_LEVELS),
          !anyNA(P[, AX]), all(as.matrix(P[, AX]) >= 0 & as.matrix(P[, AX]) <= 1))
stopifnot(identical(unname(DISEASE_GROUP[P$subtype]), P$group))
smp <- strsplit(as.character(P$samples), "+", fixed = TRUE)
stopifnot(all(lengths(smp) >= 1), !anyDuplicated(unlist(smp)), all(grepl("^[1-9][0-9]*$", unlist(smp))))
N_SAMPLES <- length(unlist(smp))
MULTI <- P$patient_id[lengths(smp) > 1]
n_sub <- table(factor(P$subtype, levels = DISEASE_LEVELS))
n_grp <- table(factor(P$group,   levels = GROUP_LEVELS))
stopifnot(all(n_sub >= 1), all(n_grp >= 3))
stopifnot(sum(P$n_tumor_cells) == sum(OLD$n_cells))
sub_cells <- tapply(P$n_tumor_cells, factor(P$subtype, levels = DISEASE_LEVELS), sum)
stopifnot(all(sub_cells == OLD$n_cells[match(DISEASE_LEVELS, OLD$subtype)]))
P$subtype <- factor(P$subtype, levels = DISEASE_LEVELS)
P$group   <- factor(P$group,   levels = GROUP_LEVELS)
P <- P[order(P$subtype), ]
wilcox_one <- function(v, g) {
  x <- v[g == "Aggressive"]; y <- v[g == "Indolent"]
  ties <- anyDuplicated(c(x, y)) > 0
  res <- stats::wilcox.test(x, y, alternative = "two.sided", exact = !ties, correct = TRUE)
  data.frame(med_agg = median(x), med_ind = median(y), n_agg = length(x), n_ind = length(y),
             W = unname(res$statistic), AUC = unname(res$statistic) / (length(x) * length(y)),
             p = res$p.value,
             method = ifelse(ties,
                           "Wilcoxon rank-sum, normal approximation with continuity correction (two-sided; tied values => exact not available)",
                           "Wilcoxon rank-sum, exact (two-sided)"),
             has_ties = ties, stringsAsFactors = FALSE)
}
ST <- do.call(rbind, lapply(AX_D, function(d) cbind(dimension = d, wilcox_one(P[[d]], P$group))))
ST$sig <- ST$p < 0.05
stopifnot(all(ST$p >= 0 & ST$p <= 1), all(ST$med_agg[ST$sig] > ST$med_ind[ST$sig]))
readr::write_csv(data.frame(
  dimension = ST$dimension, aggressive_median = round(ST$med_agg, 4), indolent_median = round(ST$med_ind, 4),
  n_aggressive_patients = ST$n_agg, n_indolent_patients = ST$n_ind, W_statistic = ST$W, AUC = round(ST$AUC, 4),
  P = ST$p, P_text = fmtP(ST$p), stars = star(ST$p), test_method = ST$method, has_ties = ST$has_ties,
  check.names = FALSE), OUT_("aggressive_vs_indolent_patient_tests.csv"))
cmp <- do.call(rbind, lapply(DISEASE_LEVELS, function(s) {
  pp <- P[P$subtype == s, ]
  data.frame(subtype = s, dimension = AX_D,
             cell_median = unlist(OLD[OLD$subtype == s, AX_D]),
             patient_median = vapply(AX_D, function(d) median(pp[[d]]), numeric(1)),
             n_patients = nrow(pp), n_tumor_cells = sum(pp$n_tumor_cells), row.names = NULL)
}))
cmp$patient_minus_cell_median <- cmp$patient_median - cmp$cell_median
readr::write_csv(transform(cmp, cell_median = round(cell_median, 4) + 0,
                           patient_median = round(patient_median, 4) + 0,
                           patient_minus_cell_median = round(patient_minus_cell_median, 4) + 0),
                 OUT_("subtype_cell_vs_patient_medians.csv"))
PAIRS <- utils::combn(DISEASE_LEVELS, 2, simplify = FALSE)
PW <- do.call(rbind, lapply(AX_D, function(d) do.call(rbind, lapply(PAIRS, function(q) {
  x <- P[[d]][P$subtype == q[1]]; y <- P[[d]][P$subtype == q[2]]
  ties <- anyDuplicated(c(x, y)) > 0
  r <- stats::wilcox.test(x, y, alternative = "two.sided", exact = !ties, correct = TRUE)
  data.frame(dimension = d, group1 = q[1], group2 = q[2], n1 = length(x), n2 = length(y), median1 = median(x), median2 = median(y),
             W = unname(r$statistic), P = r$p.value, method = ifelse(ties, "Wilcoxon, normal approximation with continuity correction (ties)", "Wilcoxon, exact"),
             exact_P_lower_bound = 2 / choose(length(x) + length(y), min(length(x), length(y))))
}))))
PW$P_BH_30 <- stats::p.adjust(PW$P, "BH"); PW$display <- PW$P < 0.05
readr::write_csv(transform(PW, median1 = round(median1, 4), median2 = round(median2, 4), P_text = fmtP(P), BH_text = fmtP(P_BH_30),
                           stars = star(P)), OUT_("subtype_pairwise_comparisons.csv"))
set.seed(JITTER_SEED, kind = "Mersenne-Twister", sample.kind = "Rejection")
off <- numeric(nrow(P))
for (g in DISEASE_LEVELS) { i <- which(P$subtype == g); off[i] <- sample(seq(-JITTER_HALF, JITTER_HALF, length.out = length(i))) }
P$xg <- as.integer(P$subtype); P$xj <- P$xg + off
stopifnot(all(abs(P$xj - P$xg) <= JITTER_HALF + 1e-9), min(P$xj) - 0.1 > X_LIM[1], max(P$xj) + 0.1 < X_LIM[2])
long <- P |>
  tidyr::pivot_longer(dplyr::all_of(AX_D), names_to = "dim", values_to = "v") |>
  dplyr::mutate(dim = factor(dim, levels = AX_D))
med <- long |> dplyr::group_by(dim, subtype, xg) |> dplyr::summarise(m = median(v), .groups = "drop")
BR0 <- 1.05; BRD <- 0.085; TICK <- 0.022
br <- do.call(rbind, lapply(AX_D, function(d) {
  b <- PW[PW$dimension == d & PW$display, ]
  if (!nrow(b)) return(NULL)
  b$x1 <- match(b$group1, DISEASE_LEVELS); b$x2 <- match(b$group2, DISEASE_LEVELS)
  b <- b[order(b$x2 - b$x1, b$x1), ]
  b$y <- BR0 + (seq_len(nrow(b)) - 1) * BRD
  data.frame(dim = factor(d, levels = AX_D), x1 = b$x1 + 0.06, x2 = b$x2 - 0.06, y = b$y, lab = star(b$P))
}))
Y_TOP <- max(br$y) + 0.06
AGG <- GROUP_LEVELS[1]; stopifnot(AGG == "Aggressive")
AGG_S <- names(DISEASE_GROUP)[DISEASE_GROUP == "Aggressive"]; IND_S <- names(DISEASE_GROUP)[DISEASE_GROUP == "Indolent"]
leg_lab <- setNames(sprintf("%s (%d patients)", DISEASE_LEVELS, as.integer(n_sub)), DISEASE_LEVELS)
n_pad <- 5L
if (requireNamespace("systemfonts", quietly = TRUE)) {
  w <- systemfonts::string_width(c("Aggressive", "Indolent", " "), family = FIG_FONT, size = 1000)
  if (all(w > 0)) n_pad <- as.integer(round((w[1] - w[2]) / w[3]))
}
leg_df <- data.frame(x = 1, y = 0, subtype = factor(DISEASE_LEVELS, levels = DISEASE_LEVELS), dim = factor(AX_D[1], levels = AX_D))
pD <- ggplot() +
  geom_segment(data = med, aes(x = xg - MED_HALF, xend = xg + MED_HALF, y = m, yend = m), linewidth = 0.7,
               colour = "black", lineend = "butt") +
  geom_point(data = long, aes(x = xj, y = v), fill = DISEASE_COLS[as.character(long$subtype)], shape = 21,
             size = PT_SIZE, stroke = 0.25, colour = "grey15") +
  geom_segment(data = br, aes(x = x1, xend = x2, y = y, yend = y), linewidth = 0.35) +
  geom_segment(data = br, aes(x = x1, xend = x1, y = y, yend = y - TICK), linewidth = 0.35) +
  geom_segment(data = br, aes(x = x2, xend = x2, y = y, yend = y - TICK), linewidth = 0.35) +
  geom_text(data = br, aes(x = (x1 + x2) / 2, y = y - 0.012, label = lab), vjust = 0, size = SZ_TICK / PT, family = FIG_FONT) +
  geom_point(data = leg_df[leg_df$subtype %in% AGG_S, ], aes(x, y, colour = subtype), size = 1, alpha = 0) +
  geom_point(data = leg_df[leg_df$subtype %in% IND_S, ], aes(x, y, fill = subtype), size = 1, alpha = 0, shape = 21) +
  scale_colour_manual(values = DISEASE_COLS, breaks = AGG_S, labels = leg_lab[AGG_S], name = "Aggressive") +
  scale_fill_manual(values = DISEASE_COLS, breaks = IND_S, labels = leg_lab[IND_S], name = paste0("Indolent", strrep(" ", max(n_pad, 0L)))) +
  guides(colour = guide_legend(order = 1, nrow = 1, title.position = "left", title.vjust = 0.5,
                               override.aes = list(alpha = 1, size = 2.4)),
         fill   = guide_legend(order = 2, nrow = 1, title.position = "left", title.vjust = 0.5,
                               override.aes = list(alpha = 1, size = 2.4, shape = 21, colour = "grey15"))) +
  scale_x_continuous(breaks = seq_along(DISEASE_LEVELS), labels = DISEASE_LEVELS, limits = X_LIM, expand = c(0, 0)) +
  scale_y_continuous("Patient median (0\u20131)", breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), expand = c(0, 0)) +
  coord_cartesian(ylim = c(-0.03, Y_TOP), clip = "off") +
  facet_wrap(~dim, nrow = 1, strip.position = "top") +
  labs(x = NULL) +
  theme_figS1() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = SZ_TITLE, colour = "black", margin = ggplot2::margin(t = 0, b = 2, unit = "pt")),
        panel.spacing.x = grid::unit(1.6, "mm"),
        axis.text.x = element_text(size = SZ_TICK, angle = 90, hjust = 1, vjust = 0.5, colour = "black"),
        legend.position = "bottom", legend.box = "vertical", legend.box.just = "left",
        legend.spacing.y = grid::unit(0, "pt"),
        legend.title = element_text(size = SZ_TICK, margin = ggplot2::margin(r = SZ_TICK * 0.5, unit = "pt")),
        legend.text = element_text(size = SZ_TICK, margin = ggplot2::margin(l = 2, r = SZ_TICK * 0.7, unit = "pt")),
        legend.key.size = grid::unit(8, "pt"), legend.key.height = grid::unit(SZ_TICK * 1.35, "pt"),
        legend.margin = ggplot2::margin(0, 0, 0, 0), legend.box.margin = ggplot2::margin(SZ_TICK * 0.2, 0, 0, 0),
        plot.margin = mar_pt(MAR))
pD <- tag_plot(pD, "D", W_MM, H_MM, MAR)
pD <- dataset_line(pD, sprintf("%s \u00b7 %d patients \u00b7 %s tumor cells", DATASETS_TXT,
                               nrow(P), format(sum(P$n_tumor_cells), big.mark = ",", trim = TRUE)))
save_panel(pD, "panelD_radar_by_subtype", W_MM, H_MM)
print(transform(PW[, c("dimension", "group1", "group2", "n1", "n2", "P", "P_BH_30")], P = fmtP(P), P_BH_30 = fmtP(P_BH_30)), row.names = FALSE)
print(data.frame(dimension = ST$dimension, aggressive = round(ST$med_agg, 3), indolent = round(ST$med_ind, 3), W = ST$W,
                 AUC = round(ST$AUC, 3), P = fmtP(ST$p), stars = star(ST$p),
                 method = ifelse(ST$has_ties, "normal approximation (ties)", "exact")), row.names = FALSE)
sep <- ST$dimension[ST$AUC == 1]
if (length(sep)) {
  invisible(NULL)
}
print(head(cmp[order(-abs(cmp$patient_minus_cell_median)), c("subtype", "dimension", "cell_median", "patient_median", "patient_minus_cell_median")], 5),
      row.names = FALSE, digits = 3)
