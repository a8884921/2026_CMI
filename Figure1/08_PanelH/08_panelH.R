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
W_MM <- 90; H_MM <- 68
AXES     <- c("Ig", "CNV", "MTP", "Prol", "Diff")
AX_ANGLE <- c(Ig = 90, CNV = 162, MTP = 234, Prol = 306, Diff = 18)
R_MAX    <- 0.80
R_TICKS  <- c(0.2, 0.4, 0.6, 0.8)
pt <- validate_anonymous_data(readr::read_csv(P("panelH_patient_five_score_medians.csv"), show_col_types = FALSE))
stopifnot(all(c("patient_id", "subtype", "group", "n_tumor_cells", AXES) %in% names(pt)),
          !anyDuplicated(pt$patient_id), sum(pt$n_tumor_cells) == 47968,
          all(pt$group == ifelse(pt$subtype %in% c("BL", "DLBCL"), "Aggressive", "Indolent")))
pt$subtype <- factor(pt$subtype, levels = DISEASE_LEVELS)
n_by <- table(pt$subtype); n_grp <- table(pt$group)
tst <- do.call(rbind, lapply(AXES, function(ax) {
  a <- pt[[ax]][pt$group == "Aggressive"]; b <- pt[[ax]][pt$group == "Indolent"]
  ties <- anyDuplicated(c(a, b)) > 0
  w <- stats::wilcox.test(a, b, exact = !ties, correct = TRUE)
  data.frame(dimension = ax, aggressive_median = round(median(a), 3), indolent_median = round(median(b), 3),
             n_aggressive_patients = length(a), n_indolent_patients = length(b), W = unname(w$statistic), P = w$p.value,
             method = if (ties) "Wilcoxon rank-sum, normal approximation with continuity correction" else "Wilcoxon rank-sum, exact")
}))
tst$P_text <- fmtP(tst$P); tst$stars <- star(tst$P)
tst$P_BH <- fmtP(stats::p.adjust(tst$P, method = "BH")); tst$stars_BH <- star(stats::p.adjust(tst$P, method = "BH"))
readr::write_csv(tst, O("aggressive_vs_indolent_tests.csv"))
stopifnot(identical(tst$stars, tst$stars_BH))
sub_test <- function(sub) do.call(rbind, lapply(AXES, function(ax) {
  a <- pt[[ax]][pt$subtype %in% sub]; b <- pt[[ax]][pt$group == "Indolent"]
  ties <- anyDuplicated(c(a, b)) > 0
  data.frame(comparison = sprintf("%s (%d patients) vs indolent (%d patients)", paste(sub, collapse = "+"), length(a), length(b)),
             dimension = ax, P = stats::wilcox.test(a, b, exact = !ties, correct = TRUE)$p.value)
}))
rob <- rbind(sub_test("DLBCL"), sub_test("BL"))
rob$P_text <- fmtP(rob$P); rob$stars <- star(rob$P)
readr::write_csv(rob, O("single_aggressive_subtype_vs_indolent_robustness.csv"))
print(rob[, c("comparison", "dimension", "P_text", "stars")], row.names = FALSE)
stopifnot(all(rob$P[rob$dimension %in% c("CNV", "MTP", "Prol", "Diff") & grepl("^DLBCL", rob$comparison)] < 0.05))
med <- pt |> dplyr::group_by(subtype) |>
  dplyr::summarise(n_patients = dplyr::n(), n_tumor_cells = sum(n_tumor_cells),
                   dplyr::across(dplyr::all_of(AXES), median), .groups = "drop")
readr::write_csv(med, O("subtype_patient_five_score_medians.csv"))
v1 <- validate_anonymous_data(readr::read_csv(P("panelH_radar_medians.csv"), show_col_types = FALSE))
stopifnot(sum(v1$n_cells) == sum(med$n_tumor_cells))
cmp <- dplyr::bind_rows(
  dplyr::mutate(as.data.frame(v1[, c("subtype", AXES)]), metric = "cell_level_median"),
  dplyr::mutate(as.data.frame(med[, c("subtype", AXES)]), subtype = as.character(subtype), metric = "patient_level_median"))
readr::write_csv(cmp, O("subtype_cell_vs_patient_profiles.csv"))
print(as.data.frame(med), digits = 3); print(tst[, c("dimension", "aggressive_median", "indolent_median", "P_text", "stars", "method")], row.names = FALSE)
stopifnot(max(as.matrix(med[, AXES])) <= R_MAX)
to_xy <- function(axis, r) { a <- AX_ANGLE[axis] * pi / 180; data.frame(x = r * cos(a), y = r * sin(a)) }
poly_df <- do.call(rbind, lapply(seq_len(nrow(med)), function(i) {
  d <- do.call(rbind, lapply(AXES, function(ax) { p <- to_xy(ax, med[[ax]][i]); p$axis <- ax; p }))
  d$subtype <- med$subtype[i]; rbind(d, d[1, ])
}))
poly_df$subtype <- factor(poly_df$subtype, levels = rev(DISEASE_LEVELS))
web <- do.call(rbind, lapply(R_TICKS, function(r) {
  d <- do.call(rbind, lapply(AXES, function(ax) { p <- to_xy(ax, r); p$axis <- ax; p }))
  d$r <- r; rbind(d, d[1, ]) }))
spokes  <- do.call(rbind, lapply(AXES, function(ax) { p <- to_xy(ax, R_MAX); p$axis <- ax; p }))
lab_pos <- do.call(rbind, lapply(AXES, function(ax) { p <- to_xy(ax, R_MAX * 1.10); p$axis <- ax; p }))
lab_pos$label <- paste0(lab_pos$axis, tst$stars[match(lab_pos$axis, tst$dimension)])
lab_pos$hjust <- c(Ig = 0.5, CNV = 1, MTP = 0.85, Prol = 0.15, Diff = 0)[lab_pos$axis]
lab_pos$vjust <- c(Ig = 0, CNV = 0.5, MTP = 1, Prol = 1, Diff = 0.5)[lab_pos$axis]
AGG <- c("BL", "DLBCL"); IND <- c("FL", "MZL")
stopifnot(all(DISEASE_GROUP[AGG] == "Aggressive"), all(DISEASE_GROUP[IND] == "Indolent"))
leg_lab <- with(med, setNames(sprintf("%s (%d patients)", subtype, n_patients), subtype))
n_pad <- 5L
if (requireNamespace("systemfonts", quietly = TRUE)) {
  w <- systemfonts::string_width(c("Aggressive", "Indolent", " "), family = FIG_FONT, size = 1000)
  if (all(w > 0)) n_pad <- as.integer(round((w[1] - w[2]) / w[3]))
}
TITLE_AGG <- "Aggressive"; TITLE_IND <- paste0("Indolent", strrep(" ", max(n_pad, 0L)))
leg_df <- data.frame(x = 0, y = 0, subtype = factor(DISEASE_LEVELS, levels = DISEASE_LEVELS))
pH <- ggplot() +
  geom_polygon(data = web, aes(x, y, group = r), fill = NA, colour = "#DDDDDD", linewidth = LW_AXIS) +
  geom_segment(data = spokes, aes(x = 0, y = 0, xend = x, yend = y), colour = "#DDDDDD", linewidth = LW_AXIS) +
  geom_polygon(data = poly_df, aes(x, y, fill = subtype, group = subtype),
               colour = NA, alpha = 0.12, show.legend = FALSE) +
  geom_polygon(data = poly_df, aes(x, y, colour = subtype, group = subtype),
               fill = NA, linewidth = LW_DATA, show.legend = FALSE) +
  geom_point(data = poly_df, aes(x, y, colour = subtype), size = 0.9, show.legend = FALSE) +
  geom_point(data = leg_df[leg_df$subtype %in% AGG, ], aes(x, y, colour = subtype), size = 0.9, alpha = 0) +
  geom_point(data = leg_df[leg_df$subtype %in% IND, ], aes(x, y, fill = subtype), size = 0.9, alpha = 0, shape = 21) +
  geom_text(data = lab_pos, aes(x, y, label = label, hjust = hjust, vjust = vjust),
            size = SZ_TICK / PT, family = FIG_FONT) +
  annotate("text", x = 0, y = -R_TICKS[c(2, 4)] * cos(36 * pi / 180), label = sprintf("%.1f", R_TICKS[c(2, 4)]),
           size = SZ_TICK / PT, family = FIG_FONT, colour = SEM[["neutral"]], hjust = 0.5, vjust = 1.2) +
  scale_colour_manual(values = DISEASE_COLS, breaks = AGG, labels = leg_lab[AGG], name = TITLE_AGG) +
  scale_fill_manual(values = DISEASE_COLS, breaks = IND, labels = leg_lab[IND], name = TITLE_IND) +
  guides(colour = guide_legend(order = 1, nrow = 1, title.position = "left", title.vjust = 0.5,
                               override.aes = list(alpha = 1, size = 2.4)),
         fill   = guide_legend(order = 2, nrow = 1, title.position = "left", title.vjust = 0.5,
                               override.aes = list(alpha = 1, size = 2.4, shape = 21, colour = unname(DISEASE_COLS[IND])))) +
  coord_equal(xlim = c(-1.45, 1.45), ylim = c(-0.80, 0.98), clip = "off") +
  theme_void(base_size = SZ_TICK, base_family = FIG_FONT) +
  theme(plot.margin = ggplot2::margin(SZ_TICK * 1.2, SZ_TICK * 0.4, 0, SZ_TICK * 0.4),
        legend.position = "bottom", legend.box = "vertical", legend.box.just = "left",
        legend.spacing.y = grid::unit(0, "pt"),
        legend.title = element_text(size = SZ_TICK, margin = ggplot2::margin(r = SZ_TICK * 0.5, unit = "pt")),
        legend.text = element_text(size = SZ_TICK, margin = ggplot2::margin(l = 2, r = SZ_TICK * 0.7, unit = "pt")),
        legend.key.size = grid::unit(8, "pt"), legend.key.height = grid::unit(SZ_TICK * 1.35, "pt"), legend.margin = ggplot2::margin(0, 0, 0, 0),
        legend.box.margin = ggplot2::margin(-SZ_TICK * 0.5, 0, 0, 0))
pH <- PANEL_TAG(pH, "H")
pH <- dataset_line(pH, sprintf("Datasets: 6 studies \u00b7 %d patients \u00b7 %s tumor cells",
                               nrow(pt), format(sum(pt$n_tumor_cells), big.mark = ",")))
save_panel(pH, W_MM, H_MM)
sig4 <- tst$P[tst$dimension %in% c("CNV", "MTP", "Prol", "Diff")]
stopifnot(all(sig4 < 0.01),
          all(tst$aggressive_median[tst$dimension %in% c("CNV", "MTP", "Prol", "Diff")] >
              tst$indolent_median[tst$dimension %in% c("CNV", "MTP", "Prol", "Diff")]),
          all(med$CNV[med$subtype %in% c("BL", "DLBCL")] > max(med$CNV[med$subtype %in% c("FL", "MZL")])),
          all(med$Diff[med$subtype %in% c("BL", "DLBCL")] > max(med$Diff[med$subtype %in% c("FL", "MZL")])))
