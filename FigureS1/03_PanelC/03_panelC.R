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
W_MM <- 68; H_MM <- 84
MAR <- c(SZ_TICK * 0.35, SZ_TICK * 1.0, SZ_TICK * 0.3, SZ_TICK * 0.3)
N_CELLS <- 48745L
TILE_EDGE <- "#BDBDBD"
sp <- validate_anonymous_data(readr::read_csv(IN_("panelC_spearman_matrix.csv"), show_col_types = FALSE))
stopifnot(all(c("component_1", "component_2", "spearman_rho", "p_value") %in% names(sp)),
          nrow(sp) == 25L, setequal(unique(sp$component_1), AX), setequal(unique(sp$component_2), AX),
          !anyDuplicated(sp[, c("component_1", "component_2")]), !anyNA(sp$spearman_rho), !anyNA(sp$p_value),
          all(sp$p_value >= 0 & sp$p_value <= 1))
to_mat <- function(col) {
  w <- sp |> dplyr::select(component_1, component_2, dplyr::all_of(col)) |>
    tidyr::pivot_wider(names_from = component_2, values_from = dplyr::all_of(col))
  m <- as.matrix(w[, AX]); rownames(m) <- w$component_1; m[AX, AX, drop = FALSE]
}
M <- to_mat("spearman_rho"); PM <- to_mat("p_value")
stopifnot(isTRUE(all.equal(M, t(M), tolerance = 1e-12)),
          isTRUE(all.equal(PM, t(PM), tolerance = 1e-12)),
          all(abs(diag(M) - 1) < 1e-12))
off <- sp[sp$component_1 != sp$component_2, ]
stopifnot(max(abs(off$spearman_rho)) < 0.65,
          all(star(off$p_value) == "****"))
AXC <- AX[-length(AX)]
AXR <- AX[-1]
cells <- expand.grid(j = seq_along(AXC), i = seq_along(AXR))
cells$col <- AXC[cells$j]; cells$row <- AXR[cells$i]
cells$ic <- match(cells$col, AX); cells$ir <- match(cells$row, AX)
cells <- cells[cells$ir > cells$ic, ]
cells <- cells[order(cells$ic, cells$ir), ]
rownames(cells) <- NULL
stopifnot(nrow(cells) == choose(length(AX), 2), !anyDuplicated(cells[, c("col", "row")]),
          all(cells$col != cells$row))
cells$rho <- M[cbind(cells$ir, cells$ic)]; cells$p <- PM[cbind(cells$ir, cells$ic)]
cells$y   <- length(AXR) + 1 - cells$i
.pal <- scales::seq_gradient_pal(COR_LOW, COR_HIGH, "Lab")
cells$fill_hex <- .pal(scales::rescale(scales::squish(cells$rho, COR_LIMITS), to = c(0, 1), from = COR_LIMITS))
cells$txt_col  <- text_on(cells$fill_hex)
cells$rho_txt  <- sub("^-", "\u2212", sprintf("%.2f", cells$rho))
cells$star     <- star(cells$p)
cells$p_txt    <- P_txt(cells$p)
readr::write_csv(data.frame(score1 = cells$col, score2 = cells$row, spearman_rho = cells$rho,
                            p_value_input = cells$p, P_text = ifelse(cells$p < 1e-300, "< 1.00e-300", fmtP(cells$p)),
                            stars = cells$star, n_cells = N_CELLS),
                 OUT_("five_score_Spearman_results.csv"))
pC <- ggplot(cells, aes(j, y)) +
  geom_tile(aes(fill = rho), width = 1, height = 1, colour = TILE_EDGE, linewidth = 0.3) +
  geom_text(aes(label = rho_txt), colour = cells$txt_col,
            size = SZ_TICK / PT, family = FIG_FONT) +
  scale_fill_gradient(low = COR_LOW, high = COR_HIGH, limits = COR_LIMITS, oob = scales::squish,
                      breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), name = "Spearman \u03c1",
                      guide = vector_colourbar(height_mm = 25, width_mm = 3.5, nbin = 100)) +
  scale_x_continuous(breaks = seq_along(AXC), labels = AXC, position = "bottom", expand = c(0, 0)) +
  scale_y_continuous(breaks = rev(seq_along(AXR)), labels = AXR, expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(xlim = c(0.5, length(AXC) + 0.5), ylim = c(0.5, length(AXR) + 0.5), expand = FALSE) +
  theme_figS1() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "right", legend.justification = "center",
        legend.title = element_text(size = SZ_TICK), legend.text = element_text(size = SZ_TICK),
        legend.margin = ggplot2::margin(0, 0, 0, 0),
        plot.margin = mar_pt(MAR))
pC <- tag_plot(pC, "C", W_MM, H_MM, MAR)
pC <- dataset_line(pC, sprintf("%s \u00b7 n = %s labelled cells", DATASETS_TXT, format(N_CELLS, big.mark = ",", trim = TRUE)))
save_panel(pC, "panelC_spearman", W_MM, H_MM)
print(round(M, 4))
k <- which.max(replace(sp$spearman_rho, sp$component_1 == sp$component_2, -Inf))
k <- which.min(abs(replace(sp$spearman_rho, sp$component_1 == sp$component_2, Inf)))
print(data.frame(column = cells$col, row = cells$row, rho = cells$rho_txt, cell_color = cells$fill_hex, text_color = cells$txt_col), row.names = FALSE)
