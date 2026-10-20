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
source_release_dependency(file.path(.release_root, "shared", "figure2_config.R"))
DATASETS_TXT <- "Analysis cohorts"
suppressPackageStartupMessages(library(grid))
W_MM <- 180; H_MM <- 46
info <- validate_anonymous_data(utils::read.csv(P("panelA_model_info.csv"), check.names = FALSE))
v <- setNames(info$value, info$item)
ord <- validate_anonymous_data(utils::read.csv(P("panelBC_row_order.csv"), check.names = FALSE))
n_cls <- v[paste0("n_metrics__", CATEGORY_LEVELS)]
stopifnot(v[["n_features"]] == 24, sum(n_cls) == 24, nrow(ord) == 24, setequal(ord$feature, FEATURE24),
          sum(ord$top4) == 4, all(ord$direction %in% names(DIR_COLS)),
          all(as.integer(table(factor(ord$category, levels = CATEGORY_LEVELS))) == as.integer(n_cls)))
mm  <- function(x) grid::unit(x, "mm")
tw  <- function(label, size = SZ_TICK, face = "plain")
  systemfonts::string_width(label, family = FIG_FONT, size = size, res = 720,
                            bold = identical(face, "bold")) / 10 * 25.4 / 72
txt <- function(label, x, y, just = c("left", "centre"), col = "black", size = SZ_TICK, face = "plain")
  grid::grid.text(label, x = mm(x), y = mm(y), just = just,
                  gp = grid::gpar(fontsize = size, col = col, fontfamily = FIG_FONT, fontface = face))
box <- function(x, y, w, h, col = "#CFCFCF", fill = "#F7F7F7", lwd = 0.8, r = 1.6)
  grid::grid.roundrect(mm(x), mm(y), mm(w), mm(h), r = mm(r), just = c("left", "bottom"),
                       gp = grid::gpar(col = col, fill = fill, lwd = lwd))
block <- function(x, y, w, h, fill = IG_COL)
  grid::grid.rect(mm(x), mm(y), mm(w), mm(h), just = c("left", "bottom"), gp = grid::gpar(col = NA, fill = fill))
dot <- function(x, y, r, fill) grid::grid.circle(mm(x), mm(y), mm(r), gp = grid::gpar(col = NA, fill = fill))
arrow_seg <- function(x0, x1, y, col = "#8A8A8A")
  grid::grid.segments(mm(x0), mm(y), mm(x1), mm(y), arrow = grid::arrow(length = mm(1.8), type = "closed"),
                      gp = grid::gpar(col = col, fill = col, lwd = 2))
fits <- function(label, x0, x1, size = SZ_TICK, what = label)
  if (tw(label, size) > (x1 - x0) + 1e-6) stop(sprintf("Label exceeds available width: %s (%.2f mm > %.2f mm)", what, tw(label, size), x1 - x0))
YC  <- 24
B1  <- c(x = 2.5,  y = 8, w = 54,   h = 32)
B2  <- c(x = 62.6, y = 8, w = 61,   h = 32)
B3  <- c(x = 130.5, y = 10, w = 22, h = 26)
S4  <- c(x0 = 161, x1 = 177)
GREY_DOT <- "#C4C4C4"
IGH <- c("M", "D", "G3", "G1", "A1", "G2", "G4", "E", "A2")
draw_A <- function() {
  grid::grid.newpage()
  box(B1[["x"]], B1[["y"]], B1[["w"]], B1[["h"]])
  x_row <- B1[["x"]] + 2.2
  txt("Ig transcripts", x_row, 36.6, c("left", "centre"), size = SZ_TITLE)
  txt("IGH", x_row, 30.0, c("left", "centre"))
  gx0 <- x_row + tw("IGH") + 1.6; gx1 <- B1[["x"]] + B1[["w"]] - 2
  slot <- (gx1 - gx0) / length(IGH); bw <- slot * 0.62
  grid::grid.segments(mm(gx0 - 0.9), mm(30.0), mm(gx1), mm(30.0), gp = grid::gpar(col = IG_COL, lwd = 0.9))
  for (k in seq_along(IGH)) {
    xk <- gx0 + (k - 0.5) * slot - bw / 2
    block(xk, 28.5, bw, 3.0)
    txt(IGH[k], xk + bw / 2, 25.6, c("centre", "centre"))
  }
  stopifnot(max(vapply(IGH, tw, numeric(1))) < slot - 1.0)
  lc <- list(c("IGK", 9.2), c("IGL", 9.2), c("JCHAIN", 15.4))
  xx <- x_row
  for (b in lc) {
    w <- as.numeric(b[2])
    grid::grid.rect(mm(xx), mm(16.5), mm(w), mm(3.4), just = c("left", "bottom"),
                    gp = grid::gpar(col = IG_COL, fill = "white", lwd = 0.9))
    fits(b[1], xx + 0.6, xx + w - 0.6)
    txt(b[1], xx + w / 2, 18.2, c("centre", "centre"))
    xx <- xx + w + 1.2
  }
  stopifnot(xx - 1.2 <= gx1)
  foot <- "scRNA-seq only \u00b7 no BCR needed"
  fits(foot, x_row, B1[["x"]] + B1[["w"]] - 1)
  txt(foot, x_row, 11.4, c("left", "centre"), col = GREY_DARK)
  arrow_seg(B1[["x"]] + B1[["w"]] + 0.8, B2[["x"]] - 0.8, YC)
  box(B2[["x"]], B2[["y"]], B2[["w"]], B2[["h"]])
  x_l <- B2[["x"]] + 2.2; x_r <- B2[["x"]] + B2[["w"]] - 2
  ttl <- sprintf("%d Ig metrics", v[["n_features"]])
  txt(ttl, x_l, 36.6, c("left", "centre"), size = SZ_TITLE)
  key <- "Top-4 drivers"; kx1 <- x_r; kx0 <- kx1 - tw(key)
  txt(key, kx1, 36.6, c("right", "centre"))
  dot(kx0 - 1.6, 36.6, 0.85, DIR_COLS[["lower in tumor"]]); dot(kx0 - 3.9, 36.6, 0.85, DIR_COLS[["higher in tumor"]])
  stopifnot(kx0 - 3.9 - 0.85 > x_l + tw(ttl, SZ_TITLE) + 2)
  colx <- c(x_l, B2[["x"]] + B2[["w"]] / 2 + 1.5); yn <- c(31.2, 22.8, 14.4)
  for (i in seq_along(CATEGORY_LEVELS)) {
    cat_ <- CATEGORY_LEVELS[i]; cx <- colx[(i - 1) %% 2 + 1]; ny <- yn[(i - 1) %/% 2 + 1]
    nm <- CAT_SHORT[[cat_]]
    fits(nm, cx, if ((i - 1) %% 2 == 0) colx[2] - 1.5 else x_r, what = nm)
    txt(nm, cx, ny, c("left", "centre"))
    f <- ord[ord$category == cat_, ]
    f <- f[order(!f$top4, -f$auc_strength), ]
    for (j in seq_len(nrow(f)))
      dot(cx + 0.95 + (j - 1) * 2.55, ny - 3.5, if (f$top4[j]) 1.0 else 0.85,
          if (f$top4[j]) DIR_COLS[[f$direction[j]]] else GREY_DOT)
  }
  arrow_seg(B2[["x"]] + B2[["w"]] + 0.8, B3[["x"]] - 0.8, YC)
  box(B3[["x"]], B3[["y"]], B3[["w"]], B3[["h"]], col = IG_COL, fill = "white", lwd = 1.1)
  cx3 <- B3[["x"]] + B3[["w"]] / 2
  txt("GBM", cx3, 32.4, c("centre", "centre"), size = SZ_TITLE)
  tree <- function(x0, y0, leaves) {
    nodes <- list(c(x0, y0), c(x0 - 1.35, y0 - 2.5), c(x0 + 1.35, y0 - 2.5))
    lv <- list(c(x0 - 2.0, y0 - 5.0), c(x0 - 0.65, y0 - 5.0), c(x0 + 0.65, y0 - 5.0), c(x0 + 2.0, y0 - 5.0))
    seg <- function(a, b) grid::grid.segments(mm(a[1]), mm(a[2]), mm(b[1]), mm(b[2]), gp = grid::gpar(col = GREY_DARK, lwd = 0.7))
    seg(nodes[[1]], nodes[[2]]); seg(nodes[[1]], nodes[[3]])
    seg(nodes[[2]], lv[[1]]); seg(nodes[[2]], lv[[2]]); seg(nodes[[3]], lv[[3]]); seg(nodes[[3]], lv[[4]])
    for (n in nodes) dot(n[1], n[2], 0.5, IG_COL)
    for (k in seq_along(lv)) dot(lv[[k]][1], lv[[k]][2], 0.52, DIR_COLS[[leaves[k]]])
  }
  TREE_DX <- 6.4
  stopifnot(TREE_DX - 4.0 - 2 * 0.52 >= 1.3, 0.65 * 2 - 2 * 0.52 > 0.2, TREE_DX + 2.0 + 0.52 < B3[["w"]] / 2 - 0.8)
  U <- "higher in tumor"; D <- "lower in tumor"
  tree(cx3 - TREE_DX, 27.6, c(D, U, U, D)); tree(cx3, 27.6, c(U, D, D, U)); tree(cx3 + TREE_DX, 27.6, c(D, D, U, U))
  txt(sprintf("%d trees", v[["n_trees"]]), cx3, 17.2, c("centre", "centre"), col = GREY_DARK)
  splits <- sprintf("%d splits/tree", v[["interaction_depth"]])
  txt(splits, cx3, 13.4, c("centre", "centre"), col = GREY_DARK)
  fits(sprintf("%d trees", v[["n_trees"]]), B3[["x"]] + 1, B3[["x"]] + B3[["w"]] - 1)
  fits(splits, B3[["x"]] + 0.8, B3[["x"]] + B3[["w"]] - 0.8)
  lab <- "Training labels (BCR / reference)"
  txt(lab, cx3, 43.4, c("centre", "centre"), col = GREY_DARK)
  stopifnot(cx3 + tw(lab) / 2 < W_MM - 0.5)
  dot(cx3 - 1.3, 40.4, 0.75, DIR_COLS[["higher in tumor"]]); dot(cx3 + 1.3, 40.4, 0.75, DIR_COLS[["lower in tumor"]])
  grid::grid.segments(mm(cx3), mm(39.3), mm(cx3), mm(37.6),
                      gp = grid::gpar(col = GREY_MID, lwd = 1.1, lty = "11", lineend = "butt"))
  grid::grid.segments(mm(cx3), mm(37.7), mm(cx3), mm(B3[["y"]] + B3[["h"]] + 0.45),
                      arrow = grid::arrow(length = mm(1.3), type = "closed"),
                      gp = grid::gpar(col = GREY_MID, fill = GREY_MID, lwd = 1.1, linejoin = "mitre"))
  arrow_seg(B3[["x"]] + B3[["w"]] + 0.8, S4[["x0"]] - 1.8, YC)
  cx4 <- mean(S4)
  txt("Ig-score", cx4, 32.4, c("centre", "centre"), size = SZ_TITLE)
  ramp <- grDevices::colorRampPalette(SCORE_RAMP)(80); wb <- (S4[["x1"]] - S4[["x0"]]) / length(ramp)
  for (k in seq_along(ramp))
    grid::grid.rect(mm(S4[["x0"]] + (k - 1) * wb), mm(21.6), mm(wb + 0.02), mm(4.4), just = c("left", "bottom"),
                    gp = grid::gpar(col = NA, fill = ramp[k]))
  grid::grid.rect(mm(S4[["x0"]]), mm(21.6), mm(S4[["x1"]] - S4[["x0"]]), mm(4.4), just = c("left", "bottom"),
                  gp = grid::gpar(col = GREY_DARK, fill = NA, lwd = 0.5))
  txt("0", S4[["x0"]], 18.6, c("centre", "centre"), col = DIR_COLS[["lower in tumor"]])
  txt("1", S4[["x1"]], 18.6, c("centre", "centre"), col = DIR_COLS[["higher in tumor"]])
  txt("P(tumor)", cx4, 14.6, c("centre", "centre"), col = GREY_DARK)
  stopifnot(S4[["x1"]] + tw("1") / 2 < W_MM - 0.5, cx4 + tw("Ig-score", SZ_TITLE) / 2 < W_MM - 0.5)
  grid::grid.text("A", x = mm(0.6), y = mm(H_MM - 0.6), just = c("left", "top"),
                  gp = grid::gpar(fontface = "bold", fontsize = SZ_TAG, col = "black", fontfamily = FIG_FONT))
}
assign("..dataset_txt", sprintf("%s \u00b7 %s B cells | GBM training set: %s B cells", DATASETS_TXT,
                                format(v[["n_cells_atlas"]], big.mark = ","), format(v[["n_train"]], big.mark = ",")),
       envir = .GlobalEnv)
pdf_f <- O("plot.pdf")
open_pdf(pdf_f, W_MM, H_MM); draw_A(); draw_dataset_line(); close_pdf()
if (capabilities("cairo")) {
  grDevices::png(O("plot.png"), width = W_MM / 25.4, height = H_MM / 25.4, units = "in", res = 600,
                 family = FIG_FONT, type = "cairo")
} else grDevices::png(O("plot.png"), width = W_MM / 25.4, height = H_MM / 25.4, units = "in", res = 600)
draw_A(); draw_dataset_line(); invisible(grDevices::dev.off())
rm("..dataset_txt", envir = .GlobalEnv)
