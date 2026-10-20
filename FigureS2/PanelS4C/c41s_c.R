.release_script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
.release_root <- normalizePath(file.path(.release_script_dir, "..", ".."))
source(file.path(.release_root, "privacy.R"), encoding = "UTF-8")
if (!exists("FIG4", inherits = TRUE)) {
  .d <- local({ a <- commandArgs(FALSE); f <- gsub("~+~", " ", sub("^--file=", "", a[grep("^--file=", a)]), fixed = TRUE)
                if (length(f)) dirname(normalizePath(f)) else getwd() })
  source_release_dependency(file.path(.release_root, "shared", "figure4_config.R"))
  source_release_dependency(file.path(.release_root, "shared", "figure4_plot_utils.R"))
}
if (!exists("build_only", inherits = TRUE)) build_only <- FALSE
PANEL <- "PanelS4C"
R <- validate_anonymous_data(read_panel_csv(PANEL, "FigS4C_robustness.csv"))
ORD <- c("all", "g1", "cycling", "cnvhi", "adj", "adj_nf", "naive", "dz", "mem", "ge10", "ge50", "spearman")
if (!all(setdiff(ORD, "adj_nf") %in% R$id)) warning("S4C missing analysis: ", paste(setdiff(setdiff(ORD, "adj_nf"), R$id), collapse = ", "))
stopifnot(all(R$id %in% ORD))
ORD <- ORD[ORD %in% R$id]
if (any(R$p >= 0.05)) warning("S4C contains analyses with P >= 0.05.")
GAP_AFTER <- c("all", tail(intersect(c("cycling", "cnvhi", "adj", "adj_nf"), ORD), 1), "mem", "ge50")
step <- 1 + 0.32 * (ORD %in% GAP_AFTER)
ypos <- setNames(c(0, -cumsum(step[-length(step)])), ORD) - min(c(0, -cumsum(step[-length(step)])))
PCOL <- intersect(c("MHC_I", "MHC_II"), unique(R$prog))
if (length(PCOL) != 2L) stop("S4C requires prepared MHC-I and MHC-II plot data.")
R[, `:=`(y = ypos[id], x = match(prog, PCOL))]
blend_color <- function(base, t) { m <- grDevices::colorRamp(c("#FFFFFF", base))(pmin(pmax(t, 0), 1)); grDevices::rgb(m[, 1], m[, 2], m[, 3], maxColorValue = 255) }
R[, t := abs(frac - 0.5) / 0.5]
R[, fill := fifelse(frac >= 0.5, blend_color(DIV_LOW, t), blend_color(DIV_HIGH, t))]
R[, txtcol := ifelse(t > 0.55, "white", "black")]
R[, txt := sprintf("%d/%d %s", n_consistent, n, stars)]
rowlab <- unique(R[, .(id, row, y)])
COLLAB <- paste0(R$label[match(PCOL, R$prog)], "\nM3 < M1")
p <- ggplot(R, aes(x, y)) +
  geom_tile(aes(fill = fill), width = 0.96, height = 0.90, colour = NA) +
  scale_fill_identity() +
  geom_text(aes(label = txt), colour = R$txtcol, size = SZ_SMALL / PT, family = "Arial") +
  scale_x_continuous(breaks = 1:2, labels = COLLAB, position = "top", limits = c(0.5, 2.5), expand = c(0, 0)) +
  scale_y_continuous(breaks = rowlab$y, labels = rowlab$row, limits = c(-0.5, max(ypos) + 0.5), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_fig4() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x.top = element_text(size = SZ_SMALL, lineheight = 0.95, colour = "black", margin = margin(b = SZ * 0.3)),
        axis.text.y = element_text(size = SZ_SMALL + 0.5, colour = "black"))
p <- preserve_header(p)
if (!build_only) {
  save_figure(draw_single_panel(p, "C", DS_TAG$S4C, PANEL_MM$S4C[2]), RES(PANEL), "FigureS2_PanelC_immune_escape_robustness",
         PANEL_MM$S4C[1], PANEL_MM$S4C[2],
         highlight_color = DIV_LOW)
}
