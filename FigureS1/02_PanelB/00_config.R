if (!exists(".script_dir", inherits = FALSE) || !is.character(.script_dir) ||
    length(.script_dir) != 1L || is.na(.script_dir))
  stop("Define .script_dir before loading 00_config.R.")
CODE_DIR  <- normalizePath(.script_dir)
PANEL_DIR <- normalizePath(CODE_DIR)
INDIR     <- Sys.getenv("SCI_INPUT_DIR", file.path(PANEL_DIR, "inputs"))
OUTDIR    <- Sys.getenv("SCI_OUTPUT_DIR", file.path(PANEL_DIR, "results"))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
IN_  <- function(...) file.path(INDIR,  ...)
OUT_ <- function(...) file.path(OUTDIR, ...)
FIG_FONT <- Sys.getenv("FIG_FONT", "Arial")
SZ_TICK       <- 8
SZ_TITLE      <- SZ_TICK * 1.1
SZ_TAG        <- SZ_TICK * 1.4
SZ_SMALL      <- SZ_TICK
DATA_TAG_SIZE <- SZ_TICK
PT <- 72.27 / 25.4
LW_DATA <- 1.1
LW_ERR  <- 0.65
LW_AXIS <- 0.55
SEM <- c(
  case    = "#E64B35",
  control = "#3C5488",
  up      = "#DC0000",
  down    = "#3C5488",
  ns      = "#8491B4",
  neutral = "#8A9197",
  accent  = "#00A087",
  warn    = "#F39B7F"
)
LABEL_LEVELS <- c("Tumor", "Ambiguous", "Normal")
LABEL_COLS   <- c(Tumor = "#D8503C", Ambiguous = "#BDBDBD", Normal = "#3C7DB1")
DISEASE_LEVELS <- c("BL", "DLBCL", "FL", "MZL")
DISEASE_COLS   <- c(BL = "#B35806", DLBCL = "#F1A340", FL = "#998EC3", MZL = "#542788")
DISEASE_GROUP  <- c(BL = "Aggressive", DLBCL = "Aggressive", FL = "Indolent", MZL = "Indolent")
GROUP_LEVELS   <- c("Aggressive", "Indolent")
GROUP_COLS     <- c(Aggressive = "#E08214", Indolent = "#8073AC")
COR_LOW    <- "#FFFFFF"
COR_HIGH   <- "#3A3A3A"
COR_LIMITS <- c(0, 1)
AX   <- c("Ig", "CNV", "MTP", "Prol", "Diff")
AX_D <- c("CNV", "MTP", "Prol", "Diff", "Ig")
SCORE_LONG <- c(Ig   = "Immunoglobulin composition score (Ig-score; 24 Ig indices, GBM probability)",
                CNV  = "Copy-number variation (infercnvpy, cluster-level)",
                MTP  = "Malignant transcriptional program score (UCell + GBM probability)",
                Prol = "Proliferation score (max of S and G2/M module scores)",
                Diff = "Differentiation potency (CCAT)")
theme_figS1 <- function(base_size = SZ_TICK, base_family = FIG_FONT) {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title    = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      axis.title    = ggplot2::element_text(size = SZ_TITLE, colour = "black"),
      axis.title.x  = ggplot2::element_text(margin = ggplot2::margin(t = SZ_TICK * 0.45)),
      axis.title.y  = ggplot2::element_text(margin = ggplot2::margin(r = SZ_TICK * 0.45)),
      axis.text     = ggplot2::element_text(size = SZ_TICK, colour = "black"),
      axis.text.x   = ggplot2::element_text(margin = ggplot2::margin(t = SZ_TICK * 0.30)),
      axis.text.y   = ggplot2::element_text(margin = ggplot2::margin(r = SZ_TICK * 0.30)),
      axis.line     = ggplot2::element_line(linewidth = LW_AXIS, colour = "black"),
      axis.ticks    = ggplot2::element_line(linewidth = LW_AXIS, colour = "black"),
      axis.ticks.length = grid::unit(SZ_TICK * 0.25, "pt"),
      legend.position   = "none",
      panel.grid        = ggplot2::element_blank(),
      plot.margin       = ggplot2::margin(SZ_TICK * 0.35, SZ_TICK * 0.45,
                                          SZ_TICK * 0.35, SZ_TICK * 0.35)
    )
}
fmtP <- function(p) vapply(p, function(q) {
  if (is.na(q)) "" else if (q < 1e-300) "< 1.00e-300" else sprintf("%.2e", q)
}, character(1))
P_txt <- function(p) ifelse(is.na(p), "", ifelse(p < 1e-300, "P < 1.00e-300", paste0("P = ", fmtP(p))))
star <- function(p) vapply(p, function(q) {
  if (is.na(q)) return("")
  if (q < 1e-4) "****" else if (q < 1e-3) "***" else if (q < 1e-2) "**" else if (q < .05) "*" else ""
}, character(1))
text_on <- function(fill) {
  rgb <- grDevices::col2rgb(fill) / 255
  L <- 0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ]
  ifelse(L > 0.55, "black", "white")
}
mar_pt <- function(m) ggplot2::margin(m[1], m[2], m[3], m[4], unit = "pt")
tag_plot <- function(p, letter, w_mm, h_mm, margin_pt = c(0, 0, 0, 0)) {
  m <- margin_pt * 25.4 / 72.27
  inner_w <- w_mm - m[2] - m[4]; inner_h <- h_mm - m[1] - m[3]
  p + ggplot2::labs(tag = letter) +
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = SZ_TAG, family = FIG_FONT,
                                                    hjust = 0, vjust = 1),
                   plot.tag.position = c((0.6 - m[4]) / inner_w, 1 - (0.6 - m[1]) / inner_h))
}
DATASETS_TXT <- "Datasets: 6 studies"
dataset_line <- function(p, txt) {
  assign("..dataset_txt", txt, envir = .GlobalEnv)
  p + ggplot2::labs(caption = " ") +
    ggplot2::theme(plot.caption = ggplot2::element_text(
      size = DATA_TAG_SIZE, colour = SEM[["neutral"]], family = FIG_FONT, hjust = 0),
      plot.caption.position = "plot")
}
draw_dataset_line <- function(txt, x_mm = 0.6, y_mm = 0.5) {
  gp <- grid::gpar(fontsize = DATA_TAG_SIZE, col = SEM[["neutral"]], fontfamily = FIG_FONT)
  desc <- grid::convertHeight(grid::grobDescent(grid::textGrob("gjpqy", gp = gp)), "mm")
  grid::grid.text(txt, x = grid::unit(x_mm, "mm"), y = grid::unit(y_mm, "mm") + desc,
                  just = c("left", "bottom"), gp = gp)
}
vector_colourbar <- function(height_mm = 20, width_mm = 3.5, nbin = 100) {
  ff <- names(formals(ggplot2::guide_colourbar))
  if (all(c("theme", "display") %in% ff)) {
    ggplot2::guide_colourbar(display = "rectangles", nbin = nbin,
      theme = ggplot2::theme(legend.key.height = grid::unit(height_mm, "mm"),
                             legend.key.width  = grid::unit(width_mm,  "mm")))
  } else {
    ggplot2::guide_colourbar(raster = FALSE, nbin = nbin,
      barheight = grid::unit(height_mm, "mm"), barwidth = grid::unit(width_mm, "mm"))
  }
}
SVG2PDF <- Sys.getenv("FIG_SVG2PDF", file.path(.release_root, "shared", "svg2pdf.py"))
PYTHON  <- local({
  cand <- c(Sys.getenv("FIG_PYTHON", ""),
            Sys.which("python3"), Sys.which("python"))
  for (p in cand[nzchar(cand)])
    if (suppressWarnings(system2(p, c("-c", shQuote("import reportlab, svglib")),
                                 stdout = FALSE, stderr = FALSE)) == 0) return(p)
  ""
})
PDF_BACKEND <- if (nzchar(PYTHON) && file.exists(SVG2PDF)) "svglite" else "cairo"
if (PDF_BACKEND != "svglite")
  warning("Using cairo PDF backend; optional SVG conversion support is unavailable.")
open_pdf <- function(path, width_mm, height_mm) {
  if (PDF_BACKEND == "svglite") {
    svg <- sub("\\.pdf$", ".svg", path)
    svglite::svglite(svg, width = width_mm / 25.4, height = height_mm / 25.4, bg = "white")
    assign("..svg_target", c(svg, path), envir = .GlobalEnv)
  } else if (capabilities("cairo")) {
    grDevices::cairo_pdf(path, width = width_mm / 25.4, height = height_mm / 25.4, family = FIG_FONT)
  } else {
    grDevices::pdf(path, width = width_mm / 25.4, height = height_mm / 25.4,
                   family = "Helvetica", useDingbats = FALSE)
  }
}
close_pdf <- function() {
  tgt <- if (exists("..svg_target", envir = .GlobalEnv)) get("..svg_target", envir = .GlobalEnv) else NULL
  grDevices::dev.off()
  if (!is.null(tgt)) {
    svg <- tgt[1]; pdf <- tgt[2]
    if (!file.exists(svg) || file.size(svg) == 0)
      stop("SVG output is empty: ", basename(svg), "\nPDF output was not refreshed.")
    if (!any(grepl("<text", readLines(svg, warn = FALSE), fixed = TRUE)))
      stop("SVG has no text content: ", basename(svg))
    ln <- readLines(svg, warn = FALSE)
    blank <- grepl("^<text [^>]*>[[:space:]]*</text>$", ln, useBytes = TRUE)
    if (any(blank)) writeLines(ln[!blank], svg, useBytes = TRUE)
    ok <- suppressWarnings(system2(PYTHON, c(shQuote(SVG2PDF), shQuote(svg), shQuote(pdf)),
                                   stdout = FALSE, stderr = FALSE)) == 0
    if (!ok || !file.exists(pdf) || file.size(pdf) == 0)
      stop("SVG to PDF conversion failed: ", basename(svg), " (SVG source remains available)")
    rm("..svg_target", envir = .GlobalEnv)
  }
}
check_canvas_bounds <- function(pdf, tol_pt = 0.5, raise_error = TRUE) {
  if (!nzchar(PYTHON) || !file.exists(pdf)) return(invisible(TRUE))
  py <- sprintf(paste(
    "import sys,json",
    "try:",
    "    import pymupdf as fitz",
    "except ImportError:",
    "    import fitz",
    "d=fitz.open(sys.argv[1]); bad=[]",
    "for pi,pg in enumerate(d):",
    "    R=pg.rect; P=300; mb=pg.mediabox",
    "    pg.set_mediabox(fitz.Rect(mb.x0-P,mb.y0-P,mb.x1+P,mb.y1+P))",
    "    for b in pg.get_text('dict')['blocks']:",
    "        for l in b.get('lines',[]):",
    "            for s in l['spans']:",
    "                x0,y0,x1,y1=[v-P for v in s['bbox']]",
    "                if x1>R.x1+%f or x0<R.x0-%f or y1>R.y1+%f or y0<R.y0-%f:",
    "                    bad.append((pi+1,round(x1-R.x1,2),round(R.x0-x0,2),s['text'][:60]))",
    "print(json.dumps(bad,ensure_ascii=False))", sep = "\n"), tol_pt, tol_pt, tol_pt, tol_pt)
  out <- suppressWarnings(system2(PYTHON, c("-c", shQuote(py), shQuote(pdf)), stdout = TRUE))
  bad <- tryCatch(jsonlite::fromJSON(paste(out, collapse = "")), error = function(e) NULL)
  if (is.null(bad)) { invisible(NULL)
                      return(invisible(TRUE)) }
  if (length(bad) == 0 || (is.matrix(bad) && nrow(bad) == 0)) {
    invisible(NULL)
    return(invisible(TRUE))
  }
  invisible(NULL); print(bad)
  if (raise_error) stop("Text exceeds canvas bounds: ", basename(pdf),
                 " (adjust the layout before exporting)")
  invisible(FALSE)
}
save_panel <- function(plot_obj, name, width_mm, height_mm, dpi = 600) {
  pdf_f <- OUT_(paste0(name, ".pdf"))
  png_f <- OUT_(paste0(name, ".png"))
  if (!exists("..dataset_txt", envir = .GlobalEnv))
    stop("Call dataset_line() before save_panel().")
  ds_txt <- get("..dataset_txt", envir = .GlobalEnv)
  rm("..dataset_txt", envir = .GlobalEnv)
  open_pdf(pdf_f, width_mm, height_mm); print(plot_obj); draw_dataset_line(ds_txt); close_pdf()
  if (capabilities("cairo")) {
    grDevices::png(png_f, width = width_mm / 25.4, height = height_mm / 25.4,
                   units = "in", res = dpi, family = FIG_FONT, type = "cairo")
  } else {
    grDevices::png(png_f, width = width_mm / 25.4, height = height_mm / 25.4,
                   units = "in", res = dpi)
  }
  print(plot_obj); draw_dataset_line(ds_txt); grDevices::dev.off()
  invisible(NULL)
  check_canvas_bounds(pdf_f)
  invisible(c(pdf_f, png_f))
}
