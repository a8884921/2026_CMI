#!/usr/bin/env Rscript
# =====================================================================================================
# rctd_mal_states.R - deconvolve one spatial section into Mal-1, Mal-2, Mal-3 and non-malignant cell
# types with RCTD, and write the weights table that the Spatial page of CMI Explorer reads.
# CMI Explorer 2.0.0.
# =====================================================================================================
#
# What the script does
#   1. Reads a single-cell reference (raw UMI counts) and a label for each cell. Cells labeled
#      "Malignant" (or "Tumor") are assigned a CMI grade from their CMI, at the fixed cut points
#      (the same rule as the Single-cell page):
#
#        CMI                                   reference class
#        below 0.9755057544904363              not malignant: left out (see below)
#        0.9755057544904363 to 0.98926105      Mal-1
#        above 0.98926105, up to 0.99619844    Mal-2
#        above 0.99619844                      Mal-3
#
#      A value equal to a boundary belongs to the lower grade; a value equal to the threshold is malignant.
#      Cells labeled malignant whose CMI is below the threshold are left out of the reference;
#      --below-threshold Mal-1 counts them as Mal-1 instead.
#   2. Reads one spatial section: a Space Ranger output folder, or counts plus a table of coordinates.
#   3. Runs RCTD (spacexr 2.2.1) in full mode with the settings that the local edition of CMI Explorer
#      also uses (all other arguments are spacexr defaults; see
#      run_rctd() below).
#   4. Writes the weights table (spot_id, x, y and one weight per class; weights sum to 1 per spot) and,
#      next to it, the reference classes used, the spots that were not analyzed, and a log.
#
# Install (once; R 4.3 or later)
#   install.packages(c("Matrix", "remotes"))
#   remotes::install_github("dmcable/spacexr@9f5dc33c8060f946c6072a138b70e189636e1435")   # spacexr 2.2.1
#
# Run
#   Rscript rctd_mal_states.R --ref-counts reference_counts.rds --ref-labels reference_labels.csv \
#       --spatial section1/outs --output section1_weights.csv --cores 4
#   or, in R:
#   source("rctd_mal_states.R")
#   cmi_rctd(ref_counts = "reference_counts.rds", ref_labels = "reference_labels.csv",
#            spatial = "section1/outs", output = "section1_weights.csv", cores = 4)
#
# Inputs
#   --ref-counts  raw UMI counts of the reference, genes x cells, as
#                 - a 10x folder (matrix.mtx, features.tsv or genes.tsv, barcodes.tsv; gzipped files are fine),
#                 - a ZIP file of such a folder (for example written by h5ad_to_mtx.py from an .h5ad file),
#                 - an .rds file with a matrix or sparse matrix, or with a Seurat object (the counts of
#                   its RNA assay; see --assay),
#                 - a CSV or TSV file (gene in the first column; small references only).
#   --ref-labels  CSV with one row per reference cell and the columns
#                 cell       the cell barcode, as in the counts (unique across samples);
#                 cell_type  "Malignant" for malignant B cells; any other value is used as the class name,
#                            for example CD4_T, CD8_T, Dendritic, Endothelial, Macrophage, Mast, NK,
#                            Normal_B and Stromal. Keep a class of normal B cells: without one, RCTD
#                            assigns normal B-cell signal to Mal-1.
#                 CMI        the CMI of each malignant cell: the CMI column of cells.csv from the
#                            Single-cell CMI workspace (one upload per sample) or from its command line.
#                 For a Seurat .rds the table can be left out: the metadata columns named by
#                 --label-column and --cmi-column are used instead.
#   --spatial     the Space Ranger "outs" folder of one 10x Visium section (filtered_feature_bc_matrix/ and
#                 spatial/tissue_positions.csv, or tissue_positions_list.csv of Space Ranger 1.x), or the
#                 counts of the section (genes x spots) in any format accepted by --ref-counts, with --coords
#   --coords      CSV with spot_id, x and y (needed when --spatial is not a Space Ranger folder)
#   --output      the weights table to write (.csv)
# Options
#   --below-threshold  cells labeled Malignant whose CMI is below the malignancy threshold: "exclude"
#                      (default; left out of the reference) or a class name, for example Normal_B, or
#                      Mal-1 as in version 1.2.0
#   --malignant-label  regular expression that marks malignant cells in cell_type (default: Malignant or Tumor)
#   --assay            assay of a Seurat object that holds the raw counts (default RNA when present)
#   --label-column     metadata column with the cell types (Seurat input; default cell_type)
#   --cmi-column       metadata column with the CMI (Seurat input; default CMI)
#   --cores            CPU cores for RCTD (default 1; memory use grows with the number of cores)
#   --seed             random seed (default 20260101, as in CMI Explorer)
#
# Notes
#   - Reference and section must use the same gene identifiers (10x folders: the gene symbols of
#     features.tsv, column 2) and must share at least 1,000 genes. Seurat writes "_" in gene names as
#     "-"; the gene names of the section are matched accordingly when the reference is a Seurat object.
#   - Every reference class needs at least 25 cells with at least 100 UMIs, and Mal-1, Mal-2 and Mal-3
#     must all be present; at most 10,000 cells per class are used. Spots with fewer than 100 UMIs and
#     spots outside the tissue (in_tissue = 0) are not analyzed.
#   - Coordinates of a Space Ranger folder: x = pxl_col_in_fullres and y = pxl_row_in_fullres (pixels of
#     the full-resolution image; the maps of the Spatial states workspace put y = 0 at the top, as in
#     the image). Visium HD bins are not read from Space Ranger folders; pass their counts and --coords.
#     Their tables usually exceed the 20,000 spots of the hosted site; use the local edition for them.
#   - Memory: RCTD needs a few GB for a Visium section (more for large references); use a computer with
#     at least 8 GB.
# =====================================================================================================

CMI_THRESHOLD     <- 0.9755057544904363   # malignancy threshold
CMI_MAL1_MAL2     <- 0.98926105           # Mal-1 / Mal-2 boundary
CMI_MAL2_MAL3     <- 0.99619844           # Mal-2 / Mal-3 boundary
SCRIPT_VERSION    <- "rctd_mal_states.R (CMI Explorer 2.0.0)"
MALIGNANT_LABEL   <- "^(Malignant|Tumou?r)(?:$|[ _:-])"
STATES            <- c("Mal-1", "Mal-2", "Mal-3")


# Malignant state of each CMI value (NA stays NA).
mal_state <- function(cmi) {
  cmi <- as.numeric(cmi)
  out <- ifelse(cmi < CMI_THRESHOLD, "Below threshold", ifelse(cmi <= CMI_MAL1_MAL2, "Mal-1", ifelse(cmi <= CMI_MAL2_MAL3, "Mal-2", "Mal-3")))
  out[is.na(cmi)] <- NA_character_
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a
say <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
fail <- function(...) stop(paste0(...), call. = FALSE)

# ---- reading counts ---------------------------------------------------------------------------------
find_file <- function(dir, names) {
  for (n in names) for (p in file.path(dir, c(n, paste0(n, ".gz")))) if (file.exists(p)) return(p)
  NA_character_
}
read_10x_dir <- function(dir) {
  mtx <- find_file(dir, "matrix.mtx"); feat <- find_file(dir, c("features.tsv", "genes.tsv")); bc <- find_file(dir, "barcodes.tsv")
  if (anyNA(c(mtx, feat, bc))) {
    if (length(list.files(dir, pattern = "filtered_feature_bc_matrix\\.h5$")))
      fail(dir, " has only filtered_feature_bc_matrix.h5. This script reads the filtered_feature_bc_matrix folder that Space Ranger writes next to it; copy that folder back, or read the .h5 file in R (Seurat::Read10X_h5) and pass the counts as .rds with --coords.")
    fail("No 10x matrix in ", dir, ": it needs matrix.mtx, features.tsv (or genes.tsv) and barcodes.tsv.")
  }
  m <- methods::as(Matrix::readMM(gzfile(mtx)), "CsparseMatrix")
  features <- utils::read.delim(gzfile(feat), header = FALSE, colClasses = "character", quote = "")
  barcodes <- readLines(gzfile(bc))
  genes <- if (ncol(features) >= 2) features[[2]] else features[[1]]
  if (nrow(m) != length(genes) || ncol(m) != length(barcodes))
    fail("The 10x files in ", dir, " do not match: the matrix has ", nrow(m), " rows and ", ncol(m), " columns.")
  dimnames(m) <- list(make.unique(genes), barcodes)
  m
}
# The assay of a Seurat object that holds the raw counts: the one asked for, else RNA, else the default.
seurat_assay <- function(obj, assay = NULL) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) fail("Reading a Seurat object needs the SeuratObject package.")
  have <- SeuratObject::Assays(obj)
  a <- assay %||% (if ("RNA" %in% have) "RNA" else SeuratObject::DefaultAssay(obj))
  if (!a %in% have) fail("The Seurat object has no assay ", a, " (it has ", paste(have, collapse = ", "), "); set --assay.")
  a
}
seurat_counts <- function(obj, assay = NULL) {
  a <- seurat_assay(obj, assay)
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    layers <- SeuratObject::Layers(obj[[a]], search = "counts")
    if (!length(layers)) fail("Assay ", a, " of the Seurat object has no counts layer (integrated assays hold no raw counts); set --assay to the assay with the raw counts, usually RNA.")
    if (length(layers) > 1) fail("The counts of assay ", a, " are split into ", length(layers), " layers (one per sample); join them first: obj <- SeuratObject::JoinLayers(obj, assay = \"", a, "\")")
    m <- SeuratObject::LayerData(obj, assay = a, layer = "counts")
  } else m <- SeuratObject::GetAssayData(obj, assay = a, slot = "counts")
  if (!nrow(m) || !ncol(m)) fail("Assay ", a, " of the Seurat object holds no counts; set --assay to the assay with the raw counts, usually RNA.")
  if (identical(a, "SCT")) say("Note: assay SCT holds SCTransform-corrected counts, not raw counts; use raw counts (--assay RNA).")
  methods::as(m, "CsparseMatrix")
}
counts_from_object <- function(x, what, assay = NULL) {
  if (inherits(x, "Seurat")) return(seurat_counts(x, assay))
  if (is.matrix(x) || inherits(x, "Matrix")) return(methods::as(x, "CsparseMatrix"))
  fail(what, ": the .rds file holds a ", class(x)[1], ", not a matrix or a Seurat object.")
}
read_counts <- function(path, what = "counts", obj = NULL, assay = NULL) {
  if (!is.null(obj)) return(counts_from_object(obj, what, assay))   # an .rds file that was already read
  if (!file.exists(path)) fail(what, ": ", path, " does not exist.")
  if (dir.exists(path)) {
    if (dir.exists(file.path(path, "filtered_feature_bc_matrix"))) path <- file.path(path, "filtered_feature_bc_matrix")
    return(read_10x_dir(path))
  }
  if (grepl("\\.zip$", path, ignore.case = TRUE)) {
    td <- tempfile("cmi_rctd_"); utils::unzip(path, exdir = td)
    hits <- list.files(td, pattern = "^matrix\\.mtx(\\.gz)?$", recursive = TRUE, full.names = TRUE)
    if (length(hits) != 1) fail(what, ": the ZIP file must contain exactly one matrix.mtx.")
    return(read_10x_dir(dirname(hits)))
  }
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(counts_from_object(readRDS(path), what, assay))
  if (grepl("\\.(csv|tsv|txt)(\\.gz)?$", path, ignore.case = TRUE)) {
    tab <- if (requireNamespace("data.table", quietly = TRUE)) data.table::fread(path, data.table = FALSE, check.names = FALSE)
           else utils::read.csv(path, check.names = FALSE, sep = if (grepl("\\.tsv", path, ignore.case = TRUE)) "\t" else ",")
    m <- as.matrix(tab[, -1, drop = FALSE]); rownames(m) <- make.unique(as.character(tab[[1]]))
    return(methods::as(Matrix::Matrix(m, sparse = TRUE), "CsparseMatrix"))
  }
  if (grepl("\\.h5$", path, ignore.case = TRUE)) fail(what, ": .h5 files are not read; use the matching 10x folder (Space Ranger and Cell Ranger write both).")
  fail(what, ": ", path, " is not a 10x folder, ZIP, .rds or CSV/TSV file.")
}

# ---- reading labels and coordinates -----------------------------------------------------------------
pick <- function(tab, names, what, required = TRUE) {
  key <- tolower(gsub("[^A-Za-z0-9]", "", names(tab)))
  hit <- match(tolower(gsub("[^A-Za-z0-9]", "", names)), key); hit <- hit[!is.na(hit)][1]
  if (is.na(hit) && required) fail(what, " needs a column named ", names[1], ".")
  if (is.na(hit)) NULL else tab[[hit]]
}
read_labels <- function(path, obj = NULL, label_column = "cell_type", cmi_column = "CMI") {
  if (is.null(path) || !nzchar(path)) {
    if (!inherits(obj, "Seurat")) fail("--ref-labels is required: a CSV with the columns cell, cell_type and CMI.")
    md <- obj[[]]
    if (!label_column %in% names(md)) fail("The Seurat metadata has no column ", label_column, " (set --label-column).")
    cmi <- if (cmi_column %in% names(md)) md[[cmi_column]] else rep(NA_real_, nrow(md))
    return(data.frame(cell = rownames(md), cell_type = as.character(md[[label_column]]), CMI = suppressWarnings(as.numeric(as.character(cmi))), stringsAsFactors = FALSE))
  }
  tab <- utils::read.csv(path, check.names = FALSE, colClasses = "character")
  data.frame(cell = pick(tab, c("cell", "barcode", "cell_id"), "The label table"),
             cell_type = pick(tab, unique(c(label_column, "cell_type", "celltype", "celltype_anno", "class")), "The label table"),
             CMI = suppressWarnings(as.numeric(pick(tab, unique(c(cmi_column, "CMI")), "The label table", required = FALSE) %||% NA_real_)),
             stringsAsFactors = FALSE)
}
# Spot coordinates (and whether each spot lies in the tissue, when the Space Ranger table says so).
read_coords <- function(spatial, coords = NULL) {
  if (!is.null(coords) && nzchar(coords)) {
    tab <- utils::read.csv(coords, check.names = FALSE, colClasses = "character")
    it <- pick(tab, "in_tissue", "The coordinate table", required = FALSE)
    return(data.frame(x = as.numeric(pick(tab, c("x", "pxl_col_in_fullres", "imagecol"), "The coordinate table")),
                      y = as.numeric(pick(tab, c("y", "pxl_row_in_fullres", "imagerow"), "The coordinate table")),
                      in_tissue = if (is.null(it)) TRUE else it != "0",
                      row.names = pick(tab, c("spot_id", "barcode", "spot"), "The coordinate table")))
  }
  sp_dir <- if (dir.exists(file.path(spatial, "spatial"))) file.path(spatial, "spatial") else file.path(dirname(spatial), "spatial")
  f2 <- file.path(sp_dir, "tissue_positions.csv"); f1 <- file.path(sp_dir, "tissue_positions_list.csv")
  pos <- if (file.exists(f2)) utils::read.csv(f2, colClasses = "character")
         else if (file.exists(f1)) stats::setNames(utils::read.csv(f1, header = FALSE, colClasses = "character"),
                                                   c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"))
         else fail("No spot coordinates: ", spatial, " is not a Space Ranger outs folder (spatial/tissue_positions.csv). Give --coords.")
  data.frame(x = as.numeric(pos$pxl_col_in_fullres), y = as.numeric(pos$pxl_row_in_fullres), in_tissue = pos$in_tissue != "0", row.names = pos$barcode)
}

# ---- reference classes --------------------------------------------------------------------------------
# Malignant cells take their state from the CMI; every other label is kept as its own class. Cells
# without a label are not used.
reference_classes <- function(labels, malignant_label = MALIGNANT_LABEL, below_threshold = "exclude") {
  type <- trimws(as.character(labels$cell_type)); type[is.na(type)] <- ""
  empty <- !nzchar(type)
  malignant <- !empty & grepl(malignant_label, type, ignore.case = TRUE, perl = TRUE)
  if (!any(malignant)) fail("No reference cell is labeled Malignant (column cell_type; see --malignant-label).")
  key <- tolower(gsub("[^A-Za-z0-9]", "", type))
  reserved <- !malignant & grepl("^(mal|malig|malignant|tumou?r)(state)?[123]$", key)
  if (any(reserved))
    fail("Mal-1, Mal-2 and Mal-3 are assigned by this script from the CMI, so ", type[reserved][1], " cannot be a class name. Label those cells Malignant and give their CMI.")
  cmi <- labels$CMI
  bad <- malignant & !(is.finite(cmi) & cmi >= 0 & cmi <= 1)
  if (any(bad)) fail(sum(bad), " malignant cells have no CMI between 0 and 1 (for example ", labels$cell[which(bad)[1]], ").")
  state <- mal_state(cmi)
  below <- malignant & state == "Below threshold"
  type[malignant] <- state[malignant]
  if (any(below)) type[below] <- if (identical(below_threshold, "exclude")) NA_character_ else below_threshold
  type[empty] <- NA_character_
  data.frame(cell = labels$cell, class = type, CMI = cmi, malignant = malignant, below_threshold = below, unlabeled = empty, stringsAsFactors = FALSE)
}
has_normal_b <- function(classes) any(grepl("(^|[^a-z])b($|[^a-z])|bcell|b-cell|normalb|naiveb|memoryb|plasma", tolower(setdiff(classes, STATES))))

# ---- RCTD -------------------------------------------------------------------------------------------
# The steps and settings of the local edition of CMI Explorer:
# spots with at least 100 UMIs; genes present in both data sets (at least 1,000); a reference of at most
# 10,000 cells per class, each with at least 100 UMIs; RCTD in full mode. Every other argument is the
# spacexr 2.2.1 default and is written out for the record.
run_rctd <- function(ref, types, spc, xy, cores = 1, seed = 20260101) {
  spc <- spc[, Matrix::colSums(spc) >= 100, drop = FALSE]
  if (!ncol(spc)) fail("No spot has 100 or more UMIs.")
  common <- intersect(rownames(ref), rownames(spc))
  if (length(common) < 1000L) fail("The reference and the section share ", length(common), " genes; at least 1,000 are needed. Use the same gene identifiers (symbols) in both.")
  set.seed(seed)
  reference <- spacexr::Reference(counts = ref[common, ], cell_types = types[colnames(ref)], n_max_cells = 10000, min_UMI = 100)
  puck <- spacexr::SpatialRNA(coords = xy[colnames(spc), c("x", "y"), drop = FALSE], counts = spc[common, ])
  obj <- spacexr::create.RCTD(spatialRNA = puck, reference = reference, max_cores = cores, test_mode = FALSE,
                              gene_cutoff = 0.000125, fc_cutoff = 0.5, gene_cutoff_reg = 2e-04, fc_cutoff_reg = 0.75,
                              UMI_min = 100, UMI_max = 2e+07, counts_MIN = 10, UMI_min_sigma = 300,
                              CELL_MIN_INSTANCE = 25, MAX_MULTI_TYPES = 4)
  obj <- spacexr::run.RCTD(obj, doublet_mode = "full")
  W <- as.matrix(obj@results$weights)
  if (!nrow(W) || any(!is.finite(W)) || any(W < 0) || any(rowSums(W) <= 0)) fail("RCTD returned invalid weights.")
  list(weights = W / rowSums(W), genes = length(common), used = table(reference@cell_types))
}

# ---- the whole run ----------------------------------------------------------------------------------
cmi_rctd <- function(ref_counts, ref_labels = NULL, spatial, output, coords = NULL, below_threshold = "exclude",
                     malignant_label = MALIGNANT_LABEL, assay = NULL, label_column = "cell_type", cmi_column = "CMI",
                     cores = 1, seed = 20260101) {
  if (!requireNamespace("spacexr", quietly = TRUE))
    fail("RCTD needs spacexr: remotes::install_github(\"dmcable/spacexr@9f5dc33c8060f946c6072a138b70e189636e1435\")")
  if (!identical(as.character(utils::packageVersion("spacexr")), "2.2.1"))
    say("Note: CMI Explorer was tested with spacexr 2.2.1; this is spacexr ", as.character(utils::packageVersion("spacexr")), ".")
  t0 <- Sys.time(); notes <- character()
  note <- function(...) { m <- paste0(...); say(m); notes <<- c(notes, m) }

  say("Reading the reference")
  ref_obj <- if (grepl("\\.rds$", ref_counts, ignore.case = TRUE) && file.exists(ref_counts)) readRDS(ref_counts) else NULL
  from_seurat <- inherits(ref_obj, "Seurat")
  ref <- read_counts(ref_counts, "--ref-counts", obj = ref_obj, assay = assay)
  if (from_seurat) note("Reference counts: assay ", seurat_assay(ref_obj, assay), " of the Seurat object")
  dup <- colnames(ref)[duplicated(colnames(ref))]
  if (length(dup)) fail("Cell barcodes of the reference are not unique (for example ", dup[1], "); add the sample to each barcode, as Seurat does when merging samples.")
  lab <- read_labels(ref_labels, ref_obj, label_column, cmi_column); rm(ref_obj)
  dup <- lab$cell[duplicated(lab$cell)]
  if (length(dup)) fail("The label table lists cell ", dup[1], " more than once.")
  found <- lab$cell %in% colnames(ref)
  if (!any(found)) fail("No cell of the label table is in the reference counts (compare its cell column with the barcodes).")
  if (sum(!found) > 0.01 * nrow(lab))
    fail(sprintf("%s of the %s cells of the label table are not in the reference counts (for example %s). If the counts come from a Seurat 5 object, join its layers first (JoinLayers) so that every sample is included.",
                 format(sum(!found), big.mark = ","), format(nrow(lab), big.mark = ","), lab$cell[!found][1]))
  if (any(!found)) note(sum(!found), " cells of the label table are not in the reference counts and are not used")
  lab <- lab[found, , drop = FALSE]
  all_cls <- reference_classes(lab, malignant_label, below_threshold)
  n_below <- sum(all_cls$below_threshold); n_unlabeled <- sum(all_cls$unlabeled)
  if (n_unlabeled) note(n_unlabeled, " cells without a cell_type are not used")
  if (n_below) note(sprintf("%d cells labeled malignant are below the malignancy threshold (CMI %.4f): %s", n_below, CMI_THRESHOLD,
                            if (identical(below_threshold, "exclude")) "left out of the reference" else paste("used as", below_threshold)))
  cls <- all_cls[!is.na(all_cls$class), , drop = FALSE]
  ref <- ref[, colnames(ref) %in% cls$cell, drop = FALSE]          # keep the order of the counts
  cls <- cls[match(colnames(ref), cls$cell), , drop = FALSE]
  if (any(abs(ref@x - round(ref@x)) > 1e-8)) fail("The reference must hold raw integer UMI counts, not normalized values.")
  types <- factor(cls$class); names(types) <- cls$cell
  missing_states <- setdiff(STATES, levels(types))
  if (length(missing_states)) fail("The reference has no ", paste(missing_states, collapse = " or "), " cells. Mal-1, Mal-2 and Mal-3 each need at least 25 cells (the Spatial page needs all three).")
  if (!has_normal_b(levels(types)))
    note("Warning: the reference has no class of normal B cells (for example Normal_B). Without them RCTD assigns the signal of normal B cells to Mal-1. Add normal B cells, or keep the cells below the threshold as Normal_B (--below-threshold Normal_B).")

  say("Reading the section")
  spc <- read_counts(spatial, "--spatial")
  if (anyDuplicated(colnames(spc))) fail("Spot barcodes of the section are not unique (for example ", colnames(spc)[duplicated(colnames(spc))][1], ").")
  if (any(abs(spc@x - round(spc@x)) > 1e-8)) fail("The section must hold raw integer UMI counts.")
  if (from_seurat && !any(grepl("_", rownames(ref))) && any(grepl("_", rownames(spc)))) {
    rownames(spc) <- make.unique(gsub("_", "-", rownames(spc)))
    note("Gene names of the section: \"_\" written as \"-\", as Seurat does")
  }
  xy <- read_coords(spatial, coords)
  miss <- setdiff(colnames(spc), rownames(xy))
  if (length(miss)) fail(length(miss), " spots have no coordinates (for example ", miss[1], ").")
  xy <- xy[colnames(spc), , drop = FALSE]
  if (any(!xy$in_tissue)) {
    note(sum(!xy$in_tissue), " spots outside the tissue (in_tissue = 0) are not analyzed")
    spc <- spc[, xy$in_tissue, drop = FALSE]; xy <- xy[xy$in_tissue, , drop = FALSE]
  }
  if (any(!is.finite(as.matrix(xy[, c("x", "y")])))) fail("Spot coordinates must be finite numbers.")
  common <- intersect(rownames(ref), rownames(spc))
  if (length(common) < 1000L) fail("The reference and the section share ", length(common), " genes; at least 1,000 are needed. Use the same gene identifiers (symbols) in both.")
  if (length(common) < 0.5 * nrow(spc))
    note(sprintf("Only %s of the %s genes of the section are in the reference; check that both use the same gene identifiers", format(length(common), big.mark = ","), format(nrow(spc), big.mark = ",")))
  enough <- table(factor(cls$class[Matrix::colSums(ref[common, , drop = FALSE]) >= 100], levels = levels(types)))
  if (any(enough < 25)) fail("Every reference class needs at least 25 cells with 100 or more UMIs over the genes shared with the section: ",
                             paste(sprintf("%s %d", names(enough)[enough < 25], enough[enough < 25]), collapse = ", "), ".")
  say(sprintf("Reference: %s cells in %d classes (%s); section: %s spots", format(ncol(ref), big.mark = ","), nlevels(types),
              paste(sprintf("%s %s", levels(types), format(as.integer(table(types)), big.mark = ",", trim = TRUE)), collapse = ", "), format(ncol(spc), big.mark = ",")))

  say("Running RCTD in full mode (this can take a while)")
  r <- run_rctd(ref, types, spc, xy, cores = cores, seed = seed)
  W <- r$weights
  out <- data.frame(spot_id = rownames(W), x = xy[rownames(W), "x"], y = xy[rownames(W), "y"], W, check.names = FALSE)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(out, output, row.names = FALSE)
  stem <- sub("\\.csv$", "", output, ignore.case = TRUE)
  used <- as.integer(r$used[levels(types)])
  utils::write.csv(data.frame(class = levels(types), cells_labeled = as.integer(table(types)), cells_used = used),
                   paste0(stem, "_reference_classes.csv"), row.names = FALSE)
  umi <- Matrix::colSums(spc); skipped <- setdiff(colnames(spc), rownames(W))
  utils::write.csv(data.frame(spot_id = skipped, umi = as.numeric(umi[skipped]),
                              reason = ifelse(umi[skipped] < 100, "fewer than 100 UMIs",
                                              "filtered by RCTD (fewer than 100 UMIs over the shared genes or fewer than 10 counts over the genes RCTD uses)")),
                   paste0(stem, "_spots_not_analyzed.csv"), row.names = FALSE)
  writeLines(c(SCRIPT_VERSION, sprintf("Run %s (%.0f s)", format(t0, "%Y-%m-%d %H:%M:%S %Z"), as.numeric(difftime(Sys.time(), t0, units = "secs"))),
               paste("Reference counts:", ref_counts), paste("Reference labels:", ref_labels %||% "metadata of the Seurat object"), paste("Section:", spatial),
               if (!is.null(coords)) paste("Coordinates:", coords),
               sprintf("Cut points (CMI): threshold %.16g, Mal-1/Mal-2 %.8f, Mal-2/Mal-3 %.8f", CMI_THRESHOLD, CMI_MAL1_MAL2, CMI_MAL2_MAL3),
               sprintf("Malignant cells below the threshold: %d (%s)", n_below, below_threshold),
               sprintf("Reference classes (cells used): %s", paste(sprintf("%s %d", levels(types), used), collapse = ", ")),
               sprintf("Shared genes: %d; spots analyzed: %d of %d", r$genes, nrow(W), ncol(spc)),
               "RCTD: doublet_mode full; Reference n_max_cells 10000, min_UMI 100; create.RCTD UMI_min 100, UMI_max 2e7, counts_MIN 10, UMI_min_sigma 300, CELL_MIN_INSTANCE 25, other arguments default",
               sprintf("Seed %s; cores %s", seed, cores),
               sprintf("R %s; spacexr %s; Matrix %s", paste(R.version$major, R.version$minor, sep = "."), utils::packageVersion("spacexr"), utils::packageVersion("Matrix")),
               if (length(notes)) c("Notes:", paste("  -", notes))),
             paste0(stem, "_log.txt"))
  say(sprintf("Wrote %s (%s spots, %d classes). Upload it to Spatial states.", output, format(nrow(W), big.mark = ","), ncol(W)))
  invisible(out)
}

# ---- command line -----------------------------------------------------------------------------------
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args) || any(args %in% c("-h", "--help"))) {
    cat("Usage: Rscript rctd_mal_states.R --ref-counts PATH --ref-labels labels.csv --spatial PATH --output weights.csv\n",
        "         [--coords coords.csv] [--below-threshold exclude|CLASS] [--malignant-label REGEX] [--assay RNA]\n",
        "         [--label-column cell_type] [--cmi-column CMI] [--cores 1] [--seed 20260101]\n",
        "The comments at the top of this file describe the inputs.\n", sep = "")
    return(invisible(NULL))
  }
  if (length(args) %% 2 || any(!grepl("^--", args[seq(1, length(args), 2)]))) fail("Arguments come in pairs: --name value.")
  o <- stats::setNames(as.list(args[seq(2, length(args), 2)]), sub("^--", "", args[seq(1, length(args), 2)]))
  known <- c("ref-counts", "ref-labels", "spatial", "coords", "output", "below-threshold", "malignant-label", "assay", "label-column", "cmi-column", "cores", "seed")
  if (length(setdiff(names(o), known))) fail("Unknown option: --", setdiff(names(o), known)[1])
  for (k in c("ref-counts", "spatial", "output")) if (is.null(o[[k]])) fail("--", k, " is required.")
  cmi_rctd(ref_counts = o[["ref-counts"]], ref_labels = o[["ref-labels"]], spatial = o[["spatial"]], output = o[["output"]], coords = o[["coords"]],
           below_threshold = o[["below-threshold"]] %||% "exclude", malignant_label = o[["malignant-label"]] %||% MALIGNANT_LABEL,
           assay = o[["assay"]], label_column = o[["label-column"]] %||% "cell_type", cmi_column = o[["cmi-column"]] %||% "CMI",
           cores = as.integer(o[["cores"]] %||% 1), seed = as.integer(o[["seed"]] %||% 20260101))
}
if (sys.nframe() == 0L) main()
