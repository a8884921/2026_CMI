#!/usr/bin/env Rscript
# Builds www/downloads/rctd_test_data.zip: the synthetic RCTD fixture of the original release (v1.2.0;
# examples/spatial/rctd_bundle_demo.zip) in the input formats of www/downloads/rctd_mal_states.R, with
# the weights v1.2.0 computed for it (tests/reference/v1.2.0_spatial_rctd_results.csv).
#   Rscript tools/make_rctd_test_data.R .
args <- commandArgs(trailingOnly = TRUE)
app <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
suppressPackageStartupMessages(library(Matrix))
tmp <- tempfile("rctd_test_"); dir.create(tmp)
utils::unzip(file.path(app, "examples", "spatial", "rctd_bundle_demo.zip"), exdir = file.path(tmp, "bundle"))
b <- function(f) list.files(file.path(tmp, "bundle"), pattern = paste0("^", f, "$"), recursive = TRUE, full.names = TRUE)[1]
read_mat <- function(f) { t <- utils::read.csv(b(f), check.names = FALSE); m <- as.matrix(t[, -1]); rownames(m) <- t[[1]]; Matrix(m, sparse = TRUE) }

out <- file.path(tmp, "rctd_test_data"); dir.create(out)
gz <- function(p) { x <- readBin(p, "raw", file.info(p)$size); con <- gzfile(paste0(p, ".gz"), "wb"); writeBin(x, con); close(con); unlink(p) }
write_10x <- function(m, dir) {
  dir.create(dir, recursive = TRUE)
  writeMM(methods::as(m, "CsparseMatrix"), file.path(dir, "matrix.mtx")); gz(file.path(dir, "matrix.mtx"))
  con <- gzfile(file.path(dir, "features.tsv.gz"), "w"); writeLines(paste(rownames(m), rownames(m), "Gene Expression", sep = "\t"), con); close(con)
  con <- gzfile(file.path(dir, "barcodes.tsv.gz"), "w"); writeLines(colnames(m), con); close(con)
}

# Reference: counts as a 10x folder and the label table (cell, cell_type, CMI)
ref <- read_mat("reference_counts.csv")
md <- utils::read.csv(b("reference_metadata.csv"), colClasses = "character")
write_10x(ref, file.path(out, "reference_counts"))
malignant <- grepl("^Malignant", md$celltype_anno)
utils::write.csv(data.frame(cell = md$cell, cell_type = ifelse(malignant, "Malignant", md$celltype_anno), CMI = ifelse(malignant, md$CMI, "")),
                 file.path(out, "reference_labels.csv"), row.names = FALSE, quote = FALSE)

# Section: a Space Ranger-style outs folder (filtered matrix and tissue_positions.csv)
spc <- read_mat("spatial_counts.csv")
sm <- utils::read.csv(b("spatial_metadata.csv"), colClasses = "character")
write_10x(spc, file.path(out, "section", "outs", "filtered_feature_bc_matrix"))
dir.create(file.path(out, "section", "outs", "spatial"))
utils::write.csv(data.frame(barcode = sm$spot_id, in_tissue = 1L, array_row = sm$y, array_col = sm$x, pxl_row_in_fullres = sm$y, pxl_col_in_fullres = sm$x),
                 file.path(out, "section", "outs", "spatial", "tissue_positions.csv"), row.names = FALSE, quote = FALSE)

# The weights the original release (v1.2.0) computed for this fixture
v <- utils::read.csv(file.path(app, "tests", "reference", "v1.2.0_spatial_rctd_results.csv"), check.names = FALSE)
w <- v[, grep("^RCTD_weight_", names(v))]
names(w) <- sub("^Malig_([123])$", "Mal-\\1", sub("^RCTD_weight_", "", names(w)))
utils::write.csv(data.frame(spot_id = v$spot_id, x = v$x, y = v$y, w, check.names = FALSE), file.path(out, "expected_weights.csv"), row.names = FALSE)

writeLines(c(
  "Test data for rctd_mal_states.R (CMI Explorer 2.0.0)",
  "",
  "A synthetic fixture: 300 reference cells (150 malignant, with CMI values in Mal-1, Mal-2 and Mal-3,",
  "and 50 each of CD4_T, Macrophage and Stromal), 2,400 genes, and a section of 36 spots. It is the RCTD",
  "test fixture of the original application (v1.2.0); the values are not biology. It has no class of",
  "normal B cells, so the script prints a warning about that.",
  "",
  "  reference_counts/        10x folder (matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz)",
  "  reference_labels.csv     cell, cell_type (Malignant or the class name), CMI",
  "  section/outs/            Space Ranger-style folder: filtered_feature_bc_matrix/ and spatial/tissue_positions.csv",
  "  expected_weights.csv     the weights v1.2.0 computed for this section",
  "",
  "Run (about a minute):",
  "  Rscript rctd_mal_states.R --ref-counts rctd_test_data/reference_counts \\",
  "      --ref-labels rctd_test_data/reference_labels.csv --spatial rctd_test_data/section/outs \\",
  "      --output test_weights.csv",
  "",
  "Compare (in R):",
  "  a <- read.csv(\"test_weights.csv\", check.names = FALSE); e <- read.csv(\"rctd_test_data/expected_weights.csv\", check.names = FALSE)",
  "  max(abs(as.matrix(a[, -(1:3)]) - as.matrix(e[match(a$spot_id, e$spot_id), names(a)[-(1:3)]])))   # about 1e-13",
  "",
  "test_weights.csv can then be uploaded to Spatial states (choose \"Other or unknown\" as spot technology)."),
  file.path(out, "README.txt"))

dest <- file.path(app, "www", "downloads", "rctd_test_data.zip")
unlink(dest)
files <- list.files(out, recursive = TRUE)
zip::zip(dest, files = file.path("rctd_test_data", files), root = tmp, mode = "mirror")
cat(sprintf("Wrote %s (%d files, %.0f KB)\n", dest, length(files), file.size(dest) / 1024))
