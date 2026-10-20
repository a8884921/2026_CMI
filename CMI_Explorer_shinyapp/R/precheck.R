# =============================================================================
# Upload check in the interface process, run as soon as a file arrives and before anything is queued.
# Reads only file listings, headers and the first rows, so it takes well under a second and little
# memory; the worker repeats every check on the whole file. Returns
#   list(ok, summary, problems, warnings, info)
# =============================================================================

PRECHECK_ROWS <- 200L

check_ok <- function(summary, warnings = character(), info = list()) list(ok = TRUE, summary = summary, problems = character(), warnings = warnings, info = info)
check_fail <- function(..., info = list()) list(ok = FALSE, summary = NULL, problems = paste0(...), warnings = character(), info = info)
# Run an engine check that signals input_error(); return its message or NULL.
check_msg <- function(expr) tryCatch({ force(expr); NULL }, input_error = function(e) conditionMessage(e))
fmt_count <- function(x) formatC(as.numeric(x), format = "f", digits = 0, big.mark = ",")

# gzfile() also reads uncompressed files, so one opener serves .csv and .csv.gz alike.
open_upload <- function(path) gzfile(path, "rt")
count_upload_lines <- function(path, chunk = 8 * 1024^2) {
  con <- gzfile(path, "rb"); on.exit(close(con))
  n <- 0; last <- as.raw(10)
  repeat { b <- readBin(con, "raw", chunk); if (!length(b)) break; n <- n + sum(b == as.raw(10)); last <- b[length(b)] }
  if (last != as.raw(10)) n <- n + 1
  n
}
head_lines <- function(path, n) { con <- open_upload(path); on.exit(close(con)); readLines(con, n = n, warn = FALSE) }

# Identifier type of genes in the first column (or features.tsv).
id_kind <- function(ids) {
  ids <- trimws(ids[!is.na(ids) & nzchar(trimws(ids))])
  if (!length(ids)) return("none")
  if (mean(grepl("^ENSMUSG[0-9]", ids)) > 0.5) return("mouse_ensembl")
  if (mean(grepl("^ENSG[0-9]+(\\.[0-9]+)?$", ids)) > 0.5) return("ensembl")
  if (mean(grepl("(_at|_st)$|^ILMN_[0-9]+$|^A_[0-9]+_P[0-9]+$|^[0-9]+$", ids)) > 0.5) return("probe")
  if (length(ids) >= 50 && mean(grepl("^[A-Z][a-z][a-z0-9]*$", ids)) > 0.5) return("mouse_symbol")
  "symbol"
}
ID_LABEL <- c(symbol = "gene symbols", ensembl = "Ensembl IDs")
# Genes found in almost every human single-cell or bulk B-cell profile; used to recognize human symbols.
COMMON_GENES <- c("ACTB", "B2M", "MALAT1", "CD74", "HLA-DRA", "MS4A1", "CD79A", "CD79B", "PTPRC", "RPLP0", "GAPDH", "EEF1A1",
                  "TMSB4X", "FTL", "FTH1", "HLA-A", "HLA-B", "TPT1", "RPS27A", "UBC")
id_problem <- function(kind, bulk = FALSE) switch(kind,
  mouse_ensembl = "The gene identifiers are mouse Ensembl IDs (ENSMUSG...). The CMI models use human genes.",
  mouse_symbol = "Most gene symbols are in mouse style (for example Cd19, Ms4a1). The CMI models use human gene symbols.",
  probe = "The row identifiers look like microarray probe or numeric IDs (for example 1007_s_at, ILMN_1343291 or 7157). Collapse probes to HGNC gene symbols first, for example by keeping the probe with the highest mean per gene.",
  ensembl = if (bulk) "The row identifiers are Ensembl IDs. Convert them to HGNC gene symbols before uploading." else NULL,
  none = "No gene identifiers were found in the first column.",
  NULL)

# Lines of one ZIP member, read without extracting the archive.
zip_member_lines <- function(zip, member, n = -1L) {
  con <- if (grepl("\\.gz$", member, ignore.case = TRUE)) gzcon(unz(zip, member, open = "rb")) else { k <- unz(zip, member); open(k, "r"); k }
  on.exit(close(con))
  readLines(con, n = n, warn = FALSE)
}

precheck_single_cell <- function(path, name, lim) {
  nm <- tolower(name)
  if (grepl("\\.zip$", nm)) return(precheck_mtx_zip(path, lim))
  if (!grepl("\\.(csv|tsv|txt)(\\.gz)?$", nm)) return(check_fail("Use a 10x Matrix Market ZIP (matrix.mtx, features.tsv, barcodes.tsv) or a CSV/TSV/TXT count matrix (optionally .gz)."))
  lines <- tryCatch(head_lines(path, PRECHECK_ROWS + 1L), error = function(e) character())
  if (length(lines) < 2L) return(check_fail("The file has no data rows."))
  sep <- if (grepl("\\.csv(\\.gz)?$", nm)) "," else if (grepl("\t", lines[1], fixed = TRUE)) "\t" else ","
  cols <- split_header(lines[1], sep); n_cells <- length(cols) - 1L
  if (n_cells < 1L) return(check_fail("Use one gene column followed by one column per cell."))
  if (n_cells > lim$max_cells) return(check_fail(sprintf("The matrix has %s cells; this server accepts up to %s per upload. Subsample cells or use the local edition.", fmt_count(n_cells), fmt_count(lim$max_cells))))
  tab <- tryCatch(data.table::fread(text = paste0(paste(lines[-1], collapse = "\n"), "\n"), sep = sep, header = FALSE, data.table = FALSE, showProgress = FALSE,
                                    colClasses = list(character = 1L)), error = function(e) NULL)
  if (is.null(tab)) return(check_fail("The first rows could not be parsed as a table."))
  if (ncol(tab) != n_cells + 1L) return(check_fail(sprintf("The header has %d fields but the first data row has %d.", n_cells + 1L, ncol(tab))))
  vals <- as.matrix(tab[, -1, drop = FALSE])
  if (!is.numeric(vals)) return(check_fail("The matrix contains non-numeric values in the first rows."))
  if (anyNA(vals) || any(vals < 0)) return(check_fail("Raw counts must be nonnegative numbers without missing values."))
  if (any(vals != floor(vals))) return(check_fail("Raw integer counts are required; the first rows contain non-integer (normalized or log-transformed) values. Upload the raw count matrix."))
  kind <- id_kind(tab[[1]]); pb <- id_problem(kind); if (!is.null(pb)) return(check_fail(pb))
  n_genes <- count_upload_lines(path) - 1
  if (n_genes < 2000) return(check_fail(sprintf("The matrix has %s genes. Upload the full gene matrix (whole transcriptome), not a marker panel.", fmt_count(n_genes))))
  lib <- check_msg(check_one_library(cols[-1])); if (!is.null(lib)) return(check_fail(lib))
  check_ok(sprintf("%s genes \u00D7 %s cells \u00B7 %s \u00B7 integer counts in the first %d rows", fmt_count(n_genes), fmt_count(n_cells), ID_LABEL[[kind]], nrow(tab)),
           info = list(genes = n_genes, cells = n_cells))
}

precheck_mtx_zip <- function(path, lim) {
  l <- tryCatch(utils::unzip(path, list = TRUE), error = function(e) NULL)
  if (is.null(l)) return(check_fail("The ZIP archive could not be read."))
  l <- l[!grepl("(^|/)(__MACOSX|\\._)", l$Name) & !grepl("/$", l$Name), , drop = FALSE]
  if (nrow(l) > 64) return(check_fail("The archive has too many files (limit 64)."))
  base <- basename(l$Name); pick <- function(p) l$Name[grepl(p, base, ignore.case = TRUE)]
  mtx <- pick("^matrix\\.mtx(\\.gz)?$"); fea <- pick("^(features|genes)\\.tsv(\\.gz)?$"); bc <- pick("^barcodes\\.tsv(\\.gz)?$")
  miss <- c(if (length(mtx) != 1) "matrix.mtx", if (length(fea) != 1) "features.tsv", if (length(bc) != 1) "barcodes.tsv")
  if (length(miss)) return(check_fail("The ZIP must contain exactly one matrix.mtx, features.tsv (or genes.tsv) and barcodes.tsv, each optionally gzipped; missing or repeated: ", paste(miss, collapse = ", "), "."))
  if (sum(l$Length) > lim$max_expanded_mb * 1024^2) return(check_fail(sprintf("The archive expands to %s MB; this server accepts up to %s MB. Gzip the files inside the ZIP or use the local edition.",
                                                                             fmt_count(sum(l$Length) / 1024^2), fmt_count(lim$max_expanded_mb))))
  head <- tryCatch(zip_member_lines(path, mtx, 60L), error = function(e) character())
  if (!length(head) || !grepl("^%%MatrixMarket matrix coordinate (integer|real|double)", head[1], ignore.case = TRUE))
    return(check_fail("matrix.mtx must be a MatrixMarket coordinate file with integer or real values."))
  k <- which(!startsWith(head, "%"))[1]
  dims <- suppressWarnings(as.numeric(strsplit(trimws(head[k]), "\\s+")[[1]]))
  if (length(dims) != 3L || anyNA(dims)) return(check_fail("matrix.mtx has an invalid size line."))
  if (dims[2] > lim$max_cells) return(check_fail(sprintf("The matrix has %s cells; this server accepts up to %s per upload. Subsample cells or use the local edition.", fmt_count(dims[2]), fmt_count(lim$max_cells))))
  if (dims[3] > lim$max_nnz) return(check_fail(sprintf("The matrix has %s non-zero counts; this server accepts up to %s. Subsample cells or use the local edition.", fmt_count(dims[3]), fmt_count(lim$max_nnz))))
  ent <- head[-seq_len(k)]
  if (length(ent)) {
    v <- suppressWarnings(as.numeric(vapply(strsplit(trimws(ent), "\\s+"), function(r) if (length(r) >= 3) r[3] else NA_character_, character(1))))
    if (anyNA(v) || any(v < 0) || any(v != floor(v))) return(check_fail("Raw integer counts are required; matrix.mtx contains non-integer (normalized or log-transformed) values."))
  }
  feats <- tryCatch(zip_member_lines(path, fea), error = function(e) character())
  bars <- tryCatch(zip_member_lines(path, bc), error = function(e) character())
  if (length(feats) == dims[2] && length(bars) == dims[1] && dims[1] != dims[2]) return(check_fail("matrix.mtx appears to be cells \u00D7 genes; 10x matrices are genes \u00D7 cells."))
  if (length(feats) != dims[1] || length(bars) != dims[2]) return(check_fail(sprintf("matrix.mtx declares %s genes \u00D7 %s cells, but features.tsv has %s rows and barcodes.tsv %s.",
                                                                                     fmt_count(dims[1]), fmt_count(dims[2]), fmt_count(length(feats)), fmt_count(length(bars)))))
  fs <- strsplit(feats, "\t", fixed = TRUE)
  genes <- vapply(fs, function(r) r[[min(2L, length(r))]], character(1))
  kind <- id_kind(genes); pb <- id_problem(kind); if (!is.null(pb)) return(check_fail(pb))
  if (kind == "symbol" && sum(COMMON_GENES %in% genes) < 8)
    return(check_fail("Few common human genes (for example ACTB, B2M, CD74, MS4A1) are among the gene symbols. Use human gene symbols, as in column 2 of the 10x features.tsv, or human Ensembl IDs."))
  if (dims[1] < 2000) return(check_fail(sprintf("The matrix has %s genes. Upload the full gene matrix (whole transcriptome), not a marker panel.", fmt_count(dims[1]))))
  cells <- vapply(strsplit(bars, "\t", fixed = TRUE), function(r) r[[1]], character(1))
  lib <- check_msg(check_one_library(cells)); if (!is.null(lib)) return(check_fail(lib))
  check_ok(sprintf("%s genes \u00D7 %s cells \u00B7 %s non-zero counts \u00B7 %s \u00B7 integer counts", fmt_count(dims[1]), fmt_count(dims[2]), fmt_count(dims[3]), ID_LABEL[[kind]]),
           info = list(genes = dims[1], cells = dims[2], nnz = dims[3]))
}

precheck_bulk <- function(path, name, lim) {
  nm <- tolower(name)
  if (!grepl("\\.(csv|tsv|txt)(\\.gz)?$", nm)) return(check_fail("Upload the expression matrix as CSV, TSV or TXT (optionally .gz): genes in rows, samples in columns."))
  sep <- if (grepl("\\.csv(\\.gz)?$", nm)) "," else "\t"
  lines <- tryCatch(head_lines(path, PRECHECK_ROWS + 1L), error = function(e) character())
  lines <- lines[nzchar(lines)]
  if (length(lines) < 2L) return(check_fail("The expression file has no data rows."))
  fields <- split_header(lines[1], sep)
  tab <- tryCatch(data.table::fread(text = paste0(paste(lines[-1], collapse = "\n"), "\n"), sep = sep, header = FALSE, data.table = FALSE, showProgress = FALSE,
                                    na.strings = BULK_NA_STRINGS, colClasses = list(character = 1L)), error = function(e) NULL)
  if (is.null(tab)) return(check_fail("The first rows could not be parsed as a ", if (sep == ",") "comma" else "tab", "-separated table."))
  if (length(fields) == ncol(tab) - 1L) fields <- c("gene", fields)
  k <- length(fields) - 1L
  if (k < 1L || ncol(tab) != k + 1L) return(check_fail(sprintf("Use gene symbols in the first column and one column per sample (the header has %d fields, the first row %d).", length(fields), ncol(tab))))
  if (k < 2L) return(check_fail("At least two samples are required; the CMI-high/low split uses the within-cohort median."))
  vals <- suppressWarnings(as.matrix(as.data.frame(lapply(tab[, -1, drop = FALSE], function(v) as.numeric(as.character(v))))))
  raw <- as.matrix(tab[, -1, drop = FALSE])
  if (any(is.na(vals) & !is.na(raw))) return(check_fail("The matrix contains non-numeric values in the first rows (for example: ", paste(utils::head(unique(raw[is.na(vals) & !is.na(raw)]), 3), collapse = ", "), ")."))
  if (any(vals < 0, na.rm = TRUE)) return(check_fail("Expression must be nonnegative (for example log2 intensities, TPM, FPKM or counts); negative values were found. Gene-centered, z-scored or log-ratio matrices change the order of genes within a sample, which the gene-pair features depend on."))
  kind <- id_kind(tab[[1]]); pb <- id_problem(kind, bulk = TRUE); if (!is.null(pb)) return(check_fail(pb))
  dup <- anyDuplicated(fields[-1]); if (dup) return(check_fail("Sample identifiers are repeated in the header (for example: ", fields[-1][dup], ")."))
  n_genes <- count_upload_lines(path) - 1
  if (n_genes * k > lim$max_bulk_entries) return(check_fail(sprintf("The matrix has %s values; this server accepts up to %s. Score the cohort in batches of samples (each sample's score depends only on its own genes) or use the local edition.",
                                                                    fmt_count(n_genes * k), fmt_count(lim$max_bulk_entries))))
  rng <- suppressWarnings(range(vals, na.rm = TRUE))
  warn <- if (n_genes < 5000) sprintf("Only %s genes: ranks are taken over all genes of a sample, so a reduced gene set changes the features. Upload the full measured gene set.", fmt_count(n_genes)) else character()
  check_ok(sprintf("%s genes \u00D7 %s samples \u00B7 gene symbols \u00B7 values %s to %s in the first %d rows", fmt_count(n_genes), fmt_count(k), signif(rng[1], 3), signif(rng[2], 3), nrow(tab)),
           warnings = warn, info = list(genes = n_genes, samples = k, sample_ids = fields[-1]))
}

precheck_spatial <- function(path, name, lim, mode = "weights") {
  nm <- tolower(name)
  if (identical(mode, "rctd")) {
    if (!grepl("\\.zip$", nm)) return(check_fail("An RCTD bundle is a ZIP file."))
    l <- tryCatch(utils::unzip(path, list = TRUE), error = function(e) NULL)
    if (is.null(l)) return(check_fail("The ZIP archive could not be read."))
    need <- c("reference_counts.csv", "reference_metadata.csv", "spatial_counts.csv", "spatial_metadata.csv", "manifest.csv")
    have <- basename(l$Name[!grepl("(^|/)(__MACOSX|\\._)", l$Name) & !grepl("/$", l$Name)])
    if (!setequal(have, need)) return(check_fail("The bundle must contain exactly: ", paste(need, collapse = ", "), "."))
    return(check_ok(sprintf("RCTD bundle \u00B7 %s MB expanded", fmt_count(sum(l$Length) / 1024^2))))
  }
  if (!grepl("\\.(csv|tsv|txt)(\\.gz)?$", nm)) return(check_fail("Upload the spot table as CSV, TSV or TXT."))
  tab <- tryCatch(data.table::fread(path, nrows = PRECHECK_ROWS, data.table = FALSE, check.names = FALSE, showProgress = FALSE), error = function(e) NULL)
  if (is.null(tab) || !nrow(tab)) return(check_fail("The spot table could not be parsed."))
  sc <- tryCatch(spatial_columns(tab), input_error = function(e) e)
  if (inherits(sc, "condition")) return(check_fail(conditionMessage(sc)))
  n_spots <- count_upload_lines(path) - 1
  if (n_spots > lim$max_spots) return(check_fail(sprintf("The table has %s spots; this server accepts up to %s per section.", fmt_count(n_spots), fmt_count(lim$max_spots))))
  warn <- if (length(sc$ignored)) sprintf("Not read as weights: %s.", paste(utils::head(sc$ignored, 6), collapse = ", "))
  check_ok(sprintf("%s spots \u00B7 %d classes, including Mal-1, Mal-2 and Mal-3", fmt_count(n_spots), length(sc$classes)),
           warnings = warn, info = list(spots = n_spots, classes = sc$classes))
}

precheck_upload <- function(module, path, name, lim, mode = "weights") {
  tryCatch(switch(module,
    single_cell = precheck_single_cell(path, name, lim),
    bulk = precheck_bulk(path, name, lim),
    spatial = precheck_spatial(path, name, lim, mode)),
    error = function(e) check_fail("The file could not be read: ", conditionMessage(e)))
}
