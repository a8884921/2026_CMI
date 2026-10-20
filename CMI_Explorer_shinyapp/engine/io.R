# Input readers and validation for CMI Explorer.
# Every reader returns sparse matrices and never materialises a dense genes x cells table,
# so input size, not the reader, bounds memory.

input_error <- function(...) {
  stop(structure(class = c("input_error", "error", "condition"),
                 list(message = paste0(...), call = NULL)))
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

limits_for <- function(profile = Sys.getenv("CMI_PROFILE", "hosted")) {
  num <- function(name, default) { v <- suppressWarnings(as.numeric(Sys.getenv(name, ""))); if (is.finite(v)) v else default }
  if (identical(profile, "local")) {
    list(profile = "local", max_upload_mb = num("CMI_MAX_UPLOAD_MB", 4096), max_expanded_mb = num("CMI_MAX_EXPANDED_MB", 16384), max_cells = num("CMI_MAX_CELLS", Inf),
         max_nnz = num("CMI_MAX_NNZ", Inf), max_bulk_entries = num("CMI_MAX_BULK_ENTRIES", Inf),
         max_spots = num("CMI_MAX_SPOTS", Inf), max_ref_cells = num("CMI_MAX_REF_CELLS", Inf))
  } else {
    # Sized for a 1 GB shinyapps.io instance with the interface process resident (docs/DEPLOYMENT_zh.md, docs/VALIDATION.md):
    # 15 million non-zero counts peak at about 0.67 GB in the worker, leaving room for the interface.
    list(profile = "hosted", max_upload_mb = num("CMI_MAX_UPLOAD_MB", 300), max_expanded_mb = num("CMI_MAX_EXPANDED_MB", 600), max_cells = num("CMI_MAX_CELLS", 10000),
         max_nnz = num("CMI_MAX_NNZ", 1.5e7), max_bulk_entries = num("CMI_MAX_BULK_ENTRIES", 2.5e7),
         max_spots = num("CMI_MAX_SPOTS", 20000), max_ref_cells = num("CMI_MAX_REF_CELLS", 3000))
  }
}

is_ascii <- function(x) {
  x <- as.character(x)
  !is.na(x) & !is.na(iconv(x, from = "UTF-8", to = "ASCII", sub = NA))
}

checked_ids <- function(values, label, allow_duplicates = FALSE) {
  values <- trimws(as.character(values))
  if (!length(values)) input_error(label, ": no identifiers found.")
  if (anyNA(values) || any(!nzchar(values))) input_error(label, " contain empty identifiers.")
  if (!all(is_ascii(values))) input_error(label, " must use ASCII characters only.")
  if (!allow_duplicates && anyDuplicated(values)) {
    d <- unique(values[duplicated(values)])
    input_error(label, " contain duplicate identifiers (for example: ", paste(utils::head(d, 3), collapse = ", "), ").")
  }
  values
}

open_text <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
}

detect_sep <- function(header_line, path) {
  ext <- tolower(sub("\\.gz$", "", path, ignore.case = TRUE))
  if (grepl("\\.csv$", ext)) return(",")
  if (grepl("\\.(tsv|txt)$", ext)) return("\t")
  if (grepl("\t", header_line, fixed = TRUE)) "\t" else ","
}

split_header <- function(line, sep) {
  fields <- strsplit(line, sep, fixed = TRUE)[[1]]
  fields <- sub('^"(.*)"$', "\\1", trimws(fields))
  sub("^\uFEFF", "", fields)
}

# Text matrix (genes in rows, observations in columns; CSV/TSV/TXT, optionally gzipped), read in
# two streaming passes: pass 1 validates and counts non-zeros per column, pass 2 writes them into
# preallocated CSC arrays. Peak memory is the final sparse matrix plus one block of rows.
text_block_reader <- function(path, counts, max_cols, block_entries) {
  con <- open_text(path)
  header <- readLines(con, n = 1L, warn = FALSE)
  if (!length(header) || !nzchar(header)) { close(con); input_error("The file is empty.") }
  sep <- detect_sep(header, path)
  cols <- split_header(header, sep)
  if (length(cols) < 2L) { close(con); input_error("Use one identifier column followed by one column per observation.") }
  ncol_obs <- length(cols) - 1L
  if (ncol_obs > max_cols) { close(con); input_error(sprintf("This server accepts up to %s observations per upload; the file has %s. Subsample or use the local edition.",
                                                     format(max_cols, big.mark = ","), format(ncol_obs, big.mark = ","))) }
  block_rows <- max(20L, as.integer(floor(block_entries / ncol_obs)))
  state <- new.env(); state$offset <- 0L
  next_block <- function() {
    repeat {
      lines <- readLines(con, n = block_rows, warn = FALSE)
      if (!length(lines)) return(NULL)
      lines <- lines[nzchar(trimws(lines))]
      if (length(lines)) break
    }
    tab <- tryCatch(data.table::fread(text = lines, sep = sep, header = FALSE, colClasses = list(character = 1L),
                                      na.strings = c("", "NA", "NaN", "nan"), showProgress = FALSE, quote = "\""),
                    error = function(e) input_error("Could not parse the matrix near data row ", state$offset + 1L, ": ", conditionMessage(e)))
    rm(lines)
    if (ncol(tab) != ncol_obs + 1L) input_error("Data row ", state$offset + 1L, " onwards has ", ncol(tab), " fields; the header has ", ncol_obs + 1L, ".")
    genes <- tab[[1]]
    vals <- as.matrix(tab[, -1, with = FALSE]); rm(tab)
    if (!is.numeric(vals)) {
      suppressWarnings(storage.mode(vals) <- "double")
      if (anyNA(vals)) input_error("The matrix contains non-numeric values near data row ", state$offset + 1L, ".")
    }
    if (counts) {
      if (anyNA(vals) || any(!is.finite(vals))) input_error("Raw counts must be finite; missing values were found near data row ", state$offset + 1L, ".")
      if (any(vals < 0)) input_error("Raw counts must be nonnegative.")
      if (any(vals != floor(vals))) input_error("Raw integer counts are required; normalized or log-transformed values were found. Upload the raw count matrix (for example, the 10x filtered_feature_bc_matrix).")
    } else if (any(vals < 0, na.rm = TRUE) || any(is.infinite(vals))) {
      input_error("Expression values must be nonnegative and finite.")
    }
    out <- list(genes = genes, vals = vals, offset = state$offset)
    state$offset <- state$offset + length(genes)
    out
  }
  list(header = cols, next_block = next_block, close = function() close(con))
}

read_text_matrix_sparse <- function(path, counts = TRUE, max_cols = Inf, max_nnz = Inf, block_entries = 4e6) {
  # pass 1: per-column non-zero counts and row identifiers
  rd <- text_block_reader(path, counts, max_cols, block_entries)
  col_nnz <- NULL; gene_chunks <- list(); n_rows <- 0L
  repeat {
    b <- rd$next_block(); if (is.null(b)) break
    if (is.null(col_nnz)) col_nnz <- numeric(ncol(b$vals))
    col_nnz <- col_nnz + colSums(b$vals != 0 | is.na(b$vals))
    gene_chunks[[length(gene_chunks) + 1L]] <- b$genes
    n_rows <- n_rows + length(b$genes)
    rm(b); invisible(gc(verbose = FALSE))
    if (sum(col_nnz) > max_nnz) { rd$close(); input_error(sprintf("The matrix has more than %s non-zero entries, above this server's limit. Subsample cells or use the local edition.",
                                                          format(max_nnz, big.mark = ",", scientific = FALSE))) }
  }
  rd$close()
  if (!n_rows) input_error("The matrix has no data rows.")
  ids <- checked_ids(rd$header[-1], "Observation (column) identifiers")
  rows <- checked_ids(unlist(gene_chunks, use.names = FALSE), "Gene (row) identifiers"); rm(gene_chunks)
  nc <- length(ids); nnz <- sum(col_nnz)
  p <- c(0, cumsum(col_nnz))
  ii <- integer(nnz); xx <- numeric(nnz); fill <- numeric(nc)
  # pass 2: rows arrive in increasing order, so each column's row indices stay sorted
  rd <- text_block_reader(path, counts, max_cols, block_entries); on.exit(rd$close(), add = TRUE)
  repeat {
    b <- rd$next_block(); if (is.null(b)) break
    nz <- which(b$vals != 0 | is.na(b$vals))
    if (!length(nz)) next
    nrb <- nrow(b$vals)
    cc <- (nz - 1L) %/% nrb + 1L                      # column-major order: sorted by column, then row
    k <- tabulate(cc, nbins = nc)
    first <- c(1, cumsum(k)[-nc] + 1)
    pos <- p[cc] + fill[cc] + (seq_along(nz) - first[cc]) + 1
    ii[pos] <- as.integer((nz - 1L) %% nrb + b$offset)
    xx[pos] <- b$vals[nz]
    fill <- fill + k
    rm(b, nz, cc, pos); invisible(gc(verbose = FALSE))
  }
  m <- methods::new("dgCMatrix", i = ii, p = as.integer(p), x = xx, Dim = as.integer(c(n_rows, nc)), Dimnames = list(rows, ids))
  list(matrix = m, rows = rows, cols = ids)
}

safe_extract <- function(path, target, max_files = 64L, max_expanded = 4 * 1024^3) {
  listing <- tryCatch(utils::unzip(path, list = TRUE), error = function(e) input_error("The ZIP archive could not be read."))
  if (nrow(listing) > max_files) input_error("The archive has too many files (limit ", max_files, ").")
  if (sum(listing$Length) > max_expanded) input_error("The archive expands beyond the size limit.")
  names_ <- gsub("\\\\", "/", listing$Name)
  if (any(!is_ascii(names_))) input_error("The archive contains a non-ASCII file name.")
  parts <- strsplit(names_, "/", fixed = TRUE)
  if (any(startsWith(names_, "/")) || any(vapply(parts, function(p) ".." %in% p, logical(1))) || any(grepl(":", names_, fixed = TRUE))) {
    input_error("The archive contains an unsafe path.")
  }
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  utils::unzip(path, exdir = target, junkpaths = FALSE)
  files <- list.files(target, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  for (f in files) if (nzchar(Sys.readlink(f))) input_error("Symbolic links are not allowed in archives.")
  files[!grepl("(^|/)(__MACOSX|\\._)", files)]
}

read_lines_any <- function(path) { con <- open_text(path); on.exit(close(con)); readLines(con, warn = FALSE) }

# 10x-style Matrix Market triplet (genes x cells), streamed with scan() into preallocated arrays.
# Cell Ranger writes entries sorted by column and then row, which lets the CSC matrix be filled in
# place; other orderings fall back to Matrix::sparseMatrix().
read_mtx_sparse <- function(mtx, features, barcodes, max_cols = Inf, max_nnz = Inf, chunk = 1e6) {
  con <- open_text(mtx); on.exit(close(con), add = TRUE)
  header <- readLines(con, n = 1L, warn = FALSE)
  if (!length(header) || !grepl("^%%MatrixMarket matrix coordinate (integer|real|double)", header, ignore.case = TRUE)) {
    input_error("matrix.mtx must be a MatrixMarket coordinate file with integer or real values.")
  }
  repeat { l <- readLines(con, n = 1L, warn = FALSE); if (!length(l) || !startsWith(l, "%")) break }
  dims <- suppressWarnings(as.numeric(strsplit(trimws(l), "\\s+")[[1]]))
  if (length(dims) != 3L || anyNA(dims)) input_error("matrix.mtx has an invalid size line.")
  if (dims[3] > max_nnz) input_error(sprintf("The matrix has %s non-zero entries, above this server's limit of %s. Subsample cells or use the local edition.",
                                             format(dims[3], big.mark = ","), format(max_nnz, big.mark = ",", scientific = FALSE)))
  gene_rows <- strsplit(read_lines_any(features), "\t", fixed = TRUE)
  genes <- vapply(gene_rows, function(r) r[[min(2L, length(r))]], character(1))
  cells <- vapply(strsplit(read_lines_any(barcodes), "\t", fixed = TRUE), function(r) r[[1]], character(1))
  nr <- dims[1]; nc <- dims[2]; nnz <- dims[3]
  if (nr == length(cells) && nc == length(genes) && nr != nc) input_error("matrix.mtx appears to be cells x genes; 10x matrices are genes x cells.")
  if (nr != length(genes) || nc != length(cells)) input_error("matrix.mtx dimensions do not match the feature and barcode files.")
  if (nc > max_cols) input_error(sprintf("This server accepts up to %s cells per upload; the file has %s. Subsample or use the local edition.",
                                         format(max_cols, big.mark = ","), format(nc, big.mark = ",")))
  # Gene symbols (column 2) with make.unique(), as Seurat::Read10X(gene.column = 2) and
  # scanpy.read_10x_mtx(var_names = "gene_symbols") do; the CMI was developed on gene symbols.
  genes <- checked_ids(make.unique(genes), "Gene identifiers"); cells <- checked_ids(cells, "Cell barcodes")
  ii <- integer(nnz); xx <- numeric(nnz); jj <- NULL
  col_n <- numeric(nc); done <- 0; last_i <- 0L; last_j <- 0L; sorted <- TRUE
  while (done < nnz) {
    tr <- scan(con, what = list(integer(), integer(), numeric()), nmax = min(chunk, nnz - done), quiet = TRUE)
    n <- length(tr[[1]]); if (!n) break
    i <- tr[[1]]; j <- tr[[2]]; x <- tr[[3]]; rm(tr)
    if (anyNA(i) || anyNA(j) || any(i < 1L | i > nr) || any(j < 1L | j > nc)) input_error("matrix.mtx has entries outside its declared dimensions.")
    if (any(!is.finite(x)) || any(x < 0)) input_error("Raw counts must be finite and nonnegative.")
    if (any(x != floor(x))) input_error("Raw integer counts are required; matrix.mtx contains non-integer values.")
    if (sorted) {
      dj <- diff(c(last_j, j)); di <- diff(c(last_i, i))
      if (any(dj < 0L) || any(di[dj == 0L] <= 0L)) {
        sorted <- FALSE; jj <- integer(nnz)
        if (done > 0) jj[seq_len(done)] <- rep.int(seq_len(nc), col_n)
      }
      last_i <- i[n]; last_j <- j[n]
    }
    idx <- done + seq_len(n)
    ii[idx] <- i - 1L; xx[idx] <- x
    if (!sorted) jj[idx] <- j
    col_n <- col_n + tabulate(j, nbins = nc)
    done <- done + n
    rm(i, j, x); invisible(gc(verbose = FALSE))
  }
  if (done != nnz) input_error("matrix.mtx has fewer entries than its header states.")
  if (sorted) {
    m <- methods::new("dgCMatrix", i = ii, p = as.integer(c(0, cumsum(col_n))), x = xx, Dim = as.integer(c(nr, nc)), Dimnames = list(genes, cells))
  } else {
    m <- Matrix::sparseMatrix(i = ii + 1L, j = jj, x = xx, dims = c(nr, nc), dimnames = list(genes, cells), repr = "C")
  }
  list(matrix = m, rows = genes, cols = cells)
}

# Human genes only: mouse Ensembl IDs or mouse-style symbols (Cd19, Ms4a1) are rejected early.
check_species <- function(genes) {
  if (mean(grepl("^ENSMUSG", genes)) > 0.5) input_error("The gene identifiers are mouse Ensembl IDs (ENSMUSG...). The CMI models use human genes.")
  sym <- genes[!grepl("^ENS", genes)]
  if (length(sym) >= 1000 && mean(grepl("^[A-Z][a-z0-9]+$", sym)) > 0.5)
    input_error("Most gene symbols are in mouse style (for example Cd19, Ms4a1). The CMI models use human gene symbols or human Ensembl IDs.")
  invisible(TRUE)
}
# One sample per upload: 10x barcodes from several libraries carry different suffixes (-1, -2, ...).
check_one_library <- function(cells) {
  has <- grepl("-[0-9]+$", cells)
  if (mean(has) < 0.9) return(invisible(TRUE))
  suf <- unique(sub("^.*-([0-9]+)$", "\\1", cells[has]))
  if (length(suf) > 1) input_error("This matrix appears to combine ", length(suf), " libraries (barcode suffixes -", paste(utils::head(sort(as.integer(suf)), 3), collapse = ", -"),
                                   if (length(suf) > 3) ", ..." else "", "). Upload one sample at a time: proliferation scoring uses the uploaded sample as its background.")
  invisible(TRUE)
}

# Single-cell counts from CSV/TSV/TXT(.gz) or a 10x Matrix Market ZIP.
read_single_cell_counts <- function(path, workspace, lim = limits_for()) {
  name <- tolower(basename(path))
  if (grepl("\\.zip$", name)) {
    files <- safe_extract(path, file.path(workspace, "mtx"), max_expanded = lim$max_expanded_mb * 1024^2)
    pick <- function(pattern, label) {
      found <- files[grepl(pattern, basename(files), ignore.case = TRUE)]
      if (length(found) != 1L) input_error("The ZIP must contain exactly one ", label, " (found ", length(found), ").")
      found
    }
    res <- read_mtx_sparse(pick("^matrix\\.mtx(\\.gz)?$", "matrix.mtx"),
                           pick("^(features|genes)\\.tsv(\\.gz)?$", "features.tsv or genes.tsv"),
                           pick("^barcodes\\.tsv(\\.gz)?$", "barcodes.tsv"), max_cols = lim$max_cells, max_nnz = lim$max_nnz)
  } else if (grepl("\\.(csv|tsv|txt)(\\.gz)?$", name)) {
    res <- read_text_matrix_sparse(path, counts = TRUE, max_cols = lim$max_cells, max_nnz = lim$max_nnz)
  } else {
    input_error("Use a CSV/TSV/TXT count matrix (optionally .gz) or a ZIP with matrix.mtx, features.tsv and barcodes.tsv.")
  }
  m <- res$matrix; rm(res)
  check_species(rownames(m)); check_one_library(colnames(m))
  if (length(m@x) && min(m@x) == 0) m <- Matrix::drop0(m)
  invisible(gc(verbose = FALSE))
  if (nrow(m) < 2000L) input_error("The matrix has ", nrow(m), " genes. Upload the full gene matrix (whole transcriptome), not a marker panel.")
  empty <- Matrix::colSums(m) <= 0
  if (any(empty)) input_error(sum(empty), " cell(s) have no counts. Remove empty barcodes before uploading.")
  m
}
