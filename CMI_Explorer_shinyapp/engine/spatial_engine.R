# =============================================================================
# Spatial CMI grading
# =============================================================================
# Input: one tissue section as a spot table with coordinates and deconvolution weights for the
# three malignant grades (Mal-1, Mal-2, Mal-3) and any non-malignant classes. Every spot takes the
# class of highest weight: Mal-1, Mal-2, Mal-3, or Non-malignant when a non-malignant class wins.
# =============================================================================

SPATIAL_STATES <- c("Mal-1", "Mal-2", "Mal-3")

canonical_class <- function(x) {
  k <- gsub("^RCTD_weight_", "", x, ignore.case = TRUE)
  k2 <- tolower(gsub("[^A-Za-z0-9]", "", k))
  out <- k
  m <- regmatches(k2, regexec("^(?:mal|malig|malignant|tumou?r)(?:state)?([123])$", k2, perl = TRUE))
  hit <- vapply(m, length, integer(1)) == 2L
  out[hit] <- paste0("Mal-", vapply(m[hit], `[`, character(1), 2L))
  out
}

# Which columns of a spot table hold the identifiers, the coordinates and the class weights. Shared by
# the reader below and by the quick upload check in the interface (which passes only the first rows).
spatial_columns <- function(tab) {
  key <- tolower(gsub("[^A-Za-z0-9]", "", names(tab)))
  idc <- which(key %in% c("spotid", "spot", "barcode", "cellid", "id"))[1]
  xc <- which(key %in% c("x", "imagecol", "pxlcolinfullres", "arraycol", "col", "xcoord"))[1]
  yc <- which(key %in% c("y", "imagerow", "pxlrowinfullres", "arrayrow", "row", "ycoord"))[1]
  if (is.na(idc) || is.na(xc) || is.na(yc)) input_error("The spot table needs spot_id, x and y columns.")
  ignore <- c(idc, xc, yc, which(key %in% c("sampleid", "sample", "disease", "malignantload", "malignantstate", "rctddominanttype",
                                             "rctdtop1weight", "rctdsecondtype", "rctdsecondweight", "rctdnormalizedentropy", "rctdmode")))
  # Count, QC and cluster columns that often travel with spot tables are not class weights.
  qc_like <- grepl("^(numi|ncount|nfeature|totalcounts?|ncounts?|ngenes|ngenesbycounts|pctcounts|percent|intissue|tissue|cluster|seuratclusters|leiden|louvain|ident|depth|librarysize|spotclass|firsttype|secondtype|imagerow|imagecol|arrayrow|arraycol|pxlrowinfullres|pxlcolinfullres)",
                   key)
  wcols <- setdiff(seq_along(tab), c(ignore, which(qc_like)))
  wcols <- wcols[vapply(tab[wcols], is.numeric, logical(1))]
  big <- wcols[vapply(tab[wcols], function(v) isTRUE(max(v, na.rm = TRUE) > 5), logical(1))]
  if (length(big)) input_error("These columns do not look like deconvolution weights (values above 5): ", paste(utils::head(names(tab)[big], 5), collapse = ", "),
                               ". Keep spot_id, x, y and one weight column per class, and remove count or QC columns.")
  classes <- canonical_class(names(tab)[wcols])
  if (!all(SPATIAL_STATES %in% classes)) input_error("Weight columns for Mal-1, Mal-2 and Mal-3 are required (found: ", paste(classes, collapse = ", "), "). Column names such as Mal-1, Malig_1 or RCTD_weight_Malig_1 are recognized.")
  if (anyDuplicated(classes)) input_error("Two weight columns map to the same class: ", paste(unique(classes[duplicated(classes)]), collapse = ", "))
  list(key = key, idc = idc, xc = xc, yc = yc, wcols = wcols, classes = classes, ignored = names(tab)[setdiff(seq_along(tab), c(idc, xc, yc, wcols))],
       cd8 = classes[grepl("cd8", classes, ignore.case = TRUE)][1])
}

read_spatial_weights <- function(path, lim = limits_for()) {
  tab <- tryCatch(data.table::fread(path, data.table = FALSE, check.names = FALSE, showProgress = FALSE),
                  error = function(e) input_error("The spot table could not be parsed: ", conditionMessage(e)))
  sc <- spatial_columns(tab); key <- sc$key; idc <- sc$idc; xc <- sc$xc; yc <- sc$yc; wcols <- sc$wcols; classes <- sc$classes
  spots <- checked_ids(tab[[idc]], "Spot identifiers")
  if (length(spots) > lim$max_spots) input_error(sprintf("This server accepts up to %s spots per section; the table has %s.", format(lim$max_spots, big.mark = ","), format(length(spots), big.mark = ",")))
  xy <- cbind(x = as.numeric(tab[[xc]]), y = as.numeric(tab[[yc]]))
  if (any(!is.finite(xy))) input_error("Coordinates must be finite numbers.")
  W <- as.matrix(tab[, wcols, drop = FALSE]); storage.mode(W) <- "double"; colnames(W) <- classes
  if (any(!is.finite(W)) || any(W < 0)) input_error("Weights must be finite and nonnegative.")
  rs <- rowSums(W)
  if (any(rs <= 0)) input_error(sum(rs <= 0), " spot(s) have all-zero weights.")
  sample_id <- if ("sampleid" %in% key) as.character(tab[[which(key == "sampleid")[1]]]) else rep("section", length(spots))
  if (length(unique(sample_id)) > 1L) input_error("Upload one tissue section per run (the table has ", length(unique(sample_id)), " sample_id values).")
  list(spots = spots, xy = xy, W = W / rs, sample_id = sample_id[1])
}

SPOT_GRADES <- c(SPATIAL_STATES, "Non-malignant")
grade_spots <- function(sp) {
  W <- sp$W; classes <- colnames(W); n <- nrow(W)
  dominant <- classes[max.col(W, ties.method = "first")]
  mal <- intersect(SPATIAL_STATES, classes)
  grade <- ifelse(dominant %in% mal, dominant, "Non-malignant")
  spots <- data.frame(spot_id = sp$spots, sample_id = sp$sample_id, x = sp$xy[, 1], y = sp$xy[, 2], grade = grade,
                      dominant_class = dominant, malignant_load = rowSums(W[, mal, drop = FALSE]), stringsAsFactors = FALSE)
  for (cl in classes) spots[[paste0("weight_", gsub("[^A-Za-z0-9]+", "_", cl))]] <- W[, cl]
  n_g <- as.integer(table(factor(grade, levels = SPOT_GRADES)))
  mal_w <- colSums(W[, mal, drop = FALSE])
  grades <- data.frame(grade = SPOT_GRADES, spots = n_g, spot_share = n_g / n,
                       malignant_weight_share = c(unname(mal_w[SPATIAL_STATES] / sum(mal_w)), NA_real_), stringsAsFactors = FALSE)
  pitch <- if (n > 1) stats::median(RANN::nn2(sp$xy, sp$xy, k = 2L)$nn.dists[, 2]) else NA_real_   # center-to-center spot spacing
  list(spots = spots, grades = grades, classes = classes, pitch = pitch)
}

# ---- Optional RCTD deconvolution (local edition; memory patch for spacexr 2.2.1) ----------------
# spacexr is optional and is looked up by name so that deployment tools do not try to install this
# GitHub-only package on the hosted server.
RCTD_PKG <- "spacexr"
sx <- function(f) getExportedValue(RCTD_PKG, f)

rctd_q_cache <- function(cache_dir) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  index_file <- file.path(cache_dir, "index.rds")
  if (!file.exists(index_file)) {
    idx <- character()
    for (k in 1:5) {
      q <- readRDS(system.file("extdata", sprintf("Qmat/Q_mat_%d.rds", k), package = RCTD_PKG))
      for (nm in names(q)) saveRDS(q[[nm]], file.path(cache_dir, paste0("Q_", nm, ".rds")), compress = FALSE)
      idx <- c(idx, names(q)); rm(q); invisible(gc(verbose = FALSE))
    }
    saveRDS(idx, index_file)
  }
  idx <- readRDS(index_file)
  env <- new.env(parent = emptyenv())
  for (nm in idx) local({ f <- file.path(cache_dir, paste0("Q_", nm, ".rds")); makeActiveBinding(nm, function() readRDS(f), env) })
  env
}
# spacexr 2.2.1 choose_sigma_c(), unchanged except that the 126 likelihood tables (431 MB) are read
# one at a time from a disk cache instead of all at once. Verified to give identical RCTD weights.
choose_sigma_c_lazy <- function(RCTD) {
  puck = RCTD@spatialRNA; MIN_UMI = RCTD@config$UMI_min_sigma; sigma = 100
  Q_mat_all <- rctd_q_cache(.cmi_q_cache_dir)
  X_vals <- readRDS(system.file("extdata", "Qmat/X_vals.rds", package = RCTD_PKG))
  N_fit = min(RCTD@config$N_fit, sum(puck@nUMI > MIN_UMI))
  if (N_fit == 0) stop(paste('choose_sigma_c determined a N_fit of 0! This is probably due to unusually low UMI counts per bead in your dataset. Try decreasing the parameter UMI_min_sigma. It currently is', MIN_UMI, 'but none of the beads had counts larger than that.'))
  fit_ind = sample(names(puck@nUMI[puck@nUMI > MIN_UMI]), N_fit)
  beads = t(as.matrix(puck@counts[RCTD@internal_vars$gene_list_reg, fit_ind]))
  message(paste('chooseSigma: using initial Q_mat with sigma = ', sigma/100))
  for (iter in 1:RCTD@config$N_epoch) {
    set_likelihood_vars(Q_mat_all[[as.character(sigma)]], X_vals)
    results = decompose_batch(puck@nUMI[fit_ind], RCTD@cell_type_info$renorm[[1]], beads, RCTD@internal_vars$gene_list_reg, constrain = F, max_cores = RCTD@config$max_cores)
    weights = Matrix(0, nrow = N_fit, ncol = RCTD@cell_type_info$renorm[[3]])
    rownames(weights) = fit_ind; colnames(weights) = RCTD@cell_type_info$renorm[[2]];
    for (i in 1:N_fit) weights[i, ] = results[[i]]$weights
    prediction <- sweep(as.matrix(RCTD@cell_type_info$renorm[[1]][RCTD@internal_vars$gene_list_reg, ]) %*% t(as.matrix(weights)), 2, puck@nUMI[fit_ind], '*')
    message(paste('Likelihood value:', calc_log_l_vec(as.vector(prediction), as.vector(t(beads)))))
    sigma_prev <- sigma
    sigma <- chooseSigma(prediction, t(beads), Q_mat_all, X_vals, sigma)
    message(paste('Sigma value: ', sigma/100))
    if (sigma == sigma_prev) break
  }
  RCTD@internal_vars$sigma <- sigma/100
  RCTD@internal_vars$Q_mat <- Q_mat_all[[as.character(sigma)]]
  RCTD@internal_vars$X_vals <- X_vals
  return(RCTD)
}
patch_spacexr <- function(cache_dir) {
  if (!identical(as.character(utils::packageVersion(RCTD_PKG)), "2.2.1")) return(invisible(FALSE))
  ns <- asNamespace(RCTD_PKG); f <- choose_sigma_c_lazy
  environment(f) <- list2env(list(rctd_q_cache = rctd_q_cache, .cmi_q_cache_dir = cache_dir, Matrix = Matrix::Matrix, RCTD_PKG = RCTD_PKG), parent = ns)
  utils::assignInNamespace("choose_sigma_c", f, ns = RCTD_PKG)
  invisible(TRUE)
}

# Five-file bundle (v1.2.0 contract): reference_counts.csv, reference_metadata.csv (cell, celltype_anno,
# CMI, disease), spatial_counts.csv, spatial_metadata.csv (spot_id, sample_id, disease, x, y), manifest.csv.
run_rctd_bundle <- function(zip_path, workspace, threads = 1L, seed = 20260101L, lim = limits_for(), progress = function(stage, pct, detail) NULL) {
  if (!requireNamespace(RCTD_PKG, quietly = TRUE)) input_error("RCTD needs the spacexr package (github.com/dmcable/spacexr), which this installation does not have. Run RCTD locally and upload the weights table instead.")
  files <- safe_extract(zip_path, file.path(workspace, "bundle"))
  need <- c("reference_counts.csv", "reference_metadata.csv", "spatial_counts.csv", "spatial_metadata.csv", "manifest.csv")
  if (!setequal(basename(files), need)) input_error("The bundle must contain exactly: ", paste(need, collapse = ", "), ".")
  f <- function(n) files[basename(files) == n]
  progress("read", 10, "Reading reference and spatial counts")
  ref <- read_text_matrix_sparse(f("reference_counts.csv"), counts = TRUE, max_cols = lim$max_ref_cells)$matrix
  spc <- read_text_matrix_sparse(f("spatial_counts.csv"), counts = TRUE, max_cols = lim$max_spots)$matrix
  refmeta <- utils::read.csv(f("reference_metadata.csv"), colClasses = "character", check.names = FALSE)
  sm <- utils::read.csv(f("spatial_metadata.csv"), colClasses = "character", check.names = FALSE)
  if (!all(c("cell", "celltype_anno", "CMI") %in% names(refmeta))) input_error("reference_metadata.csv needs cell, celltype_anno and CMI.")
  if (!all(c("spot_id", "x", "y") %in% names(sm))) input_error("spatial_metadata.csv needs spot_id, x and y.")
  if (!setequal(refmeta$cell, colnames(ref)) || !setequal(sm$spot_id, colnames(spc))) input_error("Counts and metadata identifiers must match exactly.")
  labels <- checked_ids(refmeta$celltype_anno, "Reference labels", allow_duplicates = TRUE)
  malignant <- grepl("^(?:Malignant|Tumou?r)(?:$|[ _:-])", labels, ignore.case = TRUE, perl = TRUE)
  if (!any(malignant)) input_error("The reference has no cells labeled Malignant or Tumor.")
  cmi <- suppressWarnings(as.numeric(refmeta$CMI))
  if (any(!is.finite(cmi[malignant])) || any(cmi[malignant] < 0 | cmi[malignant] > 1)) input_error("Malignant reference cells need CMI values in [0, 1] (from the Single-cell workspace).")
  # Cells labeled malignant whose CMI is below the
  # malignancy threshold are not malignant B cells and are left out of the reference (v1.2.0 counted them as
  # Malig-1). Same rule as www/downloads/rctd_mal_states.R.
  below <- malignant & cmi < CMI_THRESHOLD
  if (!any(malignant & !below)) input_error(sprintf("No reference cell labeled malignant reaches the malignancy threshold (CMI %.4f).", CMI_THRESHOLD))
  # v1.2.0 class names inside RCTD (identical factor order and therefore identical fits); renamed to Mal-1/2/3 afterwards
  labels[malignant] <- ifelse(cmi[malignant] <= CMI_INDEX2, "Malig-1", ifelse(cmi[malignant] <= CMI_INDEX3, "Malig-2", "Malig-3"))
  left_out <- data.frame(cell = refmeta$cell[below], CMI = cmi[below], stringsAsFactors = FALSE)
  if (any(below)) { ref <- ref[, !(colnames(ref) %in% left_out$cell), drop = FALSE]; refmeta <- refmeta[!below, , drop = FALSE]; labels <- labels[!below] }
  types <- factor(labels); names(types) <- refmeta$cell
  if (any(table(types) < 25L)) input_error("Every reference class needs at least 25 cells: ", paste(names(which(table(types) < 25L)), collapse = ", "))
  sm <- sm[match(colnames(spc), sm$spot_id), , drop = FALSE]
  coords <- data.frame(x = as.numeric(sm$x), y = as.numeric(sm$y), row.names = sm$spot_id)
  if (any(!is.finite(as.matrix(coords)))) input_error("Spot coordinates must be finite.")
  spc <- spc[, Matrix::colSums(spc) >= 100, drop = FALSE]
  common <- intersect(rownames(ref), rownames(spc))
  if (length(common) < 1000L) input_error("Reference and spatial counts share fewer than 1,000 genes.")
  progress("rctd", 25, "RCTD full mode (this can take several minutes)")
  suppressPackageStartupMessages(requireNamespace(RCTD_PKG))
  patch_spacexr(file.path(Sys.getenv("CMI_CACHE_DIR", file.path(dirname(tempdir()), "cmi_cache")), "spacexr_qmat"))
  set.seed(seed)
  reference <- sx("Reference")(counts = ref[common, ], cell_types = types[colnames(ref)], n_max_cells = 10000, min_UMI = 100)
  puck <- sx("SpatialRNA")(coords = coords[colnames(spc), , drop = FALSE], counts = spc[common, ])
  obj <- sx("create.RCTD")(spatialRNA = puck, reference = reference, max_cores = threads, UMI_min = 100, counts_MIN = 10,
                           UMI_min_sigma = 300, CELL_MIN_INSTANCE = 25, MAX_MULTI_TYPES = 4, test_mode = FALSE)
  obj <- suppressMessages(sx("run.RCTD")(obj, doublet_mode = "full"))
  W <- as.matrix(obj@results$weights)
  if (!nrow(W) || any(!is.finite(W)) || any(W < 0) || any(rowSums(W) <= 0)) stop("RCTD returned invalid weights.")
  colnames(W) <- canonical_class(colnames(W))
  rc <- as.data.frame(table(types), stringsAsFactors = FALSE); names(rc) <- c("reference_class", "cells"); rc$reference_class <- canonical_class(rc$reference_class)
  list(spots = rownames(W), xy = as.matrix(coords[rownames(W), c("x", "y")]), W = W / rowSums(W), raw_weights = W,
       sample_id = if ("sample_id" %in% names(sm)) sm$sample_id[1] else "section", reference_classes = rc, left_out = left_out)
}
