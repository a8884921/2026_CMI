# =============================================================================
# Bulk CMI-riskscore engine
# =============================================================================
# The REO/RSF kernel of v1.2.0 (engine/bulk_kernel.R) as functions, with the application settings
# it ran under (duplicate policy max_mean, HGNChelper on, normalized symbols and legacy aliases on,
# minimum pair coverage 0.70). The forest is scored by rsf_predict() (engine/models.R), the portable
# copy of the randomForestSRC model.
# =============================================================================

BULK_SETTINGS <- list(duplicate_gene_policy = "max_mean", use_hgnc = TRUE, use_normalized_symbols = TRUE,
                      use_legacy_aliases = TRUE, hgnc_species = "human", minimum_pair_coverage = 0.70)

bulk_alias_table <- list(
  "EEIG1" = c("FAM102A"), "GARRE1" = c("KIAA0355"), "H1-0" = c("H1F0"), "H1-10" = c("H1FX"), "HLA-DQA1" = c("HLA-DQA"),
  "ICOSLG" = c("ICOSL", "B7H2", "B7RP1", "CD275"), "IRAG2" = c("LRMP", "JAW1"), "JPT1" = c("HN1"), "KMT2E" = c("MLL5"),
  "MARCHF8" = c("MARCH8"), "RGCC" = c("C13orf15", "RGC32"), "SEPTIN9" = c("SEPT9", "MSF"), "CYTH1" = c("PSCD1"),
  "DCTPP1" = c("XTP3TPA", "RS21C6"), "PSMG1" = c("DSCR2", "PAC1"), "RABL6" = c("C9orf86", "RBEL1", "PARF"), "BOLA3" = c("C2orf14"),
  "SPOCK2" = c("KIAA0275"), "CENPM" = c("C22orf18", "PANE1"), "POLR3K" = c("C16orf14", "RPC11"), "RHOH" = c("ARHH", "TTF"),
  "HELLS" = c("SMARCA6", "LSH", "PASG"), "RAD51AP1" = c("PIR51"), "RUVBL1" = c("TIP49", "PONTIN"), "RUVBL2" = c("TIP48", "REPTIN"),
  "MAD2L1" = c("MAD2"), "NME1" = c("NM23", "NM23-H1"), "SQSTM1" = c("P62", "OSIL"), "TRIM28" = c("KAP1", "TIF1B"),
  "DUSP4" = c("MKP2"), "TFDP1" = c("DP1", "DRTF1"), "ABCG1" = c("ABC8", "WHITE1"), "MBNL1" = c("MBNL", "EXP"),
  "NR3C1" = c("GRL", "GCR"), "ARID4B" = c("RBP1L1", "BRCAA1"), "MEGF6" = c("EGFL3"), "SMC4" = c("SMC4L1", "CAPC"),
  "LY9" = c("SLAMF3", "CD229"), "PTPRC" = c("CD45"), "NUDT1" = c("MTH1"), "ECI1" = c("DCI"), "TIMM8B" = c("DDP2"),
  "RASA3" = c("GAP1IP4BP"), "AUTS2" = c("KIAA0442"), "ZFP36L1" = c("TIS11B"))

clean_text <- function(x) trimws(as.character(x))
symbol_key <- function(x) toupper(gsub("[^A-Za-z0-9]", "", clean_text(x)))

standardize_symbols <- function(genes, source_name, make_audit = TRUE) {
  original <- clean_text(genes)
  if (mean(grepl("^ENSG[0-9]+(\\.[0-9]+)?$", original, ignore.case = TRUE), na.rm = TRUE) >= 0.50) {
    input_error(source_name, " uses Ensembl identifiers. Convert them to HGNC gene symbols before uploading.")
  }
  checked <- suppressMessages(suppressWarnings(HGNChelper::checkGeneSymbols(original, unmapped.as.na = TRUE,
                                                                             species = BULK_SETTINGS$hgnc_species, expand.ambiguous = FALSE)))
  if (!is.data.frame(checked) || nrow(checked) != length(original) || !"Suggested.Symbol" %in% names(checked)) stop("Unexpected HGNChelper result: ", source_name)
  suggested <- clean_text(checked$Suggested.Symbol)
  suggested[is.na(checked$Suggested.Symbol) | !nzchar(suggested)] <- NA_character_
  ambiguous <- !is.na(suggested) & grepl("\\s*///\\s*", suggested)
  single_map <- !is.na(suggested) & !ambiguous
  standardized <- original; standardized[single_map] <- suggested[single_map]
  approved <- if ("Approved" %in% names(checked)) as.logical(checked$Approved) else rep(NA, length(original))
  status <- rep("unmapped_keep_original", length(original))
  status[ambiguous] <- "ambiguous_keep_original"
  status[single_map & standardized == original & approved %in% TRUE] <- "approved_unchanged"
  status[single_map & standardized == original & !(approved %in% TRUE)] <- "suggested_same"
  status[single_map & standardized != original] <- "updated_to_HGNC"
  list(symbol = standardized,
       audit = if (make_audit) data.frame(source = source_name, original_symbol = original, approved = approved, suggested_symbol = suggested,
                                          standardized_symbol = standardized, status = status, stringsAsFactors = FALSE) else NULL)
}

read_reo_pairs <- function(path) {
  obj <- data.table::fread(path, data.table = FALSE, check.names = FALSE, showProgress = FALSE, colClasses = "character")
  key <- gsub("[^a-z0-9]", "", tolower(names(obj)))
  col <- function(cands, label) { h <- which(key %in% cands); if (length(h) != 1L) stop("Cannot identify the ", label, " column in the REO pair file."); h }
  out <- data.frame(feature = clean_text(obj[[col(c("feature", "reoid", "reofeature", "pairid"), "feature")]]),
                    gene_A = clean_text(obj[[col(c("genea", "a", "gene1"), "gene_A")]]),
                    gene_B = clean_text(obj[[col(c("geneb", "b", "gene2"), "gene_B")]]), stringsAsFactors = FALSE)
  out
}

# Expression table: genes in the first column, samples in the remaining columns. Streamed in two
# passes (count the rows, then parse blocks of rows with fread into a preallocated matrix), so the
# peak is one numeric matrix plus one block; the validation contract is the v1.2.0 one.
BULK_NA_STRINGS <- c("", "NA", "NaN", "nan", "-", ".")
count_lines <- function(path, chunk = 16 * 1024^2) {             # newline bytes, read in 16 MB chunks
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rb") else file(path, "rb")
  on.exit(close(con))
  n <- 0; last <- as.raw(10)
  repeat { b <- readBin(con, "raw", chunk); if (!length(b)) break; n <- n + sum(b == as.raw(10)); last <- b[length(b)]; rm(b) }
  if (last != as.raw(10)) n <- n + 1
  as.integer(n)
}
read_bulk_expression <- function(path, max_entries = Inf, block_rows = 2000L) {
  name <- tolower(basename(path))
  if (!grepl("\\.(csv|tsv|txt)(\\.gz)?$", name)) input_error("Upload the expression matrix as CSV, TSV or TXT (genes in rows, samples in columns).")
  sep <- if (grepl("\\.csv(\\.gz)?$", name)) "," else "\t"
  parse_block <- function(lines) tryCatch(data.table::fread(text = paste0(paste(lines, collapse = "\n"), "\n"), sep = sep, header = FALSE, data.table = FALSE, showProgress = FALSE,
                                                            na.strings = BULK_NA_STRINGS, colClasses = list(character = 1L), fill = FALSE),
                                          error = function(e) input_error("The expression file could not be parsed: ", conditionMessage(e)))
  con <- open_text(path)
  head2 <- readLines(con, n = 2L, warn = FALSE)
  close(con)
  if (length(head2) < 2L || !nzchar(head2[2])) input_error("The expression file has no data rows.")
  n <- count_lines(path) - 1L                                     # upper bound: blank lines are skipped below
  fields <- split_header(head2[1], sep)
  k_data <- ncol(parse_block(head2[2])) - 1L
  if (length(fields) == k_data) fields <- c("gene", fields)            # header without a name for the gene column
  k <- length(fields) - 1L
  if (k < 1L || k != k_data) input_error("Use gene symbols in the first column and one column per sample (the header has ", length(fields), " fields, the first row ", k_data + 1L, ").")
  if (as.numeric(n) * k > max_entries) input_error(sprintf("The matrix has %s values, above this server's limit of %s. Score the cohort in batches of samples (each sample's score depends only on its own genes) or use the local edition.",
                                                          format(as.numeric(n) * k, big.mark = ",", scientific = FALSE), format(max_entries, big.mark = ",", scientific = FALSE)))
  ids <- checked_ids(fields[-1], "Sample identifiers")
  x <- matrix(NA_real_, n, k); genes <- character(n); r <- 0L
  con <- open_text(path); on.exit(close(con), add = TRUE)
  invisible(readLines(con, n = 1L, warn = FALSE))
  repeat {
    l <- readLines(con, n = block_rows, warn = FALSE); if (!length(l)) break
    l <- l[nzchar(l)]; if (!length(l)) next
    b <- parse_block(l)
    if (ncol(b) != k + 1L) input_error("Rows ", r + 1L, " to ", r + length(l), " do not have ", k + 1L, " fields like the header.")
    idx <- r + seq_len(nrow(b)); genes[idx] <- as.character(b[[1]])
    for (j in seq_len(k)) {
      v <- b[[j + 1L]]
      if (!is.numeric(v)) {
        num <- suppressWarnings(as.numeric(as.character(v)))
        if (any(is.na(num) & !is.na(v))) input_error("Sample column ", ids[j], " contains non-numeric values (for example: ", paste(utils::head(unique(v[is.na(num) & !is.na(v)]), 3), collapse = ", "), ").")
        v <- num
      }
      x[idx, j] <- v
    }
    r <- r + nrow(b); rm(b, l); invisible(gc(verbose = FALSE))         # keep the heap at one matrix plus one block
  }
  if (r < n) { x <- x[seq_len(r), , drop = FALSE]; genes <- genes[seq_len(r)] }
  genes <- checked_ids(genes, "Gene identifiers", allow_duplicates = TRUE)
  for (j in seq_len(k)) {
    v <- x[, j]
    if (any(is.infinite(v))) input_error("Expression values must be finite.")
    if (any(v < 0, na.rm = TRUE)) input_error("Expression must be nonnegative (for example log2 intensities, TPM, FPKM or counts); sample ", ids[j], " has negative values. ",
                                              "Gene-centered, z-scored or log-ratio matrices change the order of genes within a sample, which the REO features depend on. ",
                                              "If the negative values come only from log-transforming low intensities, add one constant to the whole matrix: that leaves the within-sample order unchanged.")
    if (all(is.na(v))) input_error("At least one sample has no observed values (", ids[j], ").")
  }
  invisible(gc(verbose = FALSE))
  dimnames(x) <- list(genes, ids)
  x
}

resolve_duplicate_genes <- function(m) {
  gene <- clean_text(rownames(m)); rownames(m) <- gene
  bad <- is.na(gene) | !nzchar(gene)
  if (any(bad)) { m <- m[!bad, , drop = FALSE]; gene <- rownames(m) }
  if (!nrow(m)) input_error("No genes remain after empty identifiers are removed.")
  if (!anyDuplicated(gene)) return(m)
  # For each symbol keep the row with the highest mean ("max_mean"); output in split() order, as v1.2.0.
  idx_list <- split(seq_len(nrow(m)), gene)
  pick <- vapply(idx_list, function(ix) {
    if (length(ix) == 1L) return(ix)
    mu <- rowMeans(m[ix, , drop = FALSE], na.rm = TRUE); mu[!is.finite(mu)] <- -Inf
    ix[which.max(mu)] }, integer(1))
  out <- m[pick, , drop = FALSE]; rownames(out) <- names(idx_list); storage.mode(out) <- "double"
  out
}

map_endpoints <- function(source_genes, target_genes) {
  source_genes <- clean_text(source_genes); target_genes <- clean_text(target_genes); source_key <- symbol_key(source_genes)
  out <- data.frame(target_gene = target_genes, source_gene = NA_character_, source_index = NA_integer_, match_type = "missing", stringsAsFactors = FALSE)
  exact <- match(target_genes, source_genes); ok <- !is.na(exact)
  out$source_gene[ok] <- source_genes[exact[ok]]; out$source_index[ok] <- exact[ok]; out$match_type[ok] <- "exact"
  if (BULK_SETTINGS$use_normalized_symbols) for (i in which(is.na(out$source_index))) {
    hits <- which(source_key == symbol_key(target_genes[i]))
    if (length(hits) == 1L) { out$source_gene[i] <- source_genes[hits]; out$source_index[i] <- hits; out$match_type[i] <- "normalized" }
  }
  if (BULK_SETTINGS$use_legacy_aliases) for (i in which(is.na(out$source_index))) {
    aliases <- bulk_alias_table[[target_genes[i]]]
    if (is.null(aliases)) next
    for (a in aliases) { hits <- which(source_key == symbol_key(a))
      if (length(hits) == 1L) { out$source_gene[i] <- source_genes[hits]; out$source_index[i] <- hits; out$match_type[i] <- paste0("alias:", a); break } }
  }
  used <- out$source_index[!is.na(out$source_index)]
  if (anyDuplicated(used)) {
    rep_idx <- unique(used[duplicated(used)])
    stop("Ambiguous one-to-many gene mapping: ", paste(vapply(rep_idx, function(ix) paste(out$target_gene[out$source_index %in% ix], collapse = " / "), character(1)), collapse = "; "))
  }
  out
}

# REO features: within-sample average ranks over all genes, divided by the number of genes, then
# rank(A) - rank(B) for each pair. Samples are processed in blocks so only the endpoint ranks are
# kept; the arithmetic (median fill, rowRanks, division, difference) is the v1.2.0 one.
build_reo <- function(m_gene_sample, dataset_name, gene_pairs, endpoint_genes, block = 200L) {
  n_gene <- nrow(m_gene_sample); n_smp <- ncol(m_gene_sample)
  if (n_gene < 2L) input_error(dataset_name, ": at least two genes are required for ranking.")
  map <- map_endpoints(rownames(m_gene_sample), endpoint_genes)
  endpoint_index <- stats::setNames(map$source_index, map$target_gene)
  pair_ok <- !is.na(endpoint_index[gene_pairs$gene_A]) & !is.na(endpoint_index[gene_pairs$gene_B])
  feat <- matrix(NA_real_, nrow = n_smp, ncol = nrow(gene_pairs), dimnames = list(colnames(m_gene_sample), gene_pairs$feature))
  if (any(pair_ok)) {
    ia <- as.integer(endpoint_index[gene_pairs$gene_A[pair_ok]]); ib <- as.integer(endpoint_index[gene_pairs$gene_B[pair_ok]])
    for (b in split(seq_len(n_smp), ceiling(seq_len(n_smp) / block))) {
      X <- t(m_gene_sample[, b, drop = FALSE])
      if (any(is.infinite(X))) input_error(dataset_name, ": expression contains Inf or -Inf.")
      for (i in seq_len(nrow(X))) {
        bad <- !is.finite(X[i, ])
        if (any(bad)) { fill <- stats::median(X[i, !bad], na.rm = TRUE); if (!is.finite(fill)) fill <- 0; X[i, bad] <- fill }
      }
      Rk <- matrixStats::rowRanks(X, ties.method = "average") / n_gene
      feat[b, pair_ok] <- Rk[, ia, drop = FALSE] - Rk[, ib, drop = FALSE]
      rm(X, Rk); invisible(gc(verbose = FALSE))
    }
  }
  pair_coverage <- mean(pair_ok)
  coverage_class <- if (pair_coverage >= 0.90) "Adequate" else if (pair_coverage >= BULK_SETTINGS$minimum_pair_coverage) "Moderate" else "Limited"
  qc <- data.frame(dataset = dataset_name, samples = n_smp, expression_genes = n_gene, endpoint_genes_total = length(endpoint_genes),
                   endpoint_genes_observed = sum(!is.na(map$source_index)), endpoint_coverage = mean(!is.na(map$source_index)),
                   gene_pairs_total = nrow(gene_pairs), gene_pairs_computable = sum(pair_ok), gene_pairs_missing = sum(!pair_ok),
                   pair_coverage = pair_coverage, exact_matches = sum(map$match_type == "exact"), normalized_matches = sum(map$match_type == "normalized"),
                   alias_matches = sum(grepl("^alias:", map$match_type)), coverage_class = coverage_class, stringsAsFactors = FALSE)
  list(feature = feat, mapping = map, pair_ok = pair_ok, missing_pairs = gene_pairs[!pair_ok, , drop = FALSE], qc = qc)
}

parameter_vector <- function(value, safe_names, original_names, default_value) {
  if (is.null(value) || !length(value)) return(stats::setNames(rep(default_value, length(safe_names)), safe_names))
  if (!is.null(names(value)) && all(nzchar(names(value)))) {
    out <- rep(NA_real_, length(safe_names)); names(out) <- safe_names
    hit <- match(safe_names, names(value)); ok <- !is.na(hit); out[ok] <- as.numeric(value[hit[ok]])
    todo <- which(!is.finite(out))
    if (length(todo)) { h2 <- match(original_names[todo], names(value)); ok2 <- !is.na(h2); out[todo[ok2]] <- as.numeric(value[h2[ok2]]) }
    out[!is.finite(out)] <- default_value
    return(out)
  }
  if (length(value) == length(safe_names)) { out <- as.numeric(value); names(out) <- safe_names; out[!is.finite(out)] <- default_value; return(out) }
  stats::setNames(rep(default_value, length(safe_names)), safe_names)
}

apply_model_preprocessing <- function(feature_matrix, prep) {
  original <- clean_text(prep$original_features %||% colnames(feature_matrix))
  miss <- setdiff(original, colnames(feature_matrix))
  if (length(miss)) stop("The REO matrix lacks model features: ", paste(utils::head(miss, 20), collapse = ", "))
  x <- feature_matrix[, original, drop = FALSE]
  safe <- prep$safe_features
  if (is.null(safe) && !is.null(prep$feature_map) && all(c("original", "safe") %in% names(prep$feature_map))) {
    hit <- match(original, clean_text(prep$feature_map$original)); if (all(!is.na(hit))) safe <- clean_text(prep$feature_map$safe[hit])
  }
  if (is.null(safe) || length(safe) != length(original) || anyNA(safe)) safe <- make.names(original, unique = TRUE)
  safe <- clean_text(safe); colnames(x) <- safe
  medians <- parameter_vector(prep$medians, safe, original, default_value = NA_real_)
  for (j in seq_len(ncol(x))) { bad <- !is.finite(x[, j]); if (any(bad)) x[bad, j] <- medians[j] }
  if (any(!is.finite(x))) stop("Missing REO values require finite stored medians.")
  do_standardize <- if (is.list(prep$standardization)) isTRUE(prep$standardization$standardize) else isTRUE(prep$standardization)
  if (do_standardize) {
    centers <- parameter_vector(prep$centers, safe, original, default_value = NA_real_)
    scales <- parameter_vector(prep$scales, safe, original, default_value = NA_real_)
    if (any(!is.finite(centers)) || any(!is.finite(scales)) || any(scales <= 0)) stop("Invalid stored standardization parameters.")
    x <- sweep(x, 2, centers, "-"); x <- sweep(x, 2, scales, "/")
  }
  list(matrix = x, original = original, safe = safe)
}

resolve_model_names <- function(needed, prep_result) {
  needed <- clean_text(needed); available <- colnames(prep_result$matrix)
  if (all(needed %in% available)) return(needed)
  converted <- unname(stats::setNames(prep_result$safe, prep_result$original)[needed])
  if (length(converted) == length(needed) && all(!is.na(converted)) && all(converted %in% available)) return(converted)
  made <- make.names(needed, unique = TRUE); if (all(made %in% available)) return(made)
  stop("Selected model features cannot be aligned with the preprocessed REO columns.")
}

predict_riskscore <- function(rsf, feature_matrix) {
  prep_result <- apply_model_preprocessing(feature_matrix, rsf$preprocessing)
  needed <- resolve_model_names(rsf$selector_selected_features %||% rsf$feature_order %||% rsf$xvar_names, prep_result)
  x <- prep_result$matrix[, needed, drop = FALSE]
  forest_needed <- clean_text(rsf$xvar_names)
  if (all(forest_needed %in% colnames(x))) x <- x[, forest_needed, drop = FALSE]
  else stop("Cannot align the forest features.")
  raw <- rsf_predict(rsf, as.data.frame(x, check.names = FALSE))
  orientation <- as.numeric(rsf$risk_orientation %||% NA_real_)
  if (length(orientation) != 1L || !is.finite(orientation) || orientation == 0) stop("The model must contain a finite non-zero risk orientation.")
  risk <- orientation * raw
  if (any(!is.finite(risk))) stop("The forest returned non-finite scores.")
  stats::setNames(risk, rownames(feature_matrix))
}

run_bulk <- function(expression_path, rsf, reo_pairs_path, dataset_name = "upload", lim = limits_for(), progress = function(stage, pct, detail) NULL) {
  progress("pairs", 10, "Standardizing the 100 REO gene-pair endpoints (HGNC)")
  gp <- read_reo_pairs(reo_pairs_path)
  original_model_features <- clean_text(rsf$preprocessing$original_features)
  if (!all(original_model_features %in% gp$feature)) stop("The REO pair file does not cover every model feature.")
  pair_meta <- gp
  gp <- gp[match(original_model_features, gp$feature), c("feature", "gene_A", "gene_B"), drop = FALSE]
  A <- standardize_symbols(gp$gene_A, "REO:gene_A"); B <- standardize_symbols(gp$gene_B, "REO:gene_B")
  A$audit$feature <- gp$feature; A$audit$pair_side <- "gene_A"; B$audit$feature <- gp$feature; B$audit$pair_side <- "gene_B"
  reo_audit <- rbind(A$audit, B$audit)
  gp$gene_A <- A$symbol; gp$gene_B <- B$symbol
  if (any(gp$gene_A == gp$gene_B)) stop("Self-pairs after HGNC standardization.")
  if (anyDuplicated(paste(gp$gene_A, gp$gene_B, sep = "\r"))) stop("Two features map to the same directed pair after HGNC standardization.")
  endpoint_genes <- unique(c(gp$gene_A, gp$gene_B))

  progress("read", 25, "Reading and validating the expression matrix")
  m <- read_bulk_expression(expression_path, max_entries = lim$max_bulk_entries)
  if (ncol(m) < 2L) input_error("At least two samples are required; the CMI-high/low split uses the within-cohort median.")
  progress("hgnc", 40, sprintf("Standardizing %s gene symbols (HGNChelper)", format(nrow(m), big.mark = ",")))
  h <- standardize_symbols(rownames(m), paste0(dataset_name, ":expression"))
  rownames(m) <- h$symbol
  m <- resolve_duplicate_genes(m)
  colnames(m) <- clean_text(colnames(m))

  progress("reo", 60, "Within-sample ranks and REO rank differences")
  reo <- build_reo(m, dataset_name, gp, endpoint_genes)
  n_gene <- nrow(m); rm(m); invisible(gc(verbose = FALSE))
  progress("rsf", 75, "Scoring the random survival forest (800 trees)")
  risk <- predict_riskscore(rsf, reo$feature)
  score <- data.frame(sample_id = names(risk), riskscore = as.numeric(risk), dataset = dataset_name, rank_pool_genes = n_gene,
                      endpoint_coverage = reo$qc$endpoint_coverage, gene_pairs_computable = reo$qc$gene_pairs_computable,
                      gene_pairs_total = reo$qc$gene_pairs_total, pair_coverage = reo$qc$pair_coverage, coverage_class = reo$qc$coverage_class,
                      stringsAsFactors = FALSE)
  s <- stats::sd(score$riskscore)
  score$riskscore_z <- if (nrow(score) > 1L && is.finite(s) && s > 0) as.numeric(scale(score$riskscore)) else 0
  score$percentile <- 100 * rank(score$riskscore, ties.method = "average") / nrow(score)
  med <- stats::median(score$riskscore)
  score$cmi_group <- ifelse(score$riskscore > med, "CMI-high", "CMI-low")
  list(scores = score, median = med, qc = reo$qc, mapping = reo$mapping, missing_pairs = reo$missing_pairs,
       features = reo$feature, reo_audit = reo_audit, expression_audit = h$audit, pairs = pair_meta)
}
