# =============================================================================
# Single-cell CMI engine (memory-lean, dependency-light)
# =============================================================================
# Reproduces the v1.2.0 chain build_object.R -> four_indices.R -> cnv_r.R -> cmi_r.R
# with four substitutions, each validated in tests/test_equivalence.R:
#   * Seurat NormalizeData + CellCycleScoring      -> exact ports (log1p(x / total * 1e4);
#                                                    AddModuleScore with seed 1, 24 bins)
#   * UCell::ScoreSignatures_UCell (2.10)           -> exact rank port (maxRank 1500, average ties)
#   * qlcMatrix::corSparse                          -> same formula
#   * AnnotationDbi + org.Hs.eg.db queries          -> ref/annotation_snapshot.rds (same tables)
# and processes cells in chunks so peak memory scales with one chunk, not the whole matrix.
# =============================================================================

# ---- gene-symbol helpers (four_indices.R build_symbol_query_map / resolve_map) ----------------
sc_symbol_map <- function(set, universe, ann, namespace) {
  tab <- ann$symbol_tables[[set]]
  symbols <- tab$symbols
  universe <- unique(as.character(universe))
  out <- stats::setNames(vector("list", length(symbols)), symbols)
  uu <- toupper(universe); uq <- !duplicated(uu) & !duplicated(uu, fromLast = TRUE)
  lut <- stats::setNames(universe[uq], uu[uq])
  for (g in symbols) {
    v <- character()
    if (g %in% universe) v <- g else if (toupper(g) %in% names(lut)) v <- unname(lut[toupper(g)])
    out[[g]] <- v
  }
  if (length(tab$valid)) {
    if (identical(namespace, "ENSEMBL")) {
      z <- tab$ensembl; ub <- sub("\\.[0-9]+$", "", universe)
      for (g in tab$valid) { ens <- unique(z$ENSEMBL[z$SYMBOL == g]); if (length(ens)) out[[g]] <- unique(c(out[[g]], universe[ub %in% ens])) }
    } else {
      z <- tab$alias
      for (g in tab$valid) { al <- unique(z$ALIAS[z$SYMBOL == g]); if (length(al)) out[[g]] <- unique(c(out[[g]], intersect(al, universe))) }
    }
  }
  out
}
sc_resolve <- function(mp, genes) unique(unlist(mp[intersect(as.character(genes), names(mp))], use.names = FALSE))

sc_map_entrez <- function(features, ann) {          # four_indices.R map_query_to_entrez (network-relevant part)
  features <- as.character(features)
  e <- ann$entrez
  if (mean(grepl("^ENSG[0-9]+(\\.[0-9]+)?$", features)) > 0.5) {
    key <- sub("\\.[0-9]+$", "", features)
    return(unname(e$ensembl_first[key]))
  }
  out <- unname(e$symbol_first[features])
  # a valid org.Hs.eg.db symbol is resolved at the symbol step even when its Entrez ID is outside the
  # PPI network; it must not fall through to the alias table
  miss <- is.na(out) & !(features %in% e$alias_blocked_by_symbol)
  out[miss] <- unname(e$alias_first[features[miss]])
  out
}

# ---- Ig metrics (four_indices.R compute_ig_indices, unchanged arithmetic) ---------------------
sc_ig_indices <- function(counts, IG_MAP) {
  gene_sum <- function(m, g) { ids <- intersect(sc_resolve(IG_MAP, g), rownames(m)); if (!length(ids)) return(rep(0, ncol(m))); Matrix::colSums(m[ids, , drop = FALSE]) }
  gene_vec <- function(m, g) as.numeric(gene_sum(m, g))
  ratio <- function(num, den, md = IG_MIN_UMI) ifelse(den >= md, num / den, NA_real_)
  evenness <- function(parts, md = IG_MIN_UMI) { P <- do.call(cbind, parts); tot <- rowSums(P); k <- ncol(P); out <- rep(NA_real_, length(tot)); ok <- tot >= md
    Pn <- P[ok, , drop = FALSE] / tot[ok]; Pn[Pn <= 0] <- 1; out[ok] <- -rowSums(Pn * log(Pn)) / log(k); out }
  total <- as.numeric(Matrix::colSums(counts))
  needed <- intersect(sc_resolve(IG_MAP, IG_CANON), rownames(counts))
  m <- methods::as(counts[needed, , drop = FALSE], "dgCMatrix")
  H <- gene_sum(m, IG_G_HEAVY); IgGt <- gene_sum(m, IG_G_IGG); SW <- gene_sum(m, IG_G_SWIT); IgAll <- gene_sum(m, IG_G_PANEL)
  vM <- gene_vec(m, "IGHM"); vD <- gene_vec(m, "IGHD"); vA1 <- gene_vec(m, "IGHA1"); vA2 <- gene_vec(m, "IGHA2"); vE <- gene_vec(m, "IGHE")
  vG1 <- gene_vec(m, "IGHG1"); vG2 <- gene_vec(m, "IGHG2"); vG3 <- gene_vec(m, "IGHG3"); vG4 <- gene_vec(m, "IGHG4")
  vJ <- gene_vec(m, c("JCHAIN", "IGJ")); ordnum <- Reduce(`+`, lapply(names(IG_LOCUS), function(g) IG_LOCUS[[g]] * gene_vec(m, g))); vIgA <- vA1 + vA2
  X <- data.frame(Ig_burden = IgAll / total, CSR_switched = ratio(SW, H), IgM_fraction = ratio(vM, H), IgD_fraction = ratio(vD, H), IgG_fraction = ratio(IgGt, H), IgA_fraction = ratio(vIgA, H),
    IgG_IgA_effector = ifelse(SW >= IG_MIN_UMI, ratio(IgGt, IgGt + vIgA, IG_MIN_UMI), NA_real_), J_chain_secretory = ratio(vJ, vJ + H), IgG1_fraction = ratio(vG1, IgGt), IgG2_fraction = ratio(vG2, IgGt),
    IgM_IgG_switch = ratio(vM, vM + IgGt), IgM_IgA_switch = ratio(vM, vM + vIgA), polymeric_isotype = ratio(vM + vIgA, H), mucosal_secretory = ratio(vA2 + vJ, vA2 + vJ + IgGt),
    switch_locus_order = ratio(ordnum, H), IgA1_fraction = ratio(vA1, H), IgA1_IgA2_ratio = ratio(vA1, vA1 + vA2), IgM_IgD_coexpr = ratio(2 * pmin(vM, vD), vM + vD),
    IgG3_fraction = ratio(vG3, IgGt), IgG4_fraction = ratio(vG4, IgGt), IgE_fraction = ratio(vE, H), IgD_IgM_ratio = ratio(vD, vD + vM), IgG_over_switched = ratio(IgGt, SW),
    heavy_isotype_evenness = evenness(list(vM, vD, vG1, vG2, vG3, vG4, vA1, vA2)), check.names = FALSE)
  X[] <- lapply(X, function(v) { v <- as.numeric(v); v[!is.finite(v)] <- NA; v })
  colnames(X) <- make.names(colnames(X)); rownames(X) <- colnames(counts)
  X
}

# ---- UCell 2.10 U statistic (data_to_ranks_data_table + u_stat), sparse-aware --------------------
# Ranks are descending with average ties over all genes of a cell; zeros share the average rank
# (nnz + 1 + G) / 2; ranks >= maxRank are set to maxRank; genes absent from the matrix count as
# zero-expressed rows, as UCell 2.10's check_genes() adds them.
ucell_scores <- function(m, programs, max_rank = 1500, chunk = 500L) {
  all_genes <- unique(unlist(programs, use.names = FALSE))
  extra <- length(setdiff(all_genes, rownames(m)))
  G <- nrow(m) + extra
  if (max_rank > G) max_rank <- G
  if (any(lengths(programs) > max_rank)) stop("A program has more genes than maxRank.")
  nc <- ncol(m)
  out <- matrix(NA_real_, nc, length(programs), dimnames = list(colnames(m), names(programs)))
  rows_needed <- intersect(all_genes, rownames(m))
  for (s in seq(1L, nc, by = chunk)) {
    cols <- s:min(s + chunk - 1L, nc)
    mc <- m[, cols, drop = FALSE]
    x <- mc@x; p <- mc@p; nnz <- diff(p); ncc <- length(cols)
    j <- rep.int(seq_len(ncc), nnz)
    o <- order(j, -x)
    xs <- x[o]; js <- j[o]
    pos <- seq_along(o) - p[js]
    brk <- c(TRUE, js[-1L] != js[-length(js)] | xs[-1L] != xs[-length(xs)])
    first <- which(brk); last <- c(first[-1L] - 1L, length(o))
    r <- numeric(length(x)); r[o] <- rep.int((pos[first] + pos[last]) / 2, last - first + 1L)
    zero_rank <- (nnz + 1 + G) / 2
    R <- mc; R@x <- r
    sub <- as.matrix(R[rows_needed, , drop = FALSE])
    zr <- matrix(zero_rank, nrow = nrow(sub), ncol = ncc, byrow = TRUE)
    is0 <- sub == 0; sub[is0] <- zr[is0]
    zcap <- ifelse(zero_rank >= max_rank, max_rank, zero_rank)
    for (k in seq_along(programs)) {
      genes <- programs[[k]]; present <- genes[genes %in% rows_needed]; L <- length(genes)
      rv <- sub[present, , drop = FALSE]
      rv[rv >= max_rank] <- max_rank
      rank_sum <- colSums(rv) + (L - length(present)) * zcap
      rank_sum_min <- L * (L + 1) / 2
      u <- 1 - (rank_sum - rank_sum_min) / (L * max_rank - rank_sum_min)
      u[u < 0] <- 0
      out[cols, k] <- u
    }
    rm(mc, x, j, o, xs, js, pos, brk, r, R, sub, zr, is0); invisible(gc(verbose = FALSE))
  }
  out
}

# ---- CNV (engine/cnv_r.R, infercnvpy 0.6.1 arithmetic), per cell batch ----------------------------
npy_msb <- function(unum) { d <- 0L; u <- unum; while (TRUE) { u <- u %/% 2L; if (u == 0L) break; d <- d + 1L }; d }
aheapsort_R <- function(v, tosort, pl, pr) {
  n <- pr - pl + 1L; if (n <= 1L) return(tosort); base <- pl - 1L
  al <- function(i) tosort[base + i]; aset <- function(i, val) tosort[base + i] <<- val
  l <- n %/% 2L
  while (l > 0L) { tmp <- al(l); i <- l; j <- l * 2L
    while (j <= n) { if (j < n && v[al(j)] < v[al(j + 1L)]) j <- j + 1L; if (v[tmp] < v[al(j)]) { aset(i, al(j)); i <- j; j <- j + j } else break }
    aset(i, tmp); l <- l - 1L }
  while (n > 1L) { tmp <- al(n); aset(n, al(1L)); n <- n - 1L; i <- 1L; j <- 2L
    while (j <= n) { if (j < n && v[al(j)] < v[al(j + 1L)]) j <- j + 1L; if (v[tmp] < v[al(j)]) { aset(i, al(j)); i <- j; j <- j + j } else break }
    aset(i, tmp) }
  tosort
}
np_order_quicksort <- function(v) {                     # numpy aquicksort tie order (cnv_r.R)
  n <- length(v); if (n <= 1L) return(seq_len(n)); if (!anyDuplicated(v)) return(order(v))
  tosort <- seq_len(n); cdepth <- npy_msb(n) * 2L; st_pl <- integer(0L); st_pr <- integer(0L); st_d <- integer(0L); pl <- 1L; pr <- n
  repeat {
    if (cdepth < 0L) { tosort <- aheapsort_R(v, tosort, pl, pr) } else {
      while ((pr - pl) > 15L) {
        pm <- pl + ((pr - pl) %/% 2L)
        a <- tosort[pm]; b <- tosort[pl]; if (v[a] < v[b]) { tosort[pm] <- b; tosort[pl] <- a }
        a <- tosort[pr]; b <- tosort[pm]; if (v[a] < v[b]) { tosort[pr] <- b; tosort[pm] <- a }
        a <- tosort[pm]; b <- tosort[pl]; if (v[a] < v[b]) { tosort[pm] <- b; tosort[pl] <- a }
        vp <- v[tosort[pm]]; pi <- pl; pj <- pr - 1L
        tmp <- tosort[pm]; tosort[pm] <- tosort[pj]; tosort[pj] <- tmp
        repeat { repeat { pi <- pi + 1L; if (!(v[tosort[pi]] < vp)) break }; repeat { pj <- pj - 1L; if (!(vp < v[tosort[pj]])) break }
          if (pi >= pj) break; tmp <- tosort[pi]; tosort[pi] <- tosort[pj]; tosort[pj] <- tmp }
        pk <- pr - 1L; tmp <- tosort[pi]; tosort[pi] <- tosort[pk]; tosort[pk] <- tmp
        if ((pi - pl) < (pr - pi)) { st_pl <- c(st_pl, pi + 1L); st_pr <- c(st_pr, pr); pr <- pi - 1L } else { st_pl <- c(st_pl, pl); st_pr <- c(st_pr, pi - 1L); pl <- pi + 1L }
        cdepth <- cdepth - 1L; st_d <- c(st_d, cdepth)
      }
      if (pr > pl) { for (pi in (pl + 1L):pr) { vi <- tosort[pi]; vp <- v[vi]; pj <- pi; pk <- pi - 1L
        while (pj > pl && vp < v[tosort[pk]]) { tosort[pj] <- tosort[pk]; pj <- pj - 1L; pk <- pk - 1L }; tosort[pj] <- vi } }
    }
    if (length(st_pl) == 0L) break
    k <- length(st_pl); pl <- st_pl[k]; pr <- st_pr[k]; cdepth <- st_d[k]; st_pl <- st_pl[-k]; st_pr <- st_pr[-k]; st_d <- st_d[-k]
  }
  tosort
}
natural_sort <- function(x) {
  if (length(x) <= 1L) return(x)
  key <- vapply(x, function(s) { toks <- strsplit(s, "([0-9]+)", perl = TRUE)[[1]]; toks <- toks[nzchar(toks)]
    parts <- vapply(toks, function(t) if (grepl("^[0-9]+$", t)) paste0("\001", sprintf("%012.0f", as.numeric(t))) else paste0("\002", tolower(t)), character(1), USE.NAMES = FALSE)
    paste(parts, collapse = "\003") }, character(1), USE.NAMES = FALSE)
  x[order(key)]
}
row_sd_population <- function(m) { mu <- rowMeans(m); sqrt(rowMeans((m - mu)^2)) }

cnv_gene_map <- function(query_genes, target_genes, ann) {   # cnv_r.R build_gene_map with the snapshot
  q <- trimws(as.character(query_genes)); target <- unique(trimws(as.character(target_genes))); target <- target[nzchar(target)]
  if (anyDuplicated(q)) stop("Query gene identifiers are duplicated.")
  n <- length(q); res_gene <- rep(NA_character_, n); res_method <- rep("unmapped", n)
  assign_map <- function(idx, val, method) { ok <- idx & is.na(res_gene) & !is.na(val) & nzchar(val) & (val %in% target); res_gene[ok] <<- val[ok]; res_method[ok] <<- method }
  all_idx <- rep(TRUE, n)
  m <- match(q, target); assign_map(all_idx, ifelse(is.na(m), NA_character_, target[m]), "exact")
  qb <- sub("\\.[0-9]+$", "", q); m <- match(qb, target); assign_map(all_idx, ifelse(is.na(m), NA_character_, target[m]), "ensembl_version_strip")
  tu <- toupper(target); uq <- !duplicated(tu) & !duplicated(tu, fromLast = TRUE); lut <- stats::setNames(target[uq], tu[uq])
  assign_map(all_idx, unname(lut[toupper(q)]), "casefold_unique")
  assign_map(all_idx, unname(ann$cnv$ensembl_to_symbol[sub("\\.[0-9]+$", "", q)]), "ensembl_to_symbol")
  assign_map(all_idx, unname(ann$cnv$alias_to_symbol[q]), "alias_to_symbol")
  list(gene = res_gene, method = res_method)
}

build_chrom_blocks <- function(chromosome, start, exclude) {
  ok <- !is.na(chromosome) & !(chromosome %in% exclude); chr_unique <- unique(chromosome[ok])
  chroms <- natural_sort(chr_unique[startsWith(chr_unique, "chr") & chr_unique != "chrM"])
  blocks <- vector("list", length(chroms)); names(blocks) <- chroms
  for (k in seq_along(chroms)) { rows <- which(chromosome == chroms[k]); blocks[[k]] <- rows[np_order_quicksort(start[rows])] }
  blocks
}
build_pyramid <- function(ngenes, window_size, step) {
  if (window_size < ngenes) {
    r <- seq_len(window_size); pyr <- pmin(r, rev(r)); nwin <- ngenes - window_size + 1L; sel <- seq.int(1L, nwin, by = step)
    i <- rep(seq_along(sel), times = window_size); j <- as.vector(outer(sel - 1L, seq_len(window_size), "+"))
    list(W = Matrix::sparseMatrix(i = i, j = j, x = rep(as.numeric(pyr), each = length(sel)), dims = c(length(sel), ngenes)), divisor = sum(pyr))
  } else {
    list(W = Matrix::sparseMatrix(i = rep.int(1L, ngenes), j = seq_len(ngenes), x = rep.int(1, ngenes), dims = c(1L, ngenes)), divisor = ngenes)
  }
}

sc_cnv_setup <- function(query_genes, assets) {
  cnv <- assets$cnv; P <- cnv$params
  gm <- cnv_gene_map(query_genes, cnv$universe, assets$ann)
  mapped_unique <- unique(gm$gene[!is.na(gm$gene)])
  common <- cnv$universe[cnv$universe %in% mapped_unique]          # choose_genes(): TASK_COMMON_NEW
  if (length(common) <= P$min_common_genes) input_error("Only ", length(common), " genes overlap the copy-number reference of normal B cells; at least ", P$min_common_genes + 1L, " are required. Upload the whole transcriptome. Use human gene symbols (for example MS4A1, CD79A) or human Ensembl IDs (ENSG...); mouse or custom identifiers are not supported.")
  manifest <- cnv$positions[common, , drop = FALSE]
  sel <- !is.na(gm$gene) & gm$gene %in% common
  gmap <- unique(data.frame(source = query_genes[sel], gene = gm$gene[sel], stringsAsFactors = FALSE))
  src <- unique(gmap$source)
  A <- Matrix::sparseMatrix(i = match(gmap$gene, common), j = match(gmap$source, src), x = 1, dims = c(length(common), length(src)))
  ref_mean <- cnv$ref_mean[match(common, cnv$ref_genes)]
  blocks <- build_chrom_blocks(as.character(manifest$chromosome), as.numeric(manifest$start), P$exclude_chromosomes)
  block_sizes <- vapply(blocks, length, integer(1))
  list(gm = gm, common = common, manifest = manifest, src = src, A = A, ref_mean = ref_mean, blocks = blocks,
       block_sizes = block_sizes, block_cols = unlist(blocks, use.names = FALSE), gene_offsets = c(0L, cumsum(block_sizes)),
       pyramids = lapply(block_sizes, build_pyramid, window_size = P$window_size, step = P$step),
       audit = data.frame(source = query_genes, gene = gm$gene, method = gm$method,
                          used_for_cnv = (!is.na(gm$gene) & gm$gene %in% common), stringsAsFactors = FALSE))
}

sc_cnv_batch <- function(counts_batch, totals_batch, setup, P) {
  raw <- methods::as(counts_batch[setup$src, , drop = FALSE], "dgCMatrix")
  m <- setup$A %*% raw
  sf <- ifelse(is.finite(totals_batch) & totals_batch > 0, 10000 / totals_batch, 0)
  m <- m %*% Matrix::Diagonal(x = sf)
  if (length(m@x)) m@x <- log1p(m@x)
  x <- t(as.matrix(m[setup$block_cols, , drop = FALSE]))
  x <- x - rep(setup$ref_mean[setup$block_cols], each = nrow(x))
  x[x > P$lfc_clip] <- P$lfc_clip; x[x < -P$lfc_clip] <- -P$lfc_clip
  parts <- vector("list", length(setup$blocks))
  for (k in seq_along(setup$blocks)) {
    idx <- setup$gene_offsets[k] + seq_len(setup$block_sizes[k])
    parts[[k]] <- as.matrix((x[, idx, drop = FALSE] %*% Matrix::t(setup$pyramids[[k]]$W)) / setup$pyramids[[k]]$divisor)
  }
  xs <- do.call(cbind, parts)
  xr <- xs - apply(xs, 1L, stats::median)
  thr <- P$threshold_amp * row_sd_population(xr)
  xr[abs(xr) < thr] <- 0
  rowSums(abs(xr)) / ncol(xr)
}

# ---- Seurat AddModuleScore / CellCycleScoring (5.0-5.4 algorithm) --------------------------------
cut_number_exact <- function(x, n) {             # ggplot2::cut_number(x, n, labels = FALSE, right = FALSE)
  brk <- stats::quantile(x, seq(0, 1, length.out = n + 1), na.rm = TRUE)
  if (anyDuplicated(brk)) stop("Insufficient data values to produce ", n, " bins.")
  cut(x, brk, include.lowest = TRUE, labels = FALSE, right = FALSE)
}
sc_module_controls <- function(row_means, features, nbin = 24, ctrl = 100, seed = 1) {
  if (!is.null(seed)) set.seed(seed)
  data.avg <- row_means[order(row_means)]
  data.cut <- cut_number_exact(data.avg + stats::rnorm(n = length(data.avg)) / 1e30, nbin)
  names(data.cut) <- names(data.avg)
  ctrl.use <- vector(mode = "list", length = length(features))
  for (i in seq_along(features)) {
    for (j in seq_along(features[[i]])) {
      ctrl.use[[i]] <- c(ctrl.use[[i]], names(sample(x = data.cut[which(data.cut == data.cut[features[[i]][j]])], size = ctrl, replace = FALSE)))
    }
  }
  lapply(ctrl.use, unique)
}
cc_phase <- function(S, G2M) {
  vapply(seq_along(S), function(k) { sc <- c(S[k], G2M[k]); if (all(sc < 0)) "G1" else if (length(which(sc == max(sc))) > 1) "Undecided" else c("S", "G2M")[which(sc == max(sc))] }, character(1))
}

# ---- CCAT (SCENT; four_indices.R my_CCAT with qlcMatrix::corSparse) -------------------------------
cor_sparse_y <- function(X, Y) {                # qlcMatrix 0.9.x corSparse(X, Y), unchanged formula
  X <- methods::as(X, "dgCMatrix"); n <- nrow(X); muX <- Matrix::colMeans(X)
  Y <- methods::as(Y, "dgCMatrix"); muY <- Matrix::colMeans(Y)
  covmat <- (as.matrix(Matrix::crossprod(X, Y)) - n * tcrossprod(muX, muY)) / (n - 1)
  sdvecX <- sqrt((Matrix::colSums(X^2) - n * muX^2) / (n - 1))
  sdvecY <- sqrt((Matrix::colSums(Y^2) - n * muY^2) / (n - 1))
  cormat <- covmat / tcrossprod(sdvecX, sdvecY)
  dimnames(cormat) <- NULL
  cormat
}

# ---- Full pipeline -----------------------------------------------------------------------------
# `loader` is a function returning the dgCMatrix of raw counts. The engine calls it itself so that it
# holds the only reference to the matrix and can log-normalize it in place (no second copy).
run_single_cell <- function(loader, assets, chunk_cells = 1000L, cnv_batch = 150L, progress = function(stage, pct, detail) NULL) {
  counts <- if (is.function(loader)) loader() else loader
  rm(loader)
  stopifnot(methods::is(counts, "dgCMatrix"))
  cells <- colnames(counts); genes <- rownames(counts); nc <- length(cells)
  ann <- assets$ann
  ens_frac <- mean(grepl("^ENSG[0-9]+(\\.[0-9]+)?$", genes))
  namespace <- if (ens_frac > 0.5) "ENSEMBL" else "SYMBOL_OR_MIXED"
  audit <- data.frame(metric = c("query_genes", "ensembl_fraction"), value = c(length(genes), ens_frac), detail = c(namespace, namespace), stringsAsFactors = FALSE)
  add_audit <- function(metric, value, detail) audit <<- rbind(audit, data.frame(metric = metric, value = value, detail = detail, stringsAsFactors = FALSE))
  totals <- as.numeric(Matrix::colSums(counts))
  chunks <- split(seq_len(nc), ceiling(seq_len(nc) / chunk_cells))
  # Lineage check (reported, not used for scoring): any B-cell marker transcript in each cell, and
  # cells in which T-cell, myeloid, NK, erythroid or platelet marker transcripts dominate (2 or more
  # UMIs and more than the B-cell markers, so that ambient RNA in true B cells is not flagged).
  ent_genes <- sc_map_entrez(genes, ann)
  marker_umis <- function(syms) {
    ent <- unname(ann$entrez$symbol_first[syms]); ent <- ent[!is.na(ent)]
    idx <- which((!is.na(ent_genes) & ent_genes %in% ent) | toupper(genes) %in% syms)
    if (!length(idx)) return(rep(0, nc))
    as.numeric(Matrix::colSums(counts[idx, , drop = FALSE]))
  }
  b_umis <- marker_umis(B_MARKERS); nb_umis <- marker_umis(NON_B_MARKERS)
  b_markers <- b_umis > 0
  non_b_markers <- nb_umis >= 2 & nb_umis > b_umis
  rm(b_umis, nb_umis)
  rm(ent_genes)

  # 1. Immunoglobulin composition -> Ig-score
  progress("ig", 12, "Immunoglobulin composition (24 metrics) and Ig-score")
  IG_MAP <- sc_symbol_map("ig", genes, ann, namespace)
  add_audit("IG_canonical_covered", sum(lengths(IG_MAP) > 0), paste0("/", length(IG_MAP)))
  Xig <- sc_ig_indices(counts, IG_MAP)
  for (mm in setdiff(assets$ig$features, colnames(Xig))) Xig[[mm]] <- NA_real_
  ig_score <- gbm_tumor_prob(assets$ig, Xig[, assets$ig$features, drop = FALSE])

  # 2. Malignant transcriptional programs (UCell) -> MTP score
  progress("mtp", 22, "Malignant transcriptional programs (UCell)")
  MTP_MAP <- sc_symbol_map("mtp", genes, ann, namespace)
  add_audit("MTP_canonical_covered", sum(lengths(MTP_MAP) > 0), paste0("/", length(MTP_MAP)))
  mtp_programs <- lapply(assets$mtp$programs, function(g) sc_resolve(MTP_MAP, g))
  empty <- names(mtp_programs)[lengths(mtp_programs) == 0]
  if (length(empty)) input_error("None of the genes of ", length(empty), " malignant transcriptional program(s) were found (", paste(gsub("_", " ", empty), collapse = ", "), "). Use human gene symbols (for example MS4A1, CD79A) or human Ensembl IDs (ENSG...); mouse or custom identifiers are not supported.")
  mtp_programs <- lapply(mtp_programs, function(g) g[g %in% genes])
  uc <- ucell_scores(counts, mtp_programs)
  uc <- as.data.frame(uc[, assets$mtp$features, drop = FALSE]); uc[is.na(uc)] <- 0
  mtp_score <- gbm_tumor_prob(assets$mtp, uc)

  # 3. Inferred copy-number variation (infercnvpy arithmetic against the fixed normal-B reference)
  progress("cnv", 35, "Inferred copy-number variation against the normal B-cell reference")
  cs <- sc_cnv_setup(genes, assets)
  cnv_score <- numeric(nc)
  for (b in split(seq_len(nc), ceiling(seq_len(nc) / cnv_batch))) {
    cnv_score[b] <- sc_cnv_batch(counts[, b, drop = FALSE], totals[b], cs, assets$cnv$params)
    invisible(gc(verbose = FALSE))
    progress("cnv", 35 + round(25 * max(b) / nc), sprintf("Copy-number inference: %s of %s cells", format(max(b), big.mark = ","), format(nc, big.mark = ",")))
  }
  cnv_qc <- list(manifest = data.frame(gene = cs$common, chromosome = as.character(cs$manifest$chromosome), start = cs$manifest$start, end = cs$manifest$end),
                 audit = cs$audit, used_genes = length(cs$common), smoothed_windows = sum(vapply(cs$pyramids, function(p) nrow(p$W), integer(1))))
  rm(cs); invisible(gc(verbose = FALSE))

  # 4. Log-normalization in place (Seurat LogNormalize: log1p(x / total * 1e4))
  progress("normalize", 62, "Log-normalizing counts")
  xv <- counts@x; counts@x <- numeric(0)
  p <- counts@p
  for (ch in chunks) {
    lo <- p[ch[1]] + 1L; hi <- p[ch[length(ch)] + 1L]
    if (hi >= lo) { seg <- lo:hi; xv[seg] <- log1p(xv[seg] / rep.int(totals[ch], diff(p[c(ch, ch[length(ch)] + 1L)])) * 1e4) }
    invisible(gc(verbose = FALSE))
  }
  counts@x <- xv; rm(xv); norm <- counts; rm(counts); invisible(gc(verbose = FALSE))
  row_means <- Matrix::rowMeans(norm)

  # 5. Proliferation: Seurat CellCycleScoring, Prol = max(S, G2M)
  progress("prol", 70, "Proliferation (S and G2/M module scores)")
  cc <- ann$cell_cycle_genes
  CC_MAP <- sc_symbol_map("cell_cycle", genes, ann, namespace)
  s.genes <- intersect(sc_resolve(CC_MAP, cc$s.genes), genes); g2m.genes <- intersect(sc_resolve(CC_MAP, cc$g2m.genes), genes)
  add_audit(c("S_genes_resolved", "G2M_genes_resolved"), c(length(s.genes), length(g2m.genes)), c(paste0("/", length(cc$s.genes)), paste0("/", length(cc$g2m.genes))))
  if (length(s.genes) < 3 || length(g2m.genes) < 3) input_error("Fewer than three S-phase or G2/M genes were found, so proliferation cannot be scored. Use human gene symbols (for example MS4A1, CD79A) or human Ensembl IDs (ENSG...); mouse or custom identifiers are not supported.")
  features <- list(S.Score = s.genes, G2M.Score = g2m.genes)
  ctrl.use <- sc_module_controls(row_means, features, ctrl = min(lengths(features)))
  S <- G2M <- numeric(nc)
  for (ch in chunks) {
    blk <- norm[, ch, drop = FALSE]
    S[ch] <- Matrix::colMeans(blk[features[[1]], , drop = FALSE]) - Matrix::colMeans(blk[ctrl.use[[1]], , drop = FALSE])
    G2M[ch] <- Matrix::colMeans(blk[features[[2]], , drop = FALSE]) - Matrix::colMeans(blk[ctrl.use[[2]], , drop = FALSE])
    rm(blk); invisible(gc(verbose = FALSE))
  }
  phase <- cc_phase(S, G2M); prol_score <- pmax(S, G2M)

  # 6. Differentiation potency: CCAT against the SCENT PPI network
  progress("diff", 80, "Differentiation potency (CCAT)")
  ent <- sc_map_entrez(genes, ann); ok <- which(!is.na(ent))
  add_audit("CCAT_network_genes_mapped", length(ok), paste0("/", length(genes)))
  keep <- ok
  if (anyDuplicated(ent[ok])) { ord <- order(row_means[ok], decreasing = TRUE); keep <- ok[ord]; keep <- keep[!duplicated(ent[keep])] }
  common <- intersect(rownames(assets$ppi), ent[keep])
  if (length(common) <= 100) input_error("Only ", length(common), " genes overlap the protein-interaction network used for differentiation potency; more than 100 are required. Use human gene symbols (for example MS4A1, CD79A) or human Ensembl IDs (ENSG...); mouse or custom identifiers are not supported.")
  rows_common <- keep[match(common, ent[keep])]
  k.v <- Matrix::rowSums(assets$ppi[common, common, drop = FALSE])
  Y <- matrix(k.v, ncol = 1)
  diff_score <- numeric(nc)
  for (ch in chunks) { diff_score[ch] <- as.numeric(cor_sparse_y(norm[rows_common, ch, drop = FALSE], Y)); invisible(gc(verbose = FALSE)) }
  add_audit("CCAT_network_genes_used", length(common), paste0("/", nrow(assets$ppi)))
  rm(norm); invisible(gc(verbose = FALSE))

  # 7. CMI from the five indices and its grade
  progress("cmi", 90, "Clonal malignancy index")
  idx <- data.frame(cell = cells, Ig_score = unname(ig_score), cnv_score = unname(cnv_score), mtp_score = unname(mtp_score),
                    diff_score = unname(diff_score), prol_score = unname(prol_score), check.names = FALSE, stringsAsFactors = FALSE)
  if (any(!is.finite(as.matrix(idx[, -1])))) stop("A non-finite index value was produced; CMI cannot be scored.")
  cmi <- cmi_predict(idx, assets$cmi)
  res <- cbind(idx, CMI = cmi, malignant = cmi >= CMI_THRESHOLD,
               state = as.character(assign_state(cmi)),
               S_score = S, G2M_score = G2M, phase = phase, n_counts = totals, b_markers_detected = b_markers,
               non_b_markers_dominant = non_b_markers, stringsAsFactors = FALSE)
  list(cells = res, ig_metrics = cbind(cell = cells, Xig), mtp_programs = cbind(cell = cells, uc),
       audit = audit, cnv = cnv_qc, namespace = namespace)
}
