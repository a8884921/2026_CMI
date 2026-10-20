#!/usr/bin/env Rscript
# Builds the demonstration inputs shipped in examples/. All values are simulated and labeled as such.
#   examples/spatial/spot_weights_demo.csv   simulated Visium-like section (hexagonal grid, 12 classes)
#   examples/single_cell/simulated_b_cells.zip  900 simulated B cells (normal profile plus a graded malignant program)
args <- commandArgs(trailingOnly = TRUE)
app <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
set.seed(20260101)

# ---- spatial: hexagonal grid ---------------------------------------------------------------
nr <- 56; nc <- 56
grid <- expand.grid(col = seq_len(nc), row = seq_len(nr))
grid$x <- (grid$col + ifelse(grid$row %% 2 == 0, 0.5, 0)) * 100
grid$y <- grid$row * 100 * sqrt(3) / 2
cx <- function(x, y, x0, y0, rx, ry) ((x - x0) / rx)^2 + ((y - y0) / ry)^2
r1 <- cx(grid$x, grid$y, 2300, 2500, 1500, 1350)       # tumor nodule 1
r2 <- cx(grid$x, grid$y, 4300, 3400, 900, 800)         # tumor nodule 2
tumor <- pmin(r1, r2)
classes <- c("Mal-1", "Mal-2", "Mal-3", "CD4_T", "CD8_T", "Dendritic", "Endothelial", "Macrophage", "Mast", "NK", "Normal_B", "Stromal")
alpha <- matrix(0.35, nrow(grid), length(classes), dimnames = list(NULL, classes))
core <- cx(grid$x, grid$y, 2050, 2350, 820, 700) + 0.35 * sin(grid$x / 260)
inside <- tumor < 1
alpha[inside, "Mal-1"] <- 2.0; alpha[inside, "Mal-2"] <- 1.6
alpha[inside & core < 1, "Mal-3"] <- 3.2; alpha[inside & core < 1, "Mal-2"] <- 1.8; alpha[inside & core < 1, "Mal-1"] <- 0.8
alpha[inside & core >= 1 & core < 1.8, "Mal-2"] <- 3.5
alpha[r2 < 1, "Mal-1"] <- 3.2; alpha[r2 < 1, "Mal-2"] <- 1.0
alpha[inside, "Macrophage"] <- 0.9; alpha[inside & core < 1, "Macrophage"] <- 1.4
margin <- tumor >= 0.75 & tumor < 1.6
alpha[margin, "CD8_T"] <- 1.9; alpha[margin, "CD4_T"] <- 1.3; alpha[margin, "Dendritic"] <- 0.5
alpha[inside & core >= 1, "CD8_T"] <- 1.2; alpha[inside & core < 1, "CD8_T"] <- 0.6
outside <- tumor >= 1
alpha[outside, "Normal_B"] <- 2.2; alpha[outside, "Stromal"] <- 2.0; alpha[outside, "Endothelial"] <- 1.0
alpha[outside, "NK"] <- 0.5; alpha[outside, "Mast"] <- 0.3; alpha[outside & tumor >= 1.6, "CD8_T"] <- 0.8
W <- matrix(stats::rgamma(length(alpha), shape = alpha), nrow(alpha), dimnames = dimnames(alpha))
W <- W / rowSums(W)
keep <- cx(grid$x, grid$y, 2900, 2500, 2600, 2300) < 1          # round tissue outline
spots <- data.frame(spot_id = sprintf("SPOT_%05d", seq_len(sum(keep))), sample_id = "SIMULATED_SECTION",
                    x = round(grid$x[keep], 1), y = round(grid$y[keep], 1), round(W[keep, ], 5), check.names = FALSE)
dir.create(file.path(app, "examples", "spatial"), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(spots, file.path(app, "examples", "spatial", "spot_weights_demo.csv"), row.names = FALSE)
cat("spatial demo:", nrow(spots), "spots\n")

# ---- single cell: simulated B cells with a graded malignant program -------------------------------
# 900 cells drawn from the normal B-cell reference profile (ref/cnv_reference.rds) with Poisson-gamma
# noise. For 70% of the cells a malignant program of random strength t (0.1 to 1) is added: kappa
# light-chain restriction with loss of IgD and more J chain, copy-number changes (gains of 1q, 3q and
# 18q, loss of 6q), higher malignant-program genes, more cycling cells and lower MHC class I and II
# genes. The result exercises every output of the Single-cell workspace; it is not biology.
# (The v1.2.0 test fixture, counts_matrix_market.zip / counts.csv, stays for tests/test_equivalence.R.)
local({
  suppressPackageStartupMessages(library(Matrix))
  set.seed(11)
  n_cells <- 900; k_ig <- 6; k_cnv <- 1.5; k_mtp <- 5; k_cyc <- 0.6
  r <- readRDS(file.path(app, "ref", "cnv_reference.rds")); genes <- r$ref_genes; G <- length(genes)
  base <- expm1(r$ref_mean); base <- base / sum(base)
  idx <- function(g) match(intersect(g, genes), genes)
  mtp <- readRDS(file.path(app, "ref", "mtp_gbm.rds"))$programs
  ann <- readRDS(file.path(app, "ref", "annotation_snapshot.rds"))
  s_genes <- idx(ann$cell_cycle_genes$s.genes); g2m_genes <- idx(ann$cell_cycle_genes$g2m.genes)
  pos <- r$positions
  seg <- function(chr, from, to) { p <- pos[pos$chromosome %in% chr & pos$start >= from & pos$start <= to, , drop = FALSE]; idx(rownames(p)) }
  gain1q <- seg("chr1", 1.5e8, 2.5e8); loss6q <- seg("chr6", 9e7, 1.7e8); gain18 <- seg("chr18", 4e7, 8e7); gain3 <- seg("chr3", 1e8, 2e8)
  mhc1 <- idx(c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2")); mhc2 <- idx(c("HLA-DRA", "HLA-DRB1", "CD74", "CIITA"))
  igk <- idx("IGKC"); igl <- idx(c("IGLC1", "IGLC2", "IGLC3", "IGLC7")); ighd <- idx("IGHD"); jch <- idx("JCHAIN"); ighm <- idx("IGHM")
  mtp_idx <- idx(unlist(mtp))
  n0 <- round(0.3 * n_cells)
  t <- c(rep(0, n0), sort(stats::runif(n_cells - n0, 0.1, 1)))
  L <- round(exp(stats::rnorm(n_cells, log(6000), 0.3)))
  cols <- vector("list", n_cells)
  for (i in seq_len(n_cells)) {
    f <- rep(1, G); ti <- t[i]
    f[igk] <- f[igk] * (1 + k_ig * ti); f[igl] <- f[igl] * (1 - 0.97 * ti); f[ighd] <- f[ighd] * (1 - 0.9 * ti); f[jch] <- f[jch] * (1 + 4 * k_ig * ti); f[ighm] <- f[ighm] * (1 + 0.5 * ti)
    f[gain1q] <- f[gain1q] * (1 + k_cnv * ti); f[loss6q] <- f[loss6q] * (1 - 0.5 * ti); f[gain18] <- f[gain18] * (1 + k_cnv * ti); f[gain3] <- f[gain3] * (1 + k_cnv * ti)
    f[mtp_idx] <- f[mtp_idx] * (1 + k_mtp * ti)
    f[mhc1] <- f[mhc1] * (1 - 0.6 * ti); f[mhc2] <- f[mhc2] * (1 - 0.85 * ti)
    if (stats::runif(1) < 0.04 + k_cyc * ti) { if (stats::runif(1) < 0.5) f[s_genes] <- f[s_genes] * 5 else f[g2m_genes] <- f[g2m_genes] * 5 }
    p <- base * f * stats::rgamma(G, shape = 2.5, rate = 2.5); p <- p / sum(p)
    x <- stats::rpois(G, L[i] * p); nz <- which(x > 0)
    cols[[i]] <- list(i = nz, x = x[nz])
  }
  m <- sparseMatrix(i = unlist(lapply(cols, `[[`, "i")), j = rep(seq_len(n_cells), vapply(cols, function(c) length(c$i), 1L)),
                    x = unlist(lapply(cols, `[[`, "x")), dims = c(G, n_cells), dimnames = list(genes, sprintf("SIMCELL_%04d", seq_len(n_cells))))
  tmp <- tempfile("simsc_"); dir.create(tmp)
  writeMM(m, file.path(tmp, "matrix.mtx"))
  gz <- function(lines, f) { con <- gzfile(file.path(tmp, f), "w", compression = 9); writeLines(lines, con); close(con) }
  gz(readLines(file.path(tmp, "matrix.mtx")), "matrix.mtx.gz"); unlink(file.path(tmp, "matrix.mtx"))
  gz(paste(genes, genes, "Gene Expression", sep = "\t"), "features.tsv.gz"); gz(colnames(m), "barcodes.tsv.gz")
  zipf <- file.path(app, "examples", "single_cell", "simulated_b_cells.zip"); unlink(zipf)
  zip::zip(zipf, files = c("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"), root = tmp)
  unlink(tmp, recursive = TRUE)
  cat(sprintf("single-cell demo: %d cells, %s non-zero counts, %.1f MB\n", n_cells, format(length(m@x), big.mark = ","), file.size(zipf) / 1e6))
})
