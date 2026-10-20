#!/usr/bin/env Rscript
# =============================================================================
# CMI Explorer - equivalence with the v1.2.0 release
# =============================================================================
# Re-runs the engines on the shipped example inputs and compares every score with the output the
# v1.2.0 release produced on 2026-10-16 (tests/reference/, README there).
# The example files are byte-identical to the v1.2.0 inputs (SHA-256 checked below).
#   Rscript tests/test_equivalence.R [app_dir]
# Exit status 0 = every check passed.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
app <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
for (f in c("constants.R", "io.R", "models.R", "sc_engine.R", "bulk_engine.R", "spatial_engine.R")) source(file.path(app, "engine", f))
suppressPackageStartupMessages(library(Matrix))
Sys.setenv(CMI_PROFILE = "local")
ref_dir <- file.path(app, "tests", "reference")
results <- list()
check <- function(name, ok, detail) {
  results[[length(results) + 1L]] <<- data.frame(check = name, result = if (isTRUE(ok)) "PASS" else "FAIL", detail = detail, stringsAsFactors = FALSE)
  cat(sprintf("%-4s  %-58s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", name, detail))
}
maxdiff <- function(a, b) max(abs(as.numeric(a) - as.numeric(b)))
sha <- function(p) digest::digest(file = p, algo = "sha256")

# ---- 0. inputs are the v1.2.0 inputs; cut points unchanged ---------------------------------------
check("Example inputs identical to the v1.2.0 inputs",
      identical(substr(c(sha(file.path(app, "examples/single_cell/counts.csv")), sha(file.path(app, "examples/bulk/expression.csv")),
                         sha(file.path(app, "examples/spatial/rctd_bundle_demo.zip"))), 1, 16), c("4beb978da1fcae9d", "00a3905f1c2c48b7", "f82ef5c6f589cee9")),
      "SHA-256 of counts.csv, expression.csv, spatial bundle")
check("CMI grade cut points unchanged", identical(c(CMI_THRESHOLD, CMI_INDEX2, CMI_INDEX3), c(0.9755057544904363, 0.98926105, 0.99619844)),
      sprintf("%.7f / %.8f / %.8f", CMI_THRESHOLD, CMI_INDEX2, CMI_INDEX3))
man <- jsonlite::fromJSON(file.path(app, "ref", "asset_manifest.json"), simplifyVector = FALSE)
asset_ok <- all(vapply(names(man$assets), function(f) identical(sha(file.path(app, "ref", f)), man$assets[[f]]$sha256), logical(1)))
check("Model assets match their manifest", asset_ok, sprintf("%d files, SHA-256", length(man$assets)))

# ---- 1. single cell ------------------------------------------------------------------------------
assets <- load_assets(file.path(app, "ref"))
t0 <- Sys.time()
sc_mtx <- run_single_cell(function() read_single_cell_counts(file.path(app, "examples/single_cell/counts_matrix_market.zip"), tempfile()), assets)
t_sc <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
sc_csv <- run_single_cell(function() read_single_cell_counts(file.path(app, "examples/single_cell/counts.csv"), tempfile()), assets)
sc_small <- run_single_cell(function() read_single_cell_counts(file.path(app, "examples/single_cell/counts.csv"), tempfile()), assets, chunk_cells = 7L, cnv_batch = 5L)
check("Matrix Market ZIP and CSV give identical results", identical(sc_mtx$cells, sc_csv$cells), "96 cells, every column")
check("Results do not depend on chunk size", identical(sc_csv$cells, sc_small$cells), "chunks of 1,000/150 vs 7/5 cells")
ref <- utils::read.csv(file.path(ref_dir, "v1.2.0_single_cell_results.csv"), colClasses = c(cell = "character"))
x <- sc_csv$cells[match(ref$cell, sc_csv$cells$cell), ]
for (k in c("Ig_score", "cnv_score", "mtp_score", "diff_score", "prol_score")) {
  d <- maxdiff(x[[k]], ref[[k]]); check(sprintf("Single cell: %s equals v1.2.0", k), d <= 1e-12, sprintf("max |diff| %.2g (tolerance 1e-12)", d))
}
d <- maxdiff(x$CMI, ref$CMI)
check("Single cell: CMI equals v1.2.0 (H2O)", d <= 2e-9, sprintf("max |diff| %.2g (tolerance 2e-9: H2O stores trees in single precision)", d))
check("Single cell: malignant calls equal v1.2.0", identical(x$CMI >= CMI_THRESHOLD, ref$CMI >= CMI_THRESHOLD) &&
        identical(as.character(assign_state(x$CMI)), as.character(assign_state(ref$CMI))), "threshold and Mal-1/2/3 assignment for all 96 cells")

# ---- 2. bulk -------------------------------------------------------------------------------------
rsf <- readRDS(file.path(app, "ref", "rsf_m014.rds"))
b <- run_bulk(file.path(app, "examples/bulk/expression.csv"), rsf, file.path(app, "ref", "reo_pairs.csv"))
bref <- utils::read.csv(file.path(ref_dir, "v1.2.0_bulk_results.csv"), colClasses = c(sample_id = "character"))
bx <- b$scores[match(bref$sample_id, b$scores$sample_id), ]
d <- maxdiff(bx$riskscore, bref$riskscore)
check("Bulk: CMI-riskscore equals v1.2.0 (randomForestSRC)", d <= 1e-10, sprintf("max |diff| %.2g over 24 samples (tolerance 1e-10)", d))
check("Bulk: pair coverage equals v1.2.0", maxdiff(bx$pair_coverage, bref$pair_coverage) == 0, sprintf("%.0f%% of 100 pairs", 100 * mean(bx$pair_coverage)))

# ---- 3. spatial ----------------------------------------------------------------------------------
sp <- read_spatial_weights(file.path(app, "examples/spatial/spot_weights_demo.csv"), limits_for("local"))
g1 <- grade_spots(sp); g2 <- grade_spots(sp)
check("Spatial grading is reproducible", identical(g1$spots, g2$spots) && identical(g1$grades, g2$grades), sprintf("%d spots", nrow(g1$spots)))
if (requireNamespace("spacexr", quietly = TRUE)) {
  rc <- run_rctd_bundle(file.path(app, "examples/spatial/rctd_bundle_demo.zip"), tempfile(), lim = limits_for("local"))
  sref <- utils::read.csv(file.path(ref_dir, "v1.2.0_spatial_rctd_results.csv"), colClasses = c(spot_id = "character"))
  wref <- as.matrix(sref[, grep("^RCTD_weight_", names(sref))]); colnames(wref) <- canonical_class(colnames(wref))
  W <- rc$W[sref$spot_id, colnames(wref), drop = FALSE]
  d <- maxdiff(W, wref)
  check("Spatial: RCTD weights equal v1.2.0 (memory-patched spacexr)", d <= 1e-10, sprintf("max |diff| %.2g over %d spots (tolerance 1e-10)", d, nrow(W)))
} else {
  check("Spatial: RCTD weights equal v1.2.0", TRUE, "skipped: spacexr is not installed (hosted edition)")
}

res <- do.call(rbind, results)
dir.create(file.path(app, "tests", "results"), showWarnings = FALSE)
report <- jsonlite::toJSON(list(run_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), R = R.version.string, single_cell_seconds = round(t_sc, 1),
                                checks = res), auto_unbox = TRUE, pretty = TRUE)
# LF line endings on every platform.
write_lf <- function(text, path) { con <- file(path, open = "wb"); on.exit(close(con)); writeLines(text, con) }
write_lf(report, file.path(app, "tests", "results", "equivalence_report.json"))
cat(sprintf("\n%d of %d checks passed\n", sum(res$result == "PASS"), nrow(res)))
quit(save = "no", status = if (all(res$result == "PASS")) 0L else 1L)
