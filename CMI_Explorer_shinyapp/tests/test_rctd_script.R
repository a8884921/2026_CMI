#!/usr/bin/env Rscript
# =============================================================================
# CMI Explorer - checks of the downloadable RCTD script (www/downloads/rctd_mal_states.R)
# =============================================================================
# The script must use the cut points of the engines, read every input format identically, refuse
# ambiguous inputs, reproduce the RCTD weights of the original release (v1.2.0) on its test data,
# agree with the local edition when reference cells fall below the threshold, and write a table that
# the Spatial page accepts. Needs spacexr 2.2.1 for the RCTD runs (skipped without it).
#   Rscript tests/test_rctd_script.R [app_dir]
# Exit status 0 = every check passed.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
app <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
for (f in c("constants.R", "io.R", "spatial_engine.R")) source(file.path(app, "engine", f))
source(file.path(app, "R", "precheck.R"))
suppressPackageStartupMessages(library(Matrix))
results <- list()
check <- function(name, ok, detail) {
  results[[length(results) + 1L]] <<- data.frame(check = name, result = if (isTRUE(ok)) "PASS" else "FAIL", detail = detail, stringsAsFactors = FALSE)
  cat(sprintf("%-4s  %-62s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", name, detail))
}
script <- file.path(app, "www", "downloads", "rctd_mal_states.R")
s <- new.env(); sys.source(script, envir = s)          # defines the functions; main() runs only under Rscript
err <- function(expr) tryCatch({ force(expr); "" }, error = conditionMessage)

# ---- 1. cut points and state rule -------------------------------------------------------------------
check("Cut points equal those of the engines",
      identical(c(s$CMI_THRESHOLD, s$CMI_MAL1_MAL2, s$CMI_MAL2_MAL3), c(CMI_THRESHOLD, CMI_INDEX2, CMI_INDEX3)),
      sprintf("CMI %.16g / %.8f / %.8f", s$CMI_THRESHOLD, s$CMI_MAL1_MAL2, s$CMI_MAL2_MAL3))
grid <- c(0, CMI_THRESHOLD - 1e-12, CMI_THRESHOLD, CMI_INDEX2, CMI_INDEX2 + 1e-12, CMI_INDEX3, CMI_INDEX3 + 1e-12, 1, seq(0.95, 1, by = 1e-6))
check("States equal those of Single-cell CMI", identical(s$mal_state(grid), as.character(assign_state(grid))),
      sprintf("%s CMI values including every boundary", format(length(grid), big.mark = ",")))
lab <- data.frame(cell = paste0("c", 1:7), cell_type = c("Malignant", "Malignant", "Malignant", "Tumor cell", "CD8_T", "Normal_B", ""),
                  CMI = c(0.97, CMI_THRESHOLD, 0.9950, 0.9990, NA, NA, NA), stringsAsFactors = FALSE)
r1 <- s$reference_classes(lab); r2 <- s$reference_classes(lab, below_threshold = "Normal_B")
e_res <- c(err(s$reference_classes(transform(lab, cell_type = replace(cell_type, 5, "Mal-2")))),
           err(s$reference_classes(transform(lab, cell_type = replace(cell_type, 5, "Malig_3")))))
check("Reference classes from labels and CMI",
      identical(r1$class, c(NA, "Mal-1", "Mal-2", "Mal-3", "CD8_T", "Normal_B", NA)) && identical(r2$class[1], "Normal_B") && all(grepl("assigned by this script", e_res)) &&
        isTRUE(s$has_normal_b(c("Mal-1", "Normal_B"))) && !isTRUE(s$has_normal_b(c("Mal-1", "CD8_T", "Macrophage"))),
      "below threshold left out (or relabeled); empty labels unused; Mal-1/2/3 and Malig_3 refused; normal-B class found")

# ---- 2. input formats ---------------------------------------------------------------------------------
tmp <- tempfile("rctd_script_test_"); dir.create(tmp)
utils::unzip(file.path(app, "www", "downloads", "rctd_test_data.zip"), exdir = tmp)
td <- file.path(tmp, "rctd_test_data")
m10 <- s$read_counts(file.path(td, "reference_counts"))
zip::zip(file.path(tmp, "ref.zip"), files = list.files(file.path(td, "reference_counts"), full.names = TRUE), mode = "cherry-pick")
utils::write.csv(data.frame(gene = rownames(m10), as.matrix(m10), check.names = FALSE), file.path(tmp, "ref.csv"), row.names = FALSE)
saveRDS(m10, file.path(tmp, "ref.rds"))
same <- function(a, b) identical(dim(a), dim(b)) && identical(dimnames(a), dimnames(b)) && max(abs(a - b)) == 0
fmt_ok <- c(zip = same(m10, s$read_counts(file.path(tmp, "ref.zip"))), csv = same(m10, s$read_counts(file.path(tmp, "ref.csv"))),
            rds = same(m10, s$read_counts(file.path(tmp, "ref.rds"))))
seurat_detail <- "SeuratObject not installed"
if (requireNamespace("SeuratObject", quietly = TRUE)) {
  so <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = m10))   # Seurat writes "_" in gene names as "-"
  labs <- utils::read.csv(file.path(td, "reference_labels.csv"), colClasses = "character")
  so$cell_type <- factor(labs$cell_type[match(colnames(so), labs$cell)]); so$CMI <- suppressWarnings(as.numeric(labs$CMI[match(colnames(so), labs$cell)]))
  ms <- s$read_counts("seurat.rds", obj = so); ls <- s$read_labels(NULL, so)
  lc <- s$read_labels(file.path(td, "reference_labels.csv"))
  fmt_ok["seurat"] <- identical(dim(ms), dim(m10)) && max(abs(ms - m10)) == 0 && identical(gsub("_", "-", rownames(m10)), rownames(ms)) &&
                      identical(ls$cell_type, lc$cell_type[match(ls$cell, lc$cell)]) && isTRUE(all.equal(ls$CMI, lc$CMI[match(ls$cell, lc$cell)]))
  # raw counts from the RNA assay even when another assay is the default; split layers refused
  sct <- suppressWarnings(SeuratObject::CreateAssay5Object(counts = m10 * 2L)); so[["SCT"]] <- sct; SeuratObject::DefaultAssay(so) <- "SCT"
  rna_used <- identical(s$seurat_assay(so), "RNA") && max(abs(s$read_counts("x.rds", obj = so) - m10)) == 0
  sp <- so; SeuratObject::DefaultAssay(sp) <- "RNA"; sp[["SCT"]] <- NULL
  sp[["RNA"]] <- split(sp[["RNA"]], f = rep(c("S1", "S2"), length.out = ncol(sp)))
  split_refused <- grepl("JoinLayers", err(s$read_counts("x.rds", obj = sp)))
  seurat_detail <- sprintf("RNA assay used when SCT is the default: %s; split layers refused: %s", rna_used, split_refused)
  fmt_ok["seurat_assay"] <- rna_used && split_refused
}
check("Reference formats read identically", all(fmt_ok), sprintf("10x folder = %s; %s", paste(setdiff(names(fmt_ok), "seurat_assay"), collapse = ", "), seurat_detail))
outs <- file.path(td, "section", "outs")
xy2 <- s$read_coords(outs)
dir.create(file.path(tmp, "sr1", "spatial"), recursive = TRUE)
pos <- utils::read.csv(file.path(outs, "spatial", "tissue_positions.csv"), colClasses = "character")
utils::write.table(pos, file.path(tmp, "sr1", "spatial", "tissue_positions_list.csv"), sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
utils::write.csv(data.frame(spot_id = pos$barcode, x = pos$pxl_col_in_fullres, y = pos$pxl_row_in_fullres), file.path(tmp, "coords.csv"), row.names = FALSE)
check("Spot coordinates read identically", identical(xy2, s$read_coords(file.path(tmp, "sr1"))) && identical(xy2, s$read_coords(outs, file.path(tmp, "coords.csv"))),
      "tissue_positions.csv = tissue_positions_list.csv = --coords")

# ---- 3. inputs that must be refused before RCTD runs ------------------------------------------------
labs <- utils::read.csv(file.path(td, "reference_labels.csv"), colClasses = "character")
wcsv <- function(df) { f <- tempfile(fileext = ".csv", tmpdir = tmp); utils::write.csv(df, f, row.names = FALSE); f }
run_err <- function(...) err(s$cmi_rctd(spatial = outs, output = file.path(tmp, "never.csv"), ...))
dup <- m10; colnames(dup)[2] <- colnames(dup)[1]; saveRDS(dup, f_dup <- file.path(tmp, "dup.rds"))
part <- m10[, 1:200]; saveRDS(part, f_part <- file.path(tmp, "part.rds"))
no3 <- labs; no3$CMI[no3$cell_type == "Malignant" & as.numeric(no3$CMI) > CMI_INDEX3] <- "0.993"
msgs <- c(duplicated = run_err(ref_counts = f_dup, ref_labels = file.path(td, "reference_labels.csv")),
          missing = run_err(ref_counts = f_part, ref_labels = file.path(td, "reference_labels.csv")),
          no_mal3 = run_err(ref_counts = file.path(td, "reference_counts"), ref_labels = wcsv(no3)))
check("Ambiguous references are refused", grepl("not unique", msgs[["duplicated"]]) && grepl("JoinLayers", msgs[["missing"]]) && grepl("no Mal-3", msgs[["no_mal3"]]),
      "duplicated barcodes; labeled cells missing from the counts (split Seurat layers); no Mal-3 class")

# ---- 4. RCTD runs, as a user runs them --------------------------------------------------------------
rscript <- function(out, labels, extra = character()) {
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  st <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(script), "--ref-counts", shQuote(file.path(td, "reference_counts")),
                "--ref-labels", shQuote(labels), "--spatial", shQuote(outs), "--output", shQuote(out), extra),
                stdout = paste0(out, ".run.log"), stderr = paste0(out, ".run.log"))
  if (st != 0 || !file.exists(out)) { cat(utils::tail(readLines(paste0(out, ".run.log")), 5), sep = "\n"); return(NULL) }
  utils::read.csv(out, check.names = FALSE)
}
wdiff <- function(a, e) if (is.null(a) || !setequal(a$spot_id, e$spot_id)) Inf else max(abs(as.matrix(a[, -(1:3)]) - as.matrix(e[match(a$spot_id, e$spot_id), names(a)[-(1:3)]])))
if (requireNamespace("spacexr", quietly = TRUE)) {
  e <- utils::read.csv(file.path(td, "expected_weights.csv"), check.names = FALSE)
  t0 <- Sys.time(); a <- rscript(file.path(tmp, "out", "test_weights.csv"), file.path(td, "reference_labels.csv"))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  check("RCTD weights equal v1.2.0", wdiff(a, e) <= 1e-10,
        sprintf("max |diff| %.2g over %d spots (tolerance 1e-10; %.0f s, spacexr %s)", wdiff(a, e), if (is.null(a)) 0L else nrow(a), secs, utils::packageVersion("spacexr")))
  out <- file.path(tmp, "out", "test_weights.csv")
  pc <- precheck_spatial(out, basename(out), limits_for("hosted"))
  sp <- tryCatch(grade_spots(read_spatial_weights(out, limits_for("hosted"))), error = function(e) NULL)
  check("Output accepted by Spatial CMI grades", isTRUE(pc$ok) && !is.null(sp) && all(c("Mal-1", "Mal-2", "Mal-3") %in% sp$spots$grade),
        if (isTRUE(pc$ok)) pc$summary else paste(pc$problems, collapse = " "))

  # 20 malignant cells moved below the threshold: left out by the script and by the local edition alike;
  # with --below-threshold Mal-1 they are counted as Mal-1, as in v1.2.0, which gives the v1.2.0 weights.
  low <- labs$cell[labs$cell_type == "Malignant" & abs(as.numeric(labs$CMI) - 0.985) < 1e-9][1:20]
  lab_low <- labs; lab_low$CMI[lab_low$cell %in% low] <- "0.96"; f_low <- wcsv(lab_low)
  a_ex <- rscript(file.path(tmp, "out", "low_exclude.csv"), f_low)
  a_m1 <- rscript(file.path(tmp, "out", "low_mal1.csv"), f_low, c("--below-threshold", "Mal-1"))
  bdir <- file.path(tmp, "bundle"); utils::unzip(file.path(app, "examples", "spatial", "rctd_bundle_demo.zip"), exdir = bdir)
  bf <- function(n) list.files(bdir, pattern = paste0("^", n, "$"), recursive = TRUE, full.names = TRUE)[1]
  md <- utils::read.csv(bf("reference_metadata.csv"), colClasses = "character"); md$CMI[md$cell %in% low] <- "0.96"
  utils::write.csv(md, bf("reference_metadata.csv"), row.names = FALSE, quote = FALSE)
  zip::zip(file.path(tmp, "bundle_low.zip"), files = list.files(dirname(bf("manifest.csv")), full.names = TRUE), mode = "cherry-pick")
  loc <- tryCatch(run_rctd_bundle(file.path(tmp, "bundle_low.zip"), file.path(tmp, "work"), lim = limits_for("local")), error = function(e) NULL)
  e_loc <- if (is.null(loc)) NULL else data.frame(spot_id = loc$spots, x = 0, y = 0, loc$W, check.names = FALSE)
  d_ex <- if (is.null(e_loc)) Inf else wdiff(a_ex, e_loc)
  check("Cells below the threshold: script equals the local edition", d_ex <= 1e-10 && !is.null(loc) && nrow(loc$left_out) == 20,
        sprintf("20 Mal-1 cells at CMI 0.96 left out by both; max |diff| %.2g", d_ex))
  check("--below-threshold Mal-1 reproduces version 1.2.0", wdiff(a_m1, e) <= 1e-10,
        sprintf("the same 20 cells counted as Mal-1; max |diff| from the v1.2.0 weights %.2g", wdiff(a_m1, e)))
} else {
  check("RCTD weights equal v1.2.0", TRUE, "skipped: spacexr is not installed")
}

res <- do.call(rbind, results)
dir.create(file.path(app, "tests", "results"), showWarnings = FALSE)
writeLines(jsonlite::toJSON(list(run_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), R = R.version.string, checks = res), auto_unbox = TRUE, pretty = TRUE),
           file.path(app, "tests", "results", "rctd_script_report.json"))
cat(sprintf("\n%d of %d checks passed\n", sum(res$result == "PASS"), nrow(res)))
unlink(tmp, recursive = TRUE)
quit(status = if (all(res$result == "PASS")) 0 else 1)
