#!/usr/bin/env Rscript
# =============================================================================
# CMI Explorer analysis worker: one analysis per process, started by the interface.
#   Rscript engine/worker.R --root <app> --module single_cell|bulk|spatial --input <file>
#           --output <dir> --state <state.json> [--spatial-mode weights|rctd] [--threads 1]
#           [--seed 20260101] [--input-name <uploaded file name>] [--synthetic]
# Analyses take an instance-wide slot so that only CMI_MAX_CONCURRENT_JOBS (hosted: 1) run at once;
# later requests wait first-come, first-served and report their place in the queue.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
arg <- function(key, default = NULL) { i <- match(paste0("--", key), args); if (is.na(i) || i == length(args)) default else args[i + 1L] }
root <- normalizePath(arg("root"), mustWork = TRUE)
module <- arg("module"); input <- arg("input"); out_dir <- arg("output"); state_path <- arg("state")
spatial_mode <- arg("spatial-mode", "weights")
threads <- as.integer(arg("threads", "1")); seed <- as.integer(arg("seed", "20260101")); synthetic <- "--synthetic" %in% args
stopifnot(module %in% c("single_cell", "bulk", "spatial"), spatial_mode %in% c("weights", "rctd"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# Absolute paths: relative ones would break once a library changes the working directory (zip::zip).
out_dir <- normalizePath(out_dir, mustWork = TRUE)
dir.create(dirname(state_path), recursive = TRUE, showWarnings = FALSE)
state_path <- file.path(normalizePath(dirname(state_path), mustWork = TRUE), basename(state_path))
input <- normalizePath(input, mustWork = TRUE)
input_name <- arg("input-name", basename(input))         # the file name the user uploaded (the server stages a copy)
for (d in c("figures", "qc")) dir.create(file.path(out_dir, d), showWarnings = FALSE)
started <- Sys.time()
log_path <- file.path(out_dir, "run.log")
say <- function(...) { msg <- paste0(format(Sys.time(), "%H:%M:%S"), "  ", sprintf(...)); cat(msg, "\n", sep = ""); cat(msg, "\n", sep = "", file = log_path, append = TRUE) }

for (f in c("constants.R", "io.R", "models.R")) source(file.path(root, "engine", f))
# Matrix (about 80 MB) is attached only where sparse matrices are used: single-cell scoring and RCTD.
if (module == "single_cell" || (module == "spatial" && spatial_mode == "rctd")) suppressPackageStartupMessages(library(Matrix))
source(file.path(root, "R", "figures.R")); ensure_utf8()
lim <- limits_for()

write_state <- function(status, progress, detail, stage = status, extra = list()) {
  payload <- c(list(status = status, stage = stage, progress = progress, detail = detail, module = module,
                    started_at = format(started, "%Y-%m-%dT%H:%M:%S%z"), updated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                    pid = Sys.getpid()), extra)
  tmp <- paste0(state_path, ".tmp")
  writeLines(jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null", na = "null", digits = NA), tmp)
  file.rename(tmp, state_path)
  heartbeat()
}

# ---- instance-wide first-come, first-served queue ------------------------------------------------
lock_root <- Sys.getenv("CMI_LOCK_DIR", file.path(dirname(tempdir()), "cmi_locks"))
queue_dir <- file.path(lock_root, "queue")
slots <- max(1L, as.integer(Sys.getenv("CMI_MAX_CONCURRENT_JOBS", if (identical(lim$profile, "local")) "2" else "1")))
dir.create(queue_dir, recursive = TRUE, showWarnings = FALSE)
pid_alive <- function(pid) {
  pid <- suppressWarnings(as.integer(pid)); if (is.na(pid)) return(FALSE)
  if (pid == Sys.getpid()) return(TRUE)
  isTRUE(tryCatch(ps::ps_is_running(ps::ps_handle(pid)), error = function(e) FALSE))
}
ticket <- file.path(queue_dir, sprintf("%020.0f_%d", as.numeric(Sys.time()) * 1e6, Sys.getpid()))
held <- NULL
heartbeat <- function() if (!is.null(held)) writeLines(c(as.character(Sys.getpid()), format(as.numeric(Sys.time()), digits = 15)), file.path(held, "owner"))
slot_free <- function(d) {
  if (!dir.exists(d)) return(TRUE)
  own <- tryCatch(readLines(file.path(d, "owner"), warn = FALSE), error = function(e) character())
  if (!length(own)) return((as.numeric(Sys.time()) - as.numeric(file.info(d)$mtime)) > 60)    # a slot being created
  !pid_alive(own[1])                                                                          # its owner has exited
}
acquire_slot <- function() {
  file.create(ticket); waited <- 0; last <- NA
  repeat {
    tk <- sort(list.files(queue_dir))
    alive <- vapply(tk, function(t) pid_alive(sub("^.*_", "", t)), logical(1))
    for (t in tk[!alive]) unlink(file.path(queue_dir, t))
    tk <- tk[alive]
    pos <- match(basename(ticket), tk)
    if (is.na(pos)) { file.create(ticket); next }
    if (pos <= slots) for (k in seq_len(slots)) {
      d <- file.path(lock_root, sprintf("slot%02d", k))
      if (dir.exists(d) && slot_free(d)) unlink(d, recursive = TRUE)
      if (dir.create(d, showWarnings = FALSE)) { held <<- d; heartbeat(); unlink(ticket); return(invisible(TRUE)) }
    }
    if (!identical(pos, last) || waited %% 15 == 0) {
      write_state("queued", 0, if (pos <= slots) "Another analysis is running on this server; yours starts as soon as it finishes."
                               else sprintf("Waiting in line: %d analys%s ahead of yours on this server.", pos - 1L, if (pos == 2L) "is" else "es"))
      last <- pos
    }
    Sys.sleep(1); waited <- waited + 1
    if (waited > 900) { unlink(ticket); stop("The analysis queue did not clear within 15 minutes. Please try again later.") }
  }
}
release_slot <- function() { if (!is.null(held)) { unlink(held, recursive = TRUE); held <<- NULL }; unlink(ticket) }

progress <- function(stage, pct, detail) { say("[%s] %s", stage, detail); write_state("running", pct, detail, stage) }

csv <- function(x, name) utils::write.csv(x, file.path(out_dir, name), row.names = FALSE, na = "")
json <- function(x, name) writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", na = "null", digits = NA, pretty = TRUE), file.path(out_dir, name))
legends <- character()
fig_plan <- 0L; fig_done <- 0L                       # figure sub-steps for the progress bar
plan_figures <- function(n) { fig_plan <<- as.integer(n); fig_done <<- 0L }
figure <- function(p, name, width = NULL, height = NULL) {
  if (is.null(p)) return(invisible())
  fig_done <<- fig_done + 1L
  if (fig_plan > 0L) write_state("running", min(99, 94 + floor(5 * (fig_done - 1L) / fig_plan)), sprintf("Writing figure %d of %d (SVG, PDF and PNG)", fig_done, fig_plan), "export")
  sz <- FIG_SIZES[[name]] %||% c(7.2, 4.6); w <- width %||% sz[1]; h <- height %||% sz[2]
  r <- tryCatch(save_figure_set(p, file.path(out_dir, "figures", name), w, h), error = function(e) { say("Figure %s skipped: %s", name, conditionMessage(e)); NULL })
  if (!is.null(r)) {
    lg <- attr(r, "legend")
    legends <<- c(legends, paste0(name, "\n  ", paste(c(lg$title, lg$subtitle), collapse = ". "), if (length(lg$caption)) paste0("\n  ", lg$caption) else "", "\n"))
  }
  invisible(r)
}

provenance <- function(extra) {
  pkgs <- c("Matrix", "data.table", "jsonlite", "ggplot2", "svglite", "ragg", "systemfonts", "HGNChelper", "matrixStats", "RANN", "spacexr")
  ver <- vapply(pkgs, function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_), character(1))
  man <- jsonlite::fromJSON(file.path(root, "ref", "asset_manifest.json"), simplifyVector = FALSE)
  c(list(app = "CMI Explorer", app_version = CMI_VERSION, app_release_date = CMI_RELEASE_DATE, module = module, synthetic_input = synthetic,
         input_file = input_name, input_sha256 = digest::digest(file = input, algo = "sha256"),
         parameters = list(threads = threads, seed = seed, spatial_mode = spatial_mode),
         profile = lim$profile, limits = lim[setdiff(names(lim), "profile")],
         started_at = format(started, "%Y-%m-%dT%H:%M:%S%z"), completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
         R = R.version.string, packages = as.list(ver[!is.na(ver)]),
         cmi_grades = list(malignancy_threshold = CMI_THRESHOLD, mal1_mal2_boundary = CMI_INDEX2, mal2_mal3_boundary = CMI_INDEX3),
         model_assets = man$assets), extra)
}

run <- function() {
  acquire_slot()
  on.exit(release_slot(), add = TRUE)
  write_state("running", 2, "Starting the analysis")
  register_figure_fonts(file.path(root, "engine", "fonts"))
  if (module == "single_cell") {
    source(file.path(root, "engine", "sc_engine.R"))
    progress("assets", 4, "Loading the models and references")
    assets <- load_assets(file.path(root, "ref"))
    dims <- NULL
    res <- run_single_cell(function() {
      progress("read", 6, "Reading and checking the count matrix")
      m <- read_single_cell_counts(input, file.path(out_dir, "work"), lim); dims <<- dim(m)
      say("Input: %d genes x %d cells, %s non-zero entries", nrow(m), ncol(m), format(length(m@x), big.mark = ","))
      m }, assets, progress = progress)
    unlink(file.path(out_dir, "work"), recursive = TRUE); rm(assets); invisible(gc(verbose = FALSE))
    cells <- res$cells[, c("cell", "CMI", "state", "malignant", "Ig_score", "cnv_score", "mtp_score", "diff_score", "prol_score",
                           "n_counts", "b_markers_detected", "non_b_markers_dominant")]
    progress("export", 94, "Writing tables, figures and the result bundle")
    csv(cells, "cells.csv"); csv(res$audit, "qc/gene_mapping_audit.csv"); csv(res$ig_metrics, "qc/ig_composition_metrics.csv")
    csv(res$mtp_programs, "qc/mtp_program_scores.csv"); csv(res$cnv$manifest, "qc/cnv_gene_manifest.csv"); csv(res$cnv$audit, "qc/cnv_gene_mapping_audit.csv")
    sc_summary <- list(cells = nrow(cells), genes = dims[1], malignant = sum(cells$malignant), malignant_pct = 100 * mean(cells$malignant),
                       state_counts = as.list(table(factor(cells$state, levels = STATE_LEVELS))), median_cmi = stats::median(cells$CMI),
                       lineage = list(no_b_marker_pct = 100 * mean(!cells$b_markers_detected), non_b_marker_pct = 100 * mean(cells$non_b_markers_dominant)),
                       gene_namespace = res$namespace, cnv_genes_used = res$cnv$used_genes, cnv_windows = res$cnv$smoothed_windows)
    json(sc_summary, "summary.json")
    plan_figures(2)
    figure(fig_sc_cmi(cells, synthetic), "fig1_cmi"); figure(fig_sc_grades(cells, synthetic), "fig2_grades")
    json(provenance(list(input_dimensions = list(genes = dims[1], cells = dims[2]))), "provenance.json")
    rows <- nrow(cells)
  } else if (module == "bulk") {
    source(file.path(root, "engine", "bulk_engine.R"))
    rsf <- readRDS(file.path(root, "ref", "rsf_m014.rds"))
    res <- run_bulk(input, rsf, file.path(root, "ref", "reo_pairs.csv"), lim = lim, progress = progress)
    progress("export", 94, "Writing tables, figures and the result bundle")
    csv(res$scores, "samples.csv"); csv(res$qc, "qc/dataset_qc.csv"); csv(res$mapping, "qc/reo_endpoint_mapping.csv")
    csv(res$reo_audit, "qc/hgnc_reo_audit.csv"); csv(res$expression_audit, "qc/hgnc_expression_audit.csv")
    if (nrow(res$missing_pairs)) csv(res$missing_pairs, "qc/missing_reo_pairs.csv")
    feat <- data.frame(sample_id = rownames(res$features), res$features, check.names = FALSE); csv(feat, "qc/reo_features.csv")
    b_summary <- list(samples = nrow(res$scores), genes_in_rank_pool = res$qc$expression_genes, pair_coverage = res$qc$pair_coverage,
                      coverage_class = res$qc$coverage_class, endpoint_coverage = res$qc$endpoint_coverage, median_riskscore = res$median,
                      groups = as.list(table(res$scores$cmi_group)))
    json(b_summary, "summary.json")
    plan_figures(1)
    figure(fig_bulk_ranked(res$scores, synthetic), "fig1_riskscore")
    json(provenance(list(input_dimensions = list(samples = nrow(res$scores), genes = res$qc$expression_genes))), "provenance.json")
    rows <- nrow(res$scores)
  } else {
    source(file.path(root, "engine", "spatial_engine.R"))
    if (spatial_mode == "rctd") {
      if (!identical(lim$profile, "local") && !identical(Sys.getenv("CMI_ENABLE_RCTD"), "true")) input_error("RCTD deconvolution is available in the local edition only (it needs several GB of memory). Run RCTD locally and upload the weights table.")
      sp <- run_rctd_bundle(input, file.path(out_dir, "work"), threads, seed, lim, progress)
      csv(sp$reference_classes, "qc/reference_classes.csv")
      if (nrow(sp$left_out)) csv(sp$left_out, "qc/reference_left_out_below_threshold.csv")
      csv(data.frame(spot_id = sp$spots, sp$raw_weights, check.names = FALSE), "qc/rctd_weights.csv")
      unlink(file.path(out_dir, "work"), recursive = TRUE)
    } else {
      progress("read", 8, "Reading the spot table")
      sp <- read_spatial_weights(input, lim)
    }
    progress("grade", 60, "Grading each spot")
    res <- grade_spots(sp)
    progress("export", 94, "Writing tables, figures and the result bundle")
    csv(res$spots, "spots.csv"); csv(res$grades, "grades.csv")
    json(list(spots = nrow(res$spots), classes = res$classes, grades = res$grades), "summary.json")
    plan_figures(1)
    figure(fig_sp_map(res, synthetic = synthetic), "fig1_grade_map")
    json(provenance(list(input_dimensions = list(spots = nrow(res$spots), classes = length(res$classes)))), "provenance.json")
    rows <- nrow(res$spots)
  }
  writeLines(c("Figure legends (title, subtitle and caption of each figure; the files in figures/panels/ omit them)", "", legends),
             file.path(out_dir, "figures", "legends.txt"))
  writeLines(readme_text(module), file.path(out_dir, "README.txt"))
  write_state("running", 99, "Packing the result bundle", "export")
  members <- list.files(out_dir, recursive = TRUE)
  members <- members[!grepl("^(result_bundle\\.zip|state\\.json.*)$|_web\\.png$", members)]
  zip::zip(file.path(out_dir, "result_bundle.zip"), files = members, root = out_dir, mode = "mirror")
  secs <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  write_state("success", 100, sprintf("Done: %s %s analyzed in %s.", format(rows, big.mark = ","),
                                      switch(module, single_cell = "cells", bulk = "samples", spatial = "spots"),
                                      if (secs < 90) sprintf("%.0f s", secs) else sprintf("%.1f min", secs / 60)), extra = list(rows = rows, synthetic = synthetic))
}

readme_text <- function(module) {
  common <- c("CMI Explorer result bundle", strrep("=", 26), "",
              sprintf("Module: %s   Version: %s   Created: %s", module, CMI_VERSION, format(Sys.time(), "%Y-%m-%d %H:%M %Z")), "",
              "provenance.json   input checksum, parameters, package versions and the SHA-256 of every model asset used",
              "summary.json      the numbers shown in the interface", "run.log           the steps the worker took",
              "figures/          each figure as SVG (editable text), PDF (vector, embedded fonts) and PNG (600 dpi)",
              "figures/panels/   the same figures without title, subtitle and caption, for assembling your own figures",
              "figures/legends.txt  title, subtitle and caption text of every figure", "qc/               audit and intermediate tables", "")
  body <- switch(module,
    single_cell = c("cells.csv columns", "  CMI          clonal malignancy index of the cell (0-1; a model score, not a probability)",
                    "  state        CMI grade: Below threshold (CMI < 0.97551), Mal-1 (to 0.98926), Mal-2 (to 0.99620) or Mal-3 (above)",
                    "  malignant    CMI at or above the malignancy threshold (0.97551)",
                    "  Ig_score, cnv_score, mtp_score, diff_score, prol_score  the five inputs of the CMI",
                    "  n_counts     UMIs per cell",
                    "  b_markers_detected  any UMI of CD19, MS4A1, CD79A, CD79B or PAX5 (upload check only)",
                    "  non_b_markers_dominant  2 or more UMIs of CD3E, CD3D, LYZ, CD14, NKG7, HBB or PPBP, and more than of the B-cell markers (upload check only)"),
    bulk = c("samples.csv columns", "  riskscore    CMI-riskscore of the sample (higher = higher predicted risk)",
             "  cmi_group    CMI-high above the median of this cohort, otherwise CMI-low",
             "  pair_coverage  share of the 100 REO gene pairs computable from this matrix"),
    spatial = c("spots.csv: CMI grade (Mal-1, Mal-2, Mal-3 or Non-malignant), class of highest weight, malignant load and normalized weights per spot",
                "grades.csv: spots per grade, share of spots and share of the summed Mal-1/2/3 weight"))
  c(common, body, "", "Research use only. Not for diagnosis or treatment decisions.")
}

tryCatch(run(), error = function(e) {
  msg <- conditionMessage(e)
  user <- inherits(e, "input_error")
  say("%s: %s", if (user) "Input problem" else "Error", msg)
  write_state("failed", 0, if (user) msg else paste("The analysis stopped because of an internal error:", msg), extra = list(input_problem = user))
  release_slot()
  quit(save = "no", status = 1L)
})
quit(save = "no", status = 0L)
