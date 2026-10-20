# =============================================================================
# Analysis workspaces (single-cell, spatial, bulk): one Shiny module, three configurations.
# Uploads are checked in the interface process as soon as they arrive (R/precheck.R); analyses run
# in a separate R process (R/jobs.R, engine/worker.R); results are read back from its output folder.
# =============================================================================

fmt_n <- function(x, digits = 0) if (is.null(x) || !length(x) || is.na(x[1])) "\u2013" else formatC(x, format = "f", digits = digits, big.mark = ",")
fmt_pct <- function(x, digits = 1) if (is.null(x) || !length(x) || is.na(x[1])) "\u2013" else paste0(formatC(x, format = "f", digits = digits), "%")
fmt_signed <- function(x, digits = 1) { v <- round(x, digits); v[!is.na(v) & v == 0] <- 0; out <- gsub("-", "\u2212", formatC(v, format = "f", digits = digits), fixed = TRUE); out[is.na(x)] <- "\u2013"; out }
fmt_num <- function(x, digits = 2) if (is.null(x) || !length(x) || is.na(x[1])) "\u2013" else gsub("-", "\u2212", formatC(x, format = "f", digits = digits), fixed = TRUE)
fmt_ci <- function(lo, hi, digits = 2) if (is.null(lo) || is.na(lo) || is.na(hi)) "" else sprintf(" (95%% CI %s\u2013%s)", fmt_num(lo, digits), fmt_num(hi, digits))
num0 <- function(x) if (is.null(x) || !length(x)) NA_real_ else suppressWarnings(as.numeric(x[1]))

WS <- list(
  single_cell = list(
    title = "Single-cell CMI", nav = "Single-cell",
    lede = "Compute the CMI of each B cell and its grade: Below threshold, Mal-1, Mal-2 or Mal-3.",
    file_label = "Count matrix", accept = c(".zip", ".csv", ".tsv", ".txt", ".gz"),
    file_help = "Raw UMI counts of the B cells of one sample: 10x ZIP (matrix.mtx, features.tsv, barcodes.tsv) or genes \u00D7 cells CSV.",
    example = "examples/single_cell/simulated_b_cells.zip", example_synthetic = TRUE,
    example_note = "900 simulated B cells.",
    steps = c("Read and check counts", "Five dimensions", "CMI and malignant states", "Figures and bundle"),
    step_of = function(stage) switch(stage, queued = 0, start = 1, assets = 1, read = 1, ig = 2, mtp = 2, cnv = 2, normalize = 2, prol = 2, diff = 2, cmi = 3, export = 4, success = 5, 0),
    tabs = c("CMI", "Grades", "Cells", "QC", "Downloads"),
    figures = c("fig1_cmi", "fig2_grades")),
  spatial = list(
    title = "Spatial CMI grades", nav = "Spatial",
    lede = "Grade each spot of a tissue section: Mal-1, Mal-2, Mal-3 or non-malignant.",
    file_label = "Spot table", accept = c(".csv", ".tsv", ".txt", ".gz", ".zip"),
    file_help = "One section: spot_id, x, y and one deconvolution weight column per class, including Mal-1, Mal-2 and Mal-3.",
    example = "examples/spatial/spot_weights_demo.csv", example_synthetic = TRUE,
    example_note = "A simulated section, 2,171 spots.",
    steps = c("Read the spot table", "Grade each spot", "Figures and bundle"),
    step_of = function(stage) switch(stage, queued = 0, start = 1, read = 1, rctd = 1, grade = 2, export = 3, success = 4, 0),
    tabs = c("Map", "Spots", "Downloads"),
    figures = c("fig1_grade_map")),
  bulk = list(
    title = "Bulk CMI-riskscore", nav = "Bulk",
    lede = "Compute the CMI-riskscore of each sample and grade the cohort into CMI-high and CMI-low.",
    file_label = "Expression matrix", accept = c(".csv", ".tsv", ".txt", ".gz"),
    file_help = "Genes \u00D7 samples: gene symbols in the first column, all measured genes, values on the platform's usual scale.",
    example = "examples/bulk/expression.csv", example_synthetic = TRUE,
    example_note = "24 example samples.",
    steps = c("Standardize gene symbols", "Gene-pair features", "CMI-riskscore", "Figures and bundle"),
    step_of = function(stage) switch(stage, queued = 0, start = 1, pairs = 1, read = 1, hgnc = 1, reo = 2, rsf = 3, export = 4, success = 5, 0),
    tabs = c("Risk score", "Samples", "QC", "Downloads"),
    figures = c("fig1_riskscore"))
)

# fileInput with an accessible name for its read-only text box, which is also taken out of the tab order.
# data-max-mb lets the page refuse a file above the upload limit before sending it (www/app.js).
file_input <- function(id, label, accept, describedby, max_mb) {
  fi <- shiny::fileInput(id, label, accept = accept, buttonLabel = "Choose file", placeholder = "No file chosen")
  htmltools::tagQuery(fi)$find("input")$each(function(el, i) {
    if (identical(el$attribs$type, "text")) { el$attribs$tabindex <- "-1"; el$attribs[["aria-label"]] <- paste(label, "(selected file)") }
    else if (identical(el$attribs$type, "file")) { el$attribs[["aria-describedby"]] <- describedby; el$attribs[["data-max-mb"]] <- format(max_mb, scientific = FALSE) }
  })$allTags()
}
disabled_button <- function(id, label, class) htmltools::tagAppendAttributes(shiny::actionButton(id, label, class = class), disabled = NA)

workspace_ui <- function(id, module, lim) {
  ns <- shiny::NS(id); spec <- WS[[module]]
  up <- sprintf("Max %s MB", fmt_n(lim$max_upload_mb))
  limits <- switch(module,
    single_cell = sprintf("%s \u00B7 %s cells \u00B7 about 1 min", up, fmt_n(lim$max_cells)),
    spatial = sprintf("%s \u00B7 %s spots", up, fmt_n(lim$max_spots)),
    bulk = sprintf("%s \u00B7 %s values", up, fmt_n(lim$max_bulk_entries)))
  if (identical(lim$profile, "local")) limits <- sprintf("Local edition \u00B7 max %s MB \u00B7 no other limits", fmt_n(lim$max_upload_mb))
  local_ed <- identical(lim$profile, "local") || identical(Sys.getenv("CMI_ENABLE_RCTD"), "true")
  options <- if (module == "spatial" && local_ed)
    shiny::radioButtons(ns("spatial_mode"), "Input type", choices = c("Weights table" = "weights", "RCTD bundle (counts and reference)" = "rctd"), selected = "weights")
  shiny::div(class = "ws",
    shiny::p(class = "ws-lede", spec$lede),
    shiny::div(class = "ws-grid",
      shiny::div(class = "ws-side", shiny::div(class = "panel",
        shiny::h2("Your data"),
        file_input(ns("file"), spec$file_label, spec$accept, ns("file-help"), lim$max_upload_mb),
        shiny::p(class = "field-help", id = ns("file-help"), spec$file_help),
        shiny::uiOutput(ns("file_check")),
        options,
        disabled_button(ns("run"), "Run analysis", "btn-primary"),
        shiny::div(class = "or-line", "or"),
        shiny::actionButton(ns("example"), "Try the example", class = "btn-outline-secondary"),
        shiny::p(class = "field-help example-note", spec$example_note),
        shiny::div(class = "limits", limits),
        shiny::p(class = "guide-link", go_button("Input format and examples \u2192", "help", class = NULL, subtab = paste0("h-", module))))),
      shiny::div(class = "ws-main", id = ns("main"),
        shiny::div(class = "visually-hidden", role = "status", `aria-live` = "polite", `aria-atomic` = "true", shiny::textOutput(ns("sr_status"), inline = TRUE)),
        shiny::div(class = "visually-hidden", role = "alert", `aria-atomic` = "true", shiny::textOutput(ns("sr_alert"), inline = TRUE)),
        shiny::uiOutput(ns("status")),
        shiny::uiOutput(ns("results")))))
}

# The Analysis page: one workspace per data type, switched with pills.
analysis_page <- function(lim) shiny::div(class = "page ws-page",
  shiny::h1(class = "visually-hidden", "Analysis"),
  shiny::div(class = "ws-switch",
    bslib::navset_pill(id = "ws",
      bslib::nav_panel(WS$single_cell$nav, value = "single_cell", workspace_ui("sc", "single_cell", lim)),
      bslib::nav_panel(WS$spatial$nav, value = "spatial", workspace_ui("sp", "spatial", lim)),
      bslib::nav_panel(WS$bulk$nav, value = "bulk", workspace_ui("bulk", "bulk", lim)))))

# Input format shown in the results area until an analysis runs: what rows, columns and values are,
# with the first rows of the example input as a table.
FORMAT_EXAMPLES <- list(
  single_cell = list(
    lead = "A count matrix: one row per gene, one column per cell.",
    axes = c(Rows = "genes (symbols or Ensembl IDs)", Columns = "cells (barcodes)", Values = "raw UMI counts"),
    head = c("gene", "SIMCELL_0001", "SIMCELL_0002", "SIMCELL_0003", "SIMCELL_0004"),
    rows = list(c("MS4A1", "17", "16", "31", "2"), c("CD79A", "6", "2", "3", "7"), c("CD74", "17", "98", "166", "43"),
                c("IGKC", "13", "8", "30", "14"), c("IGHM", "6", "4", "2", "7")),
    notes = c("Upload it as CSV/TSV (as above) or as a 10x ZIP with matrix.mtx, features.tsv and barcodes.tsv.",
              "B cells of one sample, all measured genes.")),
  spatial = list(
    lead = "A spot table: one row per spot, with its position and one deconvolution weight per class.",
    axes = c(Rows = "spots", Columns = "spot_id, x, y, then one column per class", Values = "weights from 0 to 1"),
    head = c("spot_id", "x", "y", "Mal-1", "Mal-2", "Mal-3", "CD8_T", "Macrophage"),
    rows = list(c("SPOT_00001", "2400", "259.8", "0.000", "0.000", "0.002", "0.039", "0.113"),
                c("SPOT_00002", "2500", "259.8", "0.009", "0.023", "0.089", "0.002", "0.000"),
                c("SPOT_00111", "1650", "692.8", "0.142", "0.039", "0.363", "0.026", "0.004")),
    notes = c("Mal-1, Mal-2 and Mal-3 are required; every other weight column is a non-malignant class.",
              "One tissue section per file. Making the weights with RCTD: see Help.")),
  bulk = list(
    lead = "An expression matrix: one row per gene, one column per sample.",
    axes = c(Rows = "genes (HGNC symbols)", Columns = "samples", Values = "expression on the platform's usual scale"),
    head = c("gene", "SAMPLE_001", "SAMPLE_002", "SAMPLE_003", "SAMPLE_004"),
    rows = list(c("NUDT1", "10.29", "15.84", "10.38", "6.93"), c("MGST1", "2.69", "11.33", "0.97", "7.77"),
                c("WDR76", "1.24", "10.75", "2.60", "6.03"), c("SMC4", "19.47", "3.06", "7.20", "13.12")),
    notes = c("Upload it as CSV/TSV with all measured genes and at least 2 samples.",
              "Log2 intensity, TPM or counts; do not z-score genes across samples.")))

format_card <- function(module) {
  f <- FORMAT_EXAMPLES[[module]]; more <- "\u2026"
  shiny::div(class = "format-card",
    shiny::div(class = "format-top", shiny::h2(class = "h3", "Input format"),
               shiny::span(class = "format-hint", "Results appear here after you run an analysis or try the example.")),
    shiny::p(class = "format-lead", f$lead),
    shiny::tags$dl(class = "format-axes", lapply(names(f$axes), function(k) shiny::div(shiny::tags$dt(k), shiny::tags$dd(f$axes[[k]])))),
    shiny::div(class = "table-scroll", tabindex = "0", role = "region", `aria-label` = "Example of the input format",
      shiny::tags$table(class = "md-table",
        shiny::tags$thead(shiny::tags$tr(lapply(c(f$head, more), function(h) shiny::tags$th(scope = "col", h)))),
        shiny::tags$tbody(
          lapply(f$rows, function(r) shiny::tags$tr(shiny::tags$th(scope = "row", r[1]), lapply(c(r[-1], more), shiny::tags$td))),
          shiny::tags$tr(class = "md-more", lapply(seq_len(length(f$head) + 1L), function(i) shiny::tags$td(more)))))),
    shiny::tags$ul(class = "format-notes", lapply(f$notes, shiny::tags$li)))
}

stat_tile <- function(label, value, context = NULL) shiny::div(class = "stat", shiny::div(class = "stat-label", label), shiny::div(class = "stat-value", value),
                                                               if (!is.null(context)) shiny::div(class = "stat-context", context))
interpret <- function(...) shiny::div(class = "interpret", ...)
banner <- function(..., kind = "warn") shiny::div(class = paste("banner", kind), ...)

check_box <- function(ck, what = "File") {
  if (is.null(ck)) return(NULL)
  if (!isTRUE(ck$ok)) return(shiny::div(class = "file-check fail", role = "alert", shiny::tags$b(sprintf("%s cannot be analyzed", what)), shiny::span(ck$problems[1])))
  shiny::div(class = "file-check ok", shiny::tags$b(shiny::span(`aria-hidden` = "true", "\u2713 "), "Ready: "), shiny::span(ck$summary),
             if (length(ck$warnings)) shiny::tags$ul(lapply(ck$warnings, shiny::tags$li)))
}

# Precomputed example results (examples/precomputed/<module>/), unpacked once per server into a shared
# folder so that "Show the example" needs no analysis slot. tools/precompute_examples.R makes them.
.example_keys <- new.env(parent = emptyenv())
example_results_dir <- function(root, module) {
  src <- file.path(root, "examples", "precomputed", module)
  zipf <- file.path(src, "result_bundle.zip")
  if (!file.exists(zipf)) return(NULL)
  # keyed by the bundle's checksum, so that recomputed examples replace an unpacked older copy
  key <- .example_keys[[module]]
  if (is.null(key)) { key <- substr(digest::digest(file = zipf, algo = "sha256"), 1, 16); assign(module, key, envir = .example_keys) }
  dest <- file.path(dirname(tempdir()), "cmi_examples", CMI_VERSION, paste(module, key, sep = "-"))
  if (file.exists(file.path(dest, "summary.json"))) return(dest)
  tmp <- tempfile("example_", tmpdir = dirname(tempdir())); dir.create(tmp, recursive = TRUE)
  ok <- tryCatch({
    zip::unzip(zipf, exdir = tmp)
    webs <- list.files(file.path(src, "figures"), pattern = "_web\\.png$", full.names = TRUE)
    dir.create(file.path(tmp, "figures"), showWarnings = FALSE)
    file.copy(webs, file.path(tmp, "figures", basename(webs)), overwrite = TRUE); file.copy(zipf, file.path(tmp, "result_bundle.zip"))
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    file.exists(file.path(dest, "summary.json")) || file.rename(tmp, dest)
  }, error = function(e) FALSE)
  if (dir.exists(tmp)) unlink(tmp, recursive = TRUE)
  if (isTRUE(ok) && file.exists(file.path(dest, "summary.json"))) dest else NULL
}

workspace_server <- function(id, module, root, lim) shiny::moduleServer(id, function(input, output, session) {
  ns <- session$ns; spec <- WS[[module]]
  job <- shiny::reactiveVal(NULL); state <- shiny::reactiveVal(list(status = "idle")); res <- shiny::reactiveVal(NULL)
  check <- shiny::reactiveVal(NULL); last_run <- shiny::reactiveVal(NULL)
  phase <- shiny::reactiveVal("idle"); actions <- shiny::reactiveVal("none"); step <- shiny::reactiveVal(0)
  sr_status <- shiny::reactiveVal(""); sr_alert <- shiny::reactiveVal("")
  dirs <- character()
  send <- function(type, ...) session$sendCustomMessage(type, list(...))

  # ---- upload checks --------------------------------------------------------------------------
  mode <- function() input$spatial_mode %||% "weights"
  # A file above the upload limit is not sent (www/app.js). Its refusal replaces the previous check and
  # stands until another file arrives: the earlier upload (prev) is not checked or used again.
  refused <- function(r, prev) list(ok = FALSE, refused = TRUE, datapath = NA_character_, prev = prev %||% NA_character_,
    problems = sprintf("%s is %s MB; uploads to this server are limited to %s MB. Choose a smaller file, or run larger data sets in the local edition (see Help).",
                       as.character(r$name %||% "The file")[1], fmt_n(round(as.numeric(r$mb %||% NA)[1])), fmt_n(lim$max_upload_mb)))
  superseded <- function(old, f) isTRUE(old$refused) && identical(old$prev, f$datapath)
  run_check <- function() {
    f <- input$file
    if (is.null(f) || !nrow(f)) { check(NULL); return() }
    if (superseded(shiny::isolate(check()), f)) return()
    ck <- precheck_upload(module, f$datapath, f$name, lim, mode = mode()); ck$name <- f$name; ck$datapath <- f$datapath
    check(ck)
  }
  shiny::observeEvent(input$file, run_check())
  if (module == "spatial") shiny::observeEvent(input$spatial_mode, run_check(), ignoreInit = TRUE)
  shiny::observeEvent(input$file_refused, check(refused(input$file_refused, input$file$datapath)))
  output$file_check <- shiny::renderUI(check_box(check(), "This file"))
  can_run <- shiny::reactive(isTRUE(check()$ok) && is.null(job()))
  shiny::observe(send("cmi-enable", id = ns("run"), enabled = can_run()))
  shiny::observe(send("cmi-enable", id = ns("example"), enabled = is.null(job())))

  # ---- starting and following an analysis ------------------------------------------------------
  stage_input <- function(path, jd, name) {
    ext <- sub("^.*?((\\.[A-Za-z0-9]+)+)$", "\\1", basename(name)); if (identical(ext, basename(name))) ext <- ""
    staged <- file.path(jd, paste0("input", tolower(ext)))
    if (isTRUE(suppressWarnings(file.link(path, staged))) || isTRUE(file.copy(path, staged))) staged else NULL
  }
  launch <- function(path, name, synthetic) {
    if (!is.null(job())) return(invisible())
    res(NULL)
    jd <- tempfile(paste0(module, "_"), tmpdir = jobs_root()); dir.create(jd); dirs <<- c(dirs, jd)
    err_id <- substr(basename(jd), nchar(module) + 2L, nchar(basename(jd)))
    staged <- stage_input(path, jd, name)
    if (is.null(staged)) { state(list(status = "failed", detail = "The uploaded file could not be staged on the server.", error_id = err_id)); return(invisible()) }
    extra <- if (module == "spatial") c("--spatial-mode", mode()) else character()
    shown <- sub("^-", "_", gsub("[^A-Za-z0-9._ ()-]", "_", basename(name)))
    extra <- c(extra, "--input-name", shown)
    j <- tryCatch(start_worker(root, module, staged, jd, synthetic = synthetic, extra = extra),
                  error = function(e) { message(sprintf("[CMI %s] worker start failed: %s", err_id, conditionMessage(e))); NULL })
    if (is.null(j)) { state(list(status = "failed", detail = "The analysis could not be started on this server.", error_id = err_id)); return(invisible()) }
    j$input_name <- name; j$error_id <- err_id
    job(j); last_run(list(path = path, name = name, synthetic = synthetic))
    state(list(status = "queued", stage = "queued", progress = 0, detail = "Waiting for a free analysis slot on this server."))
    send("cmi-scroll", id = ns("main")); send("cmi-busy", id = id, busy = TRUE)
  }
  shiny::observeEvent(input$run, {
    f <- input$file; shiny::req(f, isTRUE(check()$ok), identical(check()$datapath, f$datapath))   # the checked file is the current one
    launch(f$datapath, f$name, FALSE)
  })
  live_example <- function() {
    launch(file.path(root, spec$example), basename(spec$example), isTRUE(spec$example_synthetic))
  }
  shiny::observeEvent(input$example, {
    if (!is.null(job())) return()
    d <- example_results_dir(root, module)
    r <- if (is.null(d)) NULL else tryCatch(load_results(module, d, synthetic = isTRUE(spec$example_synthetic), input_name = basename(spec$example), precomputed = TRUE),
                                            error = function(e) { message("[CMI] precomputed example unreadable: ", conditionMessage(e)); NULL })
    if (is.null(r)) return(live_example())
    res(r)
    state(list(status = "success", stage = "success", progress = 100, precomputed = TRUE,
               detail = "Example results, computed in advance with this version of the software."))
    send("cmi-focus", id = ns("results-heading"))
  })
  shiny::observeEvent(input$example_live, live_example())
  shiny::observeEvent(input$retry, { a <- last_run(); shiny::req(a); launch(a$path, a$name, a$synthetic) })
  shiny::observeEvent(input$cancel, {
    stop_worker(job()); job(NULL); res(NULL)
    state(list(status = "cancelled", detail = "Analysis cancelled. No results were kept."))
    send("cmi-busy", id = id, busy = FALSE)
  })

  shiny::observe({
    j <- job(); shiny::req(j)
    shiny::invalidateLater(1000, session)
    s <- read_state(j$state)
    if (!is.null(s) && !identical(s$status, "success")) state(s)
    if (j$proc$is_alive()) return()
    s <- read_state(j$state)
    job(NULL); send("cmi-busy", id = id, busy = FALSE)
    if (is.null(s) || !identical(s$status, "success")) {
      if (!is.null(s) && identical(s$status, "failed") && isTRUE(s$input_problem)) { state(s); return() }
      if (!is.null(s) && identical(s$status, "failed") && nzchar(s$detail %||% "")) message(sprintf("[CMI %s] %s", j$error_id, s$detail))
      detail <- if (is.null(s) || !identical(s$status, "failed")) exit_reason(j)
                else if (grepl("queue did not clear", s$detail %||% "")) s$detail
                else "The analysis stopped because of an internal error. Please try again; if it happens again, report the error ID."
      state(list(status = "failed", detail = detail, error_id = j$error_id)); return()
    }
    r <- tryCatch(load_results(module, j$out, synthetic = j$synthetic, input_name = j$input_name),
                  error = function(e) { message(sprintf("[CMI %s] results unreadable: %s", j$error_id, conditionMessage(e))); NULL })
    if (is.null(r)) { state(list(status = "failed", detail = "The results could not be read. Please try again.", error_id = j$error_id)); return() }
    state(s); res(r)
    send("cmi-focus", id = ns("results-heading"))
  })

  # ---- status panel: phase-level parts re-render rarely; the detail line follows each update ---
  shiny::observe({
    s <- state(); st <- s$status %||% "idle"
    phase(st)
    actions(switch(st, queued = "cancel", running = "cancel", failed = if (isTRUE(s$input_problem)) "none" else "retry", "none"))
    ns_ <- length(spec$steps); k <- if (st == "success") ns_ + 1 else spec$step_of(s$stage %||% ""); step(k)
    sr_status(switch(st,
      queued = paste("Queued.", s$detail %||% ""),
      running = if (k >= 1 && k <= ns_) sprintf("Step %d of %d: %s.", k, ns_, tolower(spec$steps[k])) else "Running.",
      success = paste(s$detail %||% "Done.", "Results below."),
      cancelled = "Analysis cancelled.", ""))
    sr_alert(if (st == "failed") paste("The analysis stopped.", s$detail %||% "") else "")
  })
  output$sr_status <- shiny::renderText(sr_status())
  output$sr_alert <- shiny::renderText(sr_alert())
  output$status <- shiny::renderUI({
    st <- phase()
    if (st == "idle") return(format_card(module))
    if (st == "failed" && isTRUE(shiny::isolate(state())$input_problem))
      return(shiny::div(class = "alert-input", shiny::tags$b("The input needs attention"), shiny::span(shiny::isolate(state())$detail)))
    shiny::div(class = paste("panel status", st),
      shiny::uiOutput(ns("status_line")),
      if (st %in% c("queued", "running")) shiny::tagList(shiny::uiOutput(ns("stepper")), shiny::uiOutput(ns("progress"))),
      shiny::uiOutput(ns("actions")))
  })
  output$status_line <- shiny::renderUI({
    s <- state(); st <- s$status %||% "idle"
    label <- switch(st, queued = "Queued", running = "Running", success = if (isTRUE(s$precomputed)) "Example" else "Done", failed = "Stopped", cancelled = "Cancelled", st)
    detail <- s$detail %||% ""
    if (st == "queued" && grepl("ahead of yours", detail)) detail <- paste(detail, "Each analysis takes about 10 seconds to a minute.")
    shiny::div(class = "status-line", shiny::span(class = paste("status-tag", st), label), shiny::span(class = "status-detail", detail),
               if (st == "failed" && !is.null(s$error_id)) shiny::span(class = "status-id", sprintf("Error ID %s", s$error_id)))
  })
  output$stepper <- shiny::renderUI({
    k <- step(); st <- phase()
    shiny::tags$ol(class = "stepper", `aria-label` = "Analysis steps", lapply(seq_along(spec$steps), function(i) {
      done <- i < k; active <- i == k && st == "running"
      shiny::tags$li(class = if (done) "done" else if (active) "active" else "", `aria-current` = if (active) "step" else NULL,
        shiny::span(class = "n", i), shiny::span(spec$steps[i]), shiny::span(class = "visually-hidden", if (done) " (done)" else if (active) " (in progress)" else ""))
    }))
  })
  output$progress <- shiny::renderUI({
    s <- state(); p <- max(2, min(100, num0(s$progress) %||% 0)); if (is.na(p)) p <- 2
    shiny::div(class = "progress-track", role = "progressbar", `aria-label` = "Analysis progress", `aria-valuemin` = 0, `aria-valuemax` = 100,
               `aria-valuenow` = round(p), `aria-valuetext` = s$detail %||% "", shiny::div(class = "progress-fill", style = sprintf("width:%s%%", p)))
  })
  output$actions <- shiny::renderUI({
    a <- actions()
    if (a == "cancel") return(shiny::div(class = "status-actions", shiny::actionButton(ns("cancel"), "Cancel the analysis", class = "btn-outline-secondary btn-sm")))
    if (a == "retry") return(shiny::div(class = "status-actions", shiny::actionButton(ns("retry"), "Try again", class = "btn-outline-secondary btn-sm")))
    NULL
  })

  output$results <- shiny::renderUI({ r <- res(); if (is.null(r)) return(NULL); results_ui(module, r, ns) })
  register_result_outputs(module, input, output, res, ns)

  session$onSessionEnded(function() { stop_worker(shiny::isolate(job())); for (d in dirs) unlink(d, recursive = TRUE) })
  list(dirs = function() dirs, busy = shiny::reactive(!is.null(job())))
})

load_results <- function(module, out, synthetic, input_name, precomputed = FALSE) {
  rd <- function(f, key = NULL) read_csv_safe(file.path(out, f), key)
  prov <- jsonlite::fromJSON(file.path(out, "provenance.json"), simplifyVector = TRUE)
  r <- list(module = module, out = out, synthetic = synthetic, input_name = input_name, precomputed = precomputed,
            completed = tryCatch(as.POSIXct(prov$completed_at, format = "%Y-%m-%dT%H:%M:%S%z"), error = function(e) Sys.time()),
            summary = jsonlite::fromJSON(file.path(out, "summary.json"), simplifyVector = TRUE),
            provenance = prov, figures = list.files(file.path(out, "figures"), pattern = "_web\\.png$", full.names = TRUE))
  if (module == "single_cell") { r$table <- rd("cells.csv", "cell"); r$audit <- rd("qc/gene_mapping_audit.csv") }
  if (module == "bulk") { r$table <- rd("samples.csv", "sample_id"); r$qc <- rd("qc/dataset_qc.csv"); r$mapping <- rd("qc/reo_endpoint_mapping.csv") }
  if (module == "spatial") { r$table <- rd("spots.csv", "spot_id"); r$grades <- rd("grades.csv") }
  r
}

fig_path <- function(r, stem) { f <- r$figures[grepl(paste0("/", stem, "_web\\.png$"), r$figures)]; if (length(f)) f[1] else NULL }
fig_files_links <- function(ns, stem) shiny::div(class = "figure-files", shiny::span("Figure files:"),
  shiny::downloadLink(ns(paste0("fs_", stem)), "SVG"), shiny::downloadLink(ns(paste0("fp_", stem)), "PDF"), shiny::downloadLink(ns(paste0("fg_", stem)), "PNG, 600 dpi"))
fig_out <- function(ns, id, stem) shiny::div(class = "figure result-figure", shiny::imageOutput(ns(id), height = "auto"), fig_files_links(ns, stem))

# A horizontally scrollable box that keyboard users can focus and scroll; each needs a unique label.
scroll_region <- function(..., label = NULL, labelledby = NULL)
  shiny::div(class = "table-scroll", tabindex = "0", role = "region", `aria-label` = label, `aria-labelledby` = labelledby, ...)

# Plain HTML tables (UTF-8 safe in any locale; no xtable).
html_table <- function(df, align = NULL, caption = NULL, row_header = TRUE, label = NULL, code_first = FALSE, class = NULL) {
  al <- if (is.null(align)) rep("l", ncol(df)) else strsplit(align, "")[[1]]
  cell <- function(j, v, head = FALSE) {
    cls <- if (al[j] == "r") "num" else NULL
    if (code_first && j == 1L) v <- shiny::tags$code(v)          # column or file names
    if (head) shiny::tags$th(scope = "row", class = cls, v) else shiny::tags$td(class = cls, v)
  }
  scroll_region(label = label %||% caption %||% paste("Table:", paste(names(df), collapse = ", ")),
    shiny::tags$table(class = paste(c("table-plain table-compact", if (code_first) "col-code", class), collapse = " "),
    if (!is.null(caption)) shiny::tags$caption(caption),
    shiny::tags$thead(shiny::tags$tr(lapply(seq_along(df), function(j) shiny::tags$th(scope = "col", class = if (al[j] == "r") "num" else NULL, names(df)[j])))),
    shiny::tags$tbody(lapply(seq_len(nrow(df)), function(i) shiny::tags$tr(lapply(seq_along(df), function(j) cell(j, as.character(df[[j]][i]), row_header && j == 1L)))))))
}

state_table <- function(s) {
  n <- unlist(s$state_counts)
  data.frame(Grade = names(n), Cells = fmt_vec(n), Share = sprintf("%.1f%%", 100 * n / max(1, sum(n))), check.names = FALSE)
}
grade_table <- function(g) data.frame(Grade = g$grade, Spots = fmt_vec(g$spots), `Share of spots` = sprintf("%.1f%%", 100 * g$spot_share),
  `Share of Mal-1/2/3 weight` = ifelse(is.na(g$malignant_weight_share), "\u2013", sprintf("%.1f%%", 100 * g$malignant_weight_share)), check.names = FALSE)
bulk_qc_table <- function(q) data.frame(Measure = c("Samples", "Genes in the rank pool", "Pair genes found", "Gene pairs computable", "Coverage", "Pair genes matched exactly / after symbol update / as alias"),
  Value = c(fmt_vec(q$samples), fmt_vec(q$expression_genes), sprintf("%d of %d", q$endpoint_genes_observed, q$endpoint_genes_total),
            sprintf("%d of %d", q$gene_pairs_computable, q$gene_pairs_total), tolower(q$coverage_class), sprintf("%d / %d / %d", q$exact_matches, q$normalized_matches, q$alias_matches)), check.names = FALSE)

results_meta <- function(r, ns) {
  when <- format(r$completed, "%Y-%m-%d %H:%M")
  if (isTRUE(r$precomputed)) return(shiny::div(class = "results-meta",
    shiny::span(sprintf("Example: %s \u00B7 computed in advance on %s with CMI Explorer %s", r$input_name, format(r$completed, "%Y-%m-%d"), r$provenance$app_version %||% CMI_VERSION)),
    shiny::actionLink(ns("example_live"), "Run the example live")))
  shiny::div(class = "results-meta", shiny::span(sprintf("Results for %s \u00B7 %s \u00B7 CMI Explorer %s", r$input_name, when, CMI_VERSION)),
             shiny::span(class = "muted", "Results are kept while this page is open; download the bundle to keep them."))
}

results_ui <- function(module, r, ns) {
  s <- r$summary; spec <- WS[[module]]; tl <- spec$tabs
  head <- shiny::tagList(shiny::h2(id = ns("results-heading"), tabindex = "-1", class = "results-heading", "Results"), results_meta(r, ns))
  demo <- if (isTRUE(r$synthetic)) banner(kind = "demo", shiny::tags$b("Example data. "), "These results come from the example input.") else NULL
  notes <- list(); stats <- NULL; tabs <- NULL
  if (module == "single_cell") {
    lin <- s$lineage; n <- unlist(s$state_counts)
    if (!is.null(lin) && isTRUE(num0(lin$no_b_marker_pct) > 10 || num0(lin$non_b_marker_pct) > 10)) notes <- c(notes, list(banner(shiny::tags$b("Check the cell lineage. "),
      sprintf("%s of cells express no B-cell marker and %s are dominated by non-B-cell markers. The CMI applies to B cells only: subset to B cells and run again.",
              fmt_pct(lin$no_b_marker_pct, 0), fmt_pct(lin$non_b_marker_pct, 0)))))
    stats <- shiny::div(class = "stats",
      stat_tile("Cells", fmt_n(s$cells), sprintf("%s at or above CMI 0.9755", fmt_pct(s$malignant_pct))),
      stat_tile("Mal-1", fmt_n(n[["Mal-1"]]), sprintf("%s of cells", fmt_pct(100 * n[["Mal-1"]] / max(1, sum(n))))),
      stat_tile("Mal-2", fmt_n(n[["Mal-2"]]), sprintf("%s of cells", fmt_pct(100 * n[["Mal-2"]] / max(1, sum(n))))),
      stat_tile("Mal-3", fmt_n(n[["Mal-3"]]), sprintf("%s of cells", fmt_pct(100 * n[["Mal-3"]] / max(1, sum(n))))))
    tabs <- list(
      bslib::nav_panel(tl[1], fig_out(ns, "f1", "fig1_cmi")),
      bslib::nav_panel(tl[2], fig_out(ns, "f2", "fig2_grades"), html_table(state_table(s), "lrr", label = "Cells per grade")),
      bslib::nav_panel(tl[3], shiny::p(class = "small muted", "Select a column name to sort. The full table is in the result bundle."), DT::DTOutput(ns("table"))),
      bslib::nav_panel(tl[4], shiny::h3("Gene matching"), html_table(audit_table(r$audit), "ll", label = "Gene matching"),
        shiny::h3("Cell lineage"),
        html_table(data.frame(Check = c("Cells with no B-cell marker transcript (CD19, MS4A1, CD79A, CD79B, PAX5)", "Cells dominated by T-cell, myeloid, NK, erythroid or platelet markers"),
                              Result = c(fmt_pct(lin$no_b_marker_pct), fmt_pct(lin$non_b_marker_pct)), check.names = FALSE), "lr", label = "Cell lineage"),
        shiny::p(class = "small muted", "Upload check only; not part of the CMI.")),
      bslib::nav_panel(tl[5], downloads_ui(ns, module, r)))
  }
  if (module == "bulk") {
    n_s <- num0(s$samples)
    if (n_s < 20) notes <- c(notes, list(banner(kind = "info", shiny::tags$b("Few samples. "), sprintf("The CMI-high/low split at the median of %s samples is unstable.", fmt_n(n_s)))))
    stats <- shiny::div(class = "stats",
      stat_tile("Samples", fmt_n(s$samples), sprintf("%s genes in the rank pool", fmt_n(s$genes_in_rank_pool))),
      stat_tile("Gene pairs computable", fmt_pct(100 * num0(s$pair_coverage), 0), sprintf("Coverage: %s", tolower(s$coverage_class))),
      stat_tile("CMI-high / CMI-low", sprintf("%s / %s", fmt_n(s$groups[["CMI-high"]] %||% 0), fmt_n(s$groups[["CMI-low"]] %||% 0)), "samples"),
      stat_tile("Cohort median", fmt_num(s$median_riskscore), "CMI-riskscore"))
    tabs <- list(
      bslib::nav_panel(tl[1], fig_out(ns, "f1", "fig1_riskscore"),
        interpret("Higher score = higher predicted risk. CMI-high: above the cohort median.")),
      bslib::nav_panel(tl[2], DT::DTOutput(ns("table"))),
      bslib::nav_panel(tl[3], html_table(bulk_qc_table(r$qc), "ll", label = "Gene-pair coverage"),
        shiny::p(class = "small muted", "Coverage: share of the 100 gene pairs with both genes present; below 70% the score is labeled limited."),
        shiny::h3("Gene-pair genes"), DT::DTOutput(ns("mapping"))),
      bslib::nav_panel(tl[4], downloads_ui(ns, module, r)))
  }
  if (module == "spatial") {
    g <- r$grades; sh <- function(k) { v <- g$spot_share[g$grade == k]; if (length(v)) fmt_pct(100 * v) else "\u2013" }
    stats <- shiny::div(class = "stats",
      stat_tile("Spots", fmt_n(s$spots), sprintf("%d classes", length(s$classes))),
      stat_tile("Mal-1", sh("Mal-1"), "of spots"), stat_tile("Mal-2", sh("Mal-2"), "of spots"), stat_tile("Mal-3", sh("Mal-3"), "of spots"))
    tabs <- list(
      bslib::nav_panel(tl[1], fig_out(ns, "f1", "fig1_grade_map"), html_table(grade_table(g), "lrrr", label = "Spots per grade")),
      bslib::nav_panel(tl[2], DT::DTOutput(ns("table"))),
      bslib::nav_panel(tl[3], downloads_ui(ns, module, r)))
  }
  shiny::tagList(head, demo, notes, stats, shiny::div(class = "result-tabs", do.call(bslib::navset_tab, c(tabs, list(id = ns("tabs"))))))
}

downloads_ui <- function(ns, module, r) {
  size <- function(f) { p <- file.path(r$out, f); if (!file.exists(p)) return(""); b <- file.size(p); if (b > 1024^2) sprintf(" (%.1f MB)", b / 1024^2) else sprintf(" (%.0f KB)", b / 1024) }
  figs <- list.files(file.path(r$out, "figures"), recursive = TRUE, full.names = TRUE)
  figs_bytes <- sum(file.size(figs[!grepl("_web\\.png$", figs)]), na.rm = TRUE)
  item <- function(id, label, desc, primary = FALSE) shiny::div(class = "dl-item",
    shiny::downloadButton(ns(id), label, class = if (primary) "btn-primary" else "btn-outline-secondary"), shiny::span(class = "dl-desc", desc))
  tbl <- switch(module, single_cell = "cells.csv", bulk = "samples.csv", spatial = "spots.csv")
  shiny::tagList(shiny::div(class = "dl-list",
    item("dl_bundle", "Result bundle (.zip)", paste0("Every table, each figure as SVG, PDF and 600 dpi PNG with legends, quality-control files, a README and provenance.json", size("result_bundle.zip"), "."), TRUE),
    item("dl_table", sprintf("Main table (.csv)"), paste0("One row per ", switch(module, single_cell = "cell", bulk = "sample", spatial = "spot"), size(tbl), ".")),
    item("dl_figs", "Figures (.zip)", sprintf("SVG (editable text), PDF (vector, embedded fonts), 600 dpi PNG, panels without titles and legends.txt (about %.1f MB).", figs_bytes / 1024^2)),
    item("dl_prov", "Provenance (.json)", "Input checksum, parameters, software versions and the SHA-256 of every model file used.")),
    shiny::p(class = "small muted how-to-report", "Keep provenance.json with your results: it records the input checksum, parameters and software versions."))
}

AUDIT_LABELS <- c(query_genes = "Genes in the upload", ensembl_fraction = "Share of Ensembl identifiers", IG_canonical_covered = "Immunoglobulin genes found",
                  MTP_canonical_covered = "Malignant-program genes found", S_genes_resolved = "S-phase genes matched", G2M_genes_resolved = "G2/M genes matched",
                  CCAT_network_genes_mapped = "Genes mapped to the CCAT network", CCAT_network_genes_used = "Network genes used for CCAT")
audit_table <- function(a) {
  if (is.null(a) || !nrow(a)) return(data.frame(Check = "Gene matching", Result = "not available", check.names = FALSE))
  nm <- as.character(a[[1]]); val <- suppressWarnings(as.numeric(a[[2]])); det <- as.character(a[[3]])
  tot <- suppressWarnings(as.numeric(sub("^/", "", det)))
  shown <- ifelse(nm == "ensembl_fraction", sprintf("%.0f%%", 100 * val),
           ifelse(nm %in% c("S_genes_resolved", "G2M_genes_resolved") & !is.na(tot), sprintf("%s features for %s genes", fmt_vec(val), fmt_vec(tot)),
           ifelse(!is.na(tot), sprintf("%s of %s", fmt_vec(val), fmt_vec(tot)), fmt_vec(val))))
  note <- ifelse(nm == "query_genes", ifelse(det == "ENSEMBL", "read as Ensembl IDs", "read as gene symbols"), "")
  data.frame(Check = ifelse(nm %in% names(AUDIT_LABELS), AUDIT_LABELS[nm], nm), Result = paste0(shown, ifelse(nzchar(note), paste0(" (", note, ")"), "")), check.names = FALSE)
}
fmt_vec <- function(x) ifelse(is.na(x), "\u2013", formatC(x, format = "d", big.mark = ","))

# Columns of the interactive tables: order, labels and decimals.
TABLE_SPEC <- list(
  single_cell = list(cols = c(cell = "Cell", CMI = "CMI", state = "Grade", Ig_score = "Ig-score", cnv_score = "CNV", mtp_score = "MTP",
                              diff_score = "Differentiation", prol_score = "Proliferation", n_counts = "UMIs", b_markers_detected = "B-cell marker",
                              non_b_markers_dominant = "Non-B markers dominate"),
                     digits = c(CMI = 5, Ig_score = 3, cnv_score = 4, mtp_score = 3, diff_score = 3, prol_score = 3, n_counts = 0)),
  bulk = list(cols = c(sample_id = "Sample", riskscore = "CMI-riskscore", cmi_group = "Group", percentile = "Percentile in cohort", pair_coverage = "Pair coverage", coverage_class = "Coverage"),
              digits = c(riskscore = 2, percentile = 1, pair_coverage = 2)),
  spatial = list(cols = c(spot_id = "Spot", grade = "Grade", dominant_class = "Class of highest weight", malignant_load = "Mal-1/2/3 weight", x = "x", y = "y"),
                 digits = c(malignant_load = 3, x = 1, y = 1)))
MAPPING_LABELS <- c(target_gene = "Pair gene", source_gene = "Gene in upload", source_index = "Row in upload", match_type = "Match")

datatable_for <- function(d, module) {
  spec <- TABLE_SPEC[[module]]
  keep <- c(intersect(names(spec$cols), names(d)), if (module == "spatial") grep("^weight_", names(d), value = TRUE))
  d <- d[, keep, drop = FALSE]
  lab <- ifelse(names(d) %in% names(spec$cols), spec$cols[names(d)], paste("Weight:", pretty_class(canonical_class(sub("^weight_", "", names(d))))))
  if ("dominant_class" %in% names(d)) d$dominant_class <- pretty_class(d$dominant_class)
  lg <- vapply(d, is.logical, logical(1)); d[lg] <- lapply(d[lg], function(v) ifelse(v, "Yes", "No"))
  dg <- spec$digits; wcol <- grep("^weight_", names(d), value = TRUE)
  dg <- c(dg[intersect(names(dg), names(d))], stats::setNames(rep(3, length(wcol)), wcol))
  defs <- lapply(names(dg), function(k) list(targets = match(k, names(d)) - 1L, render = num_render(dg[[k]])))
  DT::datatable(d, rownames = FALSE, colnames = unname(lab), class = "compact stripe hover",
                options = list(pageLength = 15, scrollX = TRUE, dom = "ftip", autoWidth = FALSE, pagingType = "simple_numbers",
                               columnDefs = defs, language = list(search = "Search the table:")))
}
# Display rounding with thousands separators and a true minus sign; sorting and search use the raw values.
num_render <- function(digits) DT::JS(sprintf(paste0("function(d, t) { if (t !== 'display' || d === null || d === '' || isNaN(d)) return d; ",
  "return Number(d).toLocaleString('en-US', {minimumFractionDigits: %d, maximumFractionDigits: %d}).replace('-', '\\u2212'); }"), digits, digits))

# Alternative text of the result figures, written from the numbers they show.
fig_alt <- function(module, stem, r) {
  s <- r$summary
  if (module == "single_cell") {
    n <- unlist(s$state_counts)
    return(switch(stem,
      fig1_cmi = sprintf("Histogram of the CMI of %s cells: %s below the threshold, %s Mal-1, %s Mal-2 and %s Mal-3.",
                         fmt_n(s$cells), fmt_n(n[["Below threshold"]]), fmt_n(n[["Mal-1"]]), fmt_n(n[["Mal-2"]]), fmt_n(n[["Mal-3"]])),
      fig2_grades = sprintf("Bar chart of cells per CMI grade: %s below the threshold, %s Mal-1, %s Mal-2 and %s Mal-3.",
                            fmt_n(n[["Below threshold"]]), fmt_n(n[["Mal-1"]]), fmt_n(n[["Mal-2"]]), fmt_n(n[["Mal-3"]]))))
  }
  if (module == "bulk") return(sprintf("CMI-riskscores of %s samples in rank order; %s above the cohort median of %s are CMI-high.",
                                       fmt_n(s$samples), fmt_n(s$groups[["CMI-high"]] %||% 0), fmt_num(s$median_riskscore)))
  sprintf("Map of %s spots colored by CMI grade: Mal-1, Mal-2, Mal-3 or non-malignant.", fmt_n(s$spots))
}

safe_stem <- function(name) { b <- sub("(\\.[A-Za-z0-9]+)+$", "", basename(name %||% "input")); b <- gsub("[^A-Za-z0-9_-]+", "_", b); substr(if (nzchar(b)) b else "input", 1, 40) }

register_result_outputs <- function(module, input, output, res, ns) {
  spec <- WS[[module]]
  img <- function(stem) shiny::renderImage({
    r <- res(); shiny::req(r); f <- fig_path(r, stem); shiny::req(f)
    list(src = f, contentType = "image/png", alt = fig_alt(module, stem, r), class = "result-img")
  }, deleteFile = FALSE)
  ids <- switch(module, single_cell = c(f1 = "fig1_cmi", f2 = "fig2_grades"), bulk = c(f1 = "fig1_riskscore"), spatial = c(f1 = "fig1_grade_map"))
  for (k in names(ids)) local({ kk <- k; output[[kk]] <- img(ids[[kk]]) })
  stamp <- function() { r <- res(); sprintf("CMI_%s_%s_%s", module, safe_stem(r$input_name), format(r$completed %||% Sys.time(), "%Y%m%d-%H%M")) }
  out_file <- function(rel) { r <- res(); shiny::req(r); p <- file.path(r$out, rel); shiny::validate(shiny::need(file.exists(p), "These results are no longer on the server. Run the analysis again.")); p }
  # Figure files in each format.
  for (stem in spec$figures) local({
    st <- stem
    output[[paste0("fs_", st)]] <- shiny::downloadHandler(filename = function() paste0(stamp(), "_", st, ".svg"), content = function(file) file.copy(out_file(file.path("figures", paste0(st, ".svg"))), file), contentType = "image/svg+xml")
    output[[paste0("fp_", st)]] <- shiny::downloadHandler(filename = function() paste0(stamp(), "_", st, ".pdf"), content = function(file) file.copy(out_file(file.path("figures", paste0(st, ".pdf"))), file), contentType = "application/pdf")
    output[[paste0("fg_", st)]] <- shiny::downloadHandler(filename = function() paste0(stamp(), "_", st, "_600dpi.png"), content = function(file) file.copy(out_file(file.path("figures", paste0(st, "_600dpi.png"))), file), contentType = "image/png")
  })
  output$table <- DT::renderDT({ r <- res(); shiny::req(r$table); datatable_for(r$table, module) }, server = TRUE)
  if (module == "bulk") output$mapping <- DT::renderDT({ r <- res(); shiny::req(r$mapping); m <- r$mapping
    names(m) <- ifelse(names(m) %in% names(MAPPING_LABELS), MAPPING_LABELS[names(m)], names(m))
    DT::datatable(m, rownames = FALSE, class = "compact stripe", options = list(pageLength = 10, dom = "tip", scrollX = TRUE, pagingType = "simple_numbers")) })
  main_csv <- switch(module, single_cell = "cells.csv", bulk = "samples.csv", spatial = "spots.csv")
  output$dl_bundle <- shiny::downloadHandler(filename = function() paste0(stamp(), "_results.zip"), content = function(file) file.copy(out_file("result_bundle.zip"), file), contentType = "application/zip")
  output$dl_table <- shiny::downloadHandler(filename = function() paste0(stamp(), "_", main_csv), content = function(file) file.copy(out_file(main_csv), file), contentType = "text/csv")
  output$dl_figs <- shiny::downloadHandler(filename = function() paste0(stamp(), "_figures.zip"), content = function(file) {
    r <- res(); shiny::req(r); f <- list.files(file.path(r$out, "figures"), recursive = TRUE); f <- f[!grepl("_web\\.png$", f)]
    shiny::validate(shiny::need(length(f) > 0, "These results are no longer on the server. Run the analysis again."))
    zip::zip(file, files = file.path("figures", f), root = r$out, mode = "mirror") }, contentType = "application/zip")
  output$dl_prov <- shiny::downloadHandler(filename = function() paste0(stamp(), "_provenance.json"), content = function(file) file.copy(out_file("provenance.json"), file), contentType = "application/json")
}
