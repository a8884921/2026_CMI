# =============================================================================
# Help and About pages. The Analysis page is analysis_page() in R/workspace.R.
# =============================================================================

# Links that switch the page (data-go), optionally scrolling to an element (data-subtab) or pressing the
# example button of a workspace (data-example). Handled in www/app.js.
go_button <- function(label, tab, class = "btn btn-primary", subtab = NULL, example = NULL, .noWS = NULL)
  tags$a(href = paste0("#", tab), `data-go` = tab, `data-subtab` = subtab, `data-example` = example, class = class, label, .noWS = .noWS)
code_block <- function(...) tags$pre(tabindex = "0", tags$code(paste(c(...), collapse = "\n"), .noWS = c("after-begin", "before-end")), .noWS = c("after-begin", "before-end"))
size_label <- function(path) { b <- file.size(path); if (is.na(b)) "" else if (b >= 1024^2) sprintf("%.1f MB", b / 1024^2) else sprintf("%.0f KB", max(1, b / 1024)) }
# Term / description pairs.
facts <- function(...) { x <- list(...); tags$dl(class = "facts", lapply(seq(1, length(x), by = 2), function(i) div(class = "fact", tags$dt(x[[i]]), tags$dd(x[[i + 1]])))) }
# A download link to a file under www/ (served at the root) or examples/.
file_link <- function(path, root = get0("APP_ROOT", ifnotfound = ".")) {
  href <- sub("^www/", "", path)
  tags$a(class = "file-dl", href = href, download = basename(path), tags$code(basename(path), .noWS = "outside"),
         span(class = "muted", sprintf(" %s", size_label(file.path(root, path))), .noWS = "outside"), .noWS = "outside")
}
example_row <- function(module, id, files) div(class = "help-example",
  span(class = "help-label", "Example"), lapply(files, file_link),
  go_button("Try it \u2192", module, class = "btn btn-outline-secondary btn-sm", example = paste0(id, "-example")))
fold <- function(summary, ...) tags$details(class = "fold", tags$summary(summary), div(class = "fold-body", ...))

SEURAT_EXPORT <- c("library(Seurat); library(Matrix)",
                   "b <- subset(obj, subset = sample_id == \"S01\" & cell_type == \"B cell\")   # one sample, B cells",
                   "b <- JoinLayers(b, assay = \"RNA\")                                      # Seurat 5",
                   "counts <- LayerData(b, assay = \"RNA\", layer = \"counts\")                # raw UMI counts",
                   "dir.create(\"S01\")",
                   "writeMM(counts, \"S01/matrix.mtx\")",
                   "writeLines(paste(rownames(counts), rownames(counts), \"Gene Expression\", sep = \"\\t\"), \"S01/features.tsv\")",
                   "writeLines(colnames(counts), \"S01/barcodes.tsv\")",
                   "zip::zip(file.path(getwd(), \"S01.zip\"), root = \"S01\", files = c(\"matrix.mtx\", \"features.tsv\", \"barcodes.tsv\"))")
RCTD_LABELS <- c("# reference_labels.csv: one row per reference cell (cell, cell_type, CMI)",
                 "cmi <- do.call(rbind, lapply(Sys.glob(\"results/*/cells.csv\"), read.csv))   # Single-cell results",
                 "lab <- data.frame(cell = colnames(ref), cell_type = as.character(ref$cell_type),",
                 "                  CMI = cmi$CMI[match(colnames(ref), cmi$cell)])",
                 "lab$cell_type[lab$cell_type == \"B cell\" & ref$sample_type == \"lymphoma\"] <- \"Malignant\"",
                 "lab$cell_type[lab$cell_type == \"B cell\"] <- \"Normal_B\"",
                 "write.csv(lab, \"reference_labels.csv\", row.names = FALSE)",
                 "saveRDS(LayerData(ref, assay = \"RNA\", layer = \"counts\"), \"reference_counts.rds\")")
RCTD_CALL <- c("Rscript rctd_mal_states.R --ref-counts reference_counts.rds --ref-labels reference_labels.csv \\",
               "  --spatial section1/outs --output section1_weights.csv --cores 4")

state_cut_table <- function() {
  row <- function(state, cmi) tags$tr(tags$th(scope = "row", state), tags$td(cmi))
  div(class = "table-scroll", tabindex = "0", role = "region", `aria-label` = "CMI grades",
    tags$table(class = "table-plain table-compact cut-table",
      tags$thead(tags$tr(tags$th(scope = "col", "Grade"), tags$th(scope = "col", "CMI"))),
      tags$tbody(
        row("Below threshold", "< 0.9755"),
        row("Mal-1", "0.9755\u20130.98926"),
        row("Mal-2", "0.98926\u20130.99620"),
        row("Mal-3", "> 0.99620"))))
}

# ---- Help ----------------------------------------------------------------------------------------
help_page <- function(lim) {
  local_note <- "Larger data: use the local edition (see below)."
  div(class = "page", div(class = "reading help",
    div(class = "page-head", h1("Help")),
    tags$ol(class = "quick-steps",
      tags$li(tags$b("Choose a data type"), span("Single-cell, Spatial or Bulk, on the Analysis page.")),
      tags$li(tags$b("Upload a file"), span("Or select Try the example.")),
      tags$li(tags$b("Run and download"), span("Figures and tables appear on the page; download the result bundle to keep them."))),
    tags$nav(class = "page-toc", `aria-label` = "Help contents", tags$ol(
      tags$li(tags$a(href = "#h-single_cell", "Single-cell")), tags$li(tags$a(href = "#h-spatial", "Spatial")), tags$li(tags$a(href = "#h-bulk", "Bulk")),
      tags$li(tags$a(href = "#h-states", "CMI grades")), tags$li(tags$a(href = "#h-local", "Local edition")))),

    tags$section(id = "h-single_cell", class = "help-sec",
      h2("Single-cell CMI"),
      facts("File", "10x Matrix Market ZIP (matrix.mtx, features.tsv, barcodes.tsv; gzipped is fine) or a genes \u00D7 cells CSV/TSV",
            "Cells", "B cells of one sample, after quality control",
            "Values", "Raw UMI counts; human gene symbols or Ensembl IDs",
            "Limits", sprintf("%s MB, %s cells, %s non-zero counts. %s", fmt_n(lim$max_upload_mb), fmt_n(lim$max_cells), fmt_n(lim$max_nnz), local_note)),
      example_row("single_cell", "sc", "examples/single_cell/simulated_b_cells.zip"),
      fold("Export from Seurat", code_block(SEURAT_EXPORT)),
      fold("Export from Scanpy (.h5ad)",
        p("Download ", file_link("www/downloads/h5ad_to_mtx.py"), " (needs anndata, scipy, numpy), then:"),
        code_block("python3 h5ad_to_mtx.py data.h5ad S01.zip --layer counts --where sample_id=S01 --where cell_type=B_cell")),
      p(class = "help-out", tags$b("Output: "), "the CMI and CMI grade of every cell (cells.csv), with figures.")),

    tags$section(id = "h-spatial", class = "help-sec",
      h2("Spatial CMI grades"),
      facts("File", "CSV/TSV of one tissue section",
            "Columns", "spot_id, x, y and one deconvolution weight (0\u20131) per class",
            "Classes", "Mal-1, Mal-2 and Mal-3 are required; other classes count as non-malignant",
            "Limits", sprintf("%s MB, %s spots. %s", fmt_n(lim$max_upload_mb), fmt_n(lim$max_spots), local_note)),
      example_row("spatial", "sp", "examples/spatial/spot_weights_demo.csv"),
      fold("Make the weights table with RCTD",
        tags$ol(class = "help-steps",
          tags$li("Score the B cells of your single-cell reference in Single-cell (one sample per run)."),
          tags$li("Label the reference: malignant B cells as \"Malignant\" with their CMI; other cells keep their cell type. The script splits the malignant cells into Mal-1, Mal-2 and Mal-3 at the CMI cut points.", code_block(RCTD_LABELS)),
          tags$li("Run the script on each section (R \u2265 4.3 with spacexr), then upload the table it writes.", code_block(RCTD_CALL))),
        div(class = "help-files", file_link("www/downloads/rctd_mal_states.R"), file_link("www/downloads/rctd_test_data.zip"), file_link("www/downloads/state_cutpoints.csv"))),
      p(class = "help-out", tags$b("Output: "), "the CMI grade of every spot \u2013 Mal-1, Mal-2, Mal-3 or non-malignant (spots.csv, grades.csv) \u2013 with a map.")),

    tags$section(id = "h-bulk", class = "help-sec",
      h2("Bulk CMI-riskscore"),
      facts("Expression", "Genes \u00D7 samples CSV: HGNC symbols in the first column, all measured genes, at least 2 samples",
            "Values", "The platform's usual scale (log2 intensity, TPM or counts), probes collapsed to genes; do not z-score genes across samples",
            "Limits", sprintf("%s MB, %s expression values. %s", fmt_n(lim$max_upload_mb), fmt_n(lim$max_bulk_entries), local_note)),
      example_row("bulk", "bulk", "examples/bulk/expression.csv"),
      fold("Collapse probes to genes (R)",
        code_block("expr <- expr[order(-rowMeans(expr)), ]                     # probes x samples",
                   "expr <- expr[!duplicated(probe_symbol[rownames(expr)]), ]",
                   "rownames(expr) <- probe_symbol[rownames(expr)]",
                   "write.csv(data.frame(gene = rownames(expr), expr, check.names = FALSE), \"expression.csv\", row.names = FALSE)")),
      p(class = "help-out", tags$b("Output: "), "the CMI-riskscore and CMI-high/low grade of every sample (samples.csv; CMI-high = above the cohort median), with a figure.")),

    tags$section(id = "h-states", class = "help-sec",
      h2("CMI grades"),
      state_cut_table(),
      p(class = "small muted", "The CMI is a model score from 0 to 1, not a probability. A value on a boundary belongs to the lower grade.")),

    tags$section(id = "h-local", class = "help-sec",
      h2("Local edition"),
      p("Runs the same analyses on your computer, without upload limits, and adds RCTD deconvolution:"),
      code_block("Rscript tools/install_packages.R     # once", "Rscript run_local.R"),
      fold("Score many samples from the command line",
        code_block("for f in samples/*.zip; do",
                   "  s=$(basename \"$f\" .zip)",
                   "  CMI_PROFILE=local Rscript engine/worker.R --root . --module single_cell \\",
                   "    --input \"$f\" --output \"results/$s\" --state \"results/$s/state.json\"",
                   "done"),
        p(class = "small muted", "Other modules: --module bulk, --module spatial.")))))
}

# ---- About ---------------------------------------------------------------------------------------
about_page <- function() div(class = "page", div(class = "reading help",
  div(class = "page-head", h1("About"),
      p(class = "lede", "CMI Explorer computes the clonal malignancy index (CMI) and its grades for single-cell, spatial and bulk B-cell lymphoma data, from fixed model parameters.")),
  facts("Your data", "Used only for your analysis and deleted when you close the page (at the latest after 6 hours). Upload de-identified data only. No tracking cookies.",
        "Use", "Research use only; not for diagnosis or treatment decisions. Provided as is, without warranty.",
        "Run locally", tagList(tags$code("Rscript run_local.R", .noWS = "outside"), ": the same analyses without upload limits."),
        "Version", sprintf("%s (%s)", CMI_VERSION, CMI_RELEASE_DATE))))
