#!/usr/bin/env Rscript
# =============================================================================
# CMI Explorer - integrity check for work on a local copy of the package
# =============================================================================
# 1. Locked files are byte-identical to this release (SHA-256 in tools/locked_files.sha256, which must
#    list every file named in LOCKED below): the model parameters and reference files in ref/, the
#    analysis engines, the downloadable scripts and the RCTD test data, the example files the tests
#    read, the v1.2.0 reference outputs, the tests and every tool.
#    engine/constants.R is compared without CMI_VERSION and CMI_RELEASE_DATE, so the cut points and
#    every other frozen constant are covered.
# 2. R sources are ASCII (non-ASCII characters written as \uXXXX escapes) and parse.
# 3. ref/SHA256SUMS lists every file in ref/ with its current checksum.
# 4. No access token or secret is written into a file of the package, and the
#    package's .Renviron, which the deployment uploads, holds only CMI_* settings.
#   Rscript tools/check_integrity.R [app_dir]
#   Rscript tools/check_integrity.R [app_dir] --write-locks   (authors only, after an approved change)
# Exit status 0 = every check passed.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
write_locks <- "--write-locks" %in% args; args <- setdiff(args, "--write-locks")
app <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
if (!file.exists(file.path(app, "app.R"))) stop("Run this from the CMI Explorer folder (the one that contains app.R).", call. = FALSE)
lock_file <- file.path(app, "tools", "locked_files.sha256")
sha_text <- function(x) digest::digest(paste(x, collapse = "\n"), algo = "sha256", serialize = FALSE)
sha_file <- function(p) digest::digest(file = p, algo = "sha256")

LOCKED <- c(
  file.path("ref", c("annotation_snapshot.rds", "asset_manifest.json", "ccat_ppi.rds", "cmi_gbm.rds",
                     "cnv_reference.rds", "gene_sets.csv", "ig_gbm.rds", "mtp_gbm.rds", "reo_pairs.csv",
                     "rsf_m014.rds", "state_cutpoints.csv")),
  file.path("engine", c("bulk_engine.R", "io.R", "models.R", "sc_engine.R", "spatial_engine.R", "worker.R")),
  file.path("www", "downloads", c("rctd_mal_states.R", "rctd_test_data.zip", "h5ad_to_mtx.py")),
  file.path("examples", c("single_cell/counts.csv", "single_cell/counts_matrix_market.zip", "bulk/expression.csv", "spatial/rctd_bundle_demo.zip",
                          "spatial/spot_weights_demo.csv")),
  file.path("tests", c("test_equivalence.R", "test_rctd_script.R", "reference/README.md", "reference/v1.2.0_bulk_results.csv",
                       "reference/v1.2.0_single_cell_results.csv", "reference/v1.2.0_spatial_rctd_results.csv")),
  file.path("tools", c("check_integrity.R", "deploy_shinyapps.R", "export_gene_sets.R", "install_packages.R", "make_examples.R",
                       "make_rctd_test_data.R", "precompute_examples.R")))
CONSTANTS <- "engine/constants.R#frozen"

# engine/constants.R without the version and release date.
frozen_constants <- function() {
  x <- readLines(file.path(app, "engine", "constants.R"), warn = FALSE)
  x[!grepl("^CMI_(VERSION|RELEASE_DATE) <- ", x)]
}
current <- function() c(vapply(LOCKED, function(f) if (file.exists(file.path(app, f))) sha_file(file.path(app, f)) else "missing", character(1)),
                        setNames(sha_text(frozen_constants()), CONSTANTS))
write_lf <- function(text, path) { con <- file(path, open = "wb"); on.exit(close(con)); writeLines(text, con) }

if (write_locks) {
  cur <- current()
  if (any(cur == "missing")) stop("Missing: ", paste(names(cur)[cur == "missing"], collapse = ", "), call. = FALSE)
  write_lf(c("# SHA-256 of the files that local work must not change (tools/check_integrity.R). Do not edit by hand.",
             "# Rewrite only after the authors have approved a change: Rscript tools/check_integrity.R . --write-locks",
             paste0(cur, "  ", names(cur))), lock_file)
  cat("Wrote", length(cur), "checksums to tools/locked_files.sha256\n"); quit(save = "no", status = 0L)
}

problems <- character(); say <- function(ok, label, detail = "") {
  cat(sprintf("%-4s  %-34s %s\n", if (ok) "OK" else "FAIL", label, detail))
  if (!ok) problems <<- c(problems, label)
}

# ---- 1. locked files ---------------------------------------------------------------------------------
if (!file.exists(lock_file)) stop("tools/locked_files.sha256 is missing.", call. = FALSE)
lk <- grep("^#", readLines(lock_file, warn = FALSE), value = TRUE, invert = TRUE)
lk <- lk[nzchar(lk)]
locked <- setNames(sub("  .*$", "", lk), sub("^[0-9a-f]+  ", "", lk))
cur <- current()
changed <- character()
for (f in setdiff(names(cur), names(locked))) changed <- c(changed, sprintf("%s: no entry in tools/locked_files.sha256 (the lock list was edited)", f))
for (f in names(locked)) {
  now <- if (f %in% names(cur)) cur[[f]] else if (file.exists(file.path(app, f))) sha_file(file.path(app, f)) else "missing"
  if (identical(now, locked[[f]])) next
  why <- if (now == "missing") "missing" else "changed"
  p <- file.path(app, f)
  if (why == "changed" && file.exists(p) && grepl("\\.(R|py|csv|json|md|txt)$", f)) {
    lf <- gsub("\r\n", "\n", readChar(p, file.size(p), useBytes = TRUE), fixed = TRUE)
    if (identical(digest::digest(lf, algo = "sha256", serialize = FALSE), locked[[f]])) why <- "line endings changed to CRLF (set git core.autocrlf false and restore the file)"
  }
  if (f == CONSTANTS) why <- "a frozen constant or definition in engine/constants.R changed (only CMI_VERSION and CMI_RELEASE_DATE may be edited)"
  changed <- c(changed, sprintf("%s: %s", f, why))
}
say(!length(changed), "Locked files unchanged", sprintf("%d files and the frozen part of engine/constants.R", length(LOCKED)))
if (length(changed)) cat(paste0("        ", changed, collapse = "\n"), "\n")

# ---- 2. R sources: ASCII and parseable -----------------------------------------------------------------
src <- unique(c("app.R", "run_local.R",
                file.path("R", list.files(file.path(app, "R"), "\\.R$")), file.path("engine", list.files(file.path(app, "engine"), "\\.R$")),
                file.path("tools", list.files(file.path(app, "tools"), "\\.R$")), file.path("tests", list.files(file.path(app, "tests"), "\\.R$")),
                file.path("www", "downloads", list.files(file.path(app, "www", "downloads"), "\\.R$"))))
bad_ascii <- character(); bad_parse <- character()
for (f in src) {
  p <- file.path(app, f)
  raw <- readBin(p, "raw", file.size(p))
  if (any(raw > as.raw(0x7F))) {
    lines <- readLines(p, warn = FALSE, encoding = "bytes")
    bad_ascii <- c(bad_ascii, sprintf("%s line %s", f, paste(head(which(grepl("[^\\x01-\\x7F]", lines, perl = TRUE, useBytes = TRUE)), 5), collapse = ", ")))
  }
  e <- tryCatch({ parse(p, keep.source = FALSE); "" }, error = conditionMessage)
  if (nzchar(e)) bad_parse <- c(bad_parse, sprintf("%s: %s", f, e))
}
say(!length(bad_ascii), "R sources ASCII only", if (length(bad_ascii)) "write non-ASCII characters as \\uXXXX escapes" else sprintf("%d files", length(src)))
if (length(bad_ascii)) cat(paste0("        ", bad_ascii, collapse = "\n"), "\n")
say(!length(bad_parse), "R sources parse")
if (length(bad_parse)) cat(paste0("        ", bad_parse, collapse = "\n"), "\n")

# ---- 3. ref/SHA256SUMS -----------------------------------------------------------------------------------
sums_path <- file.path(app, "ref", "SHA256SUMS")
sums <- if (file.exists(sums_path)) utils::read.table(sums_path, col.names = c("sha256", "file"), stringsAsFactors = FALSE) else data.frame(sha256 = character(), file = character())
ref_files <- sort(setdiff(list.files(file.path(app, "ref")), "SHA256SUMS"))
stale <- c(setdiff(ref_files, sums$file),
           sums$file[!vapply(seq_len(nrow(sums)), function(i) file.exists(file.path(app, "ref", sums$file[i])) &&
                                                    identical(sha_file(file.path(app, "ref", sums$file[i])), sums$sha256[i]), logical(1))])
crlf <- file.exists(sums_path) && any(readBin(sums_path, "raw", file.size(sums_path)) == as.raw(0x0D))
say(!length(stale) && !crlf, "ref/SHA256SUMS up to date",
    if (length(stale) || crlf) paste0(if (length(stale)) paste("stale:", paste(unique(stale), collapse = ", ")) else "CRLF line endings",
                                      " - run: Rscript tools/export_gene_sets.R .")
    else sprintf("%d files", length(ref_files)))

# ---- 4. no secrets in files; .Renviron holds only CMI_* settings ----------------------------------------
txt <- list.files(app, pattern = "\\.(R|r|py|md|txt|json|csv|tsv|sh|yml|yaml|dcf|log|env|Renviron|Rprofile)$|^\\.(Renviron|Rprofile|env)$",
                  recursive = TRUE, all.files = TRUE, full.names = FALSE)
txt <- txt[!grepl("^(renv|rsconnect)/|^\\.git/", txt) & file.size(file.path(app, txt)) < 5e6]
# quoted assignments and known token formats (any case); env-style lines whose upper-case name ends in
# TOKEN, SECRET, _PAT, PASSWORD or API_KEY, such as GITHUB_PAT=... (DATA_PATH=... or out_path = "..." are not tokens)
pat <- c(quoted = "(?i)(token|secret|password|api[_-]?key)[\"']?\\s*[=:]\\s*[\"'][A-Za-z0-9+/_=-]{16,}[\"']",
         github = "gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{20,}",
         openai = "sk-[A-Za-z0-9]{20,}",
         env = "^\\s*(export\\s+)?[A-Z0-9_]*(TOKEN|SECRET|_PAT|PASSWORD|API_?KEY)\\s*=\\s*[\"']?[A-Za-z0-9+/_=-]{20,}[\"']?\\s*$")
leaks <- character()
for (f in txt) {
  lines <- readLines(file.path(app, f), warn = FALSE, encoding = "bytes")
  hit <- sort(unique(unlist(lapply(pat, function(p) grep(p, lines, perl = TRUE, useBytes = TRUE)))))
  if (length(hit)) leaks <- c(leaks, sprintf("%s line %s", f, paste(head(hit, 3), collapse = ", ")))
}
renv_file <- file.path(app, ".Renviron")
if (file.exists(renv_file)) {
  rl <- trimws(readLines(renv_file, warn = FALSE)); rl <- rl[nzchar(rl) & !startsWith(rl, "#")]
  if (any(!grepl("^CMI_[A-Z0-9_]+\\s*=", rl))) leaks <- c(leaks, ".Renviron: only CMI_* settings belong here (the deployment uploads this file); keep GITHUB_PAT and other credentials in the .Renviron of your home folder")
}
say(!length(leaks), "No tokens or secrets in files", if (length(leaks)) "remove them (rsconnect and git keep their credentials outside the package)" else sprintf("%d text files", length(txt)))
if (length(leaks)) cat(paste0("        ", leaks, collapse = "\n"), "\n")

cat(if (length(problems)) sprintf("\n%d problem(s): %s\n", length(problems), paste(problems, collapse = "; ")) else "\nAll integrity checks passed\n")
quit(save = "no", status = if (length(problems)) 1L else 0L)
