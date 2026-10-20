#!/usr/bin/env Rscript
# =============================================================================
# Deploy CMI Explorer to shinyapps.io
# =============================================================================
# One-time setup (token and secret from shinyapps.io > Account > Tokens):
#   install.packages("rsconnect")
#   rsconnect::setAccountInfo(name = "<account>", token = "<token>", secret = "<secret>")
# Then, from the app folder:
#   Rscript tools/deploy_shinyapps.R <account> [app-name]
#   Rscript tools/deploy_shinyapps.R --dry-run      (all checks and the dependency scan, nothing uploaded)
#
# Only the files the running site needs are uploaded (about 38 MB): app.R, R/, engine/, ref/,
# www/ and examples/ (inputs and precomputed results). tools/, tests/ and docs/ stay on your computer.
# The script refuses to deploy while the precomputed examples are stale.
# Deploy from a library of CRAN packages (not distribution packages such as r-cran-* on Linux):
# shinyapps.io reinstalls every package from the repository recorded in its DESCRIPTION.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
args <- setdiff(args, c("--allow-placeholders", "--dry-run"))
if (!length(args) && !dry_run) stop("Usage: Rscript tools/deploy_shinyapps.R <account> [app-name] | --dry-run", call. = FALSE)
account <- if (length(args)) args[1] else "dry-run"
app_name <- if (length(args) > 1) args[2] else "cmi-explorer"
if (!requireNamespace("rsconnect", quietly = TRUE)) stop("Install rsconnect first: install.packages(\"rsconnect\")", call. = FALSE)

app_dir <- normalizePath(".", mustWork = TRUE)
if (!file.exists(file.path(app_dir, "app.R"))) stop("Run this script from the CMI Explorer folder (the one that contains app.R).", call. = FALSE)

# Runtime packages (the analysis worker runs in a separate R process, so list them explicitly).
needed <- c("shiny", "bslib", "DT", "htmltools", "jsonlite", "processx", "Matrix", "data.table", "matrixStats", "HGNChelper",
            "RANN", "ggplot2", "svglite", "ragg", "systemfonts", "patchwork", "zip", "digest", "ps")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install these packages before deploying: ", paste(missing, collapse = ", "),
                          "\n  Rscript tools/install_packages.R", call. = FALSE)

source(file.path(app_dir, "engine", "constants.R"))
# The precomputed examples must come from this version of the software.
pre <- list.files(file.path(app_dir, "examples", "precomputed"), pattern = "^precomputed\\.json$", recursive = TRUE, full.names = TRUE)
stale <- pre[vapply(pre, function(f) !identical(jsonlite::fromJSON(f)$app_version, CMI_VERSION), logical(1))]
if (length(pre) < 3 || length(stale)) stop("The precomputed examples are missing or from another version; run: Rscript tools/precompute_examples.R .", call. = FALSE)

# Check the model files against their manifest before anything leaves the computer.
man <- jsonlite::fromJSON(file.path(app_dir, "ref", "asset_manifest.json"), simplifyVector = FALSE)
bad <- names(man$assets)[!vapply(names(man$assets), function(f) {
  p <- file.path(app_dir, "ref", f); file.exists(p) && identical(digest::digest(file = p, algo = "sha256"), man$assets[[f]]$sha256) }, logical(1))]
if (length(bad)) stop("These model files do not match ref/asset_manifest.json: ", paste(bad, collapse = ", "), call. = FALSE)
# Locked files unchanged, R sources ASCII, ref/SHA256SUMS current.
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
if (!identical(system2(rscript, c(shQuote(file.path(app_dir, "tools", "check_integrity.R")), shQuote(app_dir))), 0L))
  stop("tools/check_integrity.R reported problems (see above); fix them before deploying.", call. = FALSE)

within <- function(d) list.files(file.path(app_dir, d), recursive = TRUE, all.files = FALSE, no.. = TRUE)
files <- c("app.R",
           file.path("R", within("R")), file.path("engine", within("engine")), file.path("ref", within("ref")),
           file.path("www", within("www")), file.path("examples", within("examples")))
files <- files[!grepl("(^|/)\\.DS_Store$|~$", files)]
if (file.exists(file.path(app_dir, ".Renviron"))) files <- c(files, ".Renviron")   # optional upload limits (docs/DEPLOYMENT_zh.md)
# www/downloads/ holds scripts that users run on their own computer (the RCTD script needs spacexr);
# .renvignore keeps them out of the dependency scan, so the server installs only the site's packages.
if (file.exists(file.path(app_dir, ".renvignore"))) files <- c(files, ".renvignore")
size_mb <- sum(file.size(file.path(app_dir, files))) / 1024^2
cat(sprintf("%s %d files (%.1f MB)\n", if (dry_run) "Would deploy" else "Checking", length(files), size_mb))
# The packages the server will install. Analysis packages used on the authors' computers (an R script
# that loads them in a deployed folder makes them dependencies) and packages from outside CRAN would
# make the deployment fail or exceed the instance's memory.
deps <- rsconnect::appDependencies(appDir = app_dir, appFiles = files)
cat(sprintf("Dependencies: %d packages; sources: %s\n", nrow(deps), paste(unique(deps$Source), collapse = ", ")))
cat("Runtime packages detected:", sum(needed %in% deps$Package), "of", length(needed), "\n")
not_site <- c("spacexr", "SeuratObject", "Seurat", "h2o", "GEOquery", "Biobase", "BiocGenerics", "BiocManager", "celda", "scDblFinder",
              "SingleCellExperiment", "SummarizedExperiment", "hdf5r", "reticulate", "anndata", "zellkonverter", "Azimuth", "UCell",
              "caret", "gbm", "randomForestSRC", "remotes", "devtools")
extra <- intersect(not_site, deps$Package)
if (length(extra)) stop("These packages must not be dependencies of the site: ", paste(extra, collapse = ", "),
                        "\n  An R script that loads them is in a deployed folder (R/, engine/, ref/, www/, examples/); move it out of the deployed folders.", call. = FALSE)
odd <- !tolower(deps$Source) %in% c("cran", "rspm")
if (any(odd)) stop("These dependencies were not installed from CRAN: ", paste(sprintf("%s (%s)", deps$Package[odd], deps$Source[odd]), collapse = ", "),
                   "\n  Reinstall them from CRAN (Rscript tools/install_packages.R) before deploying.", call. = FALSE)
if (dry_run) quit(status = 0)
cat(sprintf("Deploying %d files (%.1f MB) as '%s' to account '%s'\n", length(files), size_mb, app_name, account))

rsconnect::deployApp(appDir = app_dir, appFiles = files, appName = app_name, appTitle = "CMI Explorer",
                     account = account, server = "shinyapps.io", forceUpdate = TRUE, launch.browser = FALSE)

# 1 GB instance ("large") is the size the hosted limits are designed for; larger paid sizes also work.
ok <- !inherits(try(rsconnect::configureApp(app_name, account = account, server = "shinyapps.io", size = "large"), silent = TRUE), "try-error")
cat(if (ok) "Instance size set to large (1 GB).\n" else "Set the instance size to Large (1 GB) in the shinyapps.io dashboard (Settings > General).\n")
cat("In the shinyapps.io dashboard, Settings > Advanced, set:\n",
    "  Max Worker Processes = 1   (one interface process per instance keeps one analysis plus the interface below 1 GB)\n",
    "  Connection Timeout = 900 s, Read Timeout = 3600 s   (analyses may take a few minutes when queued)\n",
    "  Instance Idle Timeout = 15 min or more\n", sep = "")
cat("Done. Open https://", account, ".shinyapps.io/", app_name, "/\n", sep = "")
