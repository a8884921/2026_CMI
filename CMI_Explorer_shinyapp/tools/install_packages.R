#!/usr/bin/env Rscript
# Installs the R packages CMI Explorer needs (R 4.3 or later), from CRAN.
#   Rscript tools/install_packages.R            # hosted edition (what shinyapps.io runs)
#   Rscript tools/install_packages.R --rctd     # also RCTD (spacexr) for the local edition, and SeuratObject so that
#                                                # the RCTD script can read Seurat references (tests/test_rctd_script.R)
repos <- getOption("repos"); if (is.null(repos) || identical(unname(repos["CRAN"]), "@CRAN@")) repos <- c(CRAN = "https://cloud.r-project.org")
pkgs <- c("shiny", "bslib", "DT", "htmltools", "jsonlite", "processx", "Matrix", "data.table", "matrixStats", "HGNChelper",
          "RANN", "ggplot2", "svglite", "ragg", "systemfonts", "patchwork", "zip", "digest", "ps", "rsconnect")
todo <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(todo)) install.packages(todo, repos = repos) else message("All packages are already installed.")
if ("--rctd" %in% commandArgs(trailingOnly = TRUE)) {
  # RCTD needs several GB of memory; it is used only by the local edition (Spatial > RCTD bundle).
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = repos)
  # The commit the original release (v1.2.0) used (spacexr 2.2.1); engine/spatial_engine.R patches its memory use.
  remotes::install_github("dmcable/spacexr@9f5dc33c8060f946c6072a138b70e189636e1435", build_vignettes = FALSE)
  if (!requireNamespace("SeuratObject", quietly = TRUE)) install.packages("SeuratObject", repos = repos)
  for (p in c("spacexr", "SeuratObject")) message(sprintf("%-12s %s", p, tryCatch(as.character(utils::packageVersion(p)), error = function(e) "MISSING")))
}
for (p in pkgs) message(sprintf("%-12s %s", p, tryCatch(as.character(utils::packageVersion(p)), error = function(e) "MISSING")))
