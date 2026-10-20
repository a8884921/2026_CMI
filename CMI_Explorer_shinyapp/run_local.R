# Local edition: files up to 4 GB, no cell, count, spot or value limits, two analyses at a time, and RCTD
# deconvolution in the Spatial workspace when the spacexr package is installed.
#   Rscript run_local.R          (or source("run_local.R") in RStudio)
# The app folder is found from this file's own location, so the command works from any working directory.
Sys.setenv(CMI_PROFILE = "local")

# Rscript writes a space in the --file argument as "~+~" (R's own quoting), so it has to be decoded.
app_dir_of_file <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", f[1]), fixed = TRUE), mustWork = FALSE)))
  for (i in rev(seq_len(sys.nframe()))) {                    # sourced (R / RStudio): the source file is in the call stack
    o <- sys.frame(i)$ofile
    if (!is.null(o) && length(o) && nzchar(o[1])) return(dirname(normalizePath(o[1], mustWork = FALSE)))
  }
  getwd()
}

app_dir <- app_dir_of_file()
if (!file.exists(file.path(app_dir, "app.R"))) app_dir <- getwd()
if (!file.exists(file.path(app_dir, "app.R"))) stop("Run this from the CMI Explorer folder (the one that contains app.R).", call. = FALSE)
shiny::runApp(app_dir, launch.browser = TRUE)
