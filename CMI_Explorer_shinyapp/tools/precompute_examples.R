#!/usr/bin/env Rscript
# Runs the analysis worker on the three example inputs and stores the results in
# examples/precomputed/<module>/ (the result bundle plus the screen-resolution figures), so that
# "Show the example" answers at once without taking an analysis slot. "Run the example live" in the
# interface still computes them on the server. Run after any change to engine/, R/figures.R or the
# example inputs:
#   Rscript tools/precompute_examples.R .
# The example of each workspace (input file, whether it is synthetic) is
# defined once, in WS of R/workspace.R.
args <- commandArgs(trailingOnly = TRUE)
app <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
source(file.path(app, "engine", "constants.R"))
source(file.path(app, "R", "workspace.R"))
# Rscript of this R installation; paths are quoted (folders may contain spaces) and the environment is set
# here rather than with system2(env =), which Windows does not support.
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
work <- tempfile("precompute_"); dir.create(work)
ex <- function(m, extra = character()) list(input = WS[[m]]$example, synthetic = isTRUE(WS[[m]]$example_synthetic), extra = extra)
runs <- list(single_cell = ex("single_cell"), spatial = ex("spatial", c("--spatial-mode", "weights")), bulk = ex("bulk"))
env <- c(CMI_PROFILE = "hosted", CMI_LOCK_DIR = file.path(work, "locks"), LANGUAGE = "en", OMP_NUM_THREADS = "1")
if (!isTRUE(l10n_info()$`UTF-8`) && .Platform$OS.type != "windows") env <- c(env, LANG = "C.UTF-8", LC_ALL = "C.UTF-8")
saved <- Sys.getenv(names(env), unset = NA)
do.call(Sys.setenv, as.list(env))
for (module in names(runs)) {
  r <- runs[[module]]; out <- file.path(work, module)
  indir <- file.path(work, paste0(module, "_input")); dir.create(indir)
  inp <- file.path(indir, basename(r$input)); file.copy(file.path(app, r$input), inp)
  t0 <- Sys.time()
  st <- system2(rscript, c(shQuote(file.path(app, "engine", "worker.R")), "--root", shQuote(app), "--module", module, "--input", shQuote(inp),
                           "--output", shQuote(out), "--state", shQuote(file.path(work, paste0(module, "_state.json"))), r$extra, if (r$synthetic) "--synthetic"),
                stdout = file.path(work, paste0(module, ".log")), stderr = file.path(work, paste0(module, ".log")))
  if (!identical(st, 0L) || !file.exists(file.path(out, "result_bundle.zip"))) stop("The ", module, " example failed; see ", file.path(work, paste0(module, ".log")))
  dest <- file.path(app, "examples", "precomputed", module)
  unlink(dest, recursive = TRUE); dir.create(file.path(dest, "figures"), recursive = TRUE)
  file.copy(file.path(out, "result_bundle.zip"), dest)
  file.copy(list.files(file.path(out, "figures"), pattern = "_web\\.png$", full.names = TRUE), file.path(dest, "figures"))
  writeLines(jsonlite::toJSON(list(app_version = CMI_VERSION, input = r$input, synthetic = r$synthetic, computed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                                   seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)), auto_unbox = TRUE, pretty = TRUE),
             file.path(dest, "precomputed.json"))
  cat(sprintf("%-12s %5.1f s  %5.1f MB\n", module, as.numeric(difftime(Sys.time(), t0, units = "secs")),
              sum(file.size(list.files(dest, recursive = TRUE, full.names = TRUE))) / 1024^2))
}
unlink(work, recursive = TRUE)
for (k in names(saved)) if (is.na(saved[[k]])) Sys.unsetenv(k) else do.call(Sys.setenv, setNames(list(saved[[k]]), k))
