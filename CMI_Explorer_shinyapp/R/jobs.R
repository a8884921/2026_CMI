# Starting, polling and stopping analysis workers (engine/worker.R) from the interface process.

rscript_path <- function() {
  cand <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  if (file.exists(cand)) return(cand)
  found <- Sys.which("Rscript")
  if (!nzchar(found)) stop("Rscript is not available on this server.")
  unname(found)
}

# Job folders live in one directory shared by all interface processes of an instance, so that the
# sweeper below can remove folders left behind by sessions that ended without cleaning up.
jobs_root <- function() {
  d <- Sys.getenv("CMI_JOBS_DIR", file.path(dirname(tempdir()), "cmi_jobs"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
lock_root <- function() Sys.getenv("CMI_LOCK_DIR", file.path(dirname(tempdir()), "cmi_locks"))

# Each session touches "alive" in its job folders every few minutes (keepalive from the page).
touch_alive <- function(dirs) for (d in dirs) if (dir.exists(d)) writeLines(format(Sys.time()), file.path(d, "alive"))
# Remove job folders whose session has not been seen for `idle_min` minutes, and any older than `max_h`.
sweep_jobs <- function(idle_min = 30, max_h = 24) {
  now <- Sys.time()
  for (d in list.dirs(jobs_root(), recursive = FALSE)) {
    a <- file.path(d, "alive"); seen <- if (file.exists(a)) file.info(a)$mtime else file.info(d)$mtime
    age_h <- as.numeric(difftime(now, file.info(d)$ctime, units = "hours"))
    if (as.numeric(difftime(now, seen, units = "mins")) > idle_min || isTRUE(age_h > max_h)) unlink(d, recursive = TRUE)
  }
  invisible()
}

start_worker <- function(root, module, input_path, job_dir, synthetic = FALSE, extra = character()) {
  state <- file.path(job_dir, "state.json")
  out <- file.path(job_dir, "output")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  args <- c(file.path(root, "engine", "worker.R"), "--root", root, "--module", module, "--input", input_path,
            "--output", out, "--state", state, extra, if (synthetic) "--synthetic")
  env <- c("current", R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep), LANGUAGE = "en",
           OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
  if (!isTRUE(l10n_info()$`UTF-8`)) env <- c(env, LANG = "C.UTF-8", LC_ALL = "C.UTF-8")
  proc <- processx::process$new(rscript_path(), args = args, wd = root, env = env,
                                stdout = file.path(job_dir, "worker.out"), stderr = "2>&1", cleanup_tree = TRUE)
  touch_alive(job_dir)
  list(proc = proc, dir = job_dir, out = out, state = state, module = module, synthetic = synthetic, started = Sys.time())
}

read_state <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(jsonlite::fromJSON(path, simplifyVector = TRUE), error = function(e) NULL)
}

# Why a worker ended without a success state: killed for memory (SIGKILL), or another error.
exit_reason <- function(job) {
  st <- tryCatch(job$proc$get_exit_status(), error = function(e) NA_integer_)
  if (isTRUE(st == -9L)) return("The analysis ran out of memory on this shared server (about 0.7 GB is available per analysis). Upload fewer cells or samples, or use the local edition.")
  "The analysis stopped unexpectedly. Please try again; if it happens again, report the error ID below."
}

# A cancelled worker cannot release its queue slot or ticket; free anything owned by that process.
release_slots_of <- function(pid) {
  lr <- lock_root()
  for (d in list.dirs(lr, recursive = FALSE)) {
    if (basename(d) == "queue") next
    own <- tryCatch(readLines(file.path(d, "owner"), warn = FALSE), error = function(e) character())
    if (length(own) && identical(own[1], as.character(pid))) unlink(d, recursive = TRUE)
  }
  unlink(list.files(file.path(lr, "queue"), pattern = paste0("_", pid, "$"), full.names = TRUE))
}

stop_worker <- function(job) {
  if (is.null(job)) return(invisible())
  pid <- tryCatch(job$proc$get_pid(), error = function(e) NA)
  if (job$proc$is_alive()) { try(job$proc$kill_tree(), silent = TRUE); try(job$proc$wait(3000), silent = TRUE) }
  if (!is.na(pid)) release_slots_of(pid)
  invisible()
}

read_csv_safe <- function(path, key = NULL) {
  if (!file.exists(path)) return(NULL)
  cc <- if (!is.null(key)) stats::setNames("character", key) else NA
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = cc, encoding = "UTF-8")
}
