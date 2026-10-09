# Core execution machinery: every workflow step builds a command line, runs
# it via system2(), checks the exit status and returns a `nanoflow_step`.

#' Construct a workflow step result
#'
#' Every exported step wrapper returns an object of class `nanoflow_step`
#' holding the output paths, the tool and its version, the exact command
#' line, the parameters and the runtime. These objects are collected by
#' [run_pipeline()] into the provenance log.
#'
#' @param step Name of the workflow step (e.g. `"align"`).
#' @param tool Logical tool name (e.g. `"minimap2"`).
#' @param command Full command line(s) that were run.
#' @param outputs Named list of output paths.
#' @param params Named list of user-visible parameters.
#' @param runtime Elapsed wall-clock seconds.
#' @param log Path to the step log file, or `NA`.
#' @param status `"ok"`, `"skipped"` (outputs already existed or step not
#'   applicable) or `"failed"`.
#' @param message Optional human-readable note (e.g. why a step was skipped).
#' @return A `nanoflow_step` object.
#' @keywords internal
new_step <- function(step, tool, command = character(), outputs = list(),
                     params = list(), runtime = NA_real_, log = NA_character_,
                     status = "ok", message = NULL) {
  structure(list(
    step = step,
    tool = tool,
    tool_version = tool_version(tool),
    command = command,
    outputs = outputs,
    params = params,
    runtime_sec = runtime,
    log = log,
    status = status,
    message = message,
    finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  ), class = "nanoflow_step")
}

#' @export
print.nanoflow_step <- function(x, ...) {
  mark <- switch(x$status, ok = "\u2713", skipped = "\u21b7", "\u2717")
  cat(sprintf("<nanoflow step: %s> %s [%s, %s]\n", x$step, mark, x$status,
              paste(x$tool, x$tool_version %||% "")))
  if (!is.null(x$message)) cat("  note:   ", x$message, "\n", sep = "")
  if (length(x$outputs)) {
    cat("  outputs:\n")
    for (nm in names(x$outputs)) {
      cat(sprintf("    %-12s %s\n", nm, paste(x$outputs[[nm]], collapse = ", ")))
    }
  }
  if (!is.na(x$runtime_sec)) {
    cat(sprintf("  runtime: %.1fs", x$runtime_sec))
    if (!is.na(x$log)) cat("   log: ", x$log, sep = "")
    cat("\n")
  }
  invisible(x)
}

#' Resolve the binary for a logical tool name
#'
#' Looks first at `config$tools$<name>` (absolute-path override from the
#' config file), then on the `PATH`. This is the only place binaries are
#' resolved, so no path is ever hard-coded.
#'
#' @keywords internal
nf_bin <- function(tool, config = NULL) {
  override <- config$tools[[tool]]
  if (!is.null(override) && nzchar(override)) {
    if (file.exists(override)) return(override)
    stop(sprintf("configured binary for '%s' not found: %s", tool, override),
         call. = FALSE)
  }
  entry <- nf_tool_registry()[[tool]]
  bin <- Sys.which(entry$bin %||% tool)
  if (!nzchar(bin)) {
    stop(sprintf(
      "tool '%s' (binary '%s') not found on PATH and not set in config$tools.\nSee check_tools() and inst/conda/environment.yml.",
      tool, entry$bin %||% tool), call. = FALSE)
  }
  unname(bin)
}

#' Run an external command, check its exit status, time it
#'
#' @param bin Resolved binary path or name.
#' @param args Character vector of arguments (passed to [system2()], so no
#'   shell quoting is needed).
#' @param log File to capture stdout+stderr, `NULL` to pass through.
#' @param stdout_file Redirect stdout to this file instead of the log
#'   (for tools that write results to stdout).
#' @return List with `status`, `runtime`, `command`.
#' @keywords internal
nf_run <- function(bin, args, log = NULL, stdout_file = NULL) {
  cmd <- paste(c(bin, args), collapse = " ")
  message("[nanoflow] $ ", cmd)
  t0 <- Sys.time()
  stdout <- if (!is.null(stdout_file)) stdout_file else (log %||% "")
  stderr <- log %||% ""
  status <- suppressWarnings(system2(bin, args, stdout = stdout, stderr = stderr))
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (!identical(status, 0L)) {
    tail_log <- if (!is.null(log) && file.exists(log)) {
      paste(utils::tail(readLines(log, warn = FALSE), 15), collapse = "\n")
    } else ""
    stop(sprintf("command failed (exit %s):\n  %s\n%s", status, cmd, tail_log),
         call. = FALSE)
  }
  list(status = status, runtime = runtime, command = cmd)
}

#' Run a shell pipeline (e.g. `minimap2 ... | samtools sort ...`)
#'
#' Executed with `bash -c 'set -euo pipefail; ...'` so a failure anywhere in
#' the pipe fails the step. Paths interpolated into `cmd` must already be
#' shell-quoted (use [shQuote()]).
#'
#' @keywords internal
nf_run_shell <- function(cmd, log = NULL) {
  message("[nanoflow] $ ", cmd)
  t0 <- Sys.time()
  status <- suppressWarnings(system2(
    "bash", c("-c", shQuote(paste("set -euo pipefail;", cmd))),
    stdout = log %||% "", stderr = log %||% ""))
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (!identical(status, 0L)) {
    tail_log <- if (!is.null(log) && file.exists(log)) {
      paste(utils::tail(readLines(log, warn = FALSE), 15), collapse = "\n")
    } else ""
    stop(sprintf("command failed (exit %s):\n  %s\n%s", status, cmd, tail_log),
         call. = FALSE)
  }
  list(status = status, runtime = runtime, command = cmd)
}

# Standard resume prologue shared by all wrappers: if every expected output
# already exists (and is non-empty) and overwrite = FALSE, return a
# "skipped" step so run_pipeline() can resume cheaply.
skip_if_done <- function(step, tool, outputs, params, overwrite) {
  if (!overwrite && outputs_exist(outputs)) {
    return(new_step(step, tool, outputs = outputs, params = params,
                    status = "skipped",
                    message = "outputs already exist (resume); use overwrite = TRUE to rerun"))
  }
  NULL
}
