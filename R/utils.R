`%||%` <- function(a, b) if (is.null(a)) b else a

dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  invisible(path)
}

#' @keywords internal
outputs_exist <- function(paths) {
  paths <- unlist(paths, use.names = FALSE)
  length(paths) > 0 && all(file.exists(paths)) &&
    all(file.info(paths)$size > 0, na.rm = FALSE)
}

assert_file <- function(path, what = "file") {
  if (is.null(path) || is.na(path) || !nzchar(path)) {
    stop(sprintf("no %s given", what), call. = FALSE)
  }
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", what, path), call. = FALSE)
  }
  invisible(normalizePath(path))
}

is_fastq <- function(path) {
  grepl("\\.(fastq|fq)(\\.gz)?$", path, ignore.case = TRUE)
}

is_signal_input <- function(path) {
  dir.exists(path) || grepl("\\.(pod5|fast5)$", path, ignore.case = TRUE)
}

#' Is a CUDA GPU visible on this machine?
#'
#' Used to decide whether GPU-dependent steps (basecalling) can run.
#'
#' @return `TRUE` if `nvidia-smi` is on the `PATH` and exits cleanly.
#' @export
has_gpu <- function() {
  smi <- Sys.which("nvidia-smi")
  if (!nzchar(smi)) return(FALSE)
  status <- suppressWarnings(
    system2(smi, stdout = FALSE, stderr = FALSE, timeout = 20))
  identical(status, 0L)
}
