# Tool registry and environment checks. The registry maps the logical tool
# name used throughout nanoflow to the binary name, a version command and a
# regex extracting the version string.

#' @keywords internal
nf_tool_registry <- function() {
  v <- "[0-9]+\\.[0-9][0-9a-zA-Z.-]*"
  list(
    dorado      = list(bin = "dorado", args = "--version", pattern = v,
                       step = "basecall"),
    guppy       = list(bin = "guppy_basecaller", args = "--version",
                       pattern = v, step = "basecall"),
    nanoplot    = list(bin = "NanoPlot", args = "--version", pattern = v,
                       step = "qc"),
    pycoqc      = list(bin = "pycoQC", args = "--version", pattern = v,
                       step = "qc"),
    fastqc      = list(bin = "fastqc", args = "--version", pattern = v,
                       step = "qc"),
    porechop    = list(bin = "porechop", args = "--version", pattern = v,
                       step = "qc"),
    filtlong    = list(bin = "filtlong", args = "--version", pattern = v,
                       step = "qc"),
    multiqc     = list(bin = "multiqc", args = "--version", pattern = v,
                       step = "qc"),
    minimap2    = list(bin = "minimap2", args = "--version", pattern = v,
                       step = "align"),
    bwa         = list(bin = "bwa", args = character(), pattern = v,
                       step = "align"),
    samtools    = list(bin = "samtools", args = "--version", pattern = v,
                       step = "align"),
    qualimap    = list(bin = "qualimap", args = "--help", pattern = v,
                       step = "align"),
    flye        = list(bin = "flye", args = "--version", pattern = v,
                       step = "assembly"),
    canu        = list(bin = "canu", args = "--version", pattern = v,
                       step = "assembly"),
    wengan      = list(bin = "wengan.pl", args = character(), pattern = v,
                       step = "assembly"),
    racon       = list(bin = "racon", args = "--version", pattern = v,
                       step = "assembly"),
    medaka      = list(bin = "medaka", args = "--version", pattern = v,
                       step = "assembly"),
    quast       = list(bin = "quast.py", args = "--version", pattern = v,
                       step = "assembly"),
    clair3      = list(bin = "run_clair3.sh", args = "--version", pattern = v,
                       step = "small_variants"),
    sniffles    = list(bin = "sniffles", args = "--version", pattern = v,
                       step = "sv"),
    svim        = list(bin = "svim", args = "--version", pattern = v,
                       step = "sv"),
    nanovar     = list(bin = "nanovar", args = "--version", pattern = v,
                       step = "sv"),
    survivor    = list(bin = "SURVIVOR", args = character(), pattern = v,
                       step = "sv"),
    spectre     = list(bin = "spectre", args = "--version", pattern = v,
                       step = "sv"),
    whatshap    = list(bin = "whatshap", args = "--version", pattern = v,
                       step = "phase"),
    hapcut2     = list(bin = "HAPCUT2", args = character(), pattern = v,
                       step = "phase"),
    snpeff      = list(bin = "snpEff", args = "-version", pattern = v,
                       step = "annotate"),
    annotsv     = list(bin = "AnnotSV", args = "--version", pattern = v,
                       step = "annotate"),
    happy       = list(bin = "hap.py", args = "--version", pattern = v,
                       step = "benchmark"),
    truvari     = list(bin = "truvari", args = "version", pattern = v,
                       step = "benchmark"),
    igv_reports = list(bin = "create_report", args = "--version", pattern = v,
                       step = "report"),
    bgzip       = list(bin = "bgzip", args = "--version", pattern = v,
                       step = "core"),
    tabix       = list(bin = "tabix", args = "--version", pattern = v,
                       step = "core")
  )
}

#' Version string of an external tool, or NA if unavailable
#'
#' @param tool Logical tool name from the registry (see [check_tools()]).
#' @param config Optional nanoflow config (for `config$tools` binary
#'   overrides).
#' @return Character version (e.g. `"2.26"`) or `NA_character_`.
#' @export
tool_version <- function(tool, config = NULL) {
  entry <- nf_tool_registry()[[tool]]
  if (is.null(entry)) return(NA_character_)
  bin <- config$tools[[tool]] %||% Sys.which(entry$bin)
  if (!nzchar(bin) || !file.exists(bin)) return(NA_character_)
  out <- tryCatch(
    suppressWarnings(system2(bin, entry$args, stdout = TRUE, stderr = TRUE,
                             timeout = 30)),
    error = function(e) character())
  out <- iconv(out, sub = "")  # some tools emit non-UTF-8 in version banners
  m <- regmatches(out, regexpr(entry$pattern, out))
  m <- m[nzchar(m)]
  if (length(m)) m[1] else NA_character_
}

#' Check that every external tool is available before a run
#'
#' Verifies each binary (from `config$tools` overrides or the `PATH`) and
#' queries its version, reporting what is missing. Called by
#' [run_pipeline()] for the tools needed by the enabled steps; can also be
#' run standalone after installing `inst/conda/environment.yml`.
#'
#' @param config A nanoflow config (see [read_config()]); optional.
#' @param tools Character vector of logical tool names to check; default all.
#' @param stop_on_missing Stop (rather than warn) if any requested tool is
#'   missing.
#' @return Invisibly, a data.frame with columns `tool`, `binary`, `found`,
#'   `path`, `version`, printed as a checklist.
#' @examples
#' \donttest{
#' check_tools(tools = c("minimap2", "samtools"))
#' }
#' @export
check_tools <- function(config = NULL, tools = NULL, stop_on_missing = FALSE) {
  reg <- nf_tool_registry()
  tools <- tools %||% names(reg)
  unknown <- setdiff(tools, names(reg))
  if (length(unknown)) {
    stop("unknown tool(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  rows <- lapply(tools, function(nm) {
    entry <- reg[[nm]]
    bin <- config$tools[[nm]] %||% unname(Sys.which(entry$bin))
    found <- nzchar(bin) && file.exists(bin)
    data.frame(tool = nm, binary = entry$bin, found = found,
               path = if (found) bin else NA_character_,
               version = if (found) tool_version(nm, config) else NA_character_,
               stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, rows)
  class(res) <- c("nanoflow_tools", "data.frame")
  print(res)
  missing <- res$tool[!res$found]
  if (length(missing)) {
    msg <- paste0("missing tool(s): ", paste(missing, collapse = ", "),
                  "\nInstall them with the shipped conda environment ",
                  "(inst/conda/environment.yml) or point config$tools$<name> ",
                  "at the binary.")
    if (stop_on_missing) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  invisible(res)
}

#' @export
print.nanoflow_tools <- function(x, ...) {
  cat("nanoflow external tools:\n")
  for (i in seq_len(nrow(x))) {
    cat(sprintf("  %s %-12s %-18s %s\n",
                if (x$found[i]) "\u2713" else "\u2717",
                x$tool[i],
                if (x$found[i]) paste0("v", x$version[i] %||% "?") else "NOT FOUND",
                if (x$found[i]) x$path[i] else ""))
  }
  invisible(x)
}

# Tools required for the steps enabled in a config.
needed_tools <- function(config) {
  reg <- nf_tool_registry()
  steps <- c("core", "align",
             if (isTRUE(config$steps$qc$run)) "qc",
             if (isTRUE(config$steps$assembly$run)) "assembly",
             if (isTRUE(config$steps$small_variants$run)) "small_variants",
             if (isTRUE(config$steps$sv$run)) "sv",
             if (isTRUE(config$steps$phase$run)) "phase",
             if (isTRUE(config$steps$annotate$run)) "annotate")
  all_in_steps <- names(reg)[vapply(reg, function(e) e$step, "") %in% steps]
  # Only the *selected* tool per step plus always-needed helpers.
  selected <- c("samtools", "bgzip", "tabix",
                config$steps$align$tool %||% "minimap2",
                if (isTRUE(config$steps$qc$run))
                  c("nanoplot", "porechop", "filtlong", "multiqc"),
                if (isTRUE(config$steps$assembly$run))
                  config$steps$assembly$tool %||% "flye",
                if (isTRUE(config$steps$small_variants$run))
                  config$steps$small_variants$tool %||% "clair3",
                if (isTRUE(config$steps$sv$run))
                  config$steps$sv$tool %||% "sniffles",
                if (isTRUE(config$steps$phase$run))
                  config$steps$phase$tool %||% "whatshap",
                if (isTRUE(config$steps$annotate$run))
                  config$steps$annotate$tool %||% "snpeff")
  intersect(unique(selected), all_in_steps)
}
