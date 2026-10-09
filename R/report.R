# Visualization and reporting: per-run HTML report (R Markdown), an IGV
# session file, and igv-reports snapshots of user-listed regions.
# (EPI2ME Labs reports are a layout reference, not a dependency.)

#' Render the per-run HTML report
#'
#' Gathers QC, alignment, assembly, phasing and benchmark results across
#' all samples of a run into one self-contained HTML file.
#'
#' @param run A `nanoflow_run` object returned by [run_pipeline()], or the
#'   path to the `provenance.rds` a run writes into its output directory.
#' @param out_html Output HTML path; default `report.html` inside the run's
#'   output directory.
#' @param title Report title.
#' @return A `nanoflow_step` with `outputs$html`.
#' @export
render_report <- function(run, out_html = NULL, title = "nanoflow run report") {
  if (is.character(run)) {
    assert_file(run, "provenance RDS")
    run <- readRDS(run)
  }
  stopifnot(inherits(run, "nanoflow_run"))
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("rendering the report requires the 'rmarkdown' package", call. = FALSE)
  }
  out_html <- out_html %||% file.path(run$output_dir, "report.html")
  template <- system.file("report", "report_template.Rmd",
                          package = "nanoflow", mustWork = TRUE)
  t0 <- Sys.time()
  # Render from a copy inside the output dir so intermediate files land there.
  tmp_rmd <- file.path(dirname(out_html), ".report_template.Rmd")
  file.copy(template, tmp_rmd, overwrite = TRUE)
  on.exit(unlink(tmp_rmd), add = TRUE)
  rmarkdown::render(tmp_rmd,
                    output_file = basename(out_html),
                    output_dir = dirname(out_html),
                    params = list(run = run, title = title),
                    quiet = TRUE, envir = new.env(parent = globalenv()))
  new_step("report", "rmarkdown", "rmarkdown::render(report_template.Rmd)",
           list(html = out_html), list(title = title),
           as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

#' Write an IGV session file (XML) for a sample
#'
#' Produces a session that loads the reference plus the given tracks
#' (BAMs, VCFs, GFF3, ...) with relative paths, so the session stays valid
#' when the whole results directory is copied to another machine.
#'
#' @param reference Reference genome FASTA.
#' @param tracks Character vector of track files (BAM/VCF/GFF3/BED).
#' @param out_xml Output session path; track paths are stored relative to
#'   its directory when possible.
#' @param locus Optional initial locus, e.g. `"chr1:10000-20000"`.
#' @return A `nanoflow_step` with `outputs$session`.
#' @export
write_igv_session <- function(reference, tracks, out_xml,
                              locus = NULL) {
  assert_file(reference, "reference FASTA")
  dir_create(dirname(out_xml))
  base <- normalizePath(dirname(out_xml))
  relify <- function(p) {
    p <- normalizePath(p, mustWork = FALSE)
    if (startsWith(p, paste0(base, "/"))) sub(paste0(base, "/"), "", p,
                                              fixed = TRUE) else p
  }
  resources <- paste(sprintf('    <Resource path="%s"/>',
                             vapply(tracks, relify, "")), collapse = "\n")
  xml <- c(
    '<?xml version="1.0" encoding="UTF-8" standalone="no"?>',
    sprintf('<Session genome="%s"%s version="8">',
            normalizePath(reference),
            if (!is.null(locus)) sprintf(' locus="%s"', locus) else ""),
    "  <Resources>", resources, "  </Resources>",
    "</Session>")
  writeLines(xml, out_xml)
  new_step("igv_session", "igv", "write_igv_session()",
           list(session = out_xml), list(), 0)
}

#' Static IGV snapshots of selected regions with igv-reports
#'
#' Wraps `create_report` (igv-reports) to produce a self-contained HTML with
#' embedded browser views of the listed regions.
#'
#' @param regions Character vector of regions (`"chr1:100-200"`) or a BED
#'   file path.
#' @param reference Reference genome FASTA.
#' @param tracks Character vector of track files (BAM/VCF/...).
#' @param out_html Output HTML path.
#' @param flanking Bases of context around each region.
#' @param extra_args Extra command-line arguments appended verbatim.
#' @param config Optional nanoflow config (binary overrides).
#' @return A `nanoflow_step` with `outputs$html`.
#' @export
igv_snapshots <- function(regions, reference, tracks, out_html,
                          flanking = 1000, extra_args = character(),
                          config = NULL) {
  assert_file(reference, "reference FASTA")
  dir_create(dirname(out_html))
  sites <- regions
  if (length(regions) > 1 || !file.exists(regions[1])) {
    # turn "chr:start-end" strings into a temporary BED
    sites <- file.path(dirname(out_html), "igv_regions.bed")
    parts <- regmatches(regions,
                        regexec("^([^:]+):([0-9]+)-([0-9]+)$", regions))
    bad <- vapply(parts, length, 1L) != 4
    if (any(bad)) {
      stop("malformed region(s): ", paste(regions[bad], collapse = ", "),
           call. = FALSE)
    }
    writeLines(vapply(parts, function(p) {
      sprintf("%s\t%d\t%s", p[2], max(0, as.integer(p[3]) - 1), p[4])
    }, ""), sites)
  }
  log <- paste0(out_html, ".log")
  res <- nf_run(nf_bin("igv_reports", config), c(
    sites, "--fasta", reference,
    "--flanking", flanking,
    "--tracks", tracks,
    "--output", out_html, extra_args), log = log)
  new_step("igv_snapshots", "igv_reports", res$command,
           list(html = out_html), list(flanking = flanking),
           res$runtime, log)
}
