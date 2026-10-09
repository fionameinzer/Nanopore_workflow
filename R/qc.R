# Read QC and trimming: NanoPlot, pycoQC, FastQC, Porechop, Filtlong, MultiQC.

#' Nanopore read QC with NanoPlot
#'
#' @param fastq Nanopore FASTQ(.gz).
#' @param out_dir Output directory for plots and `NanoStats.txt`.
#' @param threads CPU threads.
#' @param extra_args Extra command-line arguments appended verbatim.
#' @param config Optional nanoflow config (binary overrides).
#' @param overwrite Rerun even if outputs already exist.
#' @return A `nanoflow_step` with `outputs$stats` and `outputs$dir`.
#' @export
qc_nanoplot <- function(fastq, out_dir, threads = 4,
                        extra_args = character(), config = NULL,
                        overwrite = FALSE) {
  assert_file(fastq, "FASTQ")
  stats <- file.path(out_dir, "NanoStats.txt")
  skip <- skip_if_done("qc_nanoplot", "nanoplot",
                       list(stats = stats, dir = out_dir), list(), overwrite)
  if (!is.null(skip)) return(skip)
  dir_create(out_dir)
  log <- file.path(out_dir, "nanoplot.log")
  res <- nf_run(nf_bin("nanoplot", config),
                c("--fastq", fastq, "-o", out_dir, "-t", threads, extra_args),
                log = log)
  new_step("qc_nanoplot", "nanoplot", res$command,
           list(stats = stats, dir = out_dir), list(threads = threads),
           res$runtime, log)
}

#' Nanopore run QC with pycoQC (needs a sequencing summary file)
#'
#' Only applicable when the sequencing run's `sequencing_summary.txt` is
#' available (it is produced by the basecaller, not derivable from FASTQ).
#'
#' @param sequencing_summary ONT `sequencing_summary.txt`.
#' @inheritParams qc_nanoplot
#' @return A `nanoflow_step` with `outputs$html`.
#' @export
qc_pycoqc <- function(sequencing_summary, out_dir,
                      extra_args = character(), config = NULL,
                      overwrite = FALSE) {
  assert_file(sequencing_summary, "sequencing summary")
  dir_create(out_dir)
  html <- file.path(out_dir, "pycoqc.html")
  skip <- skip_if_done("qc_pycoqc", "pycoqc", list(html = html), list(),
                       overwrite)
  if (!is.null(skip)) return(skip)
  log <- file.path(out_dir, "pycoqc.log")
  res <- nf_run(nf_bin("pycoqc", config),
                c("-f", sequencing_summary, "-o", html, extra_args),
                log = log)
  new_step("qc_pycoqc", "pycoqc", res$command, list(html = html), list(),
           res$runtime, log)
}

#' Short/long read QC with FastQC
#'
#' Used on Illumina reads when present, and optionally on Nanopore FASTQs.
#'
#' @param fastqs Character vector of FASTQ files.
#' @inheritParams qc_nanoplot
#' @return A `nanoflow_step` with `outputs$dir`.
#' @export
qc_fastqc <- function(fastqs, out_dir, threads = 4,
                      extra_args = character(), config = NULL,
                      overwrite = FALSE) {
  for (f in fastqs) assert_file(f, "FASTQ")
  expected <- file.path(out_dir, paste0(
    sub("\\.(fastq|fq)(\\.gz)?$", "", basename(fastqs)), "_fastqc.html"))
  skip <- skip_if_done("qc_fastqc", "fastqc", list(html = expected), list(),
                       overwrite)
  if (!is.null(skip)) return(skip)
  dir_create(out_dir)
  log <- file.path(out_dir, "fastqc.log")
  res <- nf_run(nf_bin("fastqc", config),
                c("-o", out_dir, "-t", threads, extra_args, fastqs),
                log = log)
  new_step("qc_fastqc", "fastqc", res$command,
           list(html = expected, dir = out_dir), list(), res$runtime, log)
}

#' Adapter trimming with Porechop
#'
#' @param fastq Input Nanopore FASTQ(.gz).
#' @param out_fastq Output trimmed FASTQ.gz.
#' @inheritParams qc_nanoplot
#' @return A `nanoflow_step` with `outputs$fastq`.
#' @export
trim_porechop <- function(fastq, out_fastq, threads = 4,
                          extra_args = character(), config = NULL,
                          overwrite = FALSE) {
  assert_file(fastq, "FASTQ")
  skip <- skip_if_done("trim_porechop", "porechop",
                       list(fastq = out_fastq), list(), overwrite)
  if (!is.null(skip)) return(skip)
  dir_create(dirname(out_fastq))
  log <- paste0(out_fastq, ".porechop.log")
  res <- nf_run(nf_bin("porechop", config),
                c("-i", fastq, "-o", out_fastq, "-t", threads, extra_args),
                log = log)
  new_step("trim_porechop", "porechop", res$command,
           list(fastq = out_fastq), list(threads = threads), res$runtime, log)
}

#' Length/quality filtering with Filtlong
#'
#' @param fastq Input FASTQ(.gz).
#' @param out_fastq Output filtered FASTQ.gz.
#' @param min_length Minimum read length to keep.
#' @param keep_percent Keep this percentage of the best bases.
#' @param target_bases Optionally downsample to this many bases.
#' @inheritParams qc_nanoplot
#' @return A `nanoflow_step` with `outputs$fastq`.
#' @export
filter_filtlong <- function(fastq, out_fastq, min_length = 500,
                            keep_percent = 90, target_bases = NULL,
                            extra_args = character(), config = NULL,
                            overwrite = FALSE) {
  assert_file(fastq, "FASTQ")
  params <- list(min_length = min_length, keep_percent = keep_percent,
                 target_bases = target_bases)
  skip <- skip_if_done("filter_filtlong", "filtlong",
                       list(fastq = out_fastq), params, overwrite)
  if (!is.null(skip)) return(skip)
  dir_create(dirname(out_fastq))
  log <- paste0(out_fastq, ".filtlong.log")
  res <- nf_run_shell(paste(
    shQuote(nf_bin("filtlong", config)),
    "--min_length", min_length, "--keep_percent", keep_percent,
    if (!is.null(target_bases)) paste("--target_bases", target_bases),
    paste(extra_args, collapse = " "),
    shQuote(fastq), "| gzip >", shQuote(out_fastq)), log = log)
  new_step("filter_filtlong", "filtlong", res$command,
           list(fastq = out_fastq), params, res$runtime, log)
}

#' Aggregate QC reports with MultiQC
#'
#' @param input_dirs Directories to scan for tool reports.
#' @param out_dir Output directory.
#' @param name Report name (file stem).
#' @inheritParams qc_nanoplot
#' @return A `nanoflow_step` with `outputs$html`.
#' @export
qc_multiqc <- function(input_dirs, out_dir, name = "multiqc_report",
                       extra_args = character(), config = NULL,
                       overwrite = FALSE) {
  dir_create(out_dir)
  html <- file.path(out_dir, paste0(name, ".html"))
  log <- file.path(out_dir, "multiqc.log")
  res <- nf_run(nf_bin("multiqc", config),
                c(input_dirs, "-o", out_dir, "-n", name, "--force",
                  extra_args), log = log)
  new_step("qc_multiqc", "multiqc", res$command, list(html = html), list(),
           res$runtime, log)
}
