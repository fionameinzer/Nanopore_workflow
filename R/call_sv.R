# Structural variant calling: Sniffles2 (default), SVIM, NanoVar, optional
# SURVIVOR merge, and read-depth CNV calling with Spectre.

#' Call structural variants with Sniffles2
#'
#' @param bam Sorted, indexed Nanopore BAM (minimap2 `--MD` recommended).
#' @param reference Reference genome FASTA.
#' @param out_dir Output directory.
#' @param sample Sample ID.
#' @param all_contigs Pass `--all-contigs`. Sniffles >= 2.6 silently skips
#'   contigs shorter than 1 Mb by default; set this for small genomes,
#'   custom references and the synthetic test fixture.
#' @param tandem_repeats Optional tandem-repeat BED (recommended for human;
#'   improves calls in repeats).
#' @param threads CPU threads.
#' @param extra_args Extra command-line arguments appended verbatim.
#' @param config Optional nanoflow config (binary overrides).
#' @param overwrite Rerun even if outputs already exist.
#' @return A `nanoflow_step` with `outputs$vcf`.
#' @export
call_sniffles <- function(bam, reference, out_dir, sample = "sample",
                          all_contigs = FALSE, tandem_repeats = NULL,
                          threads = 4, extra_args = character(),
                          config = NULL, overwrite = FALSE) {
  assert_file(bam, "BAM")
  assert_file(reference, "reference FASTA")
  vcf <- file.path(out_dir, paste0(sample, ".sniffles.vcf"))
  params <- list(all_contigs = all_contigs)
  skip <- skip_if_done("sv", "sniffles", list(vcf = vcf), params, overwrite)
  if (!is.null(skip)) return(skip)
  dir_create(out_dir)
  log <- file.path(out_dir, "sniffles.log")
  res <- nf_run(nf_bin("sniffles", config), c(
    "--input", bam, "--vcf", vcf, "--reference", reference,
    "--threads", threads, "--sample-id", sample,
    if (all_contigs) "--all-contigs",
    if (!is.null(tandem_repeats)) c("--tandem-repeats", tandem_repeats),
    extra_args), log = log)
  new_step("sv", "sniffles", res$command, list(vcf = vcf), params,
           res$runtime, log)
}

#' Call structural variants with SVIM (alternative)
#'
#' @inheritParams call_sniffles
#' @param min_qual Filter the SVIM output to calls with `QUAL >=` this value
#'   (SVIM reports everything; 10 is a common working point).
#' @return A `nanoflow_step` with `outputs$vcf` (filtered) and
#'   `outputs$raw_vcf`.
#' @export
call_svim <- function(bam, reference, out_dir, sample = "sample",
                      min_qual = 10, extra_args = character(),
                      config = NULL, overwrite = FALSE) {
  assert_file(bam, "BAM")
  assert_file(reference, "reference FASTA")
  raw <- file.path(out_dir, "variants.vcf")
  vcf <- file.path(out_dir, paste0(sample, ".svim.q", min_qual, ".vcf"))
  params <- list(min_qual = min_qual)
  skip <- skip_if_done("sv", "svim", list(vcf = vcf, raw_vcf = raw), params,
                       overwrite)
  if (!is.null(skip)) return(skip)
  dir_create(out_dir)
  log <- file.path(out_dir, "svim.log")
  res <- nf_run(nf_bin("svim", config),
                c("alignment", out_dir, bam, reference,
                  "--sample", sample, extra_args), log = log)
  # QUAL filter; '.' QUAL lines (e.g. BNDs) are kept out by the comparison.
  nf_run_shell(paste(
    "awk -F'\t'", shQuote(sprintf(
      "/^#/ || ($6 != \".\" && $6 + 0 >= %d)", min_qual)),
    shQuote(raw), ">", shQuote(vcf)), log = log)
  new_step("sv", "svim", res$command, list(vcf = vcf, raw_vcf = raw), params,
           res$runtime, log)
}

#' Call structural variants with NanoVar (alternative)
#'
#' @inheritParams call_sniffles
#' @return A `nanoflow_step` with `outputs$vcf`.
#' @export
call_nanovar <- function(bam, reference, out_dir, sample = "sample",
                         threads = 4, extra_args = character(),
                         config = NULL, overwrite = FALSE) {
  assert_file(bam, "BAM")
  assert_file(reference, "reference FASTA")
  vcf <- file.path(out_dir, paste0(
    sub("\\.bam$", "", basename(bam)), ".nanovar.pass.vcf"))
  skip <- skip_if_done("sv", "nanovar", list(vcf = vcf), list(), overwrite)
  if (!is.null(skip)) return(skip)
  dir_create(out_dir)
  log <- file.path(out_dir, "nanovar.log")
  res <- nf_run(nf_bin("nanovar", config),
                c(bam, reference, out_dir, "-t", threads, extra_args),
                log = log)
  new_step("sv", "nanovar", res$command, list(vcf = vcf), list(),
           res$runtime, log)
}

#' Merge SV call sets from multiple callers with SURVIVOR
#'
#' @param vcfs Character vector of SV VCFs (from different callers, same
#'   sample).
#' @param out_vcf Output merged VCF path.
#' @param max_dist Maximum breakpoint distance for two calls to be merged.
#' @param min_callers Minimum number of callers supporting a merged call.
#' @param min_size Minimum SV size.
#' @inheritParams call_sniffles
#' @return A `nanoflow_step` with `outputs$vcf`.
#' @export
merge_survivor <- function(vcfs, out_vcf, max_dist = 1000, min_callers = 2,
                           min_size = 50, extra_args = character(),
                           config = NULL, overwrite = FALSE) {
  for (f in vcfs) assert_file(f, "SV VCF")
  params <- list(max_dist = max_dist, min_callers = min_callers,
                 min_size = min_size)
  skip <- skip_if_done("sv_merge", "survivor", list(vcf = out_vcf), params,
                       overwrite)
  if (!is.null(skip)) return(skip)
  dir_create(dirname(out_vcf))
  list_file <- paste0(out_vcf, ".input_list.txt")
  writeLines(vcfs, list_file)
  log <- paste0(out_vcf, ".survivor.log")
  # SURVIVOR merge <list> <max_dist> <min_callers> <type> <strand> <dup> <min_size> <out>
  res <- nf_run(nf_bin("survivor", config),
                c("merge", list_file, max_dist, min_callers, 1, 1, 0,
                  min_size, out_vcf), log = log)
  new_step("sv_merge", "survivor", res$command, list(vcf = out_vcf), params,
           res$runtime, log)
}

#' Read-depth CNV calling with Spectre
#'
#' Spectre consumes windowed coverage produced by `mosdepth`; this wrapper
#' runs mosdepth first (1 kb windows) and then Spectre's `CNVCaller`.
#'
#' @inheritParams call_sniffles
#' @param blacklist Optional BED of regions to exclude.
#' @return A `nanoflow_step` with `outputs$dir`.
#' @export
cnv_spectre <- function(bam, reference, out_dir, sample = "sample",
                        blacklist = NULL, threads = 4,
                        extra_args = character(), config = NULL,
                        overwrite = FALSE) {
  assert_file(bam, "BAM")
  assert_file(reference, "reference FASTA")
  dir_create(out_dir)
  cov_dir <- file.path(out_dir, "mosdepth")
  dir_create(cov_dir)
  log <- file.path(out_dir, "spectre.log")
  mosdepth <- Sys.which("mosdepth")
  if (!nzchar(mosdepth)) {
    stop("mosdepth (required by Spectre) not found on PATH", call. = FALSE)
  }
  nf_run(mosdepth, c("-t", threads, "-x", "-b", "1000", "-n",
                     file.path(cov_dir, sample), bam), log = log)
  res <- nf_run(nf_bin("spectre", config), c(
    "CNVCaller",
    "--coverage", cov_dir,
    "--sample-id", sample,
    "--output-dir", out_dir,
    "--reference", reference,
    if (!is.null(blacklist)) c("--blacklist", blacklist),
    extra_args), log = log)
  new_step("cnv", "spectre", res$command, list(dir = out_dir),
           list(windows = 1000), res$runtime, log)
}
