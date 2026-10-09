# Haplotype phasing. Runs after small-variant calling because it needs a
# VCF plus a BAM. Produces a phased VCF, a haplotagged BAM, phasing stats,
# and -- when a truth VCF is available -- switch-error rates.

# bgzip + tabix a VCF in place; returns the .gz path.
bgzip_tabix <- function(vcf, config = NULL) {
  gz <- if (grepl("\\.gz$", vcf)) vcf else paste0(vcf, ".gz")
  if (!grepl("\\.gz$", vcf)) {
    nf_run(nf_bin("bgzip", config), c("-f", vcf))
  }
  nf_run(nf_bin("tabix", config), c("-f", "-p", "vcf", gz))
  gz
}

#' Phase variants with WhatsHap
#'
#' Runs `whatshap phase`, computes phasing statistics (`whatshap stats`),
#' haplotags the BAM (`whatshap haplotag`), and -- if a truth VCF is given --
#' computes switch error rates with `whatshap compare`.
#'
#' @param vcf Small-variant VCF(.gz) to phase.
#' @param bam Sorted, indexed Nanopore BAM.
#' @param reference Reference genome FASTA.
#' @param out_dir Output directory.
#' @param sample Sample ID.
#' @param truth_vcf Optional phased truth VCF for switch-error evaluation.
#' @param ignore_read_groups Pass `--ignore-read-groups` (default TRUE; the
#'   usual single-sample ONT setting, and robust to caller/BAM sample-name
#'   mismatches).
#' @param extra_args Extra arguments for `whatshap phase`.
#' @param config Optional nanoflow config (binary overrides).
#' @param overwrite Rerun even if outputs already exist.
#' @return A `nanoflow_step` with `outputs$phased_vcf`,
#'   `outputs$haplotagged_bam`, `outputs$stats_tsv` and, when truth is
#'   given, `outputs$compare_tsv`.
#' @export
phase_whatshap <- function(vcf, bam, reference, out_dir, sample = "sample",
                           truth_vcf = NULL, ignore_read_groups = TRUE,
                           extra_args = character(), config = NULL,
                           overwrite = FALSE) {
  assert_file(vcf, "VCF")
  assert_file(bam, "BAM")
  assert_file(reference, "reference FASTA")
  phased <- file.path(out_dir, paste0(sample, ".phased.vcf.gz"))
  tagged <- file.path(out_dir, paste0(sample, ".haplotagged.bam"))
  stats <- file.path(out_dir, paste0(sample, ".phasing_stats.tsv"))
  outs <- list(phased_vcf = phased, haplotagged_bam = tagged,
               stats_tsv = stats)
  params <- list(ignore_read_groups = ignore_read_groups)
  skip <- skip_if_done("phase", "whatshap", outs, params, overwrite)
  if (!is.null(skip)) return(skip)
  dir_create(out_dir)
  ensure_fai(reference, config)
  bin <- nf_bin("whatshap", config)
  log <- file.path(out_dir, "whatshap.log")

  plain <- sub("\\.gz$", "", phased)
  res <- nf_run(bin, c(
    "phase", "-o", plain, "--reference", reference,
    if (ignore_read_groups) "--ignore-read-groups",
    extra_args, vcf, bam), log = log)
  bgzip_tabix(plain, config)

  nf_run(bin, c("stats", "--tsv", stats, phased), log = log)
  nf_run(bin, c("haplotag", "-o", tagged, "--reference", reference,
                if (ignore_read_groups) "--ignore-read-groups",
                "--output-threads", 2, phased, bam), log = log)
  nf_run(nf_bin("samtools", config), c("index", tagged))

  if (!is.null(truth_vcf) && !is.na(truth_vcf)) {
    assert_file(truth_vcf, "truth VCF")
    compare <- file.path(out_dir, paste0(sample, ".switch_errors.tsv"))
    nf_run(bin, c("compare", "--tsv-pairwise", compare,
                  "--ignore-sample-name", truth_vcf, phased), log = log)
    outs$compare_tsv <- compare
  }
  new_step("phase", "whatshap", res$command, outs, params, res$runtime, log)
}

#' Phase variants with HapCUT2 (alternative)
#'
#' Runs `extractHAIRS` then `HAPCUT2`. Output is a haplotype block file plus
#' a phased VCF (HapCUT2's `--outvcf`).
#'
#' @inheritParams phase_whatshap
#' @return A `nanoflow_step` with `outputs$blocks` and `outputs$phased_vcf`.
#' @export
phase_hapcut2 <- function(vcf, bam, reference, out_dir, sample = "sample",
                          extra_args = character(), config = NULL,
                          overwrite = FALSE) {
  assert_file(vcf, "VCF")
  assert_file(bam, "BAM")
  blocks <- file.path(out_dir, paste0(sample, ".hapcut2.blocks"))
  phased <- paste0(blocks, ".phased.VCF")
  outs <- list(blocks = blocks, phased_vcf = phased)
  skip <- skip_if_done("phase", "hapcut2", outs, list(), overwrite)
  if (!is.null(skip)) return(skip)
  dir_create(out_dir)
  log <- file.path(out_dir, "hapcut2.log")
  extract <- Sys.which("extractHAIRS")
  if (!nzchar(extract)) {
    stop("extractHAIRS (part of HapCUT2) not found on PATH", call. = FALSE)
  }
  frags <- file.path(out_dir, paste0(sample, ".fragments"))
  # VCF must be uncompressed for extractHAIRS.
  plain_vcf <- vcf
  if (grepl("\\.gz$", vcf)) {
    plain_vcf <- file.path(out_dir, sub("\\.gz$", "", basename(vcf)))
    nf_run_shell(paste("gunzip -c", shQuote(vcf), ">", shQuote(plain_vcf)))
  }
  nf_run(extract, c("--bam", bam, "--VCF", plain_vcf, "--ont", "1",
                    "--ref", reference, "--out", frags), log = log)
  res <- nf_run(nf_bin("hapcut2", config),
                c("--fragments", frags, "--VCF", plain_vcf,
                  "--output", blocks, "--outvcf", "1", extra_args),
                log = log)
  new_step("phase", "hapcut2", res$command, outs, list(), res$runtime, log)
}
