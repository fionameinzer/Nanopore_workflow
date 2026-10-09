# Alignment: minimap2 (Nanopore), BWA-MEM (Illumina), samtools sort/index,
# Qualimap BAM QC.

# Ensure a .fai index exists for a reference FASTA.
ensure_fai <- function(reference, config = NULL) {
  fai <- paste0(reference, ".fai")
  if (!file.exists(fai)) {
    nf_run(nf_bin("samtools", config), c("faidx", reference))
  }
  fai
}

#' Align Nanopore reads with minimap2
#'
#' Runs `minimap2 -ax <preset> | samtools sort` and indexes the result. A
#' read group with the sample ID is added so downstream callers and phasers
#' see a proper sample name.
#'
#' @param fastq Nanopore FASTQ(.gz).
#' @param reference Reference genome FASTA.
#' @param out_dir Output directory.
#' @param sample Sample ID (used for file naming and the BAM read group).
#' @param preset minimap2 preset; default `"map-ont"` (use `"lr:hq"` for
#'   high-accuracy simplex reads if preferred).
#' @param threads CPU threads.
#' @param extra_args Extra minimap2 arguments appended verbatim.
#' @param config Optional nanoflow config (binary overrides).
#' @param overwrite Rerun even if outputs already exist.
#' @return A `nanoflow_step` with `outputs$bam` and `outputs$bai`.
#' @export
align_minimap2 <- function(fastq, reference, out_dir, sample = "sample",
                           preset = "map-ont", threads = 4,
                           extra_args = character(), config = NULL,
                           overwrite = FALSE) {
  assert_file(fastq, "FASTQ")
  assert_file(reference, "reference FASTA")
  bam <- file.path(out_dir, paste0(sample, ".ont.bam"))
  outs <- list(bam = bam, bai = paste0(bam, ".bai"))
  params <- list(preset = preset, threads = threads)
  skip <- skip_if_done("align", "minimap2", outs, params, overwrite)
  if (!is.null(skip)) return(skip)
  dir_create(out_dir)
  ensure_fai(reference, config)
  log <- file.path(out_dir, paste0(sample, ".minimap2.log"))
  rg <- sprintf("@RG\\tID:%s\\tSM:%s\\tPL:ONT", sample, sample)
  samtools <- shQuote(nf_bin("samtools", config))
  res <- nf_run_shell(paste(
    shQuote(nf_bin("minimap2", config)),
    "-ax", preset, "-t", threads, "--MD",
    "-R", shQuote(rg),
    paste(extra_args, collapse = " "),
    shQuote(reference), shQuote(fastq),
    "|", samtools, "sort", "-@", threads, "-o", shQuote(bam), "-"),
    log = log)
  nf_run(nf_bin("samtools", config), c("index", bam))
  new_step("align", "minimap2", res$command, outs, params, res$runtime, log)
}

#' Align Illumina paired-end reads with BWA-MEM
#'
#' Used when a sample has optional Illumina reads. Builds the BWA index on
#' first use if it does not exist next to the reference.
#'
#' @param r1,r2 Paired-end FASTQ(.gz) files.
#' @inheritParams align_minimap2
#' @return A `nanoflow_step` with `outputs$bam` and `outputs$bai`.
#' @export
align_bwa <- function(r1, r2, reference, out_dir, sample = "sample",
                      threads = 4, extra_args = character(), config = NULL,
                      overwrite = FALSE) {
  assert_file(r1, "R1 FASTQ")
  assert_file(r2, "R2 FASTQ")
  assert_file(reference, "reference FASTA")
  bam <- file.path(out_dir, paste0(sample, ".ill.bam"))
  outs <- list(bam = bam, bai = paste0(bam, ".bai"))
  params <- list(threads = threads)
  skip <- skip_if_done("align_illumina", "bwa", outs, params, overwrite)
  if (!is.null(skip)) return(skip)
  dir_create(out_dir)
  ensure_fai(reference, config)
  bwa <- nf_bin("bwa", config)
  if (!file.exists(paste0(reference, ".bwt"))) {
    nf_run(bwa, c("index", reference),
           log = file.path(out_dir, "bwa_index.log"))
  }
  log <- file.path(out_dir, paste0(sample, ".bwa.log"))
  rg <- sprintf("@RG\\tID:%s\\tSM:%s\\tPL:ILLUMINA", sample, sample)
  res <- nf_run_shell(paste(
    shQuote(bwa), "mem", "-t", threads, "-R", shQuote(rg),
    paste(extra_args, collapse = " "),
    shQuote(reference), shQuote(r1), shQuote(r2),
    "|", shQuote(nf_bin("samtools", config)), "sort", "-@", threads,
    "-o", shQuote(bam), "-"), log = log)
  nf_run(nf_bin("samtools", config), c("index", bam))
  new_step("align_illumina", "bwa", res$command, outs, params, res$runtime,
           log)
}

#' BAM quality control with Qualimap
#'
#' @param bam Sorted, indexed BAM.
#' @param out_dir Output directory for the Qualimap report.
#' @param memory_gb Java heap for Qualimap.
#' @inheritParams align_minimap2
#' @return A `nanoflow_step` with `outputs$report`.
#' @export
bam_qc_qualimap <- function(bam, out_dir, threads = 4, memory_gb = 8,
                            extra_args = character(), config = NULL,
                            overwrite = FALSE) {
  assert_file(bam, "BAM")
  report <- file.path(out_dir, "qualimapReport.html")
  skip <- skip_if_done("bam_qc", "qualimap", list(report = report), list(),
                       overwrite)
  if (!is.null(skip)) return(skip)
  dir_create(out_dir)
  log <- file.path(out_dir, "qualimap.log")
  res <- nf_run(nf_bin("qualimap", config),
                c("bamqc", "-bam", bam, "-outdir", out_dir, "-nt", threads,
                  sprintf("--java-mem-size=%dG", memory_gb), extra_args),
                log = log)
  new_step("bam_qc", "qualimap", res$command,
           list(report = report, dir = out_dir), list(), res$runtime, log)
}
