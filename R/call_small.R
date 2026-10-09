# Small variant calling (SNVs and indels) from Nanopore alignments.
# (The design diagram labels this box "SV"; that is a typo -- this step
# calls small variants. Structural variants are handled in call_sv.R.)

#' Call SNVs and indels with Clair3
#'
#' @param bam Sorted, indexed Nanopore BAM.
#' @param reference Reference genome FASTA (a `.fai` index is created if
#'   missing).
#' @param out_dir Output directory.
#' @param sample Sample ID.
#' @param model_dir Path to the Clair3 model directory. If `NULL`, the
#'   wrapper looks for `models/<model>` next to the `run_clair3.sh` binary
#'   (the layout of the bioconda install).
#' @param model Clair3 model name, e.g. `"r1041_e82_400bps_sup_v420"`.
#'   Configurable for other chemistries/organisms.
#' @param include_all_ctgs Pass `--include_all_ctgs` so contigs outside
#'   chr1-22/X/Y are called too (needed for non-human references and the
#'   synthetic fixture).
#' @param threads CPU threads.
#' @param extra_args Extra command-line arguments appended verbatim.
#' @param config Optional nanoflow config (binary overrides).
#' @param overwrite Rerun even if outputs already exist.
#' @return A `nanoflow_step` with `outputs$vcf` (bgzipped, indexed).
#' @export
call_clair3 <- function(bam, reference, out_dir, sample = "sample",
                        model_dir = NULL, model = "r1041_e82_400bps_sup_v420",
                        include_all_ctgs = FALSE, threads = 4,
                        extra_args = character(), config = NULL,
                        overwrite = FALSE) {
  assert_file(bam, "BAM")
  assert_file(reference, "reference FASTA")
  vcf <- file.path(out_dir, "merge_output.vcf.gz")
  params <- list(model = model, include_all_ctgs = include_all_ctgs)
  skip <- skip_if_done("small_variants", "clair3", list(vcf = vcf), params,
                       overwrite)
  if (!is.null(skip)) return(skip)
  dir_create(out_dir)
  ensure_fai(reference, config)
  bin <- nf_bin("clair3", config)
  if (is.null(model_dir)) {
    model_dir <- config$steps$small_variants$clair3_model_dir %||%
      file.path(dirname(bin), "models", model)
  }
  if (!dir.exists(model_dir)) {
    stop("Clair3 model directory not found: ", model_dir,
         "\nSet model_dir= or config$steps$small_variants$clair3_model_dir.",
         call. = FALSE)
  }
  log <- file.path(out_dir, "clair3.log")
  res <- nf_run(bin, c(
    paste0("--bam_fn=", bam),
    paste0("--ref_fn=", reference),
    paste0("--output=", out_dir),
    paste0("--threads=", threads),
    "--platform=ont",
    paste0("--model_path=", model_dir),
    paste0("--sample_name=", sample),
    if (include_all_ctgs) "--include_all_ctgs",
    extra_args), log = log)
  new_step("small_variants", "clair3", res$command, list(vcf = vcf), params,
           res$runtime, log)
}

#' Call SNVs and indels with Medaka (alternative to Clair3)
#'
#' Uses `medaka_variant` on an existing alignment.
#'
#' @inheritParams call_clair3
#' @param model Medaka model; `NULL` lets Medaka pick its default.
#' @return A `nanoflow_step` with `outputs$vcf`.
#' @export
call_medaka <- function(bam, reference, out_dir, sample = "sample",
                        model = NULL, threads = 4,
                        extra_args = character(), config = NULL,
                        overwrite = FALSE) {
  assert_file(bam, "BAM")
  assert_file(reference, "reference FASTA")
  vcf <- file.path(out_dir, "medaka.annotated.vcf")
  params <- list(model = model)
  skip <- skip_if_done("small_variants", "medaka", list(vcf = vcf), params,
                       overwrite)
  if (!is.null(skip)) return(skip)
  dir_create(out_dir)
  ensure_fai(reference, config)
  log <- file.path(out_dir, "medaka_variant.log")
  res <- nf_run(Sys.which("medaka_variant") %||% "medaka_variant", c(
    "-i", bam, "-f", reference, "-o", out_dir, "-t", threads,
    if (!is.null(model)) c("-m", model),
    extra_args), log = log)
  new_step("small_variants", "medaka", res$command, list(vcf = vcf), params,
           res$runtime, log)
}
