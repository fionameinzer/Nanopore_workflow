# Basecalling: POD5/FAST5 signal data -> FASTQ. GPU-dependent; skipped with
# a clear message when no GPU is present or when the input is already FASTQ.

#' Basecall Nanopore signal data with Dorado
#'
#' Converts POD5 (or FAST5) raw signal to FASTQ. Requires a CUDA GPU; on a
#' machine without one the step is skipped with a clear message rather than
#' failing, so CPU-only pipelines (and the test suite) run end to end on
#' FASTQ input.
#'
#' @param input POD5/FAST5 file or directory. If a FASTQ is passed the step
#'   is skipped (nothing to do).
#' @param out_dir Output directory.
#' @param sample Sample ID used for file naming.
#' @param model Dorado model name or path (e.g. `"sup"`, `"hac"`, or a full
#'   model directory). Configurable for non-default chemistries.
#' @param device Dorado device string; default `"cuda:all"`. Use `"cpu"` to
#'   force CPU basecalling (very slow; intended for smoke tests only).
#' @param extra_args Extra command-line arguments appended verbatim.
#' @param config Optional nanoflow config (binary overrides).
#' @param overwrite Rerun even if outputs already exist.
#' @return A `nanoflow_step` with `outputs$fastq`.
#' @export
basecall_dorado <- function(input, out_dir, sample = "sample", model = "sup",
                            device = "cuda:all",
                            extra_args = character(), config = NULL,
                            overwrite = FALSE) {
  if (is_fastq(input)) {
    return(new_step("basecall", "dorado", status = "skipped",
                    outputs = list(fastq = input),
                    message = "input is already FASTQ; basecalling not needed"))
  }
  if (!identical(device, "cpu") && !has_gpu()) {
    return(new_step("basecall", "dorado", status = "skipped",
                    message = paste("no CUDA GPU detected (nvidia-smi not",
                                    "available); basecalling skipped. Run on a",
                                    "GPU node or pass device = \"cpu\".")))
  }
  dir_create(out_dir)
  fastq <- file.path(out_dir, paste0(sample, ".fastq.gz"))
  params <- list(model = model, device = device)
  skip <- skip_if_done("basecall", "dorado", list(fastq = fastq), params,
                       overwrite)
  if (!is.null(skip)) return(skip)
  log <- file.path(out_dir, "dorado.log")
  bin <- nf_bin("dorado", config)
  res <- nf_run_shell(paste(
    shQuote(bin), "basecaller", shQuote(model), shQuote(input),
    "--emit-fastq", "--device", shQuote(device),
    paste(extra_args, collapse = " "),
    "| gzip >", shQuote(fastq)), log = log)
  new_step("basecall", "dorado", res$command, list(fastq = fastq), params,
           res$runtime, log)
}

#' Basecall Nanopore signal data with Guppy (deprecated)
#'
#' Kept because Guppy appears in the workflow design diagram, but Guppy has
#' been discontinued by ONT in favour of Dorado. Prefer [basecall_dorado()].
#'
#' @inheritParams basecall_dorado
#' @param guppy_config Guppy basecalling config, e.g.
#'   `"dna_r10.4.1_e8.2_400bps_sup.cfg"`.
#' @return A `nanoflow_step` with `outputs$fastq_dir`.
#' @export
basecall_guppy <- function(input, out_dir, sample = "sample",
                           guppy_config = "dna_r10.4.1_e8.2_400bps_sup.cfg",
                           device = "cuda:all",
                           extra_args = character(), config = NULL,
                           overwrite = FALSE) {
  .Deprecated("basecall_dorado",
              msg = "Guppy is discontinued by ONT; prefer basecall_dorado().")
  if (is_fastq(input)) {
    return(new_step("basecall", "guppy", status = "skipped",
                    outputs = list(fastq = input),
                    message = "input is already FASTQ; basecalling not needed"))
  }
  if (!identical(device, "cpu") && !has_gpu()) {
    return(new_step("basecall", "guppy", status = "skipped",
                    message = "no CUDA GPU detected; basecalling skipped"))
  }
  dir_create(out_dir)
  log <- file.path(out_dir, "guppy.log")
  args <- c("-i", input, "-s", out_dir, "-c", guppy_config,
            if (!identical(device, "cpu")) c("--device", device),
            extra_args)
  res <- nf_run(nf_bin("guppy", config), args, log = log)
  new_step("basecall", "guppy", res$command, list(fastq_dir = out_dir),
           list(guppy_config = guppy_config, device = device),
           res$runtime, log)
}
