# Thin Slurm wrapper: one job per sample, each calling run_pipeline() for
# that sample only. The package works identically without Slurm; this just
# writes and submits batch scripts.

#' Submit per-sample pipeline jobs to Slurm
#'
#' Writes one `sbatch` script per sample (each running
#' `nanoflow::run_pipeline()` restricted to that sample) and submits them.
#' With `dry_run = TRUE` the scripts are written but not submitted, which
#' also serves as the CPU-only test path on machines without Slurm.
#'
#' @param sample_sheet Path to the sample sheet CSV (passed through to the
#'   jobs, so it must be readable from the compute nodes).
#' @param config_file Path to the YAML config (same requirement).
#' @param partition,account Slurm partition/account; omitted when `NULL`.
#' @param time Wall-clock limit per sample job.
#' @param cpus CPUs per job; default from the config's `threads`.
#' @param mem_gb Memory per job; default from the config's `memory_gb`.
#' @param gpus GPUs per job (only needed when basecalling POD5 input),
#'   e.g. `1` adds `--gres=gpu:1`.
#' @param setup_lines Shell lines inserted before the R call (e.g.
#'   `"module load R"` or `"source activate nanoflow"`).
#' @param dry_run Write scripts but do not submit.
#' @return Invisibly, a data.frame with `sample`, `script` and (when
#'   submitted) `job_id`.
#' @export
submit_slurm <- function(sample_sheet, config_file, partition = NULL,
                         account = NULL, time = "48:00:00", cpus = NULL,
                         mem_gb = NULL, gpus = 0,
                         setup_lines = character(), dry_run = FALSE) {
  cfg <- read_config(config_file)
  sheet <- read_sample_sheet(sample_sheet)
  cpus <- cpus %||% cfg$threads
  mem_gb <- mem_gb %||% cfg$memory_gb
  job_dir <- file.path(cfg$output_dir, "slurm")
  dir_create(job_dir)

  jobs <- lapply(sheet$sample, function(sample) {
    script <- file.path(job_dir, paste0("nanoflow_", sample, ".sbatch"))
    writeLines(c(
      "#!/bin/bash",
      sprintf("#SBATCH --job-name=nanoflow_%s", sample),
      sprintf("#SBATCH --cpus-per-task=%d", cpus),
      sprintf("#SBATCH --mem=%dG", mem_gb),
      sprintf("#SBATCH --time=%s", time),
      sprintf("#SBATCH --output=%s", file.path(job_dir,
              paste0(sample, "_%j.out"))),
      if (!is.null(partition)) sprintf("#SBATCH --partition=%s", partition),
      if (!is.null(account)) sprintf("#SBATCH --account=%s", account),
      if (gpus > 0) sprintf("#SBATCH --gres=gpu:%d", gpus),
      "",
      setup_lines,
      "",
      sprintf(
        "Rscript -e 'nanoflow::run_pipeline(%s, %s, samples = %s)'",
        shQuote(normalizePath(sample_sheet)),
        shQuote(normalizePath(config_file)),
        shQuote(sample))
    ), script)
    job_id <- NA_character_
    if (!dry_run) {
      sbatch <- Sys.which("sbatch")
      if (!nzchar(sbatch)) {
        stop("sbatch not found on PATH; use dry_run = TRUE to only write scripts",
             call. = FALSE)
      }
      out <- system2(sbatch, script, stdout = TRUE)
      job_id <- regmatches(out, regexpr("[0-9]+", out))[1]
      message("[nanoflow] submitted ", sample, " as job ", job_id)
    }
    data.frame(sample = sample, script = script, job_id = job_id,
               stringsAsFactors = FALSE)
  })
  invisible(do.call(rbind, jobs))
}
