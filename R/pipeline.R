# run_pipeline(): per-sample orchestration of the whole workflow with
# resume (steps whose outputs exist are skipped) and a provenance log.

# Run one step, catching errors into a "failed" step object so one broken
# optional step does not kill the sample.
step_try <- function(step_name, tool, expr) {
  tryCatch(expr, error = function(e) {
    warning(sprintf("[%s] failed: %s", step_name, conditionMessage(e)),
            call. = FALSE)
    new_step(step_name, tool, status = "failed",
             message = conditionMessage(e))
  })
}

#' Run the full Nanopore workflow over a sample sheet
#'
#' Processes every sample independently through: optional basecalling
#' (GPU), read QC and trimming, alignment, optional de novo assembly,
#' small-variant calling, SV calling, phasing, optional annotation, truth
#' benchmarking (when truth VCFs are in the sample sheet), and reporting.
#' Steps whose outputs already exist are skipped, so an interrupted run can
#' simply be restarted (resume). A provenance log (`provenance.json` +
#' `provenance.rds`) with every command line, tool version and runtime is
#' written into the output directory.
#'
#' @param sample_sheet CSV path or data.frame; see [read_sample_sheet()].
#' @param config YAML path or list; see [read_config()]. Must provide
#'   `reference$fasta`.
#' @param samples Optionally restrict to these sample IDs (used by
#'   [submit_slurm()] to run one sample per job).
#' @param overwrite Rerun steps even if their outputs exist.
#' @return A `nanoflow_run` object: per-sample lists of `nanoflow_step`
#'   objects plus the config and sample sheet.
#' @examples
#' \dontrun{
#' run <- run_pipeline("samples.csv", "nanoflow_config.yml")
#' render_report(run)
#' }
#' @export
run_pipeline <- function(sample_sheet, config = NULL, samples = NULL,
                         overwrite = FALSE) {
  cfg <- read_config(config)
  sheet <- read_sample_sheet(sample_sheet)
  if (!is.null(samples)) {
    missing <- setdiff(samples, sheet$sample)
    if (length(missing)) stop("samples not in sheet: ",
                              paste(missing, collapse = ", "), call. = FALSE)
    sheet <- sheet[sheet$sample %in% samples, ]
  }
  if (is.null(cfg$reference$fasta)) {
    stop("config$reference$fasta is required", call. = FALSE)
  }
  out_root <- cfg$output_dir
  dir_create(out_root)

  tools <- needed_tools(cfg)
  message("[nanoflow] checking ", length(tools), " required tools ...")
  check_tools(cfg, tools = tools, stop_on_missing = FALSE)

  results <- list()
  for (i in seq_len(nrow(sheet))) {
    row <- sheet[i, ]
    message("\n[nanoflow] ===== sample ", row$sample, " (", i, "/",
            nrow(sheet), ") =====")
    results[[row$sample]] <- tryCatch(
      run_sample(row, cfg, out_root, overwrite),
      error = function(e) {
        warning(sprintf("sample %s failed: %s", row$sample,
                        conditionMessage(e)), call. = FALSE)
        list(error = new_step("sample", "nanoflow", status = "failed",
                              message = conditionMessage(e)))
      })
  }

  run <- structure(list(
    samples = results,
    sample_sheet = sheet,
    config = cfg,
    output_dir = normalizePath(out_root),
    finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  ), class = "nanoflow_run")

  # Aggregate MultiQC over everything, then the HTML report.
  if (isTRUE(cfg$steps$report$run)) {
    run$samples[["_run"]] <- list(
      multiqc = step_try("qc_multiqc", "multiqc",
        qc_multiqc(out_root, file.path(out_root, "multiqc"),
                   config = cfg, overwrite = overwrite)),
      report = step_try("report", "rmarkdown", render_report(run))
    )
  }
  write_provenance(run)
  run
}

run_sample <- function(row, cfg, out_root, overwrite = FALSE) {
  sample <- row$sample
  sdir <- function(...) file.path(out_root, sample, ...)
  ref <- cfg$reference$fasta
  threads <- cfg$threads
  steps <- list()
  fastq <- row$reads

  # 1. Basecalling (optional, GPU) -------------------------------------
  if (row$input_type == "signal") {
    bc <- basecall_dorado(row$reads, sdir("basecall"), sample,
                          model = cfg$steps$basecall$model,
                          extra_args = cfg$steps$basecall$extra_args,
                          config = cfg, overwrite = overwrite)
    steps$basecall <- bc
    if (is.null(bc$outputs$fastq)) {
      message("[nanoflow] ", bc$message, " -- sample cannot proceed")
      return(steps)
    }
    fastq <- bc$outputs$fastq
  }

  # 2. Read QC and trimming --------------------------------------------
  if (isTRUE(cfg$steps$qc$run)) {
    steps$qc_nanoplot <- step_try("qc_nanoplot", "nanoplot",
      qc_nanoplot(fastq, sdir("qc", "nanoplot_raw"), threads = 1,
                  config = cfg, overwrite = overwrite))
    if (!is.na(row$sequencing_summary)) {
      steps$qc_pycoqc <- step_try("qc_pycoqc", "pycoqc",
        qc_pycoqc(row$sequencing_summary, sdir("qc", "pycoqc"),
                  config = cfg, overwrite = overwrite))
    }
    if (!is.na(row$illumina_r1)) {
      steps$qc_fastqc <- step_try("qc_fastqc", "fastqc",
        qc_fastqc(c(row$illumina_r1, row$illumina_r2), sdir("qc", "fastqc"),
                  threads, config = cfg, overwrite = overwrite))
    }
    trimmed <- sdir("qc", paste0(sample, ".trimmed.fastq.gz"))
    steps$trim <- trim_porechop(fastq, trimmed, threads,
                                config = cfg, overwrite = overwrite)
    filtered <- sdir("qc", paste0(sample, ".filtered.fastq.gz"))
    steps$filter <- filter_filtlong(
      trimmed, filtered,
      min_length = cfg$steps$qc$min_length,
      keep_percent = cfg$steps$qc$keep_percent,
      target_bases = cfg$steps$qc$target_bases,
      config = cfg, overwrite = overwrite)
    fastq <- filtered
    steps$qc_nanoplot_filtered <- step_try("qc_nanoplot", "nanoplot",
      qc_nanoplot(fastq, sdir("qc", "nanoplot_filtered"), threads = 1,
                  config = cfg, overwrite = overwrite))
  }

  # 3. Alignment ---------------------------------------------------------
  steps$align <- align_minimap2(
    fastq, ref, sdir("align"), sample,
    preset = cfg$steps$align$preset, threads = threads,
    extra_args = cfg$steps$align$extra_args,
    config = cfg, overwrite = overwrite)
  bam <- steps$align$outputs$bam
  steps$bam_qc <- step_try("bam_qc", "qualimap",
    bam_qc_qualimap(bam, sdir("align", "qualimap"), threads,
                    memory_gb = cfg$memory_gb, config = cfg,
                    overwrite = overwrite))
  if (!is.na(row$illumina_r1)) {
    steps$align_illumina <- step_try("align_illumina", "bwa",
      align_bwa(row$illumina_r1, row$illumina_r2, ref, sdir("align"),
                sample, threads, config = cfg, overwrite = overwrite))
  }

  # 4. De novo assembly (optional branch) --------------------------------
  if (isTRUE(cfg$steps$assembly$run)) {
    steps <- c(steps, run_assembly_branch(row, cfg, fastq, sdir, overwrite))
  }

  # 5/6. Variant calling -------------------------------------------------
  small_vcf <- NULL
  if (isTRUE(cfg$steps$small_variants$run)) {
    sv_cfg <- cfg$steps$small_variants
    steps$small_variants <- if (identical(sv_cfg$tool, "medaka")) {
      call_medaka(bam, ref, sdir("small_variants"), sample,
                  threads = threads, extra_args = sv_cfg$extra_args,
                  config = cfg, overwrite = overwrite)
    } else {
      call_clair3(bam, ref, sdir("small_variants"), sample,
                  model_dir = sv_cfg$clair3_model_dir, model = sv_cfg$model,
                  threads = threads, extra_args = sv_cfg$extra_args,
                  config = cfg, overwrite = overwrite)
    }
    small_vcf <- steps$small_variants$outputs$vcf
  }

  sv_vcf <- NULL
  if (isTRUE(cfg$steps$sv$run)) {
    svs <- cfg$steps$sv
    steps$sv <- switch(svs$tool %||% "sniffles",
      svim = call_svim(bam, ref, sdir("sv"), sample,
                       extra_args = svs$extra_args, config = cfg,
                       overwrite = overwrite),
      nanovar = call_nanovar(bam, ref, sdir("sv"), sample, threads,
                             extra_args = svs$extra_args, config = cfg,
                             overwrite = overwrite),
      call_sniffles(bam, ref, sdir("sv"), sample, threads = threads,
                    extra_args = svs$extra_args, config = cfg,
                    overwrite = overwrite))
    sv_vcf <- steps$sv$outputs$vcf
  }

  # 7. Phasing (needs VCF + BAM) -----------------------------------------
  if (isTRUE(cfg$steps$phase$run) && !is.null(small_vcf)) {
    steps$phase <- step_try("phase", cfg$steps$phase$tool %||% "whatshap",
      if (identical(cfg$steps$phase$tool, "hapcut2")) {
        phase_hapcut2(small_vcf, bam, ref, sdir("phase"), sample,
                      extra_args = cfg$steps$phase$extra_args,
                      config = cfg, overwrite = overwrite)
      } else {
        phase_whatshap(small_vcf, bam, ref, sdir("phase"), sample,
                       truth_vcf = row$truth_small_vcf,
                       extra_args = cfg$steps$phase$extra_args,
                       config = cfg, overwrite = overwrite)
      })
  }

  # 8. Annotation (optional) ----------------------------------------------
  if (isTRUE(cfg$steps$annotate$run)) {
    if (!is.null(small_vcf)) {
      steps$annotate <- step_try("annotate", "snpeff",
        annotate_snpeff(small_vcf, sdir("annotate"),
                        db = cfg$reference$snpeff_db,
                        snpeff_config = cfg$reference$snpeff_data_dir,
                        sample = sample, memory_gb = cfg$memory_gb,
                        extra_args = cfg$steps$annotate$extra_args,
                        config = cfg, overwrite = overwrite))
    }
    if (!is.null(sv_vcf)) {
      steps$annotate_sv <- step_try("annotate_sv", "annotsv",
        annotate_annotsv(sv_vcf, sdir("annotate"),
                         annotations_dir = cfg$reference$annotsv_dir,
                         sample = sample, config = cfg,
                         overwrite = overwrite))
    }
  }

  # 9. Benchmarking against truth (auto when truth given) -----------------
  bench <- cfg$steps$benchmark$run
  bench_on <- isTRUE(bench) || identical(bench, "auto")
  if (bench_on && !is.na(row$truth_small_vcf) && !is.null(small_vcf)) {
    steps$benchmark_small <- step_try("benchmark_small", "happy",
      benchmark_happy(row$truth_small_vcf, small_vcf, ref,
                      sdir("benchmark"), sample, threads = threads,
                      config = cfg, overwrite = overwrite))
  }
  if (bench_on && !is.na(row$truth_sv_vcf) && !is.null(sv_vcf)) {
    steps$benchmark_sv <- step_try("benchmark_sv", "truvari",
      benchmark_truvari(row$truth_sv_vcf, sv_vcf, ref,
                        sdir("benchmark"), sample,
                        extra_args = cfg$steps$benchmark$extra_args,
                        config = cfg, overwrite = overwrite))
  }

  # 10. Per-sample visualization ------------------------------------------
  if (isTRUE(cfg$steps$report$run)) {
    tracks <- stats::na.omit(c(
      bam,
      steps$phase$outputs$haplotagged_bam %||% NULL,
      small_vcf, sv_vcf))
    steps$igv_session <- step_try("igv_session", "igv",
      write_igv_session(ref, tracks, sdir("igv_session.xml")))
    regions <- cfg$steps$report$igv_regions
    if (length(regions)) {
      steps$igv_snapshots <- step_try("igv_snapshots", "igv_reports",
        igv_snapshots(regions, ref, tracks, sdir("igv_regions.html"),
                      config = cfg))
    }
  }

  saveRDS(steps, sdir("steps.rds"))
  steps
}

run_assembly_branch <- function(row, cfg, fastq, sdir, overwrite) {
  acfg <- cfg$steps$assembly
  threads <- cfg$threads
  ref <- cfg$reference$fasta
  steps <- list()
  hybrid <- identical(acfg$tool, "wengan") && !is.na(row$illumina_r1)
  steps$assembly <- step_try("assembly", acfg$tool %||% "flye", {
    if (hybrid) {
      assemble_wengan(fastq, row$illumina_r1, row$illumina_r2,
                      sdir("assembly"), row$sample,
                      genome_size_mb = acfg$genome_size_mb %||% 3100,
                      threads = threads, extra_args = acfg$extra_args,
                      config = cfg, overwrite = overwrite)
    } else if (identical(acfg$tool, "canu")) {
      assemble_canu(fastq, sdir("assembly"), row$sample,
                    genome_size = acfg$genome_size %||%
                      stop("canu requires steps$assembly$genome_size"),
                    threads = threads, extra_args = acfg$extra_args,
                    config = cfg, overwrite = overwrite)
    } else {
      assemble_flye(fastq, sdir("assembly"),
                    genome_size = acfg$genome_size, threads = threads,
                    extra_args = acfg$extra_args, config = cfg,
                    overwrite = overwrite)
    }
  })
  asm <- steps$assembly$outputs$assembly
  if (is.null(asm)) return(steps)

  if ("racon" %in% acfg$polish) {
    steps$polish_racon <- step_try("polish", "racon",
      polish_racon(fastq, asm, sdir("assembly", "racon"),
                   threads = threads, config = cfg, overwrite = overwrite))
    asm <- steps$polish_racon$outputs$assembly %||% asm
  }
  if ("medaka" %in% acfg$polish) {
    steps$polish_medaka <- step_try("polish", "medaka",
      polish_medaka(fastq, asm, sdir("assembly", "medaka"),
                    threads = threads, config = cfg, overwrite = overwrite))
    asm <- steps$polish_medaka$outputs$assembly %||% asm
  }
  steps$assembly_qc <- step_try("assembly_qc", "quast",
    asm_qc_quast(asm, sdir("assembly", "quast"), reference = ref,
                 gff3 = cfg$reference$gff3, threads = threads,
                 config = cfg, overwrite = overwrite))
  steps
}

write_provenance <- function(run) {
  rds <- file.path(run$output_dir, "provenance.rds")
  saveRDS(run, rds)
  slim <- lapply(run$samples, function(steps) {
    lapply(steps, function(s) {
      if (!inherits(s, "nanoflow_step")) return(NULL)
      list(step = s$step, tool = s$tool, tool_version = s$tool_version,
           status = s$status, command = s$command,
           outputs = s$outputs, runtime_sec = s$runtime_sec,
           finished = s$finished, message = s$message)
    })
  })
  jsonlite::write_json(
    list(finished = run$finished, output_dir = run$output_dir,
         samples = slim),
    file.path(run$output_dir, "provenance.json"),
    auto_unbox = TRUE, pretty = TRUE, null = "null")
  invisible(run)
}

#' @export
print.nanoflow_run <- function(x, ...) {
  cat("<nanoflow run>", x$output_dir, "\n")
  cat("finished:", x$finished, "\n\n")
  for (sample in names(x$samples)) {
    cat(if (sample == "_run") "run-level steps" else paste("sample", sample),
        ":\n", sep = "")
    for (nm in names(x$samples[[sample]])) {
      s <- x$samples[[sample]][[nm]]
      if (!inherits(s, "nanoflow_step")) next
      mark <- switch(s$status, ok = "\u2713", skipped = "\u21b7", "\u2717")
      cat(sprintf("  %s %-20s %-10s %s\n", mark, nm, s$tool,
                  if (s$status == "ok" && !is.na(s$runtime_sec))
                    sprintf("%.1fs", s$runtime_sec) else s$status))
    }
  }
  invisible(x)
}
