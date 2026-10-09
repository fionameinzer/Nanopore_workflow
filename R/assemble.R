# Optional de novo assembly branch: Flye (default) or Canu for
# Nanopore-only, Wengan for hybrid assembly when Illumina reads are
# present; polishing with Racon and Medaka; QUAST assessment.

#' De novo assembly with Flye
#'
#' @param fastq Nanopore FASTQ(.gz), ideally trimmed/filtered.
#' @param out_dir Output directory (Flye's working directory).
#' @param read_type Flye input mode: `"nano-hq"` (modern, high-accuracy
#'   reads) or `"nano-raw"`.
#' @param genome_size Optional genome size estimate (e.g. `"3g"`, `"200k"`);
#'   Flye does not require it.
#' @param threads CPU threads.
#' @param extra_args Extra command-line arguments appended verbatim.
#' @param config Optional nanoflow config (binary overrides).
#' @param overwrite Rerun even if outputs already exist.
#' @return A `nanoflow_step` with `outputs$assembly` (FASTA) and
#'   `outputs$info`.
#' @export
assemble_flye <- function(fastq, out_dir, read_type = "nano-hq",
                          genome_size = NULL, threads = 4,
                          extra_args = character(), config = NULL,
                          overwrite = FALSE) {
  assert_file(fastq, "FASTQ")
  fasta <- file.path(out_dir, "assembly.fasta")
  outs <- list(assembly = fasta, info = file.path(out_dir, "assembly_info.txt"))
  params <- list(read_type = read_type, genome_size = genome_size)
  skip <- skip_if_done("assembly", "flye", outs, params, overwrite)
  if (!is.null(skip)) return(skip)
  dir_create(out_dir)
  log <- file.path(out_dir, "flye_wrapper.log")
  res <- nf_run(nf_bin("flye", config), c(
    paste0("--", read_type), fastq,
    "--out-dir", out_dir, "--threads", threads,
    if (!is.null(genome_size)) c("--genome-size", genome_size),
    extra_args), log = log)
  new_step("assembly", "flye", res$command, outs, params, res$runtime, log)
}

#' De novo assembly with Canu (alternative)
#'
#' @inheritParams assemble_flye
#' @param sample Prefix for Canu output files.
#' @param genome_size Required by Canu (e.g. `"200k"`, `"3.1g"`).
#' @return A `nanoflow_step` with `outputs$assembly`.
#' @export
assemble_canu <- function(fastq, out_dir, sample = "asm", genome_size,
                          threads = 4, extra_args = character(),
                          config = NULL, overwrite = FALSE) {
  assert_file(fastq, "FASTQ")
  fasta <- file.path(out_dir, paste0(sample, ".contigs.fasta"))
  params <- list(genome_size = genome_size)
  skip <- skip_if_done("assembly", "canu", list(assembly = fasta), params,
                       overwrite)
  if (!is.null(skip)) return(skip)
  dir_create(out_dir)
  log <- file.path(out_dir, "canu_wrapper.log")
  res <- nf_run(nf_bin("canu", config), c(
    "-p", sample, "-d", out_dir,
    paste0("genomeSize=", genome_size),
    paste0("maxThreads=", threads),
    "useGrid=false",
    "-nanopore", fastq, extra_args), log = log)
  new_step("assembly", "canu", res$command, list(assembly = fasta), params,
           res$runtime, log)
}

#' Hybrid assembly with Wengan (Nanopore + Illumina)
#'
#' Only applicable when the sample has Illumina paired-end reads.
#'
#' @inheritParams assemble_flye
#' @param r1,r2 Illumina paired-end FASTQ(.gz).
#' @param sample Output prefix.
#' @param genome_size_mb Approximate genome size in Mb (Wengan `-g`).
#' @param mode Wengan assembler mode (`"M"` for WenganM/Minia, `"A"` for
#'   WenganA/Abyss).
#' @return A `nanoflow_step` with `outputs$assembly`.
#' @export
assemble_wengan <- function(fastq, r1, r2, out_dir, sample = "asm",
                            genome_size_mb, mode = "M", threads = 4,
                            extra_args = character(), config = NULL,
                            overwrite = FALSE) {
  assert_file(fastq, "FASTQ")
  assert_file(r1, "R1 FASTQ")
  assert_file(r2, "R2 FASTQ")
  prefix <- file.path(out_dir, sample)
  fasta <- paste0(prefix, ".SPolished.asm.wengan.fasta")
  params <- list(genome_size_mb = genome_size_mb, mode = mode)
  skip <- skip_if_done("assembly", "wengan", list(assembly = fasta), params,
                       overwrite)
  if (!is.null(skip)) return(skip)
  dir_create(out_dir)
  log <- file.path(out_dir, "wengan.log")
  res <- nf_run(nf_bin("wengan", config), c(
    "-x", "ontraw", "-a", mode,
    "-s", paste(r1, r2, sep = ","),
    "-l", fastq,
    "-p", prefix, "-t", threads, "-g", genome_size_mb,
    extra_args), log = log)
  new_step("assembly", "wengan", res$command, list(assembly = fasta), params,
           res$runtime, log)
}

#' Polish an assembly with Racon
#'
#' Runs `rounds` iterations of minimap2 overlap + Racon consensus.
#'
#' @param fastq Reads used for polishing.
#' @param assembly Assembly FASTA to polish.
#' @param out_dir Output directory.
#' @param rounds Number of Racon rounds (1-2 are typical).
#' @inheritParams assemble_flye
#' @return A `nanoflow_step` with `outputs$assembly` (polished FASTA).
#' @export
polish_racon <- function(fastq, assembly, out_dir, rounds = 1, threads = 4,
                         extra_args = character(), config = NULL,
                         overwrite = FALSE) {
  assert_file(fastq, "FASTQ")
  assert_file(assembly, "assembly FASTA")
  polished <- file.path(out_dir, sprintf("racon_round%d.fasta", rounds))
  params <- list(rounds = rounds)
  skip <- skip_if_done("polish", "racon", list(assembly = polished), params,
                       overwrite)
  if (!is.null(skip)) return(skip)
  dir_create(out_dir)
  log <- file.path(out_dir, "racon.log")
  current <- assembly
  cmds <- character()
  for (r in seq_len(rounds)) {
    paf <- file.path(out_dir, sprintf("overlap_round%d.paf", r))
    out_fa <- file.path(out_dir, sprintf("racon_round%d.fasta", r))
    nf_run(nf_bin("minimap2", config),
           c("-x", "map-ont", "-t", threads, current, fastq),
           log = log, stdout_file = paf)
    res <- nf_run(nf_bin("racon", config),
                  c("-t", threads, extra_args, fastq, paf, current),
                  log = log, stdout_file = out_fa)
    cmds <- c(cmds, res$command)
    current <- out_fa
    unlink(paf)
  }
  new_step("polish", "racon", cmds, list(assembly = polished), params,
           res$runtime, log)
}

#' Polish an assembly with Medaka
#'
#' @inheritParams polish_racon
#' @param model Medaka model; `NULL` lets Medaka pick its default.
#' @return A `nanoflow_step` with `outputs$assembly` (consensus FASTA).
#' @export
polish_medaka <- function(fastq, assembly, out_dir, model = NULL,
                          threads = 4, extra_args = character(),
                          config = NULL, overwrite = FALSE) {
  assert_file(fastq, "FASTQ")
  assert_file(assembly, "assembly FASTA")
  consensus <- file.path(out_dir, "consensus.fasta")
  params <- list(model = model)
  skip <- skip_if_done("polish", "medaka", list(assembly = consensus),
                       params, overwrite)
  if (!is.null(skip)) return(skip)
  dir_create(out_dir)
  log <- file.path(out_dir, "medaka.log")
  mc <- Sys.which("medaka_consensus")
  if (!nzchar(mc)) stop("medaka_consensus not found on PATH", call. = FALSE)
  res <- nf_run(mc, c(
    "-i", fastq, "-d", assembly, "-o", out_dir, "-t", threads,
    if (!is.null(model)) c("-m", model), extra_args), log = log)
  new_step("polish", "medaka", res$command, list(assembly = consensus),
           params, res$runtime, log)
}

#' Assembly assessment with QUAST
#'
#' @param assembly Assembly FASTA (or several).
#' @param reference Optional reference FASTA for reference-based metrics.
#' @param gff3 Optional annotation (gene finding stats).
#' @param out_dir Output directory.
#' @inheritParams assemble_flye
#' @return A `nanoflow_step` with `outputs$report_tsv` and `outputs$dir`.
#' @export
asm_qc_quast <- function(assembly, out_dir, reference = NULL, gff3 = NULL,
                         threads = 4, extra_args = character(),
                         config = NULL, overwrite = FALSE) {
  for (f in assembly) assert_file(f, "assembly FASTA")
  report <- file.path(out_dir, "report.tsv")
  skip <- skip_if_done("assembly_qc", "quast", list(report_tsv = report),
                       list(), overwrite)
  if (!is.null(skip)) return(skip)
  dir_create(out_dir)
  log <- file.path(out_dir, "quast.log")
  res <- nf_run(nf_bin("quast", config), c(
    assembly, "-o", out_dir, "-t", threads,
    if (!is.null(reference)) c("-r", reference),
    if (!is.null(gff3)) c("-g", gff3),
    extra_args), log = log)
  new_step("assembly_qc", "quast", res$command,
           list(report_tsv = report, dir = out_dir), list(), res$runtime, log)
}
