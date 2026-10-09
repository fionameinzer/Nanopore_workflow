# Functional annotation: SnpEff for small variants (prebuilt GRCh38 database
# by default, custom database built from a GFF3 otherwise) and AnnotSV for
# structural variants.

#' Build a custom SnpEff database from a FASTA + GFF3
#'
#' Needed when working with non-human organisms or custom references (and
#' used by the test suite on the synthetic fixture). Creates the SnpEff
#' `data/<db_name>` layout inside `data_dir`, writes a minimal
#' `snpEff.config`, and runs `snpEff build`.
#'
#' @param fasta Reference genome FASTA.
#' @param gff3 Annotation GFF3.
#' @param db_name Name for the database (e.g. `"mygenome1"`).
#' @param data_dir Directory to hold SnpEff databases and config.
#' @param extra_args Extra arguments for `snpEff build`.
#' @param config Optional nanoflow config (binary overrides).
#' @param overwrite Rebuild even if the database exists.
#' @return A `nanoflow_step` with `outputs$snpeff_config` and
#'   `outputs$db_dir`. Pass `snpeff_config` to [annotate_snpeff()].
#' @export
build_snpeff_db <- function(fasta, gff3, db_name, data_dir,
                            extra_args = character(), config = NULL,
                            overwrite = FALSE) {
  assert_file(fasta, "reference FASTA")
  assert_file(gff3, "GFF3")
  db_dir <- file.path(data_dir, "data", db_name)
  cfg_file <- file.path(data_dir, "snpEff.config")
  outs <- list(snpeff_config = cfg_file, db_dir = db_dir,
               bin = file.path(db_dir, "snpEffectPredictor.bin"))
  skip <- skip_if_done("snpeff_build", "snpeff", outs, list(db = db_name),
                       overwrite)
  if (!is.null(skip)) return(skip)
  dir_create(db_dir)
  file.copy(fasta, file.path(db_dir, "sequences.fa"), overwrite = TRUE)
  file.copy(gff3, file.path(db_dir, "genes.gff"), overwrite = TRUE)
  writeLines(c(
    paste0("data.dir = ", normalizePath(file.path(data_dir, "data"))),
    paste0(db_name, ".genome : ", db_name)
  ), cfg_file)
  log <- file.path(data_dir, "snpeff_build.log")
  res <- nf_run(nf_bin("snpeff", config), c(
    "build", "-gff3", "-c", cfg_file,
    "-noCheckCds", "-noCheckProtein", "-nodownload",
    extra_args, db_name), log = log)
  new_step("snpeff_build", "snpeff", res$command, outs,
           list(db = db_name), res$runtime, log)
}

#' Annotate small variants with SnpEff
#'
#' Defaults to the prebuilt human GRCh38 database; for other organisms
#' either name another prebuilt database or build one from a GFF3 with
#' [build_snpeff_db()] and pass its `snpeff_config`.
#'
#' @param vcf VCF(.gz) to annotate.
#' @param out_dir Output directory.
#' @param db SnpEff database name (default from config: `"GRCh38.105"`).
#' @param snpeff_config Optional path to a `snpEff.config` (required for
#'   custom databases; also where `data.dir` is taken from).
#' @param sample Sample ID used for file naming.
#' @param memory_gb Java heap.
#' @param extra_args Extra arguments for `snpEff ann`.
#' @param config Optional nanoflow config (binary overrides).
#' @param overwrite Rerun even if outputs already exist.
#' @return A `nanoflow_step` with `outputs$vcf` (annotated),
#'   `outputs$stats_html` and `outputs$genes_txt`.
#' @export
annotate_snpeff <- function(vcf, out_dir, db = "GRCh38.105",
                            snpeff_config = NULL, sample = "sample",
                            memory_gb = 8, extra_args = character(),
                            config = NULL, overwrite = FALSE) {
  assert_file(vcf, "VCF")
  ann_vcf <- file.path(out_dir, paste0(sample, ".snpeff.vcf"))
  stats <- file.path(out_dir, paste0(sample, ".snpeff.html"))
  outs <- list(vcf = ann_vcf, stats_html = stats)
  params <- list(db = db)
  skip <- skip_if_done("annotate", "snpeff", outs, params, overwrite)
  if (!is.null(skip)) return(skip)
  dir_create(out_dir)
  log <- file.path(out_dir, "snpeff.log")
  res <- nf_run_shell(paste(
    shQuote(nf_bin("snpeff", config)), "ann",
    sprintf("-Xmx%dg", memory_gb),
    if (!is.null(snpeff_config)) paste("-c", shQuote(snpeff_config)),
    "-stats", shQuote(stats),
    paste(extra_args, collapse = " "),
    shQuote(db), shQuote(vcf), ">", shQuote(ann_vcf)), log = log)
  new_step("annotate", "snpeff", res$command, outs, params, res$runtime, log)
}

#' Annotate structural variants with AnnotSV
#'
#' @param vcf SV VCF(.gz) to annotate.
#' @param out_dir Output directory.
#' @param genome_build e.g. `"GRCh38"` (AnnotSV naming).
#' @param annotations_dir AnnotSV annotation directory
#'   (`config$reference$annotsv_dir`); `NULL` uses AnnotSV's own default.
#' @param sample Sample ID used for file naming.
#' @param extra_args Extra command-line arguments appended verbatim.
#' @param config Optional nanoflow config (binary overrides).
#' @param overwrite Rerun even if outputs already exist.
#' @return A `nanoflow_step` with `outputs$tsv` (annotated SV table).
#' @export
annotate_annotsv <- function(vcf, out_dir, genome_build = "GRCh38",
                             annotations_dir = NULL, sample = "sample",
                             extra_args = character(), config = NULL,
                             overwrite = FALSE) {
  assert_file(vcf, "SV VCF")
  tsv <- file.path(out_dir, paste0(sample, ".annotsv.tsv"))
  params <- list(genome_build = genome_build)
  skip <- skip_if_done("annotate_sv", "annotsv", list(tsv = tsv), params,
                       overwrite)
  if (!is.null(skip)) return(skip)
  dir_create(out_dir)
  log <- file.path(out_dir, "annotsv.log")
  res <- nf_run(nf_bin("annotsv", config), c(
    "-SVinputFile", vcf,
    "-outputFile", tsv,
    "-genomeBuild", genome_build,
    if (!is.null(annotations_dir)) c("-annotationsDir", annotations_dir),
    extra_args), log = log)
  new_step("annotate_sv", "annotsv", res$command, list(tsv = tsv), params,
           res$runtime, log)
}
