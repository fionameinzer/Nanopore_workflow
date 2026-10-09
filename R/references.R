# Reference and database downloads into a user-chosen directory. Nothing in
# the package assumes where these live; the config points at them.

#' Download reference genome and annotation databases
#'
#' Fetches, into a user-chosen directory: the GRCh38 reference genome
#' (no-alt analysis set; URL overridable for other organisms), the SnpEff
#' annotation database, and AnnotSV annotations. Writes a ready-to-use
#' snippet of config YAML next to the downloads.
#'
#' @param dest_dir Target directory (created if needed).
#' @param what Any of `"genome"`, `"snpeff"`, `"annotsv"`.
#' @param genome_url URL of the (gzipped) reference FASTA. Default: GRCh38
#'   no-alt analysis set from NCBI. Point this at any other FASTA for other
#'   organisms.
#' @param snpeff_db SnpEff database name to download (default
#'   `"GRCh38.105"`).
#' @param annotsv_url URL of an AnnotSV annotations tarball. AnnotSV
#'   versions its annotation bundles; see
#'   <https://lbgi.fr/AnnotSV/downloads>. If `NULL`, instructions are
#'   printed instead of downloading.
#' @param config Optional nanoflow config (binary overrides).
#' @return Invisibly, a list of the resulting paths (`genome_fasta`,
#'   `snpeff_data_dir`, `annotsv_dir` as applicable).
#' @export
download_references <- function(dest_dir,
                                what = c("genome", "snpeff"),
                                genome_url = paste0(
                                  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/",
                                  "GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/",
                                  "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"),
                                snpeff_db = "GRCh38.105",
                                annotsv_url = NULL,
                                config = NULL) {
  dir_create(dest_dir)
  out <- list()
  old_timeout <- options(timeout = max(7200, getOption("timeout")))
  on.exit(options(old_timeout), add = TRUE)

  if ("genome" %in% what) {
    gz <- file.path(dest_dir, basename(genome_url))
    fasta <- sub("\\.gz$", "", gz)
    if (!file.exists(fasta)) {
      message("[nanoflow] downloading reference genome ...")
      utils::download.file(genome_url, gz, mode = "wb")
      message("[nanoflow] decompressing ...")
      gunzip_file(gz, fasta)
      unlink(gz)
    } else {
      message("[nanoflow] reference already present: ", fasta)
    }
    if (nzchar(Sys.which("samtools"))) ensure_fai(fasta, config)
    out$genome_fasta <- fasta
  }

  if ("snpeff" %in% what) {
    data_dir <- file.path(dest_dir, "snpeff_data")
    dir_create(data_dir)
    marker <- file.path(data_dir, snpeff_db, "snpEffectPredictor.bin")
    if (!file.exists(marker)) {
      message("[nanoflow] downloading SnpEff database ", snpeff_db, " ...")
      nf_run(nf_bin("snpeff", config),
             c("download", "-dataDir", normalizePath(data_dir), snpeff_db),
             log = file.path(dest_dir, "snpeff_download.log"))
    } else {
      message("[nanoflow] SnpEff database already present: ", snpeff_db)
    }
    out$snpeff_data_dir <- data_dir
  }

  if ("annotsv" %in% what) {
    annotsv_dir <- file.path(dest_dir, "annotsv")
    if (is.null(annotsv_url)) {
      message("[nanoflow] AnnotSV annotations must match your AnnotSV ",
              "version; download the bundle from ",
              "https://lbgi.fr/AnnotSV/downloads and untar it into ",
              annotsv_dir, ", or rerun with annotsv_url=.")
    } else {
      dir_create(annotsv_dir)
      tarball <- file.path(annotsv_dir, basename(annotsv_url))
      message("[nanoflow] downloading AnnotSV annotations ...")
      utils::download.file(annotsv_url, tarball, mode = "wb")
      utils::untar(tarball, exdir = annotsv_dir)
      unlink(tarball)
    }
    out$annotsv_dir <- annotsv_dir
  }

  cfg_snippet <- file.path(dest_dir, "nanoflow_reference_config.yml")
  writeLines(yaml::as.yaml(list(reference = list(
    fasta = out$genome_fasta %||% NULL,
    snpeff_db = if ("snpeff" %in% what) snpeff_db else NULL,
    snpeff_data_dir = out$snpeff_data_dir %||% NULL,
    annotsv_dir = out$annotsv_dir %||% NULL
  ))), cfg_snippet)
  message("[nanoflow] wrote config snippet: ", cfg_snippet)
  invisible(out)
}

# Stream-decompress a .gz file with base R (no external gunzip needed).
gunzip_file <- function(gz, dest) {
  inc <- gzfile(gz, "rb")
  on.exit(close(inc), add = TRUE)
  outc <- file(dest, "wb")
  on.exit(close(outc), add = TRUE)
  repeat {
    chunk <- readBin(inc, raw(), 64 * 1024^2)
    if (!length(chunk)) break
    writeBin(chunk, outc)
  }
  invisible(dest)
}
