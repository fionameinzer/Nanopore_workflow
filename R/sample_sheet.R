#' Read and validate a sample sheet
#'
#' The sample sheet is a CSV with one row per sample:
#'
#' | column            | required | content                                        |
#' |-------------------|----------|------------------------------------------------|
#' | `sample`          | yes      | sample identifier (used for directory names)   |
#' | `reads`           | yes      | Nanopore FASTQ(.gz), or POD5/FAST5 file/dir    |
#' | `illumina_r1/_r2` | no       | paired-end short reads (hybrid assembly only)  |
#' | `truth_small_vcf` | no       | truth SNV/indel VCF (enables hap.py benchmark) |
#' | `truth_sv_vcf`    | no       | truth SV VCF (enables Truvari benchmark)       |
#' | `sequencing_summary` | no    | ONT sequencing summary (enables pycoQC)        |
#'
#' Relative paths are resolved against the directory containing the sheet,
#' so a sheet can ship next to its data.
#'
#' @param path CSV file path, or a data.frame already in the above shape.
#' @param check_files Error if a referenced file is missing (default TRUE).
#' @return A `nanoflow_samples` data.frame with an added `input_type`
#'   column (`"fastq"` or `"signal"`).
#' @examples
#' sheet <- read_sample_sheet(system.file("extdata/sim/sample_sheet.csv",
#'                                        package = "nanoflow"))
#' sheet$sample
#' @export
read_sample_sheet <- function(path, check_files = TRUE) {
  if (is.data.frame(path)) {
    sheet <- path
    base <- getwd()
  } else {
    assert_file(path, "sample sheet")
    sheet <- utils::read.csv(path, stringsAsFactors = FALSE)
    base <- dirname(normalizePath(path))
  }
  required <- c("sample", "reads")
  missing <- setdiff(required, names(sheet))
  if (length(missing)) {
    stop("sample sheet lacks required column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(sheet$sample)) {
    stop("duplicate sample IDs in sample sheet", call. = FALSE)
  }
  optional <- c("illumina_r1", "illumina_r2", "truth_small_vcf",
                "truth_sv_vcf", "sequencing_summary")
  for (col in optional) if (!col %in% names(sheet)) sheet[[col]] <- NA
  path_cols <- c("reads", optional)
  for (col in path_cols) {
    sheet[[col]] <- vapply(sheet[[col]], function(p) {
      if (is.na(p) || !nzchar(p)) return(NA_character_)
      if (!grepl("^(/|~)", p)) p <- file.path(base, p)
      if (check_files && !file.exists(p)) {
        stop(sprintf("sample sheet: %s does not exist (%s)", col, p),
             call. = FALSE)
      }
      normalizePath(p, mustWork = FALSE)
    }, character(1))
  }
  sheet$input_type <- ifelse(is_signal_input(sheet$reads), "signal", "fastq")
  structure(sheet, class = c("nanoflow_samples", "data.frame"))
}
