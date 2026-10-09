# Truth comparison: hap.py for small variants, Truvari for SVs. Run by the
# pipeline whenever the sample sheet provides truth VCFs. Both report
# precision / recall / F1.

#' Benchmark small-variant calls against a truth set with hap.py
#'
#' @param truth_vcf Truth SNV/indel VCF(.gz).
#' @param query_vcf Called VCF(.gz).
#' @param reference Reference genome FASTA.
#' @param out_dir Output directory.
#' @param sample Sample ID (output prefix).
#' @param confident_bed Optional BED of confident regions (e.g. GIAB
#'   high-confidence regions).
#' @param threads CPU threads.
#' @param extra_args Extra command-line arguments appended verbatim.
#' @param config Optional nanoflow config (binary overrides).
#' @param overwrite Rerun even if outputs already exist.
#' @return A `nanoflow_step` with `outputs$summary_csv`; the parsed
#'   precision/recall/F1 per variant type is attached as `params$metrics`.
#' @export
benchmark_happy <- function(truth_vcf, query_vcf, reference, out_dir,
                            sample = "sample", confident_bed = NULL,
                            threads = 4, extra_args = character(),
                            config = NULL, overwrite = FALSE) {
  assert_file(truth_vcf, "truth VCF")
  assert_file(query_vcf, "query VCF")
  assert_file(reference, "reference FASTA")
  prefix <- file.path(out_dir, paste0(sample, ".happy"))
  summary_csv <- paste0(prefix, ".summary.csv")
  skip <- skip_if_done("benchmark_small", "happy",
                       list(summary_csv = summary_csv), list(), overwrite)
  if (!is.null(skip)) {
    skip$params$metrics <- parse_happy_summary(summary_csv)
    return(skip)
  }
  dir_create(out_dir)
  ensure_fai(reference, config)
  truth_gz <- prepare_vcf_gz(truth_vcf, out_dir, config)
  query_gz <- prepare_vcf_gz(query_vcf, out_dir, config)
  log <- file.path(out_dir, "happy.log")
  res <- nf_run(nf_bin("happy", config), c(
    truth_gz, query_gz,
    "-r", reference,
    "-o", prefix,
    "--threads", threads,
    if (!is.null(confident_bed)) c("-f", confident_bed),
    extra_args), log = log)
  metrics <- parse_happy_summary(summary_csv)
  new_step("benchmark_small", "happy", res$command,
           list(summary_csv = summary_csv),
           list(metrics = metrics), res$runtime, log)
}

#' Benchmark SV calls against a truth set with Truvari
#'
#' @inheritParams benchmark_happy
#' @param truth_vcf Truth SV VCF(.gz).
#' @param query_vcf Called SV VCF(.gz).
#' @param dup_to_ins Treat DUP calls as INS when matching (`--dup-to-ins`);
#'   useful because tandem duplications are frequently called as insertions.
#' @param passonly Only consider PASS calls.
#' @return A `nanoflow_step` with `outputs$summary_json` (and the parsed
#'   precision/recall/F1 as `params$metrics`).
#' @export
benchmark_truvari <- function(truth_vcf, query_vcf, reference, out_dir,
                              sample = "sample", dup_to_ins = TRUE,
                              passonly = TRUE, extra_args = character(),
                              config = NULL, overwrite = FALSE) {
  assert_file(truth_vcf, "truth VCF")
  assert_file(query_vcf, "query VCF")
  assert_file(reference, "reference FASTA")
  bench_dir <- file.path(out_dir, paste0(sample, "_truvari"))
  summary_json <- file.path(bench_dir, "summary.json")
  skip <- skip_if_done("benchmark_sv", "truvari",
                       list(summary_json = summary_json), list(), overwrite)
  if (!is.null(skip)) {
    skip$params$metrics <- parse_truvari_summary(summary_json)
    return(skip)
  }
  dir_create(out_dir)
  # truvari bench requires that its output directory does not yet exist
  unlink(bench_dir, recursive = TRUE)
  truth_gz <- prepare_vcf_gz(truth_vcf, out_dir, config)
  query_gz <- prepare_vcf_gz(query_vcf, out_dir, config)
  log <- file.path(out_dir, "truvari.log")
  res <- nf_run(nf_bin("truvari", config), c(
    "bench",
    "-b", truth_gz,
    "-c", query_gz,
    "-o", bench_dir,
    "--reference", reference,
    if (dup_to_ins) "--dup-to-ins",
    if (passonly) "--passonly",
    extra_args), log = log)
  metrics <- parse_truvari_summary(summary_json)
  new_step("benchmark_sv", "truvari", res$command,
           list(summary_json = summary_json, dir = bench_dir),
           list(metrics = metrics), res$runtime, log)
}

# Copy a VCF next to the benchmark outputs, bgzip + tabix it.
prepare_vcf_gz <- function(vcf, out_dir, config = NULL) {
  dir_create(out_dir)
  if (grepl("\\.gz$", vcf)) {
    gz <- file.path(out_dir, basename(vcf))
    if (!file.exists(gz)) file.copy(vcf, gz)
  } else {
    plain <- file.path(out_dir, basename(vcf))
    file.copy(vcf, plain, overwrite = TRUE)
    gz <- paste0(plain, ".gz")
    nf_run(nf_bin("bgzip", config), c("-f", plain))
  }
  nf_run(nf_bin("tabix", config), c("-f", "-p", "vcf", gz))
  gz
}

parse_happy_summary <- function(summary_csv) {
  if (!file.exists(summary_csv)) return(NULL)
  df <- utils::read.csv(summary_csv, stringsAsFactors = FALSE)
  cols <- intersect(c("Type", "Filter", "TRUTH.TOTAL", "TRUTH.TP", "TRUTH.FN",
                      "QUERY.FP", "METRIC.Recall", "METRIC.Precision",
                      "METRIC.F1_Score"), names(df))
  df[df$Filter == "PASS", cols]
}

parse_truvari_summary <- function(summary_json) {
  if (!file.exists(summary_json)) return(NULL)
  s <- jsonlite::read_json(summary_json)
  data.frame(
    precision = s$precision %||% NA,
    recall = s$recall %||% NA,
    f1 = s$f1 %||% NA,
    TP = s$`TP-base` %||% NA,
    FP = s$FP %||% NA,
    FN = s$FN %||% NA
  )
}
