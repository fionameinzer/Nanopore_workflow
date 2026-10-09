# End-to-end run_pipeline() on the fixture, plus resume and provenance.
# QC and assembly are switched off here to keep the suite fast; those
# wrappers have their own integration tests.

pipeline_cfg <- function(out_dir) {
  list(
    reference = list(fasta = fixture("genome.fa"),
                     gff3 = fixture("genes.gff3")),
    output_dir = out_dir,
    threads = 2,
    steps = list(
      qc = list(run = FALSE),
      assembly = list(run = FALSE),
      small_variants = list(run = FALSE),  # Clair3 needs its model dir
      sv = list(run = TRUE, tool = "sniffles",
                extra_args = "--all-contigs"),
      report = list(run = FALSE)
    )
  )
}

test_that("run_pipeline processes the fixture end to end", {
  skip_if_tool_missing("minimap2", "samtools", "sniffles", "truvari",
                       "bgzip", "tabix")
  out_dir <- file.path(nf_test_dir(), "pipeline_run")
  run <- suppressWarnings(
    run_pipeline(fixture("sample_sheet.csv"), pipeline_cfg(out_dir)))
  expect_s3_class(run, "nanoflow_run")
  steps <- run$samples$sim1

  expect_equal(steps$align$status, "ok")
  expect_true(file.exists(steps$align$outputs$bam))
  expect_equal(steps$sv$status, "ok")
  # truth VCFs in the sheet trigger benchmarking automatically
  expect_equal(steps$benchmark_sv$status, "ok")
  expect_gte(steps$benchmark_sv$params$metrics$recall, 0.8)

  # provenance log
  prov <- file.path(run$output_dir, "provenance.json")
  expect_true(file.exists(prov))
  j <- jsonlite::read_json(prov)
  expect_true("sim1" %in% names(j$samples))
  expect_match(j$samples$sim1$align$command[[1]], "minimap2")
  expect_true(file.exists(file.path(run$output_dir, "provenance.rds")))
})

test_that("a second run resumes: every step is skipped", {
  skip_if_tool_missing("minimap2", "samtools", "sniffles", "truvari",
                       "bgzip", "tabix")
  out_dir <- file.path(nf_test_dir(), "pipeline_run")
  if (!dir.exists(out_dir)) skip("first pipeline test did not run")
  run2 <- suppressWarnings(
    run_pipeline(fixture("sample_sheet.csv"), pipeline_cfg(out_dir)))
  statuses <- vapply(run2$samples$sim1[c("align", "sv", "benchmark_sv")],
                     function(s) s$status, "")
  expect_true(all(statuses == "skipped"))
})

test_that("run_pipeline restricts to requested samples", {
  expect_error(
    run_pipeline(fixture("sample_sheet.csv"),
                 pipeline_cfg(tempfile()), samples = "nope"),
    "not in sheet")
})

test_that("the HTML report renders from a run object", {
  skip_if_tool_missing("minimap2", "samtools", "sniffles", "truvari",
                       "bgzip", "tabix")
  skip_if_not_installed("rmarkdown")
  skip_if(!rmarkdown::pandoc_available("2.0"), "pandoc not available")
  prov <- file.path(nf_test_dir(), "pipeline_run", "provenance.rds")
  if (!file.exists(prov)) skip("pipeline test did not run")
  step <- render_report(prov,
                        out_html = file.path(nf_test_dir(), "report.html"))
  expect_true(file.exists(step$outputs$html))
  html <- paste(readLines(step$outputs$html, warn = FALSE), collapse = "")
  expect_match(html, "Benchmarks against truth")
  expect_match(html, "sim1")
})
