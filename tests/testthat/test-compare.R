# Phase 2 bake-off: compare_sv_callers() / compare_small_variant_callers().

test_that("unknown callers are rejected before anything runs", {
  expect_error(
    compare_sv_callers(fixture("ont_reads.fastq.gz"), fixture("genome.fa"),
                       fixture("truth_sv.vcf"), tempfile(),
                       tools = "not-a-caller"),
    "unknown SV caller")
  expect_error(
    compare_small_variant_callers(fixture("ont_reads.fastq.gz"),
                                  fixture("genome.fa"),
                                  fixture("truth_small.vcf"), tempfile(),
                                  tools = "nope"),
    "unknown small-variant caller")
})

test_that("a missing caller is omitted, not fatal, and the table still builds", {
  skip_if_tool_missing("sniffles", "truvari", "bgzip", "tabix",
                       "minimap2", "samtools")
  # pair a real caller (sniffles) with one that is almost certainly absent
  # (nanovar) so the bake-off must drop the missing one gracefully.
  skip_if(nzchar(Sys.which("nanovar")), "nanovar is installed; absence untestable")
  out <- file.path(nf_test_dir(), "bakeoff_partial")
  cmp <- suppressWarnings(compare_sv_callers(
    test_bam(), test_reference(), fixture("truth_sv.vcf"), out,
    tools = c("sniffles", "nanovar"), sample = "sim1",
    all_contigs = TRUE, merge = FALSE, threads = 2))
  expect_s3_class(cmp, "nanoflow_benchmark")
  expect_true("sniffles" %in% cmp$comparison$tool)
  expect_false("nanovar" %in% cmp$comparison$tool)
})

test_that("SV bake-off ranks callers by F1 and recovers the truth SVs", {
  skip_if_tool_missing("sniffles", "svim", "truvari", "bgzip", "tabix",
                       "minimap2", "samtools")
  out <- file.path(nf_test_dir(), "bakeoff_sv")
  cmp <- suppressWarnings(compare_sv_callers(
    test_bam(), test_reference(), fixture("truth_sv.vcf"), out,
    tools = c("sniffles", "svim"), sample = "sim1",
    all_contigs = TRUE, merge = TRUE, threads = 2))
  expect_s3_class(cmp, "nanoflow_benchmark")
  df <- as.data.frame(cmp)
  expect_true(all(c("tool", "precision", "recall", "f1") %in% names(df)))
  expect_true(all(c("sniffles", "svim") %in% df$tool))
  # table is sorted best-F1 first
  expect_true(all(diff(df$f1) <= 1e-9))
  # at least one caller does well on the clean synthetic data
  expect_gte(max(df$recall, na.rm = TRUE), 0.8)

  csv <- file.path(out, "sv_comparison.csv")
  write_benchmark_csv(cmp, csv)
  expect_true(file.exists(csv))
  expect_equal(nrow(utils::read.csv(csv)), nrow(df))

  printed <- paste(capture.output(print(cmp)), collapse = "\n")
  expect_match(printed, "structural variants")
  expect_match(printed, "best F1")
})

test_that("the SURVIVOR merge appears as its own row when available", {
  skip_if_tool_missing("sniffles", "svim", "survivor", "truvari",
                       "bgzip", "tabix", "minimap2", "samtools")
  out <- file.path(nf_test_dir(), "bakeoff_sv")  # reuse (resumes)
  cmp <- suppressWarnings(compare_sv_callers(
    test_bam(), test_reference(), fixture("truth_sv.vcf"), out,
    tools = c("sniffles", "svim"), sample = "sim1",
    all_contigs = TRUE, merge = TRUE, min_callers = 1, threads = 2))
  expect_true("survivor_merge" %in% cmp$comparison$tool)
})
