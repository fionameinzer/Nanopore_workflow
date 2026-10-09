# Integration tests on the fixture. Each test skips cleanly when its tool
# is missing, so the suite runs CPU-only anywhere; the conda environment
# gives full coverage.

test_that("minimap2 alignment produces a sorted, indexed BAM with ~30x depth", {
  bam <- test_bam()
  expect_true(file.exists(bam))
  expect_true(file.exists(paste0(bam, ".bai")))
  idx <- system2(Sys.which("samtools"), c("idxstats", bam), stdout = TRUE)
  mapped <- as.numeric(strsplit(idx[1], "\t")[[1]][3])
  expect_gt(mapped, 100)
  depth_out <- system2(Sys.which("samtools"),
                       c("depth", "-a", bam), stdout = TRUE)
  depth <- mean(as.numeric(vapply(strsplit(depth_out, "\t"), `[`, "", 3)))
  expect_gt(depth, 20)
  expect_lt(depth, 40)
})

test_that("alignment resumes instead of re-running", {
  bam <- test_bam()
  step <- align_minimap2(fixture("ont_reads.fastq.gz"), test_reference(),
                         dirname(bam), sample = "sim1", threads = 2)
  expect_equal(step$status, "skipped")
})

test_that("bwa aligns the Illumina reads", {
  skip_if_tool_missing("bwa", "samtools")
  step <- align_bwa(fixture("ill_R1.fastq.gz"), fixture("ill_R2.fastq.gz"),
                    test_reference(), file.path(nf_test_dir(), "align"),
                    sample = "sim1", threads = 2)
  expect_equal(step$status, "ok")
  expect_true(file.exists(step$outputs$bam))
})

test_that("NanoPlot runs on the fixture reads", {
  skip_if_tool_missing("nanoplot")
  step <- qc_nanoplot(fixture("ont_reads.fastq.gz"),
                      file.path(nf_test_dir(), "nanoplot"), threads = 2)
  expect_equal(step$status, "ok")
  expect_true(file.exists(step$outputs$stats))
  stats <- readLines(step$outputs$stats)
  expect_true(any(grepl("Number of reads", stats)))
})

test_that("porechop + filtlong trim and filter the reads", {
  skip_if_tool_missing("porechop", "filtlong")
  trimmed <- file.path(nf_test_dir(), "trimmed.fastq.gz")
  step <- trim_porechop(fixture("ont_reads.fastq.gz"), trimmed, threads = 2)
  expect_equal(step$status, "ok")
  expect_true(file.size(trimmed) > 1e5)
  filtered <- file.path(nf_test_dir(), "filtered.fastq.gz")
  step <- filter_filtlong(trimmed, filtered, min_length = 1000,
                          keep_percent = 90)
  expect_equal(step$status, "ok")
  con <- gzfile(filtered, "r")
  on.exit(close(con))
  lens <- nchar(readLines(con, 4000)[c(FALSE, TRUE, FALSE, FALSE)])
  expect_true(all(lens >= 1000))
})
