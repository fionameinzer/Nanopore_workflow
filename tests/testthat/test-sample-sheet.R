test_that("the fixture sample sheet loads with resolved paths", {
  sheet <- read_sample_sheet(fixture("sample_sheet.csv"))
  expect_s3_class(sheet, "nanoflow_samples")
  expect_equal(sheet$sample, "sim1")
  expect_true(file.exists(sheet$reads))
  expect_true(file.exists(sheet$truth_small_vcf))
  expect_true(file.exists(sheet$truth_sv_vcf))
  expect_equal(sheet$input_type, "fastq")
})

test_that("missing required columns and duplicate IDs are rejected", {
  path <- file.path(tempdir(), "bad_sheet.csv")
  write.csv(data.frame(sample = "a"), path, row.names = FALSE)
  expect_error(read_sample_sheet(path), "required column")

  write.csv(data.frame(sample = c("a", "a"),
                       reads = fixture("ont_reads.fastq.gz")),
            path, row.names = FALSE)
  expect_error(read_sample_sheet(path), "duplicate")
})

test_that("missing files are reported with the offending column", {
  path <- file.path(tempdir(), "missing_sheet.csv")
  write.csv(data.frame(sample = "a", reads = "/no/such/reads.fastq.gz"),
            path, row.names = FALSE)
  expect_error(read_sample_sheet(path), "reads does not exist")
  expect_s3_class(read_sample_sheet(path, check_files = FALSE),
                  "nanoflow_samples")
})

test_that("POD5 input is classified as signal", {
  sheet <- read_sample_sheet(
    data.frame(sample = "s", reads = "run/raw.pod5"), check_files = FALSE)
  expect_equal(sheet$input_type, "signal")
})
