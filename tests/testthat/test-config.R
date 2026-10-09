test_that("defaults are returned and user values merge recursively", {
  cfg <- read_config(NULL)
  expect_s3_class(cfg, "nanoflow_config")
  expect_equal(cfg$threads, 4)
  expect_equal(cfg$steps$align$preset, "map-ont")

  cfg <- read_config(list(threads = 8,
                          steps = list(sv = list(extra_args = "--all-contigs"))))
  expect_equal(cfg$threads, 8)
  expect_equal(cfg$steps$sv$extra_args, "--all-contigs")
  # untouched defaults survive a partial override
  expect_equal(cfg$steps$sv$tool, "sniffles")
  expect_equal(cfg$steps$align$preset, "map-ont")
})

test_that("YAML round-trip works", {
  path <- file.path(tempdir(), "cfg.yml")
  writeLines(yaml::as.yaml(list(threads = 2, output_dir = "out",
                                steps = list(qc = list(min_length = 123)))),
             path)
  cfg <- read_config(path)
  expect_equal(cfg$threads, 2)
  expect_equal(cfg$steps$qc$min_length, 123)
  expect_equal(cfg$steps$qc$keep_percent, 90)
})

test_that("validation rejects bad values and missing files", {
  expect_error(read_config(list(threads = 0)))
  expect_error(read_config(list(reference = list(fasta = "/no/such/file.fa"))),
               "not found")
  expect_error(read_config(42), "must be")
})

test_that("config template writes and reads back", {
  path <- file.path(tempdir(), "template.yml")
  write_config_template(path, overwrite = TRUE)
  cfg <- read_config(path)
  expect_equal(cfg$steps$small_variants$tool, "clair3")
})

test_that("unknown binary overrides are caught", {
  expect_error(
    read_config(list(tools = list(minimap2 = "/no/such/minimap2"))),
    "not found")
})
