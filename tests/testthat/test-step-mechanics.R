# Command construction, exit-status checking, resume and the GPU guard --
# none of these need bioinformatics tools installed.

test_that("nf_run reports failing commands with their exit status", {
  expect_error(nanoflow:::nf_run("false", character()), "exit 1")
  expect_error(nanoflow:::nf_run_shell("true | false"), "exit 1")
  res <- nanoflow:::nf_run("true", character())
  expect_equal(res$status, 0L)
  expect_true(res$runtime >= 0)
})

test_that("failing commands surface the log tail", {
  log <- file.path(tempdir(), "failing.log")
  expect_error(
    nanoflow:::nf_run_shell("echo some-diagnostic-output; exit 3", log = log),
    "some-diagnostic-output")
})

test_that("unknown binaries give an actionable error", {
  skip_if(nzchar(Sys.which("run_clair3.sh")), "clair3 is installed here")
  expect_error(nanoflow:::nf_bin("clair3"), "not found on PATH")
})

test_that("steps resume: existing outputs short-circuit to 'skipped'", {
  out <- file.path(tempdir(), "resume-test.txt")
  writeLines("x", out)
  skip <- nanoflow:::skip_if_done("align", "minimap2", list(bam = out),
                                  list(), overwrite = FALSE)
  expect_s3_class(skip, "nanoflow_step")
  expect_equal(skip$status, "skipped")
  # overwrite = TRUE disables the skip; missing file disables it too
  expect_null(nanoflow:::skip_if_done("align", "minimap2", list(bam = out),
                                      list(), overwrite = TRUE))
  expect_null(nanoflow:::skip_if_done("align", "minimap2",
                                      list(bam = "/no/such/file"),
                                      list(), overwrite = FALSE))
})

test_that("basecalling is skipped with a clear message without a GPU", {
  skip_if(has_gpu(), "machine has a GPU; skip-path not testable")
  step <- basecall_dorado("/data/run.pod5", tempdir(), "s1")
  expect_equal(step$status, "skipped")
  expect_match(step$message, "no CUDA GPU")
  # FASTQ input never needs basecalling, GPU or not
  step <- basecall_dorado(fixture("ont_reads.fastq.gz"), tempdir(), "s1")
  expect_equal(step$status, "skipped")
  expect_match(step$message, "already FASTQ")
})

test_that("step objects print a useful summary", {
  s <- nanoflow:::new_step("align", "minimap2", "minimap2 -ax map-ont ...",
                           list(bam = "x.bam"), list(), 1.23, status = "ok")
  out <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(out, "align")
  expect_match(out, "x\\.bam")
})

test_that("check_tools returns a complete table and flags missing tools", {
  suppressWarnings(
    res <- capture.output(tbl <- check_tools(tools = c("samtools", "dorado"))))
  expect_s3_class(tbl, "data.frame")
  expect_equal(tbl$tool, c("samtools", "dorado"))
  expect_error(check_tools(tools = "not-a-tool"), "unknown tool")
})

test_that("submit_slurm dry run writes valid per-sample scripts", {
  cfg_path <- file.path(tempdir(), "slurm_cfg.yml")
  out_dir <- file.path(tempdir(), "slurm_out")
  writeLines(yaml::as.yaml(list(output_dir = out_dir, threads = 2)), cfg_path)
  jobs <- submit_slurm(fixture("sample_sheet.csv"), cfg_path,
                       partition = "batch", gpus = 1,
                       setup_lines = "module load R", dry_run = TRUE)
  expect_equal(nrow(jobs), 1)
  script <- readLines(jobs$script)
  expect_true(any(grepl("--cpus-per-task=2", script)))
  expect_true(any(grepl("--partition=batch", script)))
  expect_true(any(grepl("--gres=gpu:1", script)))
  expect_true(any(grepl("module load R", script)))
  expect_true(any(grepl("run_pipeline.*samples = 'sim1'", script)))
})

test_that("IGV session files are written with relative track paths", {
  dir <- file.path(tempdir(), "igv")
  dir.create(dir, showWarnings = FALSE)
  bam <- file.path(dir, "a.bam")
  writeLines("x", bam)
  step <- write_igv_session(fixture("genome.fa"), bam,
                            file.path(dir, "session.xml"),
                            locus = "chr1:1-1000")
  xml <- readLines(step$outputs$session)
  expect_true(any(grepl('<Resource path="a.bam"/>', xml)))
  expect_true(any(grepl('locus="chr1:1-1000"', xml)))
})
