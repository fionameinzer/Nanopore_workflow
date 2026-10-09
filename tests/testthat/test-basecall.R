# Basecalling cannot be tested on synthetic data (Badread produces no raw
# signal), so this test uses a real POD5 sample and runs only on a GPU
# machine with dorado installed. Point NANOFLOW_POD5_DIR at a directory
# containing POD5 files (e.g. a download from ONT Open Data,
# https://labs.epi2me.io/dataindex/); everywhere else it skips.

test_that("dorado basecalls a public POD5 sample (GPU machines only)", {
  skip_if_not(has_gpu(), "no CUDA GPU on this machine")
  skip_if_tool_missing("dorado")
  pod5_dir <- Sys.getenv("NANOFLOW_POD5_DIR")
  skip_if(!nzchar(pod5_dir),
          "set NANOFLOW_POD5_DIR to a directory of POD5 files")
  step <- basecall_dorado(pod5_dir, file.path(nf_test_dir(), "basecall"),
                          sample = "pod5test", model = "fast")
  expect_equal(step$status, "ok")
  expect_true(file.exists(step$outputs$fastq))
  con <- gzfile(step$outputs$fastq, "r")
  on.exit(close(con))
  first <- readLines(con, 4)
  expect_match(first[1], "^@")
  expect_gt(nchar(first[2]), 10)
})
