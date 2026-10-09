# Shared helpers for the test suite. All integration tests run on the tiny
# synthetic fixture in inst/extdata/sim and skip cleanly when the external
# tool they exercise is not installed, so the suite passes CPU-only on any
# machine (full coverage requires the conda environment).

fixture <- function(...) {
  system.file("extdata", "sim", ..., package = "nanoflow", mustWork = TRUE)
}

skip_if_tool_missing <- function(...) {
  for (tool in c(...)) {
    entry <- nanoflow:::nf_tool_registry()[[tool]]
    if (!nzchar(Sys.which(entry$bin))) {
      testthat::skip(paste("external tool not installed:", entry$bin))
    }
  }
}

# One shared scratch dir per test run; aligning once and reusing the BAM
# keeps the suite fast.
nf_test_dir <- local({
  dir <- NULL
  function() {
    if (is.null(dir)) {
      dir <<- file.path(tempdir(), "nanoflow-tests")
      dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    }
    dir
  }
})

# Reference copied to scratch so index files (.fai/.bwt) never land in the
# installed package directory.
test_reference <- local({
  ref <- NULL
  function() {
    if (is.null(ref)) {
      ref <<- file.path(nf_test_dir(), "genome.fa")
      file.copy(fixture("genome.fa"), ref)
    }
    ref
  }
})

test_bam <- local({
  bam <- NULL
  function() {
    skip_if_tool_missing("minimap2", "samtools")
    if (is.null(bam)) {
      step <- align_minimap2(fixture("ont_reads.fastq.gz"), test_reference(),
                             file.path(nf_test_dir(), "align"),
                             sample = "sim1", threads = 2)
      bam <<- step$outputs$bam
    }
    bam
  }
})
