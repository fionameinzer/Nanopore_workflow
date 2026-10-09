# The committed synthetic fixture must stay internally consistent: every
# truth REF allele must match the reference sequence at its position, and
# the haplotype FASTAs must equal the reference with the truth variants
# applied.

read_fasta_seq <- function(path) {
  lines <- readLines(path)
  paste(lines[!startsWith(lines, ">")], collapse = "")
}

read_vcf_body <- function(path) {
  lines <- readLines(path)
  body <- lines[!startsWith(lines, "#")]
  f <- strsplit(body, "\t", fixed = TRUE)
  data.frame(
    pos = as.integer(vapply(f, `[`, "", 2)),
    id = vapply(f, `[`, "", 3),
    ref = vapply(f, `[`, "", 4),
    alt = vapply(f, `[`, "", 5),
    gt = vapply(f, `[`, "", 10),
    stringsAsFactors = FALSE
  )
}

test_that("fixture files exist and have the documented shape", {
  ref <- read_fasta_seq(fixture("genome.fa"))
  expect_equal(nchar(ref), 200000L)
  small <- read_vcf_body(fixture("truth_small.vcf"))
  svs <- read_vcf_body(fixture("truth_sv.vcf"))
  expect_equal(nrow(small), 290L)
  expect_equal(nrow(svs), 10L)
  expect_true(all(grepl("^[01]\\|[01]$", c(small$gt, svs$gt))))
  # small = SNVs + indels < 50 bp; SVs >= 50 bp
  expect_true(all(abs(nchar(small$alt) - nchar(small$ref)) < 50))
})

test_that("every truth REF allele matches the reference sequence", {
  ref <- read_fasta_seq(fixture("genome.fa"))
  vars <- rbind(read_vcf_body(fixture("truth_small.vcf")),
                read_vcf_body(fixture("truth_sv.vcf")))
  observed <- substring(ref, vars$pos, vars$pos + nchar(vars$ref) - 1L)
  expect_identical(observed, vars$ref)
})

test_that("haplotype FASTAs equal reference + applied truth variants", {
  ref <- read_fasta_seq(fixture("genome.fa"))
  vars <- rbind(read_vcf_body(fixture("truth_small.vcf")),
                read_vcf_body(fixture("truth_sv.vcf")))
  vars <- vars[order(vars$pos), ]
  apply_hap <- function(hap_idx) {
    gt <- vapply(strsplit(vars$gt, "|", fixed = TRUE),
                 function(g) as.integer(g[hap_idx]), 1L)
    v <- vars[gt == 1L, ]
    pieces <- character(0)
    cursor <- 1L
    for (i in seq_len(nrow(v))) {
      pieces <- c(pieces, substr(ref, cursor, v$pos[i] - 1L), v$alt[i])
      cursor <- v$pos[i] + nchar(v$ref[i])
    }
    paste(c(pieces, substr(ref, cursor, nchar(ref))), collapse = "")
  }
  expect_identical(apply_hap(1), read_fasta_seq(fixture("hap1.fa")))
  expect_identical(apply_hap(2), read_fasta_seq(fixture("hap2.fa")))
})

test_that("the GFF3 is well-formed and in range", {
  lines <- readLines(fixture("genes.gff3"))
  expect_identical(lines[1], "##gff-version 3")
  body <- lines[!startsWith(lines, "#")]
  f <- strsplit(body, "\t", fixed = TRUE)
  expect_true(all(lengths(f) == 9L))
  start <- as.integer(vapply(f, `[`, "", 4))
  end <- as.integer(vapply(f, `[`, "", 5))
  expect_true(all(start >= 1 & end <= 200000 & start <= end))
  types <- vapply(f, `[`, "", 3)
  expect_setequal(unique(types), c("gene", "mRNA", "exon", "CDS"))
  # CDS lengths are multiples of 3 so SnpEff can build a database
  cds_len <- (end - start + 1)[types == "CDS"]
  expect_true(all(cds_len %% 3 == 0))
})
