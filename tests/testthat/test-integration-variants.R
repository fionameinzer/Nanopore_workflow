# Variant calling, phasing and benchmarking on the fixture.

test_that("Sniffles2 recovers the truth SVs (needs --all-contigs!)", {
  skip_if_tool_missing("sniffles", "minimap2", "samtools")
  step <- call_sniffles(test_bam(), test_reference(),
                        file.path(nf_test_dir(), "sv"), sample = "sim1",
                        all_contigs = TRUE, threads = 2)
  expect_equal(step$status, "ok")
  vcf <- readLines(step$outputs$vcf)
  calls <- vcf[!startsWith(vcf, "#")]
  # 10 truth SVs; allow caller noise in both directions
  expect_gte(length(calls), 8)
  expect_lte(length(calls), 15)
})

test_that("Truvari confirms high precision/recall against truth", {
  skip_if_tool_missing("sniffles", "truvari", "bgzip", "tabix",
                       "minimap2", "samtools")
  sv_vcf <- file.path(nf_test_dir(), "sv", "sim1.sniffles.vcf")
  if (!file.exists(sv_vcf)) skip("sniffles test did not run")
  step <- benchmark_truvari(fixture("truth_sv.vcf"), sv_vcf,
                            test_reference(),
                            file.path(nf_test_dir(), "bench"),
                            sample = "sim1")
  expect_equal(step$status, "ok")
  m <- step$params$metrics
  expect_gte(m$recall, 0.8)
  expect_gte(m$precision, 0.7)
})

test_that("WhatsHap phases variants and evaluates switch errors vs truth", {
  skip_if_tool_missing("whatshap", "bgzip", "tabix", "minimap2", "samtools")
  # unphase the truth VCF -> a realistic unphased caller output
  unphased <- file.path(nf_test_dir(), "unphased.vcf")
  lines <- readLines(fixture("truth_small.vcf"))
  body <- !startsWith(lines, "#")
  lines[body] <- sub("(\\d)\\|(\\d)$", "\\1/\\2", lines[body])
  writeLines(lines, unphased)

  step <- phase_whatshap(unphased, test_bam(), test_reference(),
                         file.path(nf_test_dir(), "phase"), sample = "sim1",
                         truth_vcf = fixture("truth_small.vcf"))
  expect_equal(step$status, "ok")
  expect_true(file.exists(step$outputs$phased_vcf))
  expect_true(file.exists(step$outputs$haplotagged_bam))
  phased <- readLines(gzfile(step$outputs$phased_vcf))
  n_phased <- sum(grepl("\\d\\|\\d", phased[!startsWith(phased, "#")]))
  expect_gt(n_phased, 50)

  expect_true(file.exists(step$outputs$compare_tsv))
  cmp <- read.delim(step$outputs$compare_tsv)
  expect_true("all_switches" %in% names(cmp))
  # 30x synthetic reads on 200 kb: phasing should be near-perfect
  expect_lte(cmp$all_switches[1], 5)
})

test_that("Flye assembles the fixture genome", {
  skip_if_tool_missing("flye")
  skip_on_cran()
  step <- assemble_flye(fixture("ont_reads.fastq.gz"),
                        file.path(nf_test_dir(), "flye"),
                        genome_size = "200k", threads = 2)
  expect_equal(step$status, "ok")
  fa <- readLines(step$outputs$assembly)
  asm_len <- sum(nchar(fa[!startsWith(fa, ">")]))
  # one contig close to 200 kb
  expect_gt(asm_len, 150000)
  expect_lt(asm_len, 260000)
})

test_that("a custom SnpEff database builds from the fixture GFF3", {
  skip_if_tool_missing("snpeff")
  data_dir <- file.path(nf_test_dir(), "snpeff")
  built <- build_snpeff_db(fixture("genome.fa"), fixture("genes.gff3"),
                           "simgenome", data_dir)
  expect_equal(built$status, "ok")
  step <- annotate_snpeff(fixture("truth_small.vcf"),
                          file.path(nf_test_dir(), "snpeff_ann"),
                          db = "simgenome",
                          snpeff_config = built$outputs$snpeff_config,
                          sample = "sim1")
  expect_equal(step$status, "ok")
  ann <- readLines(step$outputs$vcf)
  expect_true(any(grepl("ANN=", ann)))
})

test_that("hap.py benchmarks small variants when available", {
  skip_if_tool_missing("happy", "bgzip", "tabix")
  skip("requires a small-variant callset; exercised in the real-data vignette")
})
