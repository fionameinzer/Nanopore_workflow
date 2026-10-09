# ---------------------------------------------------------------------------
# Synthetic diploid genome + read simulator for nanoflow development fixtures.
#
# Not part of the installed package. Sourced by data-raw/make_fixture.R.
# External tools used only here (not user dependencies of nanoflow):
#   - badread      (Nanopore read simulation)
#   - art_illumina (paired-end Illumina read simulation)
#
# Everything is driven by the R RNG, so a single set.seed() in the calling
# script makes the genome, variants and truth files fully deterministic.
# Badread/ART are seeded explicitly through their own --seed/-rs options.
# ---------------------------------------------------------------------------

BASES <- c("A", "C", "G", "T")

random_seq <- function(n, gc = 0.42) {
  prob <- c((1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2)
  paste(sample(BASES, n, replace = TRUE, prob = prob), collapse = "")
}

revcomp <- function(s) {
  chartr("ACGT", "TGCA", paste(rev(strsplit(s, "")[[1]]), collapse = ""))
}

# --- reference -------------------------------------------------------------

#' Simulate a reference contig: random sequence of given GC content with a few
#' planted repeat structures (dispersed copies of one donor segment plus one
#' tandem repeat array), so that alignment/assembly see non-trivial sequence.
sim_reference <- function(length = 2e5, gc = 0.42,
                          n_dispersed = 2, dispersed_len = 2000,
                          tandem_motif_len = 50, tandem_copies = 30) {
  stopifnot(length > 20 * dispersed_len)
  seq <- random_seq(length, gc)

  # Choose non-overlapping windows for: 1 donor + n_dispersed copies + tandem.
  tandem_len <- tandem_motif_len * tandem_copies
  win_lens <- c(dispersed_len, rep(dispersed_len, n_dispersed), tandem_len)
  starts <- place_windows(length, win_lens, buffer = 2000, margin = 5000)

  donor <- substr(seq, starts[1], starts[1] + dispersed_len - 1)
  for (i in seq_len(n_dispersed)) {
    s <- starts[1 + i]
    substr(seq, s, s + dispersed_len - 1) <- donor
  }
  t_start <- starts[length(starts)]
  substr(seq, t_start, t_start + tandem_len - 1) <-
    strrep(random_seq(tandem_motif_len, gc), tandem_copies)

  attr(seq, "repeats") <- data.frame(
    type  = c("dispersed_donor", rep("dispersed_copy", n_dispersed), "tandem"),
    start = starts,
    end   = starts + win_lens - 1
  )
  seq
}

# Rejection-sample non-overlapping window start positions.
place_windows <- function(genome_len, win_lens, buffer = 1000, margin = 1000) {
  starts <- integer(0)
  occupied <- matrix(numeric(0), ncol = 2)
  for (len in win_lens) {
    for (try in 1:10000) {
      s <- sample(seq.int(margin, genome_len - margin - len), 1)
      e <- s + len - 1
      if (!any(s - buffer <= occupied[, 2] & e + buffer >= occupied[, 1])) {
        occupied <- rbind(occupied, c(s, e))
        starts <- c(starts, s)
        break
      }
      if (try == 10000) stop("could not place repeat windows; genome too small")
    }
  }
  starts
}

# --- variants --------------------------------------------------------------

#' Plant phased variants on a reference contig.
#'
#' Returns a data.frame with one row per variant, all sequence-resolved
#' (explicit REF/ALT strings), reference coordinates, 1-based:
#'   chrom pos id ref alt class(type small/sv) svtype svlen end gt1 gt2
#'
#' SVs are placed first (large exclusion buffer), then small indels, then
#' SNVs, so nothing overlaps anything else.
sim_variants <- function(ref_seq, chrom = "chr1",
                         n_snv = 250, n_indel = 40, max_indel = 30,
                         n_del = 3, n_ins = 3, n_dup = 2, n_inv = 2,
                         sv_len_range = c(300, 1500),
                         margin = 1500) {
  L <- nchar(ref_seq)
  base_at <- function(p) substr(ref_seq, p, p)
  seg <- function(s, e) substr(ref_seq, s, e)

  occupied <- matrix(numeric(0), ncol = 2)
  claim <- function(s, e, buffer) {
    if (any(s - buffer <= occupied[, 2] & e + buffer >= occupied[, 1])) {
      return(FALSE)
    }
    occupied <<- rbind(occupied, c(s, e))
    TRUE
  }
  draw_pos <- function(span, buffer) {
    for (try in 1:20000) {
      s <- sample(seq.int(margin, L - margin - span), 1)
      if (claim(s, s + span - 1, buffer)) return(s)
    }
    stop("could not place variant; reduce counts or enlarge genome")
  }
  rand_gt <- function() {
    gt <- sample(list(c(1L, 0L), c(0L, 1L), c(1L, 1L)), 1,
                 prob = c(0.4, 0.4, 0.2))[[1]]
    gt
  }

  rows <- list()
  add <- function(pos, ref, alt, class, svtype = NA_character_) {
    gt <- rand_gt()
    svlen <- nchar(alt) - nchar(ref)
    if (!is.na(svtype) && svtype == "INV") svlen <- nchar(ref) - 1L
    rows[[length(rows) + 1L]] <<- data.frame(
      chrom = chrom, pos = pos,
      id = sprintf("%s_%d", ifelse(is.na(svtype), class, tolower(svtype)),
                   length(rows) + 1L),
      ref = ref, alt = alt, class = class, svtype = svtype,
      svlen = svlen, end = pos + nchar(ref) - 1L,
      gt1 = gt[1], gt2 = gt[2], stringsAsFactors = FALSE
    )
  }

  sv_len <- function() sample(seq(sv_len_range[1], sv_len_range[2]), 1)

  # Structural variants (anchor-base padded, fully sequence-resolved).
  for (i in seq_len(n_del)) {              # deletion of [s..e]
    len <- sv_len()
    s <- draw_pos(len + 1, buffer = 500)
    add(s, seg(s, s + len), base_at(s), "sv", "DEL")
  }
  for (i in seq_len(n_ins)) {              # novel insertion after anchor s
    len <- sv_len()
    s <- draw_pos(1, buffer = 500)
    add(s, base_at(s), paste0(base_at(s), random_seq(len)), "sv", "INS")
  }
  for (i in seq_len(n_dup)) {              # tandem duplication of [s+1..e]
    len <- sv_len()
    s <- draw_pos(len + 1, buffer = 500)
    dup <- seg(s + 1, s + len)
    add(s + len, base_at(s + len), paste0(base_at(s + len), dup), "sv", "DUP")
  }
  for (i in seq_len(n_inv)) {              # inversion of [s+1..e]
    len <- sv_len()
    s <- draw_pos(len + 1, buffer = 500)
    add(s, paste0(base_at(s), seg(s + 1, s + len)),
        paste0(base_at(s), revcomp(seg(s + 1, s + len))), "sv", "INV")
  }

  # Small indels, 1..max_indel bp.
  for (i in seq_len(n_indel)) {
    len <- sample(max_indel, 1)
    if (runif(1) < 0.5) {                  # deletion
      s <- draw_pos(len + 1, buffer = 10)
      add(s, seg(s, s + len), base_at(s), "small")
    } else {                               # insertion
      s <- draw_pos(1, buffer = 10)
      add(s, base_at(s), paste0(base_at(s), random_seq(len)), "small")
    }
  }

  # SNVs.
  for (i in seq_len(n_snv)) {
    s <- draw_pos(1, buffer = 2)
    ref <- base_at(s)
    add(s, ref, sample(setdiff(BASES, ref), 1), "small")
  }

  out <- do.call(rbind, rows)
  out <- out[order(out$pos), ]
  rownames(out) <- NULL

  # Safety: REF strings must match the reference, spans must not overlap.
  stopifnot(all(mapply(function(p, r) seg(p, p + nchar(r) - 1) == r,
                       out$pos, out$ref)))
  stopifnot(all(out$pos[-1] > out$end[-nrow(out)]))
  out
}

# --- haplotype construction ------------------------------------------------

#' Apply the variants carried by one haplotype (gt column == 1) to the
#' reference sequence.
build_haplotype <- function(ref_seq, vars, hap = 1) {
  gt <- if (hap == 1) vars$gt1 else vars$gt2
  v <- vars[gt == 1L, ]
  v <- v[order(v$pos), ]
  pieces <- character(0)
  cursor <- 1L
  for (i in seq_len(nrow(v))) {
    pieces <- c(pieces, substr(ref_seq, cursor, v$pos[i] - 1L), v$alt[i])
    cursor <- v$pos[i] + nchar(v$ref[i])
  }
  pieces <- c(pieces, substr(ref_seq, cursor, nchar(ref_seq)))
  hap_seq <- paste(pieces, collapse = "")
  expected <- nchar(ref_seq) + sum(nchar(v$alt) - nchar(v$ref))
  stopifnot(nchar(hap_seq) == expected)
  hap_seq
}

# --- writers ---------------------------------------------------------------

write_fasta <- function(seqs, path, width = 70) {
  con <- file(path, "w")
  on.exit(close(con))
  for (nm in names(seqs)) {
    writeLines(paste0(">", nm), con)
    s <- seqs[[nm]]
    starts <- seq(1, nchar(s), by = width)
    writeLines(substring(s, starts, pmin(starts + width - 1, nchar(s))), con)
  }
  invisible(path)
}

write_truth_vcf <- function(vars, contig, contig_len, path,
                            sample = "TRUTH", sv = FALSE) {
  v <- vars[vars$class == (if (sv) "sv" else "small"), ]
  header <- c(
    "##fileformat=VCFv4.2",
    "##source=nanoflow-data-raw-simulator",
    sprintf("##contig=<ID=%s,length=%d>", contig, contig_len),
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End of REF span">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotype">',
    paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
          "FORMAT", sample, sep = "\t")
  )
  info <- if (sv) {
    svlen <- ifelse(v$svtype == "DEL", -abs(v$svlen), abs(v$svlen))
    sprintf("SVTYPE=%s;SVLEN=%d;END=%d", v$svtype, svlen, v$end)
  } else {
    rep(".", nrow(v))
  }
  body <- sprintf("%s\t%d\t%s\t%s\t%s\t.\tPASS\t%s\tGT\t%d|%d",
                  v$chrom, v$pos, v$id, v$ref, v$alt, info, v$gt1, v$gt2)
  writeLines(c(header, body), path)
  invisible(path)
}

#' Minimal but SnpEff-buildable GFF3: n_genes genes, each gene -> mRNA ->
#' exons + CDS. Every CDS block length is a multiple of 3 and phase is 0.
write_gff3 <- function(contig, contig_len, path, n_genes = 5,
                       exon_len = 300, intron_len = 200, n_exons = 3) {
  gene_len <- n_exons * exon_len + (n_exons - 1) * intron_len
  gap <- (contig_len - n_genes * gene_len) %/% (n_genes + 1)
  stopifnot(gap > 0)
  lines <- c("##gff-version 3",
             sprintf("##sequence-region %s 1 %d", contig, contig_len))
  for (g in seq_len(n_genes)) {
    gstart <- gap * g + gene_len * (g - 1) + 1
    gend <- gstart + gene_len - 1
    strand <- if (g %% 2 == 0) "-" else "+"
    gid <- sprintf("gene%02d", g)
    lines <- c(lines,
      sprintf("%s\tsim\tgene\t%d\t%d\t.\t%s\t.\tID=%s;Name=%s",
              contig, gstart, gend, strand, gid, toupper(gid)),
      sprintf("%s\tsim\tmRNA\t%d\t%d\t.\t%s\t.\tID=%s.t1;Parent=%s",
              contig, gstart, gend, strand, gid, gid))
    for (e in seq_len(n_exons)) {
      estart <- gstart + (e - 1) * (exon_len + intron_len)
      eend <- estart + exon_len - 1
      lines <- c(lines,
        sprintf("%s\tsim\texon\t%d\t%d\t.\t%s\t.\tID=%s.t1.exon%d;Parent=%s.t1",
                contig, estart, eend, strand, gid, e, gid),
        sprintf("%s\tsim\tCDS\t%d\t%d\t.\t%s\t0\tID=%s.t1.cds;Parent=%s.t1",
                contig, estart, eend, strand, gid, gid))
    }
  }
  writeLines(lines, path)
  invisible(path)
}

# --- read simulation (external tools) --------------------------------------

run_tool <- function(cmd, args, stdout = "", stderr = "") {
  message("  $ ", cmd, " ", paste(args, collapse = " "))
  status <- system2(cmd, args, stdout = stdout, stderr = stderr)
  if (!identical(status, 0L)) {
    stop(sprintf("%s failed with exit status %s", cmd, status))
  }
  invisible(TRUE)
}

#' Simulate Nanopore reads from each haplotype FASTA with Badread and write a
#' single combined, deterministically gzipped FASTQ. `depth` is per haplotype,
#' so total depth over the diploid sample is 2 * depth.
simulate_nanopore <- function(hap_fastas, out_fastq_gz, depth = 15,
                              seed = 1, mean_len = 10000, sd_len = 8000) {
  tmp <- sub("\\.gz$", "", out_fastq_gz)
  unlink(c(tmp, out_fastq_gz))
  for (i in seq_along(hap_fastas)) {
    part <- paste0(tmp, ".hap", i)
    run_tool("badread", c(
      "simulate",
      "--reference", hap_fastas[i],
      "--quantity", paste0(depth, "x"),
      "--seed", seed + i - 1,
      "--length", paste0(mean_len, ",", sd_len),
      "--identity", "95,99,2.5"
    ), stdout = part, stderr = paste0(part, ".log"))
    file.append(tmp, part)
    unlink(c(part, paste0(part, ".log")))
  }
  run_tool("gzip", c("-n", "-f", tmp))
  invisible(out_fastq_gz)
}

#' Simulate paired-end Illumina reads from each haplotype FASTA with ART and
#' write combined R1/R2 gzipped FASTQs. `depth` is per haplotype. Haplotype
#' FASTA headers must differ (chr1_hap1 / chr1_hap2) so read names are unique.
simulate_illumina <- function(hap_fastas, out_r1_gz, out_r2_gz, depth = 15,
                              seed = 1, read_len = 150,
                              frag_mean = 400, frag_sd = 50) {
  r1 <- sub("\\.gz$", "", out_r1_gz)
  r2 <- sub("\\.gz$", "", out_r2_gz)
  unlink(c(r1, r2, out_r1_gz, out_r2_gz))
  for (i in seq_along(hap_fastas)) {
    prefix <- tempfile(sprintf("art_hap%d_", i))
    run_tool("art_illumina", c(
      "-ss", "HS25", "-na", "-p",
      "-i", hap_fastas[i],
      "-l", read_len, "-f", depth,
      "-m", frag_mean, "-s", frag_sd,
      "-rs", seed + i - 1,
      "-o", prefix
    ), stdout = FALSE)
    file.append(r1, paste0(prefix, "1.fq"))
    file.append(r2, paste0(prefix, "2.fq"))
    unlink(paste0(prefix, c("1.fq", "2.fq")))
  }
  run_tool("gzip", c("-n", "-f", r1))
  run_tool("gzip", c("-n", "-f", r2))
  invisible(c(out_r1_gz, out_r2_gz))
}
