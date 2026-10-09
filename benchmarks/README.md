# Benchmark results

Small, text-only result snapshots produced on the committed synthetic
fixture (`inst/extdata/sim`, a 200 kb diploid genome with phased truth
variants and ~30× simulated reads). These are the numbers; the full
browsable HTML report is published to GitHub Pages by CI (see the link in
the top-level README) and uploaded as a downloadable artifact on every run.

Regenerate locally with:

```r
library(nanoflow)
ref <- "genome.fa"
file.copy(system.file("extdata/sim/genome.fa", package = "nanoflow"), ref)
aln <- align_minimap2(
  system.file("extdata/sim/ont_reads.fastq.gz", package = "nanoflow"),
  ref, "align", sample = "sim1")
cmp <- compare_sv_callers(
  aln$outputs$bam, ref,
  system.file("extdata/sim/truth_sv.vcf", package = "nanoflow"),
  "bakeoff", tools = c("sniffles", "svim"), all_contigs = TRUE)
write_benchmark_csv(cmp, "benchmarks/sv_caller_bakeoff.csv")
```

## SV caller bake-off (`sv_caller_bakeoff.csv`)

Each caller runs on the same minimap2 alignment and is scored against the
truth SV set with Truvari. The 10 planted SVs are 3 DEL, 3 INS, 2 DUP,
2 INV.

| tool | precision | recall | F1 | TP | FP | FN | caller_sec |
|---|---|---|---|---|---|---|---|
| **sniffles** | 0.909 | **1.000** | **0.952** | 10 | 1 | 0 | 0.7 |
| **svim** | **1.000** | 0.600 | 0.750 | 6 | 0 | 4 | 10 |
| survivor_merge | NA | NA | NA | NA | NA | NA | — |

The classic precision/recall tradeoff: Sniffles recovers all 10 SVs (one
false positive); SVIM makes no false calls but misses four. `survivor_merge`
is the SURVIVOR consensus of both callers — it is reported as `NA` here
because Truvari 5.5 cannot score SVIM's breakend-style inversion ALT in the
merged record; the bake-off degrades to `NA` rather than failing, so the
rest of the table is unaffected.
