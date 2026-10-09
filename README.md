# nanoflow

[![CI](https://github.com/fionameinzer/Nanopore_workflow/actions/workflows/ci.yml/badge.svg)](https://github.com/fionameinzer/Nanopore_workflow/actions/workflows/ci.yml)

**End-to-end Oxford Nanopore long-read whole-genome analysis, orchestrated from R.**

📊 **Results:** CI rebuilds the workflow on the synthetic fixture on every push —
see the [live HTML run report](https://fionameinzer.github.io/Nanopore_workflow/)
(GitHub Pages) and the caller-comparison numbers in
[`benchmarks/`](benchmarks/).

nanoflow runs a complete ONT WGS workflow on real sequencing data —
basecalling → QC/trimming → alignment → (optional) de novo assembly →
small-variant calling → SV calling → phasing → annotation → truth
benchmarking → HTML report — by wrapping established command-line tools
via `system2()`. It reimplements no algorithm; it manages inputs, outputs,
resume and provenance between steps and returns standard formats (FASTQ,
BAM, VCF, FASTA). Nanopore reads are the primary input; Illumina short
reads are an optional per-sample extra used for hybrid assembly.

| step | default tool | alternatives |
|---|---|---|
| Basecalling (GPU, optional) | Dorado | Guppy (deprecated) |
| Read QC + trimming | NanoPlot, Porechop, Filtlong, MultiQC | pycoQC, FastQC |
| Alignment | minimap2 (map-ont) + samtools | BWA-MEM (Illumina), Qualimap |
| De novo assembly (optional) | Flye + Racon/Medaka + QUAST | Canu, Wengan (hybrid) |
| Small variants (SNV/indel) | Clair3 | Medaka |
| Structural variants | Sniffles2 | SVIM, NanoVar, SURVIVOR merge, Spectre (CNV) |
| Phasing | WhatsHap | HapCUT2 |
| Annotation | SnpEff (small), AnnotSV (SV) | custom DB from GFF3 |
| Benchmarking (auto with truth) | hap.py (small), Truvari (SV) | |
| Reporting | R Markdown HTML, IGV session, igv-reports | |

Defaults target human GRCh38; the reference genome, tool models and
databases are all configurable for other diploid organisms. **No path is
hard-coded** — everything comes from a YAML config.

## Install on a new server

1. **Clone and inspect**

   ```bash
   git clone https://github.com/fionameinzer/Nanopore_workflow.git
   cd Nanopore_workflow
   ```

2. **External tools** — pick one:

   *With conda/mamba (recommended):*

   ```bash
   mamba env create -f inst/conda/environment.yml        # all CPU tools
   mamba env create -f inst/conda/environment-happy.yml  # hap.py (needs py2.7)
   conda activate nanoflow
   ```

   *Without conda (clusters):* build the Apptainer image instead —

   ```bash
   apptainer build nanoflow.sif inst/container/nanoflow.def
   ```

   *GPU basecalling only:* install ONT's [Dorado](https://github.com/nanoporetech/dorado)
   on the GPU node(s); it is not on bioconda. Without a GPU the basecalling
   step is skipped with a clear message, and everything downstream runs on
   CPU from FASTQ input.

3. **The R package** (R ≥ 4.1):

   ```bash
   R CMD INSTALL .        # or: Rscript -e 'devtools::install()'
   ```

4. **Verify the environment** — this is the portability safety net; run it
   before every first run on a new machine:

   ```r
   library(nanoflow)
   check_tools()   # ✓/✗ per tool with versions; missing tools are listed
   ```

   Tools that live outside the conda env (e.g. dorado, hap.py) can be
   pointed at per-config: `tools: {dorado: /opt/dorado/bin/dorado}`.

5. **References and databases** (user-chosen directory, nothing is
   assumed):

   ```r
   download_references("/data/refs", what = c("genome", "snpeff"))
   ```

6. **Smoke test** on the shipped 200 kb synthetic fixture (~1 minute):

   ```r
   Rscript -e 'testthat::test_package("nanoflow")'
   ```

## Quick start

```r
library(nanoflow)
write_config_template("config.yml")   # edit: reference, output_dir, threads
# samples.csv:  sample,reads[,illumina_r1,illumina_r2,truth_small_vcf,truth_sv_vcf]
run <- run_pipeline("samples.csv", "config.yml")
render_report(run)
```

- **Resume:** rerun the same call after an interruption; steps whose
  outputs exist are skipped (`overwrite = TRUE` forces reruns).
- **Provenance:** `results/provenance.json` records every command line,
  tool version and runtime.
- **Any tool option:** every wrapper has `extra_args`, e.g.
  `sv: {extra_args: ["--all-contigs"]}` (required for contigs < 1 Mb with
  Sniffles ≥ 2.6).
- **Slurm:** `submit_slurm("samples.csv", "config.yml", partition = "...")`
  submits one job per sample; `dry_run = TRUE` just writes the scripts.

## Benchmarking bake-off

Compare callers on the same alignment against a truth set and get one
ranked precision/recall/F1 table:

```r
cmp <- compare_sv_callers(bam, reference, "truth_sv.vcf", "bakeoff",
                          tools = c("sniffles", "svim", "nanovar"),
                          all_contigs = TRUE)
cmp                               # ranked table, best F1 first
write_benchmark_csv(cmp, "sv_callers.csv")
# compare_small_variant_callers(bam, ref, "truth_small.vcf", "bakeoff",
#                               tools = c("clair3", "medaka"))  # via hap.py
```

On the shipped fixture this cleanly surfaces the classic tradeoff — Sniffles
recall 1.0 / precision 0.91 (all 10 SVs, 1 FP) vs SVIM precision 1.0 /
recall 0.6 (no false calls, misses 4). Missing callers are skipped with a
warning, so the table always reflects whatever is installed.

## Vignettes

- `vignette("nanoflow")` — the full workflow on the synthetic fixture.
- `vignette("giab-hg002")` — a documented real-data run on GIAB HG002,
  including truth-set benchmarking.

## Synthetic development data

No real data was available during development, so the package is built and
tested against a deterministic synthetic fixture
(`inst/extdata/sim`, generated by `data-raw/make_fixture.R`): a 200 kb
diploid genome (planted repeats, 250 SNVs, 40 indels, 10 SVs with phased
truth genotypes), ~30× Badread Nanopore reads and ~30× ART Illumina pairs.
Validated by realignment: Sniffles2 recovers 10/10 truth SVs, and an
independent short-read caller recovers 280/290 small variants with 0 false
positives. Basecalling cannot be tested synthetically (Badread produces no
raw signal); its test uses a small public POD5 sample and is skipped
without a GPU. Badread/ART are development-only, never user dependencies.

## License

MIT © Fiona Katharina Meinzer
