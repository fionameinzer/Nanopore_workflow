#' nanoflow: end-to-end Oxford Nanopore long-read WGS workflow
#'
#' An orchestration layer over established command-line tools: basecalling
#' (Dorado), QC/trimming (NanoPlot, Porechop, Filtlong, MultiQC), alignment
#' (minimap2), optional de novo assembly (Flye/Canu/Wengan + Racon/Medaka +
#' QUAST), small-variant calling (Clair3/Medaka), SV calling
#' (Sniffles2/SVIM/NanoVar + SURVIVOR + Spectre), phasing
#' (WhatsHap/HapCUT2), annotation (SnpEff/AnnotSV), truth benchmarking
#' (hap.py/Truvari) and HTML reporting.
#'
#' Start with [write_config_template()], [check_tools()] and
#' [run_pipeline()]; `vignette("nanoflow")` walks through the synthetic
#' fixture end to end.
#'
#' @keywords internal
"_PACKAGE"
