# Phase 2 benchmarking bake-off: run several callers over the SAME alignment,
# benchmark each against a truth set, and return one tidy comparison table.
# Each caller/benchmark reuses the existing step wrappers, so the bake-off
# inherits their resume, logging and extra_args behaviour.

# Run one caller + its truth benchmark, returning a one-row data.frame (or
# NULL if the tool is missing or the step fails -- the bake-off then simply
# omits that caller instead of aborting).
bakeoff_one <- function(label, call_fn, bench_fn, out_dir, overwrite) {
  step <- tryCatch(call_fn(out_dir), error = function(e) {
    warning(sprintf("[compare] caller '%s' skipped: %s", label,
                    conditionMessage(e)), call. = FALSE)
    NULL
  })
  if (is.null(step) || is.null(step$outputs$vcf) ||
      !file.exists(step$outputs$vcf)) {
    return(list(row = NULL, vcf = NULL, call_step = step, bench_step = NULL))
  }
  bench <- tryCatch(bench_fn(step$outputs$vcf, file.path(out_dir, "benchmark")),
                    error = function(e) {
    warning(sprintf("[compare] benchmark for '%s' failed: %s", label,
                    conditionMessage(e)), call. = FALSE)
    NULL
  })
  m <- bench$params$metrics
  row <- data.frame(
    tool = label,
    precision = as.numeric(m$precision %||% NA),
    recall = as.numeric(m$recall %||% NA),
    f1 = as.numeric(m$f1 %||% NA),
    TP = as.numeric(m$TP %||% NA),
    FP = as.numeric(m$FP %||% NA),
    FN = as.numeric(m$FN %||% NA),
    caller_sec = round(step$runtime_sec %||% NA_real_, 1),
    stringsAsFactors = FALSE)
  list(row = row, vcf = step$outputs$vcf, call_step = step, bench_step = bench)
}

finalize_bakeoff <- function(results, kind) {
  results <- Filter(function(r) !is.null(r$row), results)
  comparison <- if (length(results)) {
    df <- do.call(rbind, lapply(results, `[[`, "row"))
    df[order(-df$f1, -df$recall), , drop = FALSE]
  } else {
    data.frame()
  }
  rownames(comparison) <- NULL
  structure(list(
    kind = kind,
    comparison = comparison,
    vcfs = stats::setNames(lapply(results, `[[`, "vcf"),
                           vapply(results, function(r) r$row$tool, "")),
    steps = results
  ), class = "nanoflow_benchmark")
}

#' Compare structural-variant callers against a truth set
#'
#' Runs each requested SV caller on the same BAM, benchmarks every call set
#' against the truth VCF with Truvari, and returns one tidy table of
#' precision / recall / F1 (plus TP/FP/FN and caller runtime), sorted best
#' F1 first. Callers whose binary is missing or that error out are reported
#' as a warning and omitted, so the bake-off still produces a table for
#' whatever is installed. Optionally also benchmarks the SURVIVOR-merged
#' consensus of all successful callers, shown as an extra `survivor_merge`
#' row.
#'
#' @param bam Sorted, indexed Nanopore BAM.
#' @param reference Reference genome FASTA.
#' @param truth_vcf Truth SV VCF.
#' @param out_dir Output directory (one sub-directory per caller).
#' @param tools SV callers to run; any of `"sniffles"`, `"svim"`,
#'   `"nanovar"`.
#' @param sample Sample ID.
#' @param merge Also benchmark the SURVIVOR merge of the successful callers
#'   (requires >= 2 callers and the SURVIVOR binary).
#' @param min_callers SURVIVOR: minimum callers supporting a merged call.
#' @param all_contigs Pass `--all-contigs` to Sniffles (needed for contigs
#'   < 1 Mb, e.g. the synthetic fixture).
#' @param dup_to_ins Truvari `--dup-to-ins` (match DUP calls as INS).
#' @param threads CPU threads.
#' @param extra_args Named list of per-tool extra argument vectors, e.g.
#'   `list(sniffles = "--minsvlen 30")`.
#' @param config Optional nanoflow config (binary overrides).
#' @param overwrite Rerun steps even if their outputs exist.
#' @return A `nanoflow_benchmark` object with `$comparison` (data.frame),
#'   `$vcfs` (named call sets) and `$steps`.
#' @examples
#' \dontrun{
#' cmp <- compare_sv_callers(bam, ref, "truth_sv.vcf", "bakeoff",
#'                           tools = c("sniffles", "svim"), all_contigs = TRUE)
#' cmp                       # prints the ranked table
#' as.data.frame(cmp)
#' }
#' @export
compare_sv_callers <- function(bam, reference, truth_vcf, out_dir,
                               tools = c("sniffles", "svim", "nanovar"),
                               sample = "sample", merge = TRUE,
                               min_callers = 2, all_contigs = FALSE,
                               dup_to_ins = TRUE, threads = 4,
                               extra_args = list(), config = NULL,
                               overwrite = FALSE) {
  assert_file(bam, "BAM")
  assert_file(reference, "reference FASTA")
  assert_file(truth_vcf, "truth SV VCF")
  dir_create(out_dir)
  ea <- function(tool) extra_args[[tool]] %||% character()

  callers <- list(
    sniffles = function(od) call_sniffles(bam, reference, od, sample,
      all_contigs = all_contigs, threads = threads, extra_args = ea("sniffles"),
      config = config, overwrite = overwrite),
    svim = function(od) call_svim(bam, reference, od, sample,
      extra_args = ea("svim"), config = config, overwrite = overwrite),
    nanovar = function(od) call_nanovar(bam, reference, od, sample, threads,
      extra_args = ea("nanovar"), config = config, overwrite = overwrite)
  )
  unknown <- setdiff(tools, names(callers))
  if (length(unknown)) {
    stop("unknown SV caller(s): ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }

  bench_fn <- function(query_vcf, bdir) {
    benchmark_truvari(truth_vcf, query_vcf, reference, bdir, sample = sample,
                      dup_to_ins = dup_to_ins, config = config,
                      overwrite = overwrite)
  }

  results <- lapply(tools, function(tool) {
    bakeoff_one(tool, callers[[tool]], bench_fn,
                file.path(out_dir, tool), overwrite)
  })

  # SURVIVOR merge of the callers that produced a VCF.
  got <- Filter(function(r) !is.null(r$vcf), results)
  if (merge && length(got) >= 2 && nzchar(Sys.which(
        config$tools$survivor %||% "SURVIVOR"))) {
    merged_dir <- file.path(out_dir, "survivor_merge")
    merge_step <- tryCatch(
      merge_survivor(vapply(got, `[[`, "", "vcf"),
                     file.path(merged_dir, paste0(sample, ".merged.vcf")),
                     min_callers = min_callers, config = config,
                     overwrite = overwrite),
      error = function(e) {
        warning("[compare] SURVIVOR merge skipped: ", conditionMessage(e),
                call. = FALSE)
        NULL
      })
    if (!is.null(merge_step)) {
      results <- c(results, list(bakeoff_one(
        "survivor_merge",
        function(od) merge_step,
        bench_fn, merged_dir, overwrite)))
    }
  }

  finalize_bakeoff(results, "structural variants")
}

#' Compare small-variant callers against a truth set
#'
#' Runs each requested SNV/indel caller on the same BAM and benchmarks every
#' call set against the truth VCF with hap.py, returning one tidy table of
#' precision / recall / F1 by variant type. Like [compare_sv_callers()],
#' missing or failing callers are warned about and omitted.
#'
#' @param bam Sorted, indexed Nanopore BAM.
#' @param reference Reference genome FASTA.
#' @param truth_vcf Truth small-variant VCF.
#' @param out_dir Output directory (one sub-directory per caller).
#' @param tools Callers to run; any of `"clair3"`, `"medaka"`.
#' @param sample Sample ID.
#' @param clair3_model,clair3_model_dir Clair3 model name and directory
#'   (see [call_clair3()]).
#' @param confident_bed Optional confident-regions BED for hap.py.
#' @param threads CPU threads.
#' @param extra_args Named list of per-tool extra argument vectors.
#' @param config Optional nanoflow config (binary overrides).
#' @param overwrite Rerun steps even if their outputs exist.
#' @return A `nanoflow_benchmark` object. The `$comparison` table reports
#'   hap.py's overall (type = `"ALL"`) precision/recall/F1 per caller; the
#'   full per-type hap.py tables remain in each caller's `$steps` entry.
#' @examples
#' \dontrun{
#' compare_small_variant_callers(bam, ref, "truth_small.vcf", "bakeoff",
#'                               tools = c("clair3", "medaka"))
#' }
#' @export
compare_small_variant_callers <- function(bam, reference, truth_vcf, out_dir,
                                          tools = c("clair3", "medaka"),
                                          sample = "sample",
                                          clair3_model = "r1041_e82_400bps_sup_v420",
                                          clair3_model_dir = NULL,
                                          confident_bed = NULL, threads = 4,
                                          extra_args = list(), config = NULL,
                                          overwrite = FALSE) {
  assert_file(bam, "BAM")
  assert_file(reference, "reference FASTA")
  assert_file(truth_vcf, "truth small-variant VCF")
  dir_create(out_dir)
  ea <- function(tool) extra_args[[tool]] %||% character()

  callers <- list(
    clair3 = function(od) call_clair3(bam, reference, od, sample,
      model_dir = clair3_model_dir, model = clair3_model,
      include_all_ctgs = TRUE, threads = threads, extra_args = ea("clair3"),
      config = config, overwrite = overwrite),
    medaka = function(od) call_medaka(bam, reference, od, sample,
      threads = threads, extra_args = ea("medaka"), config = config,
      overwrite = overwrite)
  )
  unknown <- setdiff(tools, names(callers))
  if (length(unknown)) {
    stop("unknown small-variant caller(s): ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }

  bench_fn <- function(query_vcf, bdir) {
    step <- benchmark_happy(truth_vcf, query_vcf, reference, bdir,
                            sample = sample, confident_bed = confident_bed,
                            threads = threads, config = config,
                            overwrite = overwrite)
    # reduce hap.py's per-type metrics to the overall ALL/indel/SNV summary
    m <- step$params$metrics
    if (!is.null(m) && "Type" %in% names(m)) {
      overall <- m[m$Type %in% c("INDEL", "SNP"), , drop = FALSE]
      rec <- stats::weighted.mean(overall$METRIC.Recall,
                                  overall$TRUTH.TOTAL, na.rm = TRUE)
      prec <- stats::weighted.mean(overall$METRIC.Precision,
                                   overall$TRUTH.TOTAL, na.rm = TRUE)
      step$params$metrics <- data.frame(
        precision = prec, recall = rec,
        f1 = 2 * prec * rec / (prec + rec),
        TP = sum(overall$TRUTH.TP, na.rm = TRUE),
        FP = sum(overall$QUERY.FP, na.rm = TRUE),
        FN = sum(overall$TRUTH.FN, na.rm = TRUE))
    }
    step
  }

  results <- lapply(tools, function(tool) {
    bakeoff_one(tool, callers[[tool]], bench_fn,
                file.path(out_dir, tool), overwrite)
  })
  finalize_bakeoff(results, "small variants")
}

#' @export
as.data.frame.nanoflow_benchmark <- function(x, ...) x$comparison

#' @export
print.nanoflow_benchmark <- function(x, ...) {
  cat(sprintf("<nanoflow benchmark: %s> %d caller(s)\n", x$kind,
              nrow(x$comparison)))
  if (!nrow(x$comparison)) {
    cat("  (no caller produced a benchmarkable call set)\n")
    return(invisible(x))
  }
  df <- x$comparison
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], function(v) formatC(v, format = "g", digits = 4))
  print.data.frame(df, row.names = FALSE, right = FALSE)
  best <- x$comparison$tool[1]
  cat(sprintf("  best F1: %s\n", best))
  invisible(x)
}

#' Write a benchmark comparison table to CSV
#'
#' @param x A `nanoflow_benchmark` from [compare_sv_callers()] or
#'   [compare_small_variant_callers()].
#' @param path Output CSV path.
#' @return `path`, invisibly.
#' @export
write_benchmark_csv <- function(x, path) {
  stopifnot(inherits(x, "nanoflow_benchmark"))
  dir_create(dirname(path))
  utils::write.csv(x$comparison, path, row.names = FALSE)
  invisible(path)
}
