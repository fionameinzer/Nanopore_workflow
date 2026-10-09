# Configuration handling. All locations (reference, databases, binaries,
# output directory, scratch) come from here -- never from hard-coded paths.

#' Default nanoflow configuration
#'
#' Returns the full default configuration as a nested list. Every value can
#' be overridden by a YAML config file ([read_config()]) or an R list.
#' Defaults target human GRCh38 but nothing is organism-specific: the
#' reference FASTA, tool models and annotation databases are all
#' configurable for other diploid organisms.
#'
#' @return Nested named list with entries `reference`, `output_dir`,
#'   `scratch_dir`, `threads`, `memory_gb`, `tools` (binary overrides) and
#'   `steps` (per-step switches, tool selection and `extra_args`).
#' @export
default_config <- function() {
  list(
    reference = list(
      fasta = NULL,               # path to reference genome FASTA
      gff3 = NULL,                # annotation for custom SnpEff databases
      snpeff_db = "GRCh38.105",   # prebuilt SnpEff database name
      snpeff_data_dir = NULL,     # where SnpEff databases live
      annotsv_dir = NULL          # AnnotSV annotation directory
    ),
    output_dir = "nanoflow_results",
    scratch_dir = NULL,           # defaults to tempdir() when NULL
    threads = 4,
    memory_gb = 16,
    tools = list(),               # per-tool absolute binary overrides
    steps = list(
      basecall = list(run = "auto", tool = "dorado",
                      model = "sup", extra_args = character()),
      qc = list(run = TRUE, min_length = 500, keep_percent = 90,
                target_bases = NULL, extra_args = character()),
      align = list(tool = "minimap2", preset = "map-ont",
                   extra_args = character()),
      assembly = list(run = FALSE, tool = "flye",
                      polish = c("racon"), genome_size = NULL,
                      extra_args = character()),
      small_variants = list(run = TRUE, tool = "clair3",
                            model = "r1041_e82_400bps_sup_v420",
                            clair3_model_dir = NULL,
                            extra_args = character()),
      sv = list(run = TRUE, tool = "sniffles", extra_args = character()),
      phase = list(run = TRUE, tool = "whatshap", extra_args = character()),
      annotate = list(run = FALSE, tool = "snpeff",
                      extra_args = character()),
      benchmark = list(run = "auto", extra_args = character()),
      report = list(run = TRUE, igv_regions = character())
    )
  )
}

#' Read a nanoflow configuration
#'
#' Reads a YAML file (or takes an R list) and merges it recursively over
#' [default_config()], so a config file only needs the values that differ
#' from the defaults.
#'
#' @param config Path to a YAML file, a named list of overrides, or `NULL`
#'   for pure defaults.
#' @return Validated config list (class `nanoflow_config`).
#' @examples
#' cfg <- read_config(list(threads = 8, steps = list(sv = list(
#'   extra_args = "--all-contigs"))))
#' cfg$threads
#' @export
read_config <- function(config = NULL) {
  base <- default_config()
  user <- if (is.null(config)) {
    list()
  } else if (is.character(config)) {
    assert_file(config, "config file")
    yaml::read_yaml(config)
  } else if (is.list(config)) {
    config
  } else {
    stop("config must be a YAML file path, a list, or NULL", call. = FALSE)
  }
  cfg <- utils::modifyList(base, user)
  validate_config(cfg)
  structure(cfg, class = c("nanoflow_config", "list"))
}

#' @keywords internal
validate_config <- function(cfg) {
  stopifnot(is.numeric(cfg$threads), cfg$threads >= 1,
            is.numeric(cfg$memory_gb), cfg$memory_gb >= 1)
  if (!is.null(cfg$reference$fasta)) {
    assert_file(cfg$reference$fasta, "reference FASTA")
  }
  for (nm in names(cfg$tools)) {
    if (!is.null(cfg$tools[[nm]])) assert_file(cfg$tools[[nm]],
                                               paste0("tools$", nm, " binary"))
  }
  invisible(cfg)
}

#' Write a commented template config file
#'
#' @param path Where to write the YAML template.
#' @param overwrite Overwrite an existing file?
#' @return The path, invisibly.
#' @export
write_config_template <- function(path = "nanoflow_config.yml",
                                  overwrite = FALSE) {
  if (file.exists(path) && !overwrite) {
    stop(path, " exists; use overwrite = TRUE", call. = FALSE)
  }
  header <- c(
    "# nanoflow configuration.",
    "# Only values that differ from the defaults are needed; everything",
    "# else falls back to default_config(). No path below is required to",
    "# exist until the step that uses it actually runs.",
    "")
  writeLines(c(header, yaml::as.yaml(default_config())), path)
  invisible(path)
}
