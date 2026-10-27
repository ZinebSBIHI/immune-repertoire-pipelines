# ---------------------------------------------------------------------------
# Configuration loading and validation
#
# Every pipeline starts by calling load_config(). Nothing downstream reads a
# path, a sample label, a contrast or a threshold from anywhere else, so a
# pipeline can be moved to a new cohort by editing one YAML file.
#
# Validation is deliberately strict and happens before any data is read. A
# mistyped condition label that is caught here costs a second; the same
# mistake caught by a downstream join produces an NA group that survives
# into a figure.
# ---------------------------------------------------------------------------

#' Load and validate the project configuration
#'
#' @param path Path to the YAML configuration file.
#' @param strict If TRUE (default), unknown top-level keys raise an error.
#'   Set FALSE when extending the configuration with project-specific blocks.
#' @param root Repository root against which relative paths are resolved; by default the parent of the configuration file's directory.
#' @return A validated configuration list, with `.meta` holding the resolved
#'   path and the load time.
#' @export
load_config <- function(path = "config/config.yml", strict = TRUE,
                        root = NULL) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read the configuration. ",
         "Install it with install.packages('yaml').", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("Configuration file not found: ", path, "\n",
         "Copy the template first:  cp config/config.example.yml ", path,
         call. = FALSE)
  }

  # Read as UTF-8 explicitly rather than relying on the session locale. A
  # repository is cloned onto machines whose locale is anything at all, and
  # under a C/POSIX locale a non-ASCII character anywhere in the file -- in a
  # comment, in a display label -- otherwise makes the whole configuration
  # unreadable.
  lines <- readLines(path, encoding = "UTF-8", warn = FALSE)
  cfg <- yaml::yaml.load(paste(lines, collapse = "\n"))
  if (!is.list(cfg)) {
    stop("Configuration did not parse as a mapping: ", path,
         "\nCheck the YAML indentation.", call. = FALSE)
  }

  # Relative paths in the configuration are resolved against the repository
  # root, not the working directory, so a pipeline behaves the same however
  # it was launched.
  if (is.null(root)) root <- dirname(dirname(normalizePath(path)))
  for (p in c("data_dir", "output_dir", "reference_dir", "sample_sheet",
              "seurat_object")) {
    v <- cfg$paths[[p]]
    if (!is.null(v) && !grepl("^(/|~|[A-Za-z]:)", v)) {
      cfg$paths[[p]] <- file.path(root, v)
    }
  }

  cfg <- validate_config(cfg, strict = strict)

  cfg$.meta <- list(
    config_path = normalizePath(path, mustWork = TRUE),
    loaded_at   = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
  cfg
}

#' Validate a configuration list
#'
#' Checks structure, types and internal consistency. Returns the
#' configuration with defaults filled in.
#' @param cfg Configuration list, as parsed from the YAML file.
#' @param strict If TRUE, unknown top-level blocks are an error.
#' @export
validate_config <- function(cfg, strict = TRUE) {
  known_blocks <- c("project", "paths", "cohort", "repertoire",
                    "differential_abundance", "convergence", "bcr", "gex",
                    "specificity", "statistics", "output")

  if (strict) {
    unknown <- setdiff(names(cfg), c(known_blocks, ".meta"))
    if (length(unknown)) {
      stop("Unknown configuration block(s): ", paste(unknown, collapse = ", "),
           "\nKnown blocks: ", paste(known_blocks, collapse = ", "),
           call. = FALSE)
    }
  }

  required <- c("project", "paths", "cohort")
  missing_blocks <- setdiff(required, names(cfg))
  if (length(missing_blocks)) {
    stop("Configuration is missing required block(s): ",
         paste(missing_blocks, collapse = ", "), call. = FALSE)
  }

  # --- project -------------------------------------------------------------
  cfg$project$name <- cfg$project$name %||% "unnamed-project"
  cfg$project$seed <- as.integer(cfg$project$seed %||% 1234L)

  # --- paths ---------------------------------------------------------------
  for (p in c("data_dir", "output_dir")) {
    if (is.null(cfg$paths[[p]])) {
      stop("paths$", p, " must be set in the configuration.", call. = FALSE)
    }
  }
  cfg$paths$reference_dir <- cfg$paths$reference_dir %||% "reference"

  if (!dir.exists(cfg$paths$data_dir)) {
    stop("paths$data_dir does not exist: ", cfg$paths$data_dir,
         "\nInputs are never created by the pipelines; point this at the ",
         "directory that holds them.", call. = FALSE)
  }

  # --- cohort --------------------------------------------------------------
  if (!length(cfg$cohort$subjects %||% character(0))) {
    stop("cohort$subjects must list at least one subject.", call. = FALSE)
  }
  cfg$cohort$subjects   <- as.character(cfg$cohort$subjects)
  cfg$cohort$conditions <- as.character(cfg$cohort$conditions %||% character(0))
  cfg$cohort$timepoints <- as.character(cfg$cohort$timepoints %||% character(0))

  if (anyDuplicated(cfg$cohort$subjects)) {
    stop("cohort$subjects contains duplicates.", call. = FALSE)
  }
  if (anyDuplicated(cfg$cohort$conditions)) {
    stop("cohort$conditions contains duplicates. The order of this vector ",
         "defines factor level order everywhere.", call. = FALSE)
  }

  # Display labels must refer to declared conditions, or a legend will
  # silently fall back to the raw label for the one that was mistyped.
  lab <- cfg$cohort$condition_labels
  if (length(lab)) {
    unknown_lab <- setdiff(names(lab), cfg$cohort$conditions)
    if (length(unknown_lab)) {
      stop("cohort$condition_labels names not found in cohort$conditions: ",
           paste(unknown_lab, collapse = ", "), call. = FALSE)
    }
  }

  # --- repertoire ----------------------------------------------------------
  rep <- cfg$repertoire %||% list()
  rep$locus                <- rep$locus %||% "TRB"
  rep$clonotype_definition <- rep$clonotype_definition %||% "strict"
  rep$abundance_column     <- rep$abundance_column %||% "umi_count"
  rep$min_abundance        <- as.numeric(rep$min_abundance %||% 1)
  rep$min_abundance_sweep  <- as.numeric(rep$min_abundance_sweep %||% c(1, 2, 3, 5))
  rep$diversity_indices    <- as.character(rep$diversity_indices %||%
                                             c("chao1", "shannon", "inv_simpson"))
  rep$rarefaction          <- rep$rarefaction %||% list(enabled = TRUE,
                                                        n_bootstrap = 200)
  rep$rarefaction$n_bootstrap <- as.integer(rep$rarefaction$n_bootstrap %||% 200L)

  valid_def <- c("aa", "nt", "gene", "strict")
  if (!rep$clonotype_definition %in% valid_def) {
    stop("repertoire$clonotype_definition must be one of: ",
         paste(valid_def, collapse = ", "),
         " (got '", rep$clonotype_definition, "')", call. = FALSE)
  }
  if (rep$min_abundance < 1) {
    stop("repertoire$min_abundance must be at least 1.", call. = FALSE)
  }
  cfg$repertoire <- rep

  # --- differential abundance ---------------------------------------------
  da <- cfg$differential_abundance %||% list()
  da$fold_change_threshold <- as.numeric(da$fold_change_threshold %||% 2)
  da$fdr_threshold         <- as.numeric(da$fdr_threshold %||% 0.05)
  da$fdr_method            <- da$fdr_method %||% "BH"
  da$detection_floor       <- as.numeric(da$detection_floor %||% rep$min_abundance)

  if (!is.null(da$reference_timepoint) && length(cfg$cohort$timepoints) &&
      !da$reference_timepoint %in% cfg$cohort$timepoints) {
    stop("differential_abundance$reference_timepoint ('",
         da$reference_timepoint, "') is not among cohort$timepoints.",
         call. = FALSE)
  }
  if (da$fdr_threshold <= 0 || da$fdr_threshold >= 1) {
    stop("differential_abundance$fdr_threshold must lie strictly between ",
         "0 and 1.", call. = FALSE)
  }
  if (da$fold_change_threshold <= 1) {
    stop("differential_abundance$fold_change_threshold must be greater ",
         "than 1 (got ", da$fold_change_threshold, ").\nAt 1 every ",
         "clonotype satisfies both the expansion and the contraction ",
         "criterion, so the categories stop being exclusive.",
         call. = FALSE)
  }
  if (da$detection_floor != rep$min_abundance) {
    message("Note: differential_abundance$detection_floor (", da$detection_floor,
            ") differs from repertoire$min_abundance (", rep$min_abundance,
            "). De novo calls and repertoire counts will therefore use ",
            "different floors; state both in Methods.")
  }
  cfg$differential_abundance <- da

  # --- convergence ---------------------------------------------------------
  cv <- cfg$convergence %||% list()
  cv$cdr3_identity  <- as.numeric(cv$cdr3_identity %||% 0.75)
  cv$min_subjects   <- as.integer(cv$min_subjects %||% 2L)
  cv$n_permutations <- as.integer(cv$n_permutations %||% 200L)
  cv$distance_normalisation <- cv$distance_normalisation %||% "shorter"
  cv$require_mutated_or_switched <- isTRUE(cv$require_mutated_or_switched %||% TRUE)
  cv$shm_threshold  <- as.numeric(cv$shm_threshold %||% 0.02)

  if (cv$cdr3_identity <= 0 || cv$cdr3_identity > 1) {
    stop("convergence$cdr3_identity must lie in (0, 1].", call. = FALSE)
  }
  if (cv$min_subjects < 2) {
    stop("convergence$min_subjects must be at least 2: a group confined to ",
         "one subject is clonal relatedness, not convergence.", call. = FALSE)
  }
  cfg$convergence <- cv

  # --- bcr -----------------------------------------------------------------
  bcr <- cfg$bcr %||% list()
  bcr$clone_threshold_mode  <- bcr$clone_threshold_mode %||% "estimate"
  bcr$clone_threshold_fixed <- as.numeric(bcr$clone_threshold_fixed %||% 0.15)
  bcr$clone_method  <- bcr$clone_method %||% "nt"
  bcr$clone_linkage <- bcr$clone_linkage %||% "single"
  bcr$clone_fields  <- as.character(bcr$clone_fields %||% "subject_id")
  if (!bcr$clone_threshold_mode %in% c("estimate", "fixed")) {
    stop("bcr$clone_threshold_mode must be 'estimate' or 'fixed'.",
         call. = FALSE)
  }
  bcr$cell_resolution <- bcr$cell_resolution %||% list()
  bcr$cell_resolution$prefer_switched_isotype <-
    isTRUE(bcr$cell_resolution$prefer_switched_isotype %||% TRUE)
  bcr$cell_resolution$aggregate_umi <-
    isTRUE(bcr$cell_resolution$aggregate_umi %||% TRUE)
  bcr$cell_resolution$doublet_min_umi <-
    as.numeric(bcr$cell_resolution$doublet_min_umi %||% 10)
  bcr$shm <- bcr$shm %||% list()
  bcr$shm$regions <- as.character(bcr$shm$regions %||% c("IMGT_V", "IMGT_VDJ"))
  bcr$shm$expansion_min_cells_resting <-
    as.integer(bcr$shm$expansion_min_cells_resting %||% 2L)
  bcr$shm$expansion_min_cells_stimulated <-
    as.integer(bcr$shm$expansion_min_cells_stimulated %||% 4L)
  bcr$report_unswitched_combined <-
    isTRUE(bcr$report_unswitched_combined %||% TRUE)
  cfg$bcr <- bcr

  # --- gex -----------------------------------------------------------------
  gex <- cfg$gex %||% list()
  gex$mito_pattern <- gex$mito_pattern %||% "^MT-"
  gex$qc <- gex$qc %||% list()
  if (is.null(gex$qc$default)) {
    gex$qc$default <- list(min_features = 500, max_features = 4000,
                           max_counts = 10000, max_percent_mt = 5)
  }
  gex$normalisation <- gex$normalisation %||% list()
  gex$normalisation$method <- gex$normalisation$method %||% "LogNormalize"
  gex$normalisation$scale_factor <-
    as.numeric(gex$normalisation$scale_factor %||% 1e4)
  gex$normalisation$n_variable_features <-
    as.integer(gex$normalisation$n_variable_features %||% 2000L)
  gex$integration <- gex$integration %||% list(method = "none")
  gex$integration$method <- gex$integration$method %||% "none"
  if (!gex$integration$method %in% c("none", "harmony")) {
    stop("gex$integration$method must be 'none' or 'harmony'.", call. = FALSE)
  }
  if (identical(gex$integration$method, "harmony") &&
      is.null(gex$integration$batch_variable)) {
    stop("gex$integration$method is 'harmony' but ",
         "gex$integration$batch_variable is not set. Name the variable the ",
         "correction is applied over: it has to be a decision, since ",
         "correcting on a variable that coincides with the condition ",
         "removes the effect under study.", call. = FALSE)
  }
  gex$clustering <- gex$clustering %||% list()
  gex$clustering$n_pcs      <- as.integer(gex$clustering$n_pcs %||% 30L)
  gex$clustering$resolution <- as.numeric(gex$clustering$resolution %||% 0.5)
  gex$clustering$sweep_n_pcs <-
    as.integer(gex$clustering$sweep_n_pcs %||% c(10L, 20L, 30L))
  gex$clustering$sweep_resolution <-
    as.numeric(gex$clustering$sweep_resolution %||% c(0.1, 0.3, 0.5))
  gex$max_single_subject_fraction <-
    as.numeric(gex$max_single_subject_fraction %||% 0.80)
  cfg$gex <- gex

  # --- statistics ----------------------------------------------------------
  st <- cfg$statistics %||% list()
  st$paired_test        <- st$paired_test %||% "wilcoxon"
  st$alpha              <- as.numeric(st$alpha %||% 0.05)
  st$report_attainable_p <- isTRUE(st$report_attainable_p %||% TRUE)
  cfg$statistics <- st

  # --- output --------------------------------------------------------------
  out <- cfg$output %||% list()
  out$formats <- as.character(out$formats %||% c("pdf", "png"))
  out$dpi     <- as.numeric(out$dpi %||% 300)
  out$write_session_info     <- isTRUE(out$write_session_info %||% TRUE)
  out$write_resolved_config  <- isTRUE(out$write_resolved_config %||% TRUE)
  bad_fmt <- setdiff(out$formats, c("pdf", "png", "svg", "tiff"))
  if (length(bad_fmt)) {
    stop("output$formats contains unsupported format(s): ",
         paste(bad_fmt, collapse = ", "), call. = FALSE)
  }
  cfg$output <- out

  cfg
}

#' Null-coalescing operator
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Resolve display labels for a factor
#'
#' @param values Character vector of raw labels.
#' @param labels Named list mapping raw label to display label.
#' @return Character vector of display labels; raw label where no mapping
#'   exists.
#' @export
display_labels <- function(values, labels = NULL) {
  values <- as.character(values)
  if (!length(labels)) return(values)
  vapply(values, function(v) {
    lab <- labels[[v]]
    if (is.null(lab) || !length(lab) || is.na(lab)) v else as.character(lab)[1]
  }, character(1), USE.NAMES = FALSE)
}

#' Turn a column into an ordered factor following the configured order
#'
#' Any value absent from `levels` is an error rather than an NA: an
#' unexpected group label is a data or configuration problem, and silently
#' dropping it from a figure hides it.
#' @param x Values to convert.
#' @param levels Declared levels, in display order.
#' @param what What the values are, named in the error message.
#' @export
as_configured_factor <- function(x, levels, what = "value") {
  x <- as.character(x)
  extra <- setdiff(unique(x[!is.na(x)]), levels)
  if (length(extra)) {
    stop("Unexpected ", what, "(s) not declared in the configuration: ",
         paste(sort(extra), collapse = ", "),
         "\nDeclared: ", paste(levels, collapse = ", "), call. = FALSE)
  }
  factor(x, levels = levels)
}

#' Parse metadata out of input file names
#'
#' @param files Character vector of file paths.
#' @param patterns Named list of single-capture-group regexes, as given by
#'   `cohort$filename_pattern`.
#' @param required Fields that must match for every file.
#' @return A data frame with one row per file and one column per pattern.
#' @export
parse_filename_metadata <- function(files, patterns, required = names(patterns)) {
  stopifnot(length(files) > 0L)
  base <- basename(files)

  out <- data.frame(file = files, basename = base, stringsAsFactors = FALSE)
  for (field in names(patterns)) {
    value <- rep(NA_character_, length(base))
    hit <- regexpr(patterns[[field]], base, perl = TRUE) > 0L
    if (any(hit)) {
      # Take the first capture group where the pattern defines one, so that
      # separators written into the pattern (e.g. "_(t[0-9]+)") do not end up
      # in the parsed value.
      grp <- regmatches(base[hit],
                        regexec(patterns[[field]], base[hit], perl = TRUE))
      value[hit] <- vapply(grp,
                           function(g) if (length(g) >= 2L) g[[2]] else g[[1]],
                           character(1))
    }
    out[[field]] <- value
  }

  for (field in intersect(required, names(patterns))) {
    bad <- out$basename[is.na(out[[field]])]
    if (length(bad)) {
      stop("Could not parse '", field, "' from ", length(bad), " file name(s):\n  ",
           paste(utils::head(bad, 10), collapse = "\n  "),
           if (length(bad) > 10) "\n  ..." else "",
           "\nPattern used: ", patterns[[field]],
           "\nFix cohort$filename_pattern in the configuration, or rename the ",
           "inputs.", call. = FALSE)
    }
  }
  out
}
