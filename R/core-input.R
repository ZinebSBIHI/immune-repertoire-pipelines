# ---------------------------------------------------------------------------
# Reading the inputs of the three pipelines
#
# Sample metadata come from two sources, combined in this order:
#
#   1. the input file names, parsed with cohort$filename_pattern (one regular
#      expression with one capture group per field: subject, timepoint,
#      condition, or any other name);
#   2. an optional sample sheet (paths$sample_sheet), a CSV with a sample_id
#      column and any other columns -- compartment, stimulation, batch --
#      which override and extend what the file names give.
#
# Either is enough on its own. The sample sheet is the place for metadata
# that does not belong in a file name, and keeps it out of the code.
# ---------------------------------------------------------------------------

#' Read a sample sheet
#'
#' @param path CSV file with a `sample_id` column; every other column is
#'   sample metadata. A `subject` column is accepted for `subject_id`.
#' @return Data frame, one row per sample, all columns as character.
#' @export
read_sample_sheet <- function(path) {
  if (!file.exists(path)) stop("Sample sheet not found: ", path, call. = FALSE)
  x <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                       colClasses = "character", na.strings = c("", "NA"),
                       strip.white = TRUE)
  if ("subject" %in% names(x) && !"subject_id" %in% names(x)) {
    names(x)[names(x) == "subject"] <- "subject_id"
  }
  if (!"sample_id" %in% names(x)) {
    stop("The sample sheet must have a 'sample_id' column: ", path,
         call. = FALSE)
  }
  dup <- unique(x$sample_id[duplicated(x$sample_id)])
  if (length(dup)) {
    stop("Duplicated sample_id in the sample sheet: ",
         paste(dup, collapse = ", "), call. = FALSE)
  }
  x
}

#' Metadata of a set of samples, from file names and the sample sheet
#'
#' @param sample_ids Sample identifiers, one per input.
#' @param cfg Validated configuration.
#' @param names_to_parse Strings the filename patterns are applied to (the
#'   file names by default).
#' @param required Metadata fields that must be known for every sample.
#' @return Data frame: `sample_id` and one column per metadata field, with
#'   the filename field `subject` renamed `subject_id`.
#' @export
sample_metadata <- function(sample_ids, cfg, names_to_parse = sample_ids,
                            required = "subject_id") {
  meta <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
  pat <- cfg$cohort$filename_pattern
  if (length(pat)) {
    parsed <- parse_filename_metadata(names_to_parse, pat, required = character(0))
    for (f in names(pat)) {
      col <- if (identical(f, "subject")) "subject_id" else f
      meta[[col]] <- parsed[[f]]
    }
  }
  sheet_path <- cfg$paths$sample_sheet
  if (!is.null(sheet_path)) {
    sheet <- read_sample_sheet(sheet_path)
    missing_ids <- setdiff(sample_ids, sheet$sample_id)
    if (length(missing_ids)) {
      stop("Samples absent from the sample sheet: ",
           paste(utils::head(missing_ids, 10), collapse = ", "),
           if (length(missing_ids) > 10) " ..." else "", call. = FALSE)
    }
    idx <- match(meta$sample_id, sheet$sample_id)
    for (col in setdiff(names(sheet), "sample_id")) {
      meta[[col]] <- sheet[[col]][idx]
    }
  }
  for (col in required) {
    if (!col %in% names(meta) || anyNA(meta[[col]])) {
      bad <- if (col %in% names(meta)) meta$sample_id[is.na(meta[[col]])] else
        meta$sample_id
      stop("Metadata field '", col, "' is unknown for ", length(bad),
           " sample(s): ", paste(utils::head(bad, 10), collapse = ", "),
           "\nDeclare it in cohort$filename_pattern or in the sample sheet ",
           "(paths$sample_sheet).", call. = FALSE)
    }
  }
  meta
}

#' Apply the configured factor levels to the cohort columns
#'
#' @param df Data frame with some of `subject_id`, `timepoint`, `condition`.
#' @param cfg Validated configuration.
#' @return `df` with those columns as factors in the configured order.
#' @export
apply_cohort_levels <- function(df, cfg) {
  df$subject_id <- as_configured_factor(df$subject_id, cfg$cohort$subjects,
                                        "subject")
  for (col in c("timepoint", "condition")) {
    lv <- cfg$cohort[[paste0(col, "s")]]
    if (col %in% names(df) && length(lv) && !all(is.na(df[[col]]))) {
      df[[col]] <- as_configured_factor(df[[col]], lv, col)
    }
  }
  df
}

#' Input file stem: the name without directory and AIRR extensions
#' @noRd
airr_stem <- function(path) {
  sub("[.](airr[.])?tsv([.]gz)?$", "", basename(path))
}

#' Read the inputs of a bulk repertoire run
#'
#' One AIRR Rearrangement TSV per sample (MiXCR, IgBLAST + MakeDb, or any
#' AIRR-compliant tool), optionally gzipped. The sample is the file stem.
#'
#' @param cfg Validated configuration.
#' @return List: `db` (all samples, metadata joined), `abundance_col`,
#'   `sample_meta` (one row per sample).
#' @export
read_bulk_inputs <- function(cfg) {
  files <- list.files(cfg$paths$data_dir, pattern = "[.](airr[.])?tsv([.]gz)?$",
                      full.names = TRUE, recursive = TRUE)
  if (!length(files)) {
    stop("No AIRR .tsv files under ", cfg$paths$data_dir,
         "\nSee data/README.md, or run the example script to create a ",
         "synthetic dataset.", call. = FALSE)
  }
  ids <- airr_stem(files)
  if (anyDuplicated(ids)) {
    stop("Two input files share the sample name '",
         ids[duplicated(ids)][1], "'.", call. = FALSE)
  }
  meta <- sample_metadata(ids, cfg, names_to_parse = basename(files))
  db <- do.call(rbind, lapply(seq_along(files), function(i) {
    x <- utils::read.delim(files[i], stringsAsFactors = FALSE,
                           check.names = FALSE, na.strings = c("", "NA"))
    x$sample_id <- ids[i]
    for (col in setdiff(names(meta), "sample_id")) x[[col]] <- meta[[col]][i]
    x
  }))
  abundance <- cfg$repertoire$abundance_column
  if (!abundance %in% names(db)) {
    stop("repertoire$abundance_column is '", abundance, "' but no such ",
         "column is present. Available: ", paste(names(db), collapse = ", "),
         call. = FALSE)
  }
  db <- apply_cohort_levels(db, cfg)
  sample_meta <- apply_cohort_levels(meta, cfg)
  for (col in c("timepoint", "condition")) {
    if (!col %in% names(sample_meta)) sample_meta[[col]] <- NA_character_
    if (!col %in% names(db)) db[[col]] <- NA_character_
  }
  list(db = db, abundance_col = abundance, sample_meta = sample_meta)
}

#' Read the inputs of a single-cell repertoire run
#'
#' Either Cell Ranger V(D)J `filtered_contig_annotations.csv`, one per
#' library directory (the sample is the directory name), or AIRR TSV files
#' carrying a per-cell barcode (`cell_id` or `barcode`; the sample is the
#' file stem). Cell identifiers are made unique across libraries as
#' `<sample>_<barcode>`.
#'
#' @param cfg Validated configuration.
#' @return List: `db` (one row per contig), `abundance_col`, `sample_meta`.
#' @export
read_single_cell_inputs <- function(cfg) {
  csv <- list.files(cfg$paths$data_dir, recursive = TRUE, full.names = TRUE,
                    pattern = "filtered_contig_annotations[.]csv$")
  tsv <- list.files(cfg$paths$data_dir, recursive = TRUE, full.names = TRUE,
                    pattern = "[.](airr[.])?tsv([.]gz)?$")
  if (length(csv)) {
    ids <- basename(dirname(csv))
    if (anyDuplicated(ids)) {
      stop("Two contig files sit in directories of the same name: ",
           ids[duplicated(ids)][1], call. = FALSE)
    }
    meta <- sample_metadata(ids, cfg)
    db <- do.call(rbind, lapply(seq_along(csv), function(i) {
      x <- utils::read.csv(csv[i], stringsAsFactors = FALSE,
                           na.strings = c("", "NA", "None"))
      opt <- function(col) if (col %in% names(x)) x[[col]] else NA_character_
      out <- data.frame(
        sample_id = ids[i],
        cell_id = paste(ids[i], x$barcode, sep = "_"),
        locus = x$chain, v_call = x$v_gene, d_call = opt("d_gene"),
        j_call = x$j_gene, c_call = opt("c_gene"),
        junction = x$cdr3_nt, junction_aa = x$cdr3,
        umi_count = x$umis, productive = x$productive,
        stringsAsFactors = FALSE)
      for (col in setdiff(names(meta), "sample_id")) {
        out[[col]] <- meta[[col]][i]
      }
      out
    }))
    abundance <- "umi_count"
  } else if (length(tsv)) {
    ids <- airr_stem(tsv)
    meta <- sample_metadata(ids, cfg, names_to_parse = basename(tsv))
    db <- do.call(rbind, lapply(seq_along(tsv), function(i) {
      x <- utils::read.delim(tsv[i], stringsAsFactors = FALSE,
                             check.names = FALSE, na.strings = c("", "NA"))
      # A bulk file pointed at this pipeline by mistake has no barcode, and
      # would otherwise have every molecule counted as a cell.
      bc <- intersect(c("cell_id", "barcode"), names(x))
      if (!length(bc)) {
        stop("File has no per-cell barcode column (cell_id or barcode): ",
             basename(tsv[i]), "\nBulk repertoire data belong to the ",
             "bulk-repertoire-pipeline.", call. = FALSE)
      }
      x$sample_id <- ids[i]
      x$cell_id <- paste(ids[i], x[[bc[1]]], sep = "_")
      for (col in setdiff(names(meta), "sample_id")) x[[col]] <- meta[[col]][i]
      x
    }))
    abundance <- cfg$repertoire$abundance_column
    if (!abundance %in% names(db)) {
      stop("repertoire$abundance_column is '", abundance, "' but no such ",
           "column is present. Available: ",
           paste(names(db), collapse = ", "), call. = FALSE)
    }
  } else {
    stop("No input found under ", cfg$paths$data_dir,
         "\nExpected filtered_contig_annotations.csv per library, or AIRR ",
         ".tsv files with a cell barcode. See data/README.md.",
         call. = FALSE)
  }
  db <- apply_cohort_levels(db, cfg)
  sample_meta <- apply_cohort_levels(meta, cfg)
  for (col in c("timepoint", "condition")) {
    if (!col %in% names(sample_meta)) sample_meta[[col]] <- NA_character_
    if (!col %in% names(db)) db[[col]] <- NA_character_
  }
  list(db = db, abundance_col = abundance, sample_meta = sample_meta)
}

#' Read single-cell expression libraries
#'
#' Cell Ranger filtered feature-barcode matrices, either as `.h5` files or
#' as Matrix Market directories (`matrix.mtx(.gz)`, `features.tsv(.gz)` or
#' `genes.tsv`, `barcodes.tsv(.gz)`). The library name is the file stem, or
#' the enclosing directory when Cell Ranger's default name is kept.
#'
#' @param cfg Validated configuration.
#' @param min_cells,min_features Passed to `Seurat::CreateSeuratObject()`.
#' @return List: `objects` (named list of Seurat objects, metadata set) and
#'   `meta` (one row per library).
#' @export
read_gex_inputs <- function(cfg, min_cells = 3, min_features = 200) {
  need(c("Seurat", "SeuratObject"), "read_gex_inputs()")
  dd <- cfg$paths$data_dir
  h5 <- list.files(dd, pattern = "[.]h5$", recursive = TRUE, full.names = TRUE)
  mtx <- list.files(dd, pattern = "^matrix[.]mtx([.]gz)?$", recursive = TRUE,
                    full.names = TRUE)
  mtx_dirs <- unique(dirname(mtx))
  inputs <- c(h5, mtx_dirs)
  if (!length(inputs)) {
    stop("No expression matrix under ", dd, "\nExpected Cell Ranger ",
         ".h5 files or filtered_feature_bc_matrix directories, one per ",
         "library. See data/README.md.", call. = FALSE)
  }
  generic <- c("filtered_feature_bc_matrix", "raw_feature_bc_matrix",
               "sample_filtered_feature_bc_matrix", "outs")
  name_of <- function(p) {
    nm <- sub("[.]h5$", "", basename(p))
    while (nm %in% generic) {
      p <- dirname(p)
      nm <- basename(p)
    }
    nm
  }
  ids <- vapply(inputs, name_of, character(1), USE.NAMES = FALSE)
  if (anyDuplicated(ids)) {
    stop("Two or more matrices resolve to the same library name: ",
         paste(unique(ids[duplicated(ids)]), collapse = ", "),
         "\nLibrary names key the QC thresholds and the cell-name prefixes ",
         "and must be unique.", call. = FALSE)
  }
  meta <- sample_metadata(ids, cfg)
  objects <- stats::setNames(vector("list", length(inputs)), ids)
  for (i in seq_along(inputs)) {
    counts <- if (grepl("[.]h5$", inputs[i])) {
      Seurat::Read10X_h5(inputs[i])
    } else {
      Seurat::Read10X(inputs[i])
    }
    if (is.list(counts)) counts <- counts[["Gene Expression"]]
    o <- Seurat::CreateSeuratObject(counts, project = ids[i],
                                    min.cells = min_cells,
                                    min.features = min_features)
    o$sample_id <- ids[i]
    for (col in setdiff(names(meta), "sample_id")) o[[col]] <- meta[[col]][i]
    objects[[i]] <- o
  }
  list(objects = objects, meta = apply_cohort_levels(meta, cfg))
}
