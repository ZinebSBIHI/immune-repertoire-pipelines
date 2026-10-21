#' immunerepertoire: shared core for immune repertoire and single-cell pipelines
#'
#' The functions shared by three companion pipelines:
#'
#' * bulk-repertoire-pipeline: bulk TCR or BCR repertoires;
#' * single-cell-repertoire-pipeline: single-cell TCR or BCR repertoires;
#' * single-cell-gex-pipeline: single-cell gene expression.
#'
#' Function families are grouped by file prefix: `core-` (configuration,
#' run directories, figures), `rep-` (clonotypes, diversity, rarefaction,
#' shared analysis sections), `ana-` (differential abundance, convergence,
#' statistics), `bcr-` (B-cell receptor processing).
#'
#' Every analytical decision that is easy to get wrong in repertoire work is
#' argued in `docs/design_notes.md` of the source repository, with the code
#' that implements it.
#'
#' @keywords internal
#' @importFrom rlang .data
"_PACKAGE"

#' Arguments shared by the section functions
#'
#' Every `section_*()` function writes one part of a pipeline's report: it
#' adds its tables and figures to the run directory and its numbers to the
#' report, and returns its results invisibly so the pipeline can pass them
#' on.
#'
#' @param db AIRR Rearrangement data frame, one row per rearrangement (per
#'   chain for single-cell data).
#' @param clono Clonotype table: `sample_id`, `subject_id`, `clonotype`,
#'   `abundance`, `frequency` and the sample metadata (output of
#'   `collapse_clonotypes()` joined to the sample table).
#' @param sample_meta One row per sample: `sample_id`, `subject_id` and the
#'   timepoint or condition columns the design has.
#' @param abundance_col Abundance column of `db` (UMI or read counts, named
#'   explicitly because the two are not interchangeable).
#' @param umi_col UMI count column of `db`.
#' @param locus Analysed locus or loci (e.g. "TRB", "IGH").
#' @param run_dir Run directory from `init_run()`.
#' @param cfg Validated configuration from `load_config()`.
#' @param report Report writer from `open_report()`.
#' @param unit Unit of observation shown in labels ("molecules" for bulk,
#'   "cells" for single-cell data).
#' @param within Grouping inside which a common depth or reference is taken
#'   (usually "subject_id"), or NULL for the whole cohort.
#' @param prefix File-name prefix of the section's outputs, which orders
#'   them in the run directory.
#' @name section_arguments
#' @keywords internal
NULL
