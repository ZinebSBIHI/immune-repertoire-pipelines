# ---------------------------------------------------------------------------
# Single-cell BCR: one heavy chain per cell, isotype, clonal QC
#
# A 10x V(D)J library routinely reports more than one heavy chain contig per
# cell barcode. The reasons differ and so should the treatment:
#
#   - ambient mRNA and index hopping give a dominant contig and a minor one,
#     and the cell should be kept with its dominant chain;
#   - a genuine doublet gives two well-supported contigs, and the barcode
#     should be dropped rather than arbitrarily resolved;
#   - IGHM and IGHD are co-transcribed in naive and marginal-zone B cells,
#     so a cell can legitimately carry both with good support.
#
# Keeping the highest-UMI contig without distinguishing these cases silently
# converts doublets into cells and makes the isotype call a function of
# transcript abundance.
#
# ISOTYPE REPORTING. Because one isotype is retained per cell, an IgM
# versus IgD split across naive or marginal-zone cells describes which
# transcript happened to dominate, not two populations. The combined
# unswitched fraction (IGHM or IGHD) is the quantity that means something.
# `summarise_isotypes()` reports it that way by default.
# ---------------------------------------------------------------------------

UNSWITCHED_ISOTYPES <- c("IGHM", "IGHD")

#' Is this isotype call class-switched?
#'
#' An absent call is NOT switched. Testing `!x %in% UNSWITCHED_ISOTYPES`
#' directly would return TRUE for an empty string, so a contig with no
#' C-gene assignment would outrank a well-supported IGHM contig under the
#' switched-isotype preference, and would survive the germline-proximity
#' filter in the convergence analysis. Both are the opposite of intended.
#'
#' @param c_gene Character vector of C-gene calls, alleles already stripped.
#' @return Logical vector; FALSE where the call is missing or empty.
is_switched <- function(c_gene) {
  g <- as.character(c_gene)
  assigned <- !is.na(g) & nzchar(trimws(g))
  assigned & !g %in% UNSWITCHED_ISOTYPES
}

#' Resolve single-cell BCR data to one heavy chain per cell
#'
#' @param db AIRR data frame with `cell_id`, `locus`, `c_call` and a UMI
#'   count column.
#' @param umi_col Name of the UMI count column.
#' @param prefer_switched Place class-switched isotypes ahead of IGHM/IGHD
#'   when choosing which chain represents the cell.
#' @param aggregate_umi Sum UMIs of contigs sharing a junction within a cell
#'   before choosing, so a sequence split across contigs is not
#'   under-weighted.
#' @param doublet_min_umi A cell with two or more heavy chains each at or
#'   above this count is called a doublet and removed.
#' @return List with `db` (resolved, one heavy chain per cell, plus light
#'   chains), `doublets` (removed cell ids) and `filter_log`.
resolve_one_heavy_per_cell <- function(db,
                                       umi_col = "umi_count",
                                       prefer_switched = TRUE,
                                       aggregate_umi = TRUE,
                                       doublet_min_umi = 10) {
  need("dplyr", "resolve_one_heavy_per_cell()")
  for (cl in c("cell_id", "locus", "c_call", umi_col)) {
    if (!cl %in% names(db)) {
      stop("resolve_one_heavy_per_cell() needs column '", cl, "'.",
           call. = FALSE)
    }
  }

  log <- list(step0_input_rows = nrow(db),
              step0_input_cells = dplyr::n_distinct(db$cell_id))

  # Productive only. A non-productive heavy chain cannot be the cell's BCR.
  if ("productive" %in% names(db)) {
    prod <- db$productive
    if (is.character(prod)) prod <- toupper(prod) %in% c("T", "TRUE")
    db <- db[which(as.logical(prod)), , drop = FALSE]
  }
  log$step1_productive_rows <- nrow(db)

  db$.umi <- as.numeric(db[[umi_col]])
  db$.c_gene <- sub("\\*.*$", "", sub(",.*$", "", as.character(db$c_call)))

  heavy <- db[db$locus == "IGH", , drop = FALSE]
  light <- db[db$locus %in% c("IGK", "IGL"), , drop = FALSE]
  log$step2_heavy_rows <- nrow(heavy)
  log$step2_light_rows <- nrow(light)

  if (!nrow(heavy)) {
    stop("No productive IGH rearrangements found.", call. = FALSE)
  }

  # Aggregate UMIs across contigs that carry the same junction in the same
  # cell: assembly can split one molecule population across contigs, and
  # treating them separately understates its support.
  #
  # The representative row for each (cell, junction) is chosen by an
  # explicit ordering, never by position in the file. Several contigs of one
  # junction can carry different C-gene calls, so taking whichever row came
  # first would make the isotype a function of the order the assembler
  # happened to emit, which is neither reproducible nor meaningful. The
  # switched-isotype preference is applied here as well as at the chain
  # level below, so that it is not overridden by this step.
  if (aggregate_umi && "junction" %in% names(heavy)) {
    # The contig's own UMI count is kept: it, not the aggregated total, is
    # what distinguishes the contigs sharing a junction. Ordering on the
    # aggregated value would tie every one of them and fall back on row
    # order, which is the non-determinism this step exists to remove.
    heavy$.umi_contig <- heavy$.umi
    heavy <- dplyr::ungroup(dplyr::mutate(
      dplyr::group_by(heavy, .data$cell_id, .data$junction),
      .umi = sum(.data$.umi, na.rm = TRUE)))

    heavy$.switched <- is_switched(heavy$.c_gene)
    # The C-gene call is the final tiebreak, so that two contigs with equal
    # support still resolve the same way on every run and on every machine.
    ord <- if (prefer_switched) {
      order(heavy$cell_id, heavy$junction, -as.integer(heavy$.switched),
            -heavy$.umi_contig, heavy$.c_gene, na.last = TRUE)
    } else {
      order(heavy$cell_id, heavy$junction, -heavy$.umi_contig,
            heavy$.c_gene, na.last = TRUE)
    }
    heavy <- heavy[ord, , drop = FALSE]
    heavy <- heavy[!duplicated(heavy[, c("cell_id", "junction")]), , drop = FALSE]
    heavy$.umi_contig <- NULL
  }
  log$step3_heavy_rows_after_aggregation <- nrow(heavy)

  # Doublet detection before resolution. A barcode with two well-supported
  # heavy chains is two cells; picking one would make it a cell with an
  # arbitrary BCR.
  support <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(heavy, .data$cell_id),
    n_supported = sum(.data$.umi >= doublet_min_umi),
    .groups = "drop"))
  doublets <- support$cell_id[support$n_supported >= 2L]

  heavy <- heavy[!heavy$cell_id %in% doublets, , drop = FALSE]
  light <- light[!light$cell_id %in% doublets, , drop = FALSE]
  log$step4_doublet_cells_removed <- length(doublets)
  log$step4_heavy_rows <- nrow(heavy)

  # Choose the chain representing each cell: switched isotype first if
  # requested, then UMI support.
  heavy$.switched <- is_switched(heavy$.c_gene)

  ord <- if (prefer_switched) {
    order(heavy$cell_id, -as.integer(heavy$.switched), -heavy$.umi)
  } else {
    order(heavy$cell_id, -heavy$.umi)
  }
  heavy <- heavy[ord, , drop = FALSE]
  heavy <- heavy[!duplicated(heavy$cell_id), , drop = FALSE]
  log$step5_cells_one_heavy <- nrow(heavy)

  # Light chains: keep the best-supported per cell per locus.
  if (nrow(light)) {
    light <- light[order(light$cell_id, light$locus, -light$.umi), , drop = FALSE]
    light <- light[!duplicated(light[, c("cell_id", "locus")]), , drop = FALSE]
  }
  log$step5_light_rows <- nrow(light)

  out <- dplyr::bind_rows(heavy, light)
  out$.switched <- NULL

  message(sprintf(
    paste0("Heavy chain resolution: %d cells in, %d doublets removed, ",
           "%d cells retained (%.1f%%)."),
    log$step0_input_cells, length(doublets), nrow(heavy),
    100 * nrow(heavy) / max(1L, log$step0_input_cells)))

  list(db = out, doublets = doublets,
       filter_log = as.data.frame(log, stringsAsFactors = FALSE))
}

#' Summarise isotype composition
#'
#' @param db Resolved heavy chains (one per cell).
#' @param group_cols Grouping columns.
#' @param combine_unswitched Report IGHM and IGHD together as a single
#'   unswitched fraction. TRUE by default, because with one isotype per cell
#'   an IgM/IgD split is a statement about transcript dominance rather than
#'   about two populations.
#' @return Data frame of counts and fractions per isotype class.
summarise_isotypes <- function(db, group_cols = "sample_id",
                               combine_unswitched = TRUE) {
  need("dplyr", "summarise_isotypes()")
  heavy <- db[db$locus == "IGH", , drop = FALSE]
  heavy$isotype <- sub("\\*.*$", "", sub(",.*$", "",
                                         as.character(heavy$c_call)))

  unassigned <- is.na(heavy$isotype) | !nzchar(trimws(heavy$isotype))
  if (combine_unswitched) {
    heavy$isotype_class <- ifelse(
      unassigned, "unassigned",
      ifelse(heavy$isotype %in% UNSWITCHED_ISOTYPES,
             "unswitched_IGHM_or_IGHD", heavy$isotype))
  } else {
    heavy$isotype_class <- ifelse(unassigned, "unassigned", heavy$isotype)
  }

  out <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(heavy, dplyr::across(dplyr::all_of(group_cols)),
                    .data$isotype_class),
    n_cells = dplyr::n(), .groups = "drop"))

  out <- dplyr::ungroup(dplyr::mutate(
    dplyr::group_by(out, dplyr::across(dplyr::all_of(group_cols))),
    fraction = .data$n_cells / sum(.data$n_cells)))

  attr(out, "unswitched_combined") <- combine_unswitched
  out
}

#' Class-switched fraction per group
switched_fraction <- function(db, group_cols = "sample_id") {
  need("dplyr", "switched_fraction()")
  heavy <- db[db$locus == "IGH", , drop = FALSE]
  heavy$isotype <- sub("\\*.*$", "", sub(",.*$", "",
                                         as.character(heavy$c_call)))
  # Cells without an isotype call are excluded from the denominator rather
  # than counted as unswitched: the isotype is unknown, not IGHM.
  heavy <- heavy[!is.na(heavy$isotype) & nzchar(trimws(heavy$isotype)), ,
                 drop = FALSE]

  dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(heavy, dplyr::across(dplyr::all_of(group_cols))),
    n_cells_with_isotype = dplyr::n(),
    n_switched = sum(is_switched(.data$isotype)),
    pct_switched = 100 * mean(is_switched(.data$isotype)),
    .groups = "drop"))
}

#' Flag clones large enough to drive a subject's summary statistics
#'
#' A single lineage can account for most of a compartment. Any
#' repertoire-level quantity computed on that compartment then describes the
#' clone, not the repertoire: gene usage, CDR3 length and mutation load all
#' follow it.
#'
#' This does not decide anything. It produces the profile needed to decide,
#' and `verdict` records which criteria a clone meets. Exclusion is a
#' judgement to be made and stated, and any statistic affected by it should
#' be reported both ways in a sensitivity analysis.
#'
#' Two criteria separate real expansion from a technical artefact:
#'   - UMI support at the level of the rest of that subject's cells
#'     (a clone built from low-UMI contigs is suspect);
#'   - a paired light chain across its cells
#'     (a heavy-chain-only expansion is usually an assembly artefact).
#' Intraclonal diversity is NOT a criterion: a genuine expansion over a few
#' days produces identical daughter cells, so low diversity is expected.
#'
#' @param db Resolved cells with clone assignments.
#' @param clone_col,subject_col,umi_col Column names.
#' @param min_fraction Report clones above this fraction of their subject's
#'   cells.
#' @return Data frame, one row per flagged clone.
profile_dominant_clones <- function(db,
                                    clone_col = "clone_id",
                                    subject_col = "subject_id",
                                    umi_col = "umi_count",
                                    min_fraction = 0.01) {
  need("dplyr", "profile_dominant_clones()")
  heavy <- db[db$locus == "IGH", , drop = FALSE]
  heavy$.umi <- as.numeric(heavy[[umi_col]])

  per_subject <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(heavy, dplyr::across(dplyr::all_of(subject_col))),
    subject_cells = dplyr::n(),
    subject_median_umi = stats::median(.data$.umi, na.rm = TRUE),
    .groups = "drop"))

  light_by_cell <- if (any(db$locus %in% c("IGK", "IGL"))) {
    unique(db$cell_id[db$locus %in% c("IGK", "IGL")])
  } else character(0)

  clones <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(heavy, dplyr::across(dplyr::all_of(c(subject_col, clone_col)))),
    n_cells = dplyr::n(),
    median_umi = stats::median(.data$.umi, na.rm = TRUE),
    n_distinct_junction = dplyr::n_distinct(.data$junction),
    n_with_light_chain = sum(.data$cell_id %in% light_by_cell),
    v_gene = dplyr::first(sub("\\*.*$", "", sub(",.*$", "",
                                                as.character(.data$v_call)))),
    isotype = dplyr::first(sub("\\*.*$", "", sub(",.*$", "",
                                                 as.character(.data$c_call)))),
    cdr3_length = stats::median(nchar(as.character(.data$junction_aa)),
                                na.rm = TRUE),
    .groups = "drop"))

  clones <- dplyr::left_join(clones, per_subject, by = subject_col)
  clones$fraction_of_subject <- clones$n_cells / clones$subject_cells
  clones$umi_ratio <- clones$median_umi / clones$subject_median_umi
  clones$pct_with_light_chain <- 100 * clones$n_with_light_chain / clones$n_cells

  flagged <- clones[clones$fraction_of_subject >= min_fraction, , drop = FALSE]
  if (!nrow(flagged)) return(flagged)

  # Two criteria, reported rather than applied.
  flagged$umi_at_baseline   <- flagged$umi_ratio >= 0.5
  flagged$has_paired_light  <- flagged$pct_with_light_chain >= 50

  flagged$verdict <- ifelse(
    flagged$umi_at_baseline & flagged$has_paired_light, "supported",
    ifelse(!flagged$umi_at_baseline & !flagged$has_paired_light,
           "review: low UMI support and no paired light chain",
           ifelse(!flagged$umi_at_baseline,
                  "review: low UMI support",
                  "review: no paired light chain")))

  flagged[order(-flagged$fraction_of_subject), ]
}

#' Expansion threshold appropriate to the material
#'
#' Resting and stimulated material are not comparable on one threshold. A
#' multi-day stimulation culture produces identical daughter cells, so a
#' cut-off that marks genuine expansion in resting material is met trivially
#' after stimulation.
#'
#' @param stimulated Logical vector, or a single logical.
#' @param min_resting,min_stimulated Thresholds in cells per clone.
expansion_threshold <- function(stimulated, min_resting = 2L,
                                min_stimulated = 4L) {
  ifelse(as.logical(stimulated), min_stimulated, min_resting)
}
