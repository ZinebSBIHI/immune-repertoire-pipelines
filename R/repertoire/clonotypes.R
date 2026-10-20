# ---------------------------------------------------------------------------
# Clonotype definition, collapsing and gene usage
#
# Input throughout is a long data frame in AIRR Rearrangement format
# (Vander Heiden et al. 2018, PMID 30323809), one row per rearrangement,
# with at least:
#
#   sample_id, subject_id, v_call, j_call, junction, junction_aa, locus
#   and an abundance column (UMI or read count)
#
# ABUNDANCE COLUMN. UMI-based and read-based protocols are not
# interchangeable. `duplicate_count` in an AIRR file means molecules in some
# assemblers' output and reads in others'; using reads where UMIs exist
# inflates abundance by the PCR duplication factor, which differs per
# sample. The column is therefore named by the caller, never guessed.
# ---------------------------------------------------------------------------

#' Build the clonotype key
#'
#' @param db AIRR data frame.
#' @param definition One of "aa", "nt", "gene", "strict".
#' @param strip_allele Drop allele suffixes from gene calls (IGHV1-69*01 ->
#'   IGHV1-69). Allele resolution is not reliable at typical repertoire
#'   sequencing depth and allele-level keys split genuine clonotypes.
#' @return Character vector of clonotype keys.
clonotype_key <- function(db, definition = "strict", strip_allele = TRUE) {
  required <- switch(
    definition,
    "aa"     = "junction_aa",
    "nt"     = "junction",
    "gene"   = c("v_call", "j_call"),
    "strict" = c("v_call", "j_call", "junction"),
    stop("Unknown clonotype definition: ", definition, call. = FALSE)
  )
  missing_cols <- setdiff(required, names(db))
  if (length(missing_cols)) {
    stop("Clonotype definition '", definition, "' needs column(s): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  gene <- function(x) {
    x <- as.character(x)
    # A rearrangement may carry several ambiguous calls separated by commas;
    # the first is the primary assignment.
    x <- sub(",.*$", "", x)
    if (strip_allele) x <- sub("\\*.*$", "", x)
    x
  }

  switch(
    definition,
    "aa"     = as.character(db$junction_aa),
    "nt"     = as.character(db$junction),
    "gene"   = paste(gene(db$v_call), gene(db$j_call), sep = "|"),
    "strict" = paste(gene(db$v_call), gene(db$j_call),
                     as.character(db$junction), sep = "|")
  )
}

#' Collapse an AIRR table to clonotype counts per sample
#'
#' @param db AIRR data frame.
#' @param abundance_col Name of the abundance column (UMI or read count).
#' @param definition Clonotype definition.
#' @param min_abundance Clonotypes below this abundance are dropped.
#' @param productive_only Keep only productive rearrangements.
#' @return Data frame: sample_id, subject_id, clonotype, abundance, frequency.
collapse_clonotypes <- function(db, abundance_col,
                                definition = "strict",
                                min_abundance = 1,
                                productive_only = TRUE,
                                locus = NULL) {
  need("dplyr", "collapse_clonotypes()")
  if (!abundance_col %in% names(db)) {
    stop("Abundance column '", abundance_col, "' not found. Available: ",
         paste(names(db), collapse = ", "),
         "\nName the column explicitly: UMI counts and read counts are not ",
         "interchangeable.", call. = FALSE)
  }

  if (productive_only) {
    if ("productive" %in% names(db)) {
      prod <- db$productive
      if (is.character(prod)) prod <- toupper(prod) %in% c("T", "TRUE")
      db <- db[which(as.logical(prod)), , drop = FALSE]
    } else {
      warning("No 'productive' column; productive filtering skipped. ",
              "Non-productive rearrangements inflate diversity.",
              call. = FALSE)
    }
  }

  if (!is.null(locus)) {
    if ("locus" %in% names(db)) {
      db <- db[db$locus %in% locus, , drop = FALSE]
    } else {
      # Silently returning every locus would analyse a mixed TRA/TRB or
      # IGH/IGK file as a single repertoire, which inflates richness and
      # makes every diversity index meaningless.
      warning("locus = '", paste(locus, collapse = "/"), "' was requested ",
              "but the input has no 'locus' column, so NO locus filtering ",
              "was applied. If the file holds more than one locus, every ",
              "count and diversity index below is computed over all of ",
              "them together.", call. = FALSE)
    }
  }

  if (!nrow(db)) {
    stop("No rearrangements left after filtering. Check the productive flag, ",
         "the locus and the abundance column.", call. = FALSE)
  }

  db$clonotype <- clonotype_key(db, definition)
  db$.abund <- as.numeric(db[[abundance_col]])

  if (!"subject_id" %in% names(db)) db$subject_id <- NA_character_

  out <- dplyr::summarise(
    dplyr::group_by(db, sample_id = .data$sample_id,
                    subject_id = .data$subject_id,
                    clonotype = .data$clonotype),
    abundance = sum(.data$.abund, na.rm = TRUE),
    n_rearrangements = dplyr::n(),
    .groups = "drop"
  )

  out <- out[out$abundance >= min_abundance, , drop = FALSE]

  # Frequency is computed after the abundance filter, so frequencies sum to
  # one within the filtered repertoire. Computing it before would make
  # frequencies depend on material that no longer contributes to any count.
  out <- dplyr::mutate(
    dplyr::group_by(out, .data$sample_id),
    frequency = .data$abundance / sum(.data$abundance)
  )
  dplyr::ungroup(out)
}

#' V or J gene usage per sample
#'
#' @param db AIRR data frame.
#' @param gene_col "v_call" or "j_call".
#' @param level "clonotype" (one count per clonotype) or "molecule"
#'   (abundance-weighted).
#'
#' At molecule level a single large clone dominates its subject's usage
#' profile, so an apparent gene bias can be one expanded lineage. Computing
#' usage at clonotype level separates repertoire-wide bias from clonal
#' expansion; where they disagree, the clone is the explanation. Report both.
gene_usage <- function(db, gene_col = "v_call", level = c("clonotype", "molecule"),
                       abundance_col = NULL, definition = "strict") {
  level <- match.arg(level)
  need(c("dplyr", "tidyr"), "gene_usage()")

  db$gene <- sub("\\*.*$", "", sub(",.*$", "", as.character(db[[gene_col]])))
  db <- db[!is.na(db$gene) & nzchar(db$gene), , drop = FALSE]

  if (level == "clonotype") {
    db$clonotype <- clonotype_key(db, definition)
    u <- unique(db[, c("sample_id", "clonotype", "gene")])
    counts <- dplyr::summarise(
      dplyr::group_by(u, .data$sample_id, .data$gene),
      n = dplyr::n(), .groups = "drop")
  } else {
    if (is.null(abundance_col)) {
      stop("level = 'molecule' needs abundance_col.", call. = FALSE)
    }
    db$.abund <- as.numeric(db[[abundance_col]])
    counts <- dplyr::summarise(
      dplyr::group_by(db, .data$sample_id, .data$gene),
      n = sum(.data$.abund, na.rm = TRUE), .groups = "drop")
  }

  # Complete the grid over every (sample, gene) pair seen anywhere, filling
  # absences with zero. Without this, a gene missing from a sample
  # contributes no row, so a mean usage taken across samples averages only
  # the samples where the gene was observed: a gene at 30% in one sample of
  # three outranks a gene at 1% in all three.
  counts <- tidyr::complete(
    counts,
    sample_id = unique(counts$sample_id),
    gene = unique(counts$gene),
    fill = list(n = 0))

  dplyr::ungroup(dplyr::mutate(
    dplyr::group_by(counts, .data$sample_id),
    frequency = .data$n / sum(.data$n),
    level = level))
}

#' CDR3 length distribution weighted by abundance
cdr3_length_spectrum <- function(clono, junction_aa) {
  need("dplyr", "cdr3_length_spectrum()")
  clono$cdr3_length <- nchar(junction_aa)
  dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(clono, .data$sample_id, .data$cdr3_length),
    n_clonotypes = dplyr::n(),
    abundance = sum(.data$abundance),
    .groups = "drop"))
}
