# ---------------------------------------------------------------------------
# B-cell receptor: somatic hypermutation profiles and selection
#
# Mutation is described at three levels, in the order the original analyses
# used them:
#
#   load      mutation frequency of each heavy chain (one per cell in
#             single-cell data), by group, by isotype, per clone, and for
#             class-switched against unswitched cells;
#   R/S       replacement against silent mutations, CDR against FWR, summed
#             over each group. Descriptive: its expected value without
#             selection depends on codon composition and mutation hotspots;
#   BASELINe  selection strength (Sigma) estimated per clone within each
#             group, then combined per group. This is the test.
#
# Counts are always taken on heavy chains, against the germline chosen by
# section_bcr_clones() (D-masked where available).
#
# The single-cell sections (load, R/S, BASELINe, switched/unswitched) take
# the cohort's grouping columns as arguments: which compartments, which
# conditions and which sets of subjects are compared is a property of the
# study, not of the method. The bulk sections (isotype, V gene by CDR3
# length, largest clones) follow the bulk analyses step by step.
#
# References
#   Gupta NT, Vander Heiden JA, Uduman M, Gadala-Maria D, Yaari G,
#   Kleinstein SH. Change-O: a toolkit for analyzing large-scale B cell
#   immunoglobulin repertoire sequencing data. Bioinformatics.
#   2015;31(20):3356-3358. doi:10.1093/bioinformatics/btv359
#   Hershberg U, Uduman M, Shlomchik MJ, Kleinstein SH. Improved methods for
#   detecting selection by mutation analysis of Ig V region sequences.
#   Int Immunol. 2008;20(5):683-694. doi:10.1093/intimm/dxn026
#   Uduman M, Yaari G, Hershberg U, Stern JA, Shlomchik MJ, Kleinstein SH.
#   Detecting selection in immunoglobulin sequences. Nucleic Acids Res.
#   2011;39(Web Server issue):W499-W504. doi:10.1093/nar/gkr413
#   Yaari G, Uduman M, Kleinstein SH. Quantifying selection in
#   high-throughput immunoglobulin sequencing data sets. Nucleic Acids Res.
#   2012;40(17):e134. doi:10.1093/nar/gks457
# ---------------------------------------------------------------------------

#' Isotypes counted as class-switched in the switched/unswitched comparison
#'
#' The explicit list of the original analysis. Under this rule any other
#' call counts as unswitched: IGHM and IGHD, but also a call without a
#' subclass (IGHG, IGHA). section_shm_switched() reports such calls so that
#' the effect of the rule is visible.
#' @noRd
SWITCHED_ISOTYPES <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4",
                       "IGHA1", "IGHA2", "IGHE")

#' Isotype order of the bulk figures, from unswitched to IgE
#' @noRd
ISOTYPE_ORDER_IGH <- c("IGHD", "IGHM", "IGHG1", "IGHG2", "IGHG3", "IGHG4",
                       "IGHA1", "IGHA2", "IGHE")

#' ColorBrewer "Blues", for more groups than the brewer scale holds
#' @noRd
BLUES_9 <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6",
             "#4292C6", "#2171B5", "#08519C", "#08306B")

# --- settings ----------------------------------------------------------------

#' Defaults of the advanced BCR sections, overridable under `cfg$bcr`
#'
#' The values are those of the original analyses. They are read here rather
#' than in validate_config(), so that a configuration without these keys
#' behaves exactly as before.
#' @noRd
bcr_advanced_settings <- function(cfg) {
  b   <- cfg$bcr %||% list()
  sel <- b$selection %||% list()
  lin <- b$lineage %||% list()
  iso <- b$mutation_by_isotype %||% list()
  ov  <- b$overview %||% list()
  list(
    nproc                  = as.integer(b$nproc %||% 1L),
    rs_region              = as.character(sel$rs_region %||% "IMGT_V"),
    baseline_region        = as.character(sel$baseline_region %||%
                                            "IMGT_V_BY_REGIONS"),
    baseline_test          = as.character(sel$baseline_test_statistic %||%
                                            "focused"),
    collapse_method        = as.character(sel$collapse_method %||% "mostCommon"),
    baseline_min_sequences = as.integer(sel$baseline_min_sequences %||% 5L),
    baseline_min_clones    = as.integer(sel$baseline_min_clones %||% 2L),
    switched_isotypes      = as.character(sel$switched_isotypes %||%
                                            SWITCHED_ISOTYPES),
    exclude_stimulated     = isTRUE(sel$exclude_stimulated %||% TRUE),
    lineage_min_cells      = as.integer(lin$min_cells %||% 3L),
    lineage_minseq_sc      = as.integer(lin$minseq_single_cell %||% 3L),
    lineage_minseq_bulk    = as.integer(lin$minseq_bulk %||% 2L),
    lineage_build          = as.character(lin$build %||% "pml"),
    lineage_per_page       = as.integer(lin$trees_per_page %||% 6L),
    lineage_ncol           = as.integer(lin$page_columns %||% 3L),
    lineage_tip_size       = as.numeric(lin$tip_size %||% 2.5),
    lineage_max_clones     = if (is.null(lin$max_clones)) NULL else
                               as.integer(lin$max_clones),
    isotype_region         = as.character(iso$region %||% "IMGT_VDJ_BY_REGIONS"),
    isotype_levels         = as.character(iso$isotype_levels %||%
                                            ISOTYPE_ORDER_IGH),
    overview_exclude       = as.character(ov$exclude_isotypes %||% "IGHE"),
    overview_max_seq       = as.numeric(ov$max_seq_count %||% 250),
    top_n                  = as.integer(ov$top_n_clones %||% 100L)
  )
}

# --- small helpers -----------------------------------------------------------

#' Resolve a shazam region definition given by name
#' @noRd
shazam_region <- function(region) {
  if (!is.character(region)) return(region)
  if (length(region) != 1L) {
    stop("A region definition is named by a single string.", call. = FALSE)
  }
  # Region definitions are lazy-loaded *datasets* of shazam:
  # getExportedValue() -- what `::` uses -- finds them in the lazy-data
  # table, get(envir = asNamespace()) does not.
  out <- tryCatch(getExportedValue("shazam", region), error = function(e) NULL)
  if (is.null(out)) {
    stop("Unknown shazam region definition '", region, "'. Use IMGT_V, ",
         "IMGT_VDJ, IMGT_V_BY_REGIONS or IMGT_VDJ_BY_REGIONS.", call. = FALSE)
  }
  out
}

#' Column suffix naming the region of a total, e.g. "imgt_v"
#' @noRd
region_column_label <- function(region) {
  if (is.character(region)) tolower(region) else "region"
}

#' The germline column to count mutations against
#'
#' The one section_bcr_clones() recorded, else the D-masked germline when the
#' column exists, else the plain germline: the rule of the original analyses.
#' @noRd
selection_germline_col <- function(db, germline_col = NULL) {
  if (!is.null(germline_col)) {
    if (!germline_col %in% names(db)) {
      stop("Germline column '", germline_col, "' not found.", call. = FALSE)
    }
    return(germline_col)
  }
  recorded <- attr(db, "germline_col")
  if (!is.null(recorded) && recorded %in% names(db)) return(recorded)
  for (cand in c("germline_alignment_d_mask", "germline_alignment")) {
    if (cand %in% names(db)) return(cand)
  }
  stop("No germline column: expected germline_alignment_d_mask or ",
       "germline_alignment (section_bcr_clones() records the one it used).",
       call. = FALSE)
}

#' Heavy chains carrying both an observed and a germline alignment
#' @noRd
selection_heavy <- function(db, germline_col,
                            sequence_col = "sequence_alignment") {
  for (cl in c(sequence_col, germline_col)) {
    if (!cl %in% names(db)) {
      stop("Column '", cl, "' not found.", call. = FALSE)
    }
  }
  sq <- as.character(db[[sequence_col]])
  gl <- as.character(db[[germline_col]])
  keep <- !is.na(sq) & nzchar(sq) & !is.na(gl) & nzchar(gl)
  if ("locus" %in% names(db)) keep <- keep & db$locus %in% "IGH"
  db[which(keep), , drop = FALSE]
}

#' Heavy-chain rows (all of them when there is no locus column)
#' @noRd
igh_rows <- function(db) {
  if ("locus" %in% names(db)) db[which(db$locus == "IGH"), , drop = FALSE] else db
}

#' Subject-qualified clone identifier
#'
#' Clones are inferred within subject, so a clone number is an identity only
#' together with its subject. Qualifying it keeps two subjects' clone 12
#' apart wherever clones are grouped across subjects.
#' @noRd
subject_clone_key <- function(db, clone_col = "clone_id",
                              subject_col = "subject_id") {
  if (!clone_col %in% names(db)) {
    stop("Clone column '", clone_col, "' not found.", call. = FALSE)
  }
  clone <- as.character(db[[clone_col]])
  if (!is.null(subject_col) && subject_col %in% names(db)) {
    paste(db[[subject_col]], clone, sep = "_")
  } else {
    clone
  }
}

#' Silence shazam's notice that recomputed columns replace existing ones
#'
#' The replacement is intended wherever this is used; the warning would
#' otherwise be repeated once per pass.
#' @noRd
muffle_overwrite <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("exist and will be overwritten", conditionMessage(w),
              fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
}

#' Add allele-stripped v_gene and j_gene where absent
#' @noRd
with_gene_cols <- function(db) {
  if (!"v_gene" %in% names(db)) {
    db$v_gene <- if ("v_call" %in% names(db)) strip_allele(db$v_call) else
      NA_character_
  }
  if (!"j_gene" %in% names(db)) {
    db$j_gene <- if ("j_call" %in% names(db)) strip_allele(db$j_call) else
      NA_character_
  }
  db
}

#' Add the isotype column, from c_call, where absent
#' @noRd
with_isotype_col <- function(db, isotype_col) {
  if (isotype_col %in% names(db)) return(db)
  if (!"c_call" %in% names(db)) {
    stop("No isotype column '", isotype_col, "' and no c_call to derive it ",
         "from.", call. = FALSE)
  }
  db[[isotype_col]] <- strip_allele(db$c_call)
  db
}

#' Add junction_aa, translated from the junction, where absent
#' @noRd
with_junction_aa <- function(db) {
  if ("junction_aa" %in% names(db)) return(db)
  if (!"junction" %in% names(db)) {
    stop("Neither junction_aa nor junction is available.", call. = FALSE)
  }
  need("alakazam", "Translating junctions")
  db$junction_aa <- alakazam::translateDNA(as.character(db$junction))
  db
}

#' A factor following a configured order; values outside it become NA
#' @noRd
ordered_groups <- function(x, levels = NULL) {
  if (is.null(levels)) {
    levels <- if (is.factor(x)) levels(droplevels(x)) else
      sort(unique(as.character(x[!is.na(x)])))
  }
  factor(as.character(x), levels = as.character(levels))
}

#' Stop with the list of missing columns
#' @noRd
require_columns <- function(db, cols, what) {
  miss <- setdiff(cols, names(db))
  if (length(miss)) {
    stop(what, " needs column(s) not found: ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  invisible(TRUE)
}

# --- mutation profiles -------------------------------------------------------

#' Mutation counts and frequencies by region, replacement and silent
#'
#' Four passes of `shazam::observedMutations()`, as in the original
#' analyses: counts by sub-region split into replacement (R) and silent (S);
#' the total count; frequencies by sub-region; the total frequency. The
#' total frequency is computed rather than derived, because the sub-region
#' frequencies have different denominators.
#'
#' Region definitions: `IMGT_V` covers FWR1-3 and CDR1-2 and excludes CDR3
#' and FWR4, giving `cdr` and `fwr` sub-regions; `IMGT_VDJ` covers the whole
#' V(D)J with the same two sub-regions; `IMGT_V_BY_REGIONS` and
#' `IMGT_VDJ_BY_REGIONS` split them into CDR1-3 and FWR1-4. The totals of
#' `IMGT_VDJ` and `IMGT_VDJ_BY_REGIONS` are identical.
#'
#' @param db AIRR data frame.
#' @param region shazam region definition, by name.
#' @param germline_col Germline column. By default the one recorded by
#'   `section_bcr_clones()`, else `germline_alignment_d_mask` when present,
#'   else `germline_alignment`.
#' @param sequence_col Observed alignment column.
#' @param counts_only Stop after the first pass (counts by sub-region), which
#'   is all the R/S ratio needs.
#' @param nproc Cores passed to shazam.
#' @return The heavy chains with complete alignments, with
#'   `mu_count_<subregion>_<r|s>`, `mu_freq_<subregion>_<r|s>` and the totals
#'   `mu_count_<region>` and `mu_freq_<region>` (region in lower case, e.g.
#'   `mu_freq_imgt_v`). Plain `mu_count` / `mu_freq` columns already present
#'   are left as they were.
#' @export
observed_mutation_profile <- function(db, region = "IMGT_V", germline_col = NULL,
                                      sequence_col = "sequence_alignment",
                                      counts_only = FALSE, nproc = 1L) {
  need("shazam", "observed_mutation_profile()")
  germline_col <- selection_germline_col(db, germline_col)
  region_def <- shazam_region(region)
  label <- region_column_label(region)

  heavy <- selection_heavy(db, germline_col, sequence_col)
  if (!nrow(heavy)) {
    stop("No heavy chain carries both a sequence and a germline alignment.",
         call. = FALSE)
  }
  # V(D)J definitions place the CDR3/FWR4 boundaries from the junction
  # length, so rows without one are removed first, as in the original bulk
  # analysis.
  if (region_def@name %in% c("IMGT_VDJ", "IMGT_VDJ_BY_REGIONS")) {
    if (!"junction_length" %in% names(heavy)) {
      if (!"junction" %in% names(heavy)) {
        stop(region_def@name, " needs junction_length (or junction).",
             call. = FALSE)
      }
      heavy$junction_length <- nchar(as.character(heavy$junction))
    }
    heavy <- heavy[which(!is.na(heavy$junction_length)), , drop = FALSE]
    if (!nrow(heavy)) stop("No heavy chain has a junction length.", call. = FALSE)
  }

  sub_labels  <- region_def@labels
  total_count <- paste0("mu_count_", label)
  total_freq  <- paste0("mu_freq_", label)
  recomputed  <- c(paste0("mu_count_", sub_labels),
                   paste0("mu_freq_", sub_labels), total_count, total_freq)
  heavy <- heavy[, setdiff(names(heavy), recomputed), drop = FALSE]

  # The combined passes write to mu_count / mu_freq. section_bcr_clones()
  # already uses mu_freq for its primary region, so those columns are put
  # aside and restored rather than silently replaced.
  kept  <- intersect(c("mu_count", "mu_freq"), names(heavy))
  stash <- heavy[, kept, drop = FALSE]
  heavy <- heavy[, setdiff(names(heavy), kept), drop = FALSE]
  heavy$.profile_row <- seq_len(nrow(heavy))

  pass <- function(x, frequency, combine) {
    shazam::observedMutations(x, sequenceColumn = sequence_col,
                              germlineColumn = germline_col,
                              regionDefinition = region_def,
                              frequency = frequency, combine = combine,
                              nproc = nproc)
  }
  heavy <- pass(heavy, frequency = FALSE, combine = FALSE)
  if (!counts_only) {
    heavy <- pass(heavy, frequency = FALSE, combine = TRUE)
    names(heavy)[names(heavy) == "mu_count"] <- total_count
    heavy <- pass(heavy, frequency = TRUE, combine = FALSE)
    heavy <- pass(heavy, frequency = TRUE, combine = TRUE)
    names(heavy)[names(heavy) == "mu_freq"] <- total_freq
  }

  heavy <- heavy[order(heavy$.profile_row), , drop = FALSE]
  for (nm in kept) heavy[[nm]] <- stash[[nm]][heavy$.profile_row]
  heavy$.profile_row <- NULL
  attr(heavy, "germline_col") <- germline_col
  attr(heavy, "total_mutation_col") <- if (counts_only) NULL else total_freq
  heavy
}

#' Replacement-to-silent ratio in CDR and FWR, per group
#'
#' Mutations are summed over every heavy chain of a group before the ratio
#' is taken, and the silent count is floored at 1, so a group without silent
#' mutations yields its replacement count rather than an infinite ratio.
#'
#' @param db Output of `observed_mutation_profile()` with a CDR/FWR
#'   definition (`IMGT_V` or `IMGT_VDJ`), i.e. with `mu_count_cdr_r`,
#'   `mu_count_cdr_s`, `mu_count_fwr_r` and `mu_count_fwr_s`.
#' @param group_cols Grouping columns, e.g. subject and compartment.
#' @return One row per group: `n`, the four summed counts, `rs_cdr` and
#'   `rs_fwr`.
#' @export
replacement_silent_ratio <- function(db, group_cols = c("subject_id",
                                                        "compartment")) {
  counts <- c("mu_count_cdr_r", "mu_count_cdr_s",
              "mu_count_fwr_r", "mu_count_fwr_s")
  require_columns(db, c(counts, group_cols), "replacement_silent_ratio()")
  out <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(db, dplyr::across(dplyr::all_of(group_cols))),
    n     = dplyr::n(),
    cdr_r = sum(.data$mu_count_cdr_r, na.rm = TRUE),
    cdr_s = sum(.data$mu_count_cdr_s, na.rm = TRUE),
    fwr_r = sum(.data$mu_count_fwr_r, na.rm = TRUE),
    fwr_s = sum(.data$mu_count_fwr_s, na.rm = TRUE),
    .groups = "drop"))
  out$rs_cdr <- out$cdr_r / pmax(out$cdr_s, 1)
  out$rs_fwr <- out$fwr_r / pmax(out$fwr_s, 1)
  out
}

#' Selection strength per group with BASELINe
#'
#' Within each group separately: mutations are counted, each clone is
#' collapsed to its most common sequence, BASELINe Sigma is estimated per
#' clone with the focused test statistic, and the clone-level densities are
#' combined per group. A clone spanning several groups is collapsed
#' separately in each, so a group's estimate rests on that group's sequences
#' only. Sigma > 0 indicates positive selection, Sigma < 0 negative
#' (purifying) selection.
#'
#' @param db Heavy chains with clone assignments and alignments.
#' @param group_col Column defining the groups (e.g. compartment).
#' @param groups Groups to analyse, in order; all by default.
#' @param clone_col,subject_col Clone and subject columns; clones are
#'   qualified by subject.
#' @param germline_col See `observed_mutation_profile()`.
#' @param region Region definition.
#' @param test_statistic BASELINe test statistic.
#' @param targeting_model SHM targeting model giving the expected mutation
#'   pattern without selection; shazam's default, the human 5-mer model
#'   `HH_S5F`.
#' @param collapse_method Method of `shazam::collapseClones()`.
#' @param min_sequences,min_clones A group with fewer sequences or fewer
#'   clones is skipped, and the reason recorded.
#' @param nproc Cores passed to shazam.
#' @return A list: `baseline`, the grouped shazam Baseline object (NULL when
#'   nothing could be computed); `stats`, Sigma with its confidence interval
#'   per group and region; `groups`, sequences, clones and status per group.
#'   An `error` attribute carries the message if `calcBaseline()` failed.
#' @export
baseline_selection <- function(db, group_col, groups = NULL,
                               clone_col = "clone_id",
                               subject_col = "subject_id",
                               germline_col = NULL,
                               region = "IMGT_V_BY_REGIONS",
                               test_statistic = "focused",
                               targeting_model = "HH_S5F",
                               collapse_method = "mostCommon",
                               min_sequences = 5L, min_clones = 2L,
                               nproc = 1L) {
  need("shazam", "baseline_selection()")
  require_columns(db, c(group_col, clone_col), "baseline_selection()")
  germline_col <- selection_germline_col(db, germline_col)
  region_def <- shazam_region(region)
  # calcBaseline()'s default model is a lazy-loaded dataset that its default
  # argument only finds when shazam is attached, so it is passed explicitly.
  targeting <- if (is.character(targeting_model)) {
    getExportedValue("shazam", targeting_model)
  } else {
    targeting_model
  }

  heavy <- selection_heavy(db, germline_col)
  heavy <- heavy[!is.na(heavy[[clone_col]]), , drop = FALSE]
  heavy$.clone_key <- subject_clone_key(heavy, clone_col, subject_col)
  if (is.null(groups)) {
    groups <- if (is.factor(db[[group_col]])) levels(droplevels(db[[group_col]]))
              else sort(unique(stats::na.omit(as.character(heavy[[group_col]]))))
  }
  groups <- as.character(groups)
  heavy[[group_col]] <- as.character(heavy[[group_col]])
  heavy <- heavy[, setdiff(names(heavy),
                           paste0("mu_count_", region_def@labels)), drop = FALSE]

  status <- data.frame(group = groups, n_sequences = 0L, n_clones = 0L,
                       status = "", stringsAsFactors = FALSE)
  per <- list()
  for (i in seq_along(groups)) {
    sub <- heavy[which(heavy[[group_col]] == groups[i]), , drop = FALSE]
    status$n_sequences[i] <- nrow(sub)
    status$n_clones[i] <- length(unique(sub$.clone_key))
    if (nrow(sub) < min_sequences || status$n_clones[i] < min_clones) {
      status$status[i] <- sprintf("skipped: fewer than %d sequences or %d clones",
                                  min_sequences, min_clones)
      next
    }
    res <- tryCatch(muffle_overwrite({
      obs <- shazam::observedMutations(
        sub, sequenceColumn = "sequence_alignment",
        germlineColumn = germline_col, regionDefinition = region_def,
        frequency = FALSE, combine = FALSE, nproc = nproc)
      obs <- shazam::collapseClones(
        obs, cloneColumn = ".clone_key", sequenceColumn = "sequence_alignment",
        germlineColumn = germline_col, regionDefinition = region_def,
        method = collapse_method, nproc = nproc)
      obs[[group_col]] <- groups[i]
      obs
    }), error = function(e) e)
    if (inherits(res, "error")) {
      status$status[i] <- paste("failed:", conditionMessage(res))
      next
    }
    status$status[i] <- "computed"
    per[[length(per) + 1L]] <- res
  }
  names(status)[1] <- group_col

  out <- list(baseline = NULL, stats = NULL, groups = status)
  if (!length(per)) return(out)

  bl <- tryCatch(muffle_overwrite({
    b <- shazam::calcBaseline(
      dplyr::bind_rows(per), testStatistic = test_statistic,
      sequenceColumn = "clonal_sequence", germlineColumn = "clonal_germline",
      regionDefinition = region_def, targetingModel = targeting,
      nproc = nproc)
    shazam::groupBaseline(b, groupBy = group_col)
  }), error = function(e) e)
  if (inherits(bl, "error")) {
    attr(out, "error") <- conditionMessage(bl)
    return(out)
  }
  out$baseline <- bl
  out$stats <- as.data.frame(bl@stats)
  out
}

#' Mutation load of class-switched and unswitched cells
#'
#' Tells apart two explanations of a difference in mean mutation load
#' between compartments, or between groups of subjects: different maturation
#' of the class-switched cells, or a different share of unswitched cells,
#' which carry few mutations.
#'
#' The rule of the original analysis: cells without an isotype call are
#' excluded, stimulated material is excluded, and a cell is class-switched
#' when its isotype is in `switched`, unswitched otherwise.
#'
#' @param db Heavy chains, one per cell.
#' @param shm_col Mutation frequency column.
#' @param isotype_col Isotype column (alleles stripped); derived from
#'   `c_call` when absent.
#' @param compartment_col,subject_col Compartment and subject columns.
#' @param subject_group_col Optional column grouping subjects (e.g. two
#'   sampling campaigns); a single group "all" otherwise.
#' @param switched Isotypes counted as class-switched. Default IGHG1-4,
#'   IGHA1-2 and IGHE.
#' @param exclude_stimulated Drop rows whose logical `stimulated` column is
#'   TRUE.
#' @param scale Multiplier applied to `shm_col`; 100 reports percentages, as
#'   in the original tables.
#' @return A list: `by_subject` (cells and mean load per group of subjects,
#'   subject, compartment and isotype class); `summary` (number of subjects,
#'   mean, minimum and maximum of the subject means); `unswitched_share`
#'   (percentage of unswitched cells and overall load per subject and
#'   compartment); `outside_lists` (calls in neither the switched list nor
#'   IGHM/IGHD, counted as unswitched); `excluded` (cells dropped for a
#'   missing isotype or as stimulated).
#' @export
shm_by_switch_status <- function(db, shm_col = "mu_freq", isotype_col = "c_gene",
                                 compartment_col = "compartment",
                                 subject_col = "subject_id",
                                 subject_group_col = NULL,
                                 switched = SWITCHED_ISOTYPES,
                                 exclude_stimulated = TRUE, scale = 100) {
  db <- with_isotype_col(db, isotype_col)
  require_columns(db, c(shm_col, compartment_col, subject_col,
                        subject_group_col), "shm_by_switch_status()")
  d <- igh_rows(db)

  iso <- as.character(d[[isotype_col]])
  has_iso <- !is.na(iso) & nzchar(trimws(iso)) & iso != "NA"
  n_no_isotype <- sum(!has_iso)
  d <- d[which(has_iso), , drop = FALSE]
  n_stimulated <- 0L
  if (exclude_stimulated && "stimulated" %in% names(d)) {
    stim <- as.logical(d$stimulated) %in% TRUE
    n_stimulated <- sum(stim)
    d <- d[which(!stim), , drop = FALSE]
  }
  if (!nrow(d)) {
    stop("No heavy chain with an isotype call remains.", call. = FALSE)
  }

  iso <- as.character(d[[isotype_col]])
  d$subject_group <- if (is.null(subject_group_col)) "all" else
    as.character(d[[subject_group_col]])
  d$isotype_class <- ifelse(iso %in% switched, "class-switched", "unswitched")
  d$.shm <- scale * as.numeric(d[[shm_col]])

  by_subject <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(d, dplyr::across(dplyr::all_of(
      c("subject_group", subject_col, compartment_col, "isotype_class")))),
    cells = dplyr::n(),
    shm = mean(.data$.shm, na.rm = TRUE),
    .groups = "drop"))

  summary <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(by_subject, dplyr::across(dplyr::all_of(
      c("subject_group", compartment_col, "isotype_class")))),
    n_subjects = dplyr::n(),
    shm_mean = mean(.data$shm),
    shm_min = min(.data$shm),
    shm_max = max(.data$shm),
    .groups = "drop"))

  share <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(d, dplyr::across(dplyr::all_of(
      c("subject_group", subject_col, compartment_col)))),
    n_cells = dplyr::n(),
    pct_unswitched = 100 * mean(.data$isotype_class == "unswitched"),
    shm_all = mean(.data$.shm, na.rm = TRUE),
    .groups = "drop"))

  outside <- iso[!iso %in% c(switched, UNSWITCHED_ISOTYPES)]
  outside_tbl <- as.data.frame(table(isotype = outside),
                               stringsAsFactors = FALSE)
  names(outside_tbl)[names(outside_tbl) == "Freq"] <- "n_cells"

  list(by_subject = by_subject, summary = summary, unswitched_share = share,
       outside_lists = outside_tbl,
       excluded = c(no_isotype = n_no_isotype, stimulated = n_stimulated))
}

#' Mean mutation frequency per clone, group and isotype
#'
#' @param db Heavy chains with clone assignments.
#' @param shm_col Mutation frequency column.
#' @param clone_col,subject_col Clone and subject columns.
#' @param group_col Sample or condition column.
#' @param isotype_col Isotype column (derived from `c_call` when absent).
#' @return One row per subject, clone, group and isotype: `n_sequences` and
#'   `mean_mut_freq`.
#' @export
clone_mean_mutation <- function(db, shm_col = "mu_freq", clone_col = "clone_id",
                                subject_col = "subject_id",
                                group_col = "sample_id", isotype_col = "c_gene") {
  db <- with_isotype_col(db, isotype_col)
  require_columns(db, c(shm_col, clone_col, group_col),
                  "clone_mean_mutation()")
  heavy <- igh_rows(db)
  heavy <- heavy[!is.na(heavy[[clone_col]]), , drop = FALSE]
  keys <- c(intersect(subject_col, names(heavy)), clone_col, group_col,
            isotype_col)
  dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(heavy, dplyr::across(dplyr::all_of(keys))),
    n_sequences = dplyr::n(),
    mean_mut_freq = mean(.data[[shm_col]], na.rm = TRUE),
    .groups = "drop"))
}

#' The largest clones, with their mutation load and genes
#'
#' @param db Heavy chains with clone assignments.
#' @param n Number of clones returned.
#' @param shm_col Mutation frequency column.
#' @param clone_col,subject_col Clone and subject columns.
#' @param group_col Sample or condition column, listed per clone.
#' @return The `n` largest clones (ties in their original order):
#'   `n_sequences` (cells in single-cell data), mean and median mutation
#'   frequency, the V and J genes and the groups the clone was found in,
#'   each listed once.
#' @export
top_clones_summary <- function(db, n = 100L, shm_col = "mu_freq",
                               clone_col = "clone_id",
                               subject_col = "subject_id",
                               group_col = "sample_id") {
  require_columns(db, c(shm_col, clone_col, group_col), "top_clones_summary()")
  heavy <- with_gene_cols(igh_rows(db))
  heavy <- heavy[!is.na(heavy[[clone_col]]), , drop = FALSE]
  keys <- c(intersect(subject_col, names(heavy)), clone_col)
  out <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(heavy, dplyr::across(dplyr::all_of(keys))),
    n_sequences = dplyr::n(),
    mean_mu_freq = mean(.data[[shm_col]], na.rm = TRUE),
    median_mu_freq = stats::median(.data[[shm_col]], na.rm = TRUE),
    v_genes = paste(unique(.data$v_gene), collapse = "; "),
    j_genes = paste(unique(.data$j_gene), collapse = "; "),
    present_in = paste(unique(.data[[group_col]]), collapse = "; "),
    .groups = "drop"))
  out <- out[order(-out$n_sequences), , drop = FALSE]
  utils::head(out, n)
}

#' IGHV gene by CDR3 length overview, one point per clone and group
#'
#' Each point is a clone within a group, V gene, isotype and junction
#' length, carrying the median of a mutation measure (in percent) and its
#' number of sequences. The length is `nchar(junction_aa)`, the junction
#' (CDR3 with its two anchor residues), as in the original analyses.
#'
#' @param db Heavy chains with clone assignments.
#' @param value_col Mutation measure, e.g. `mu_freq_cdr3_r` (CDR3
#'   replacement frequency) or a total mutation frequency.
#' @param clone_col,subject_col Clone and subject columns.
#' @param group_col Sample or condition column.
#' @param isotype_col Isotype column (derived from `c_call` when absent).
#' @param clones Optional subject-qualified clone keys ("subject_clone") to
#'   restrict to, e.g. the largest clones.
#' @return One row per subject, clone, V gene, isotype, group and junction
#'   length: `median_mut` (percent) and `seq_count`.
#' @export
vgene_cdr3_overview <- function(db, value_col = "mu_freq_cdr3_r",
                                clone_col = "clone_id",
                                subject_col = "subject_id",
                                group_col = "sample_id", isotype_col = "c_gene",
                                clones = NULL) {
  db <- with_isotype_col(db, isotype_col)
  require_columns(db, c(value_col, clone_col, group_col),
                  "vgene_cdr3_overview()")
  heavy <- with_junction_aa(with_gene_cols(igh_rows(db)))
  iso <- as.character(heavy[[isotype_col]])
  vg <- as.character(heavy$v_gene)
  heavy <- heavy[which(!is.na(iso) & nzchar(iso) & !is.na(vg) & nzchar(vg)), ,
                 drop = FALSE]
  if (!is.null(clones)) {
    heavy <- heavy[subject_clone_key(heavy, clone_col, subject_col) %in% clones, ,
                   drop = FALSE]
  }
  heavy$junction_aa_length <- nchar(as.character(heavy$junction_aa))
  keys <- c(intersect(subject_col, names(heavy)), clone_col, "v_gene",
            isotype_col, group_col, "junction_aa_length")
  dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(heavy, dplyr::across(dplyr::all_of(keys))),
    median_mut = stats::median(.data[[value_col]], na.rm = TRUE) * 100,
    seq_count = dplyr::n(),
    .groups = "drop"))
}

# --- figures (internal) ------------------------------------------------------

#' Box plots of mutation load per group
#' @noRd
mutation_load_plot <- function(d, x_col, facet_rows = NULL, facet_cols = NULL,
                               fill_values = NULL, title = NULL,
                               y_label = "V mutation frequency (%)") {
  n_x <- nlevels(d[[x_col]])
  p <- ggplot2::ggplot(d, ggplot2::aes(.data[[x_col]], .data$.shm,
                                       fill = .data[[x_col]])) +
    ggplot2::geom_boxplot(outlier.size = 0.3)
  if (length(facet_rows) && length(facet_cols)) {
    p <- p + ggplot2::facet_grid(rows = ggplot2::vars(!!!rlang::syms(facet_rows)),
                                 cols = ggplot2::vars(!!!rlang::syms(facet_cols)))
  } else if (length(facet_cols)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(!!!rlang::syms(facet_cols)),
                                 nrow = 1)
  } else if (length(facet_rows)) {
    p <- p + ggplot2::facet_grid(rows = ggplot2::vars(!!!rlang::syms(facet_rows)))
  }
  fill <- if (!is.null(fill_values)) {
    ggplot2::scale_fill_manual(values = unlist(fill_values), guide = "none")
  } else if (n_x <= 9) {
    ggplot2::scale_fill_brewer(palette = "Blues", guide = "none")
  } else {
    ggplot2::scale_fill_manual(values = grDevices::colorRampPalette(BLUES_9)(n_x),
                               guide = "none")
  }
  p + fill +
    ggplot2::labs(title = title, x = NULL, y = y_label) +
    theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Box plots of a mutation measure per group, one panel per isotype
#' @noRd
mutation_by_isotype_plot <- function(d, value_col, group_col, isotype_col,
                                     title, y_label = "Mutation frequency") {
  ggplot2::ggplot(d, ggplot2::aes(.data[[group_col]], .data[[value_col]],
                                  fill = .data[[isotype_col]])) +
    ggplot2::geom_boxplot(outlier.alpha = 0) +
    ggplot2::labs(title = title, x = NULL, y = y_label) +
    ggplot2::scale_fill_viridis_d(option = "B") +
    ggplot2::facet_grid(cols = ggplot2::vars(.data[[isotype_col]]),
                        scales = "free_y") +
    theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5,
                                                       hjust = 1)) +
    ggplot2::guides(fill = "none")
}

#' Density of the per-clone mean mutation, isotype by group
#' @noRd
clone_mutation_density_plot <- function(cm, group_col, isotype_col) {
  ggplot2::ggplot(cm, ggplot2::aes(.data$mean_mut_freq,
                                   fill = .data[[isotype_col]])) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::facet_grid(rows = ggplot2::vars(.data[[isotype_col]]),
                        cols = ggplot2::vars(.data[[group_col]]),
                        scales = "free") +
    ggplot2::scale_fill_viridis_d(option = "B") +
    ggplot2::labs(x = "Mean mutation frequency per clone", y = "Density",
                  title = "Somatic hypermutation rate by isotype and condition",
                  fill = "Isotype") +
    theme_pub(base_size = 8) +
    ggplot2::guides(fill = "none")
}

#' The IGHV by CDR3 length overview figure, styled as in the original
#' @noRd
vgene_cdr3_overview_plot <- function(ov, group_col, isotype_col, title,
                                     colour_title, size_title,
                                     colour_values, colour_limits, size_breaks) {
  ggplot2::ggplot(ov, ggplot2::aes(.data$v_gene, .data$junction_aa_length,
                                   colour = .data$median_mut)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$seq_count)) +
    ggplot2::labs(title = title, x = "IGHV genes",
                  y = "Junction (CDR3) length (aa)",
                  colour = colour_title, size = size_title) +
    ggplot2::theme_classic() +
    ggplot2::facet_grid(rows = ggplot2::vars(.data[[group_col]]),
                        cols = ggplot2::vars(.data[[isotype_col]])) +
    ggplot2::scale_colour_gradientn(
      colours = c("green", "white", "magenta", "#6f0795"),
      values = scales::rescale(colour_values), limits = colour_limits) +
    ggplot2::scale_size(breaks = size_breaks, range = c(1, 15)) +
    ggplot2::theme(
      panel.background  = ggplot2::element_rect(fill = "black"),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key        = ggplot2::element_rect(colour = "white", fill = "white"),
      axis.text.x       = ggplot2::element_blank(),
      axis.ticks.x      = ggplot2::element_blank(),
      axis.text.y       = ggplot2::element_text(colour = "black", size = 15,
                                                face = "bold"),
      axis.title        = ggplot2::element_text(size = 15),
      legend.text       = ggplot2::element_text(size = 15),
      strip.text.x      = ggplot2::element_text(face = "bold", size = 15),
      strip.text.y      = ggplot2::element_text(face = "bold", size = 15))
}

#' Aggregate, tabulate and draw one IGHV by CDR3 length overview
#'
#' Shared by section_vgene_cdr3_overview() and section_top_clones(). Points
#' of an excluded isotype (IgE in the original) or with more sequences than
#' `max_seq_count` are tabulated but not drawn, as in the original figures.
#' @noRd
write_vgene_cdr3_overview <- function(db, value_col, clone_col, subject_col,
                                      group_col, group_levels, isotype_col,
                                      clones, title, colour_title, unit,
                                      colour_values, colour_limits, size_breaks,
                                      name, run_dir, cfg, report) {
  s <- bcr_advanced_settings(cfg)
  ov <- vgene_cdr3_overview(db, value_col = value_col, clone_col = clone_col,
                            subject_col = subject_col, group_col = group_col,
                            isotype_col = isotype_col, clones = clones)
  save_table(ov, name, run_dir,
             sprintf("median %s (percent) per clone, %s, V gene, isotype and junction length",
                     value_col, group_col))

  shown_iso <- setdiff(s$isotype_levels, s$overview_exclude)
  iso <- as.character(ov[[isotype_col]])
  keep <- iso %in% shown_iso & ov$seq_count <= s$overview_max_seq
  d <- ov[which(keep), , drop = FALSE]
  d[[isotype_col]] <- factor(as.character(d[[isotype_col]]), levels = shown_iso)
  d[[group_col]] <- ordered_groups(d[[group_col]], group_levels)
  d <- d[!is.na(d[[group_col]]), , drop = FALSE]

  report$p(sprintf(paste0(
    "%d points (clone x %s x V gene x isotype x junction length). Not drawn: ",
    "%d with an isotype outside the drawn set (%s) and %d with more than %g ",
    "%s."), nrow(ov), group_col, sum(!iso %in% shown_iso),
    paste(shown_iso, collapse = ", "),
    sum(iso %in% shown_iso & ov$seq_count > s$overview_max_seq),
    s$overview_max_seq, unit))

  if (!nrow(d)) {
    report$p("Nothing left to draw.")
    return(invisible(ov))
  }
  p <- vgene_cdr3_overview_plot(
    d, group_col = group_col, isotype_col = isotype_col, title = title,
    colour_title = colour_title,
    size_title = sprintf("%s per clone", tools_title(unit)),
    colour_values = colour_values, colour_limits = colour_limits,
    size_breaks = size_breaks)
  save_figure(p, name, run_dir, cfg, width = 18, height = 14,
              legend = paste0(
                "Each point is a clone within a ", group_col, " (rows) and ",
                "isotype (columns), placed by IGHV gene (x, labels omitted) ",
                "and junction length (y); colour: ", colour_title,
                " (values beyond the scale are grey); size: ", unit,
                " per clone."))
  invisible(ov)
}

# --- sections ------------------------------------------------------------------

#' Somatic hypermutation load per cell, by group
#'
#' Box plots of the mutation frequency (in percent) of each heavy chain, one
#' box per group, faceted by subject or by any other variable. The load
#' describes the maturation state of the cells sampled; after a short
#' culture it reflects mutation already present, not mutation acquired in
#' culture.
#'
#' The original single-cell analysis draws two such views: ex vivo
#' compartments faceted by subject, and, within stimulated material,
#' proliferating against non-proliferating cells faceted by subject and
#' condition. Restrict `db` to the cells of each view and call the section
#' once per view, with a distinct `name`.
#'
#' @param db Heavy chains with `shm_col`, one per cell in single-cell data.
#' @param x_col Column on the x axis (e.g. compartment).
#' @param facet_rows,facet_cols Faceting columns: a grid with both, one row
#'   of panels with `facet_cols` alone.
#' @param x_levels Optional order of the groups; groups outside it are left
#'   out of the figure.
#' @param shm_col Mutation frequency column (a fraction, shown in percent).
#' @param fill_values Optional named colours, one per group; shades of blue
#'   otherwise.
#' @param title Figure and report title.
#' @param name File stem, distinct for each call.
#' @param run_dir Run directory from `init_run()`.
#' @param cfg Configuration list from `load_config()`.
#' @param report Report writer from `open_report()`.
#' @param prefix File-name prefix ordering the outputs of a run.
#' @return The summary table (n, median and mean load per group), invisibly.
#' @export
section_mutation_load <- function(db, x_col, facet_rows = NULL,
                                  facet_cols = "subject_id", x_levels = NULL,
                                  shm_col = "mu_freq", fill_values = NULL,
                                  title = "VH somatic hypermutation load",
                                  name = "mutation_load",
                                  run_dir, cfg, report, prefix = "10") {
  report$h1(title)
  require_columns(db, c(x_col, facet_rows, facet_cols, shm_col),
                  "section_mutation_load()")
  d <- igh_rows(db)
  d$.shm <- 100 * as.numeric(d[[shm_col]])
  d[[x_col]] <- ordered_groups(d[[x_col]], x_levels)
  d <- d[!is.na(d$.shm) & !is.na(d[[x_col]]), , drop = FALSE]
  if (!nrow(d)) {
    report$p("No heavy chain with a mutation frequency in the groups requested.")
    return(invisible(NULL))
  }

  tab <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(d, dplyr::across(dplyr::all_of(
      unique(c(facet_rows, facet_cols, x_col))))),
    n = dplyr::n(),
    median_pct = stats::median(.data$.shm),
    mean_pct = mean(.data$.shm),
    .groups = "drop"))
  report$p(sprintf(paste0(
    "Mutation frequency (%s, percent) of each heavy chain by %s. The load ",
    "describes the maturation state of the cells sampled, not mutation ",
    "acquired during a culture."), shm_col, x_col))
  report$p("")
  report$tbl(as.data.frame(tab))
  file_stem <- paste0(prefix, "_", name)
  save_table(tab, file_stem, run_dir, "mutation load per group (percent)")

  p <- mutation_load_plot(d, x_col, facet_rows = facet_rows,
                          facet_cols = facet_cols, fill_values = fill_values,
                          title = title)
  size <- if (length(facet_rows)) c(9, 8) else c(11, 5)
  save_figure(p, file_stem, run_dir, cfg, width = size[1], height = size[2],
              legend = paste0(
                "Somatic hypermutation per heavy chain (", shm_col,
                ", percent) by ", x_col, ". It reflects the maturation state ",
                "of the cells sampled, not de novo mutation."))
  invisible(tab)
}

#' Replacement and silent mutations, CDR against FWR, per group
#'
#' R/S ratios summed over each subject and group. A higher ratio in the CDR
#' than in the FWR is consistent with antigen-driven positive selection on
#' the paratope. The ratio is descriptive: its expectation without selection
#' depends on codon composition and mutation hotspots, which is why
#' section_baseline_selection() is the test (Hershberg et al. 2008).
#'
#' Counts are computed on `IMGT_V` (CDR1-2, FWR1-3; configuration key
#' `bcr$selection$rs_region`) when `db` does not carry them already.
#'
#' @inheritParams section_mutation_load
#' @param group_col Grouping column (e.g. compartment).
#' @param subject_col Subject column; one panel per subject.
#' @param group_levels Optional order of the groups; others are left out.
#' @param region Region definition used when counts must be computed.
#' @return The R/S table, invisibly.
#' @export
section_rs_ratio <- function(db, group_col, subject_col = "subject_id",
                             group_levels = NULL, region = NULL,
                             name = "rs_ratio", run_dir, cfg, report,
                             prefix = "11") {
  s <- bcr_advanced_settings(cfg)
  report$h1("Replacement and silent mutations, CDR against FWR")
  require_columns(db, c(group_col, subject_col), "section_rs_ratio()")
  counts <- c("mu_count_cdr_r", "mu_count_cdr_s",
              "mu_count_fwr_r", "mu_count_fwr_s")
  if (!all(counts %in% names(db))) {
    region <- region %||% s$rs_region
    db <- observed_mutation_profile(db, region = region, counts_only = TRUE,
                                    nproc = s$nproc)
    report$p(sprintf("Replacement and silent mutations counted on %s.", region))
  }
  heavy <- igh_rows(db)
  heavy[[group_col]] <- ordered_groups(heavy[[group_col]], group_levels)
  heavy <- heavy[!is.na(heavy[[group_col]]), , drop = FALSE]

  rs <- replacement_silent_ratio(heavy, c(subject_col, group_col))
  report$p(paste0(
    "Mutations are summed over all heavy chains of each group before the ",
    "ratio is taken; the silent count is floored at 1. R/S is descriptive: ",
    "BASELINe is the test of selection."))
  report$p("")
  report$tbl(as.data.frame(rs))
  file_stem <- paste0(prefix, "_", name)
  save_table(rs, file_stem, run_dir,
             "summed R and S counts and R/S ratios, CDR and FWR, per group")

  long <- rbind(
    data.frame(rs[, c(subject_col, group_col)], region = "CDR",
               rs = rs$rs_cdr, check.names = FALSE),
    data.frame(rs[, c(subject_col, group_col)], region = "FWR",
               rs = rs$rs_fwr, check.names = FALSE))
  p <- ggplot2::ggplot(long, ggplot2::aes(.data[[group_col]], .data$rs,
                                          fill = .data$region)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap(ggplot2::vars(.data[[subject_col]]), nrow = 1) +
    ggplot2::scale_fill_manual(values = c(CDR = "#d95f02", FWR = "#7570b3"),
                               name = "Region") +
    ggplot2::labs(title = "Replacement/Silent (R/S) ratio, CDR vs FWR",
                  x = NULL, y = "R/S ratio") +
    theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_figure(p, file_stem, run_dir, cfg, width = 12, height = 4,
              legend = paste(
                "Ratio of replacement to silent mutations, summed over the",
                "heavy chains of each group, in CDR and in FWR. A higher R/S",
                "in CDR than in FWR is consistent with antigen-driven",
                "positive selection on the paratope."))
  invisible(rs)
}

#' Selection strength per group (BASELINe)
#'
#' BASELINe Sigma estimated per clone within each group, then combined per
#' group (Yaari et al. 2012); clones spanning several groups are treated
#' independently in each. Defaults follow the original analysis and can be
#' changed under `bcr$selection`: region `IMGT_V_BY_REGIONS`, focused test
#' statistic, clones collapsed to their most common sequence, groups with
#' fewer than 5 sequences or 2 clones skipped.
#'
#' The original single-cell analysis runs it twice: across ex vivo
#' compartments, and across proliferating and non-proliferating cells of
#' each stimulation condition. Restrict `db` and call it once per view.
#'
#' @inheritParams section_mutation_load
#' @param group_col Column defining the groups.
#' @param groups Groups, in order; all by default.
#' @param label Optional label of this view, used in titles.
#' @param clone_col,subject_col Clone and subject columns.
#' @return The result of `baseline_selection()`, invisibly.
#' @export
section_baseline_selection <- function(db, group_col, groups = NULL,
                                       label = NULL, clone_col = "clone_id",
                                       subject_col = "subject_id",
                                       name = "baseline_selection",
                                       run_dir, cfg, report, prefix = "12") {
  need("shazam", "section_baseline_selection()")
  s <- bcr_advanced_settings(cfg)
  report$h1(paste0("Selection (BASELINe)",
                   if (!is.null(label)) paste0(": ", label) else ""))
  report$p(sprintf(paste0(
    "BASELINe Sigma per clone within each %s (clones spanning several ",
    "groups are treated independently), combined per group. Region %s, %s ",
    "test statistic, clones collapsed to their %s sequence; groups with ",
    "fewer than %d sequences or %d clones are skipped. Sigma > 0: positive ",
    "selection; Sigma < 0: negative selection."),
    group_col, s$baseline_region, s$baseline_test,
    sub("^most", "most ", tolower(s$collapse_method)),
    s$baseline_min_sequences, s$baseline_min_clones))

  res <- baseline_selection(db, group_col = group_col, groups = groups,
                            clone_col = clone_col, subject_col = subject_col,
                            region = s$baseline_region,
                            test_statistic = s$baseline_test,
                            collapse_method = s$collapse_method,
                            min_sequences = s$baseline_min_sequences,
                            min_clones = s$baseline_min_clones,
                            nproc = s$nproc)
  report$h2("Groups")
  report$tbl(res$groups)
  file_stem <- paste0(prefix, "_", name)
  save_table(res$groups, paste0(file_stem, "_groups"), run_dir,
             "sequences, clones and BASELINe status per group")

  if (is.null(res$baseline)) {
    err <- attr(res, "error")
    report$p(if (!is.null(err)) paste("calcBaseline failed:", err) else
      "No group had enough sequences and clones; BASELINe not computed.")
    return(invisible(res))
  }
  report$h2("Sigma per group and region")
  report$tbl(res$stats)
  save_table(res$stats, paste0(file_stem, "_sigma"), run_dir,
             "BASELINe Sigma and confidence interval per group and region")

  p <- tryCatch(
    shazam::plotBaselineSummary(res$baseline, group_col, style = "summary",
                                silent = TRUE) +
      ggplot2::labs(title = paste0("BASELINe selection (Sigma)",
                                   if (!is.null(label)) paste0(" - ", label) else "")),
    error = function(e) e)
  if (inherits(p, "error")) {
    report$p(paste("BASELINe figure not drawn:", conditionMessage(p)))
  } else {
    n_groups <- sum(res$groups$status == "computed")
    save_figure(p, file_stem, run_dir, cfg,
                width = max(7, 1.4 * n_groups + 3), height = 5,
                legend = paste(
                  "BASELINe Sigma computed per clone within each group",
                  "(clones spanning several groups are treated",
                  "independently). Sigma > 0 positive selection, Sigma < 0",
                  "negative. Yaari et al. 2012."))
  }
  invisible(res)
}

#' Mutation load of class-switched and unswitched cells
#'
#' Mean mutation load per subject, compartment and isotype class, its range
#' across subjects, and the share of unswitched cells per subject and
#' compartment. The comparison shows whether a difference in mean mutation
#' load between compartments, or between groups of subjects, comes from the
#' class-switched cells themselves or from a different proportion of
#' unswitched cells. Missing isotype calls and stimulated material are
#' excluded; the switched list is `bcr$selection$switched_isotypes`.
#'
#' @inheritParams section_mutation_load
#' @param compartment_col,subject_col Compartment and subject columns.
#' @param subject_group_col Optional column grouping subjects.
#' @param shm_col Mutation frequency column (a fraction, reported in percent).
#' @param isotype_col Isotype column (derived from `c_call` when absent).
#' @return The result of `shm_by_switch_status()`, invisibly.
#' @export
section_shm_switched <- function(db, compartment_col = "compartment",
                                 subject_col = "subject_id",
                                 subject_group_col = NULL,
                                 shm_col = "mu_freq", isotype_col = "c_gene",
                                 name = "shm_switched", run_dir, cfg, report,
                                 prefix = "13") {
  s <- bcr_advanced_settings(cfg)
  report$h1("Mutation load of class-switched and unswitched cells")
  report$p(paste0(
    "A difference in mean mutation load can come from the class-switched ",
    "cells themselves or from a different share of unswitched cells, which ",
    "carry few mutations. Load is therefore reported per isotype class, ",
    "with the unswitched share beside it."))
  res <- shm_by_switch_status(db, shm_col = shm_col, isotype_col = isotype_col,
                              compartment_col = compartment_col,
                              subject_col = subject_col,
                              subject_group_col = subject_group_col,
                              switched = s$switched_isotypes,
                              exclude_stimulated = s$exclude_stimulated)
  report$p(sprintf(paste0(
    "Class-switched: %s. Excluded: %d heavy chains without an isotype call, ",
    "%d from stimulated material."),
    paste(s$switched_isotypes, collapse = ", "),
    res$excluded[["no_isotype"]], res$excluded[["stimulated"]]))
  if (nrow(res$outside_lists)) {
    report$p(sprintf(paste0(
      "Calls in neither the switched list nor IGHM/IGHD, counted as ",
      "unswitched under this rule: %s."),
      paste(sprintf("%s (%d)", res$outside_lists$isotype,
                    res$outside_lists$n_cells), collapse = ", ")))
  }

  summ <- as.data.frame(res$summary)
  summ$range <- sprintf("%.1f%% to %.1f%%", summ$shm_min, summ$shm_max)
  report$h2("Per group of subjects, compartment and isotype class")
  report$tbl(summ)
  report$h2("Unswitched share and overall load per subject and compartment")
  report$tbl(as.data.frame(res$unswitched_share))

  file_stem <- paste0(prefix, "_", name)
  save_table(res$by_subject, paste0(file_stem, "_by_subject"), run_dir,
             "cells and mean mutation load per subject, compartment and isotype class")
  save_table(res$summary, paste0(file_stem, "_summary"), run_dir,
             "mean, minimum and maximum of subject means per compartment and isotype class")
  save_table(res$unswitched_share, paste0(file_stem, "_unswitched_share"),
             run_dir, "unswitched share and overall load per subject and compartment")
  invisible(res)
}

#' Somatic hypermutation by isotype (bulk)
#'
#' The bulk analysis: mutation counts and frequencies by region on
#' `IMGT_VDJ_BY_REGIONS` (configuration key `bcr$mutation_by_isotype$region`),
#' then the total mutation frequency and the CDR3 replacement frequency by
#' group and isotype, and the density of the per-clone mean frequency.
#' Isotypes outside the configured order (IGHD, IGHM, IGHG1-4, IGHA1-2, IGHE
#' by default) are left out of the figures, and counted in the report.
#'
#' @inheritParams section_mutation_load
#' @param group_col Sample or condition column on the x axis.
#' @param group_levels Optional order of the groups; others are left out.
#' @param isotype_col Isotype column (derived from `c_call` when absent).
#' @param clone_col,subject_col Clone and subject columns.
#' @param region Region definition; the configured one by default.
#' @return The heavy chains with the mutation profile, invisibly; the name of
#'   the total frequency column is in attribute `total_mutation_col`.
#' @export
section_mutation_by_isotype <- function(db, group_col = "sample_id",
                                        group_levels = NULL,
                                        isotype_col = "c_gene",
                                        clone_col = "clone_id",
                                        subject_col = "subject_id",
                                        region = NULL, run_dir, cfg, report,
                                        prefix = "10") {
  need("shazam", "section_mutation_by_isotype()")
  s <- bcr_advanced_settings(cfg)
  region <- region %||% s$isotype_region
  report$h1("Somatic hypermutation by isotype")

  prof <- observed_mutation_profile(db, region = region, nproc = s$nproc)
  total_col <- attr(prof, "total_mutation_col")
  prof <- with_isotype_col(prof, isotype_col)
  report$p(sprintf(paste0(
    "Mutations counted on %s against %s: counts and frequencies by region, ",
    "replacement and silent, and in total (%s)."),
    region, attr(prof, "germline_col"), total_col))

  iso <- as.character(prof[[isotype_col]])
  in_order <- iso %in% s$isotype_levels
  if (any(!in_order)) {
    left_out <- table(ifelse(is.na(iso[!in_order]), "(none)", iso[!in_order]))
    report$p(sprintf("Isotypes outside the configured order, left out of the figures: %s.",
                     paste(sprintf("%s (%d)", names(left_out), as.integer(left_out)),
                           collapse = ", ")))
  }
  d <- prof[which(in_order), , drop = FALSE]
  d[[isotype_col]] <- factor(as.character(d[[isotype_col]]),
                             levels = s$isotype_levels)
  d[[group_col]] <- ordered_groups(d[[group_col]], group_levels)
  d <- d[!is.na(d[[group_col]]), , drop = FALSE]
  if (!nrow(d)) {
    report$p("No heavy chain left in the configured isotypes and groups.")
    return(invisible(prof))
  }

  has_cdr3 <- "mu_freq_cdr3_r" %in% names(d)
  tab <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(d, dplyr::across(dplyr::all_of(c(group_col, isotype_col)))),
    n = dplyr::n(),
    median_mu_freq = stats::median(.data[[total_col]], na.rm = TRUE),
    median_mu_freq_cdr3_r = if (has_cdr3)
      stats::median(.data$mu_freq_cdr3_r, na.rm = TRUE) else NA_real_,
    .groups = "drop"))
  report$tbl(as.data.frame(tab))
  save_table(tab, paste0(prefix, "a_mutation_by_isotype"), run_dir,
             "median mutation frequencies per group and isotype")

  p1 <- mutation_by_isotype_plot(d, total_col, group_col, isotype_col,
                                 title = "VH mutation frequency")
  save_figure(p1, paste0(prefix, "a_vh_mutation_by_isotype"), run_dir, cfg,
              width = 12, height = 7,
              legend = sprintf(paste(
                "Mutation frequency of each heavy chain over %s,",
                "by %s (x) and isotype (panels). Outliers not drawn."),
                region, group_col))
  if (has_cdr3) {
    p2 <- mutation_by_isotype_plot(d, "mu_freq_cdr3_r", group_col, isotype_col,
                                   title = "CDR3 replacement mutations")
    save_figure(p2, paste0(prefix, "b_cdr3_replacement_by_isotype"), run_dir,
                cfg, width = 12, height = 7,
                legend = sprintf(paste(
                  "Frequency of replacement mutations in CDR3 of each heavy",
                  "chain, by %s (x) and isotype (panels). Outliers not drawn."),
                  group_col))
  } else {
    report$p(sprintf("%s has no CDR3 region; the CDR3 replacement figure is skipped.",
                     region))
  }

  cm <- clone_mean_mutation(d, shm_col = total_col, clone_col = clone_col,
                            subject_col = subject_col, group_col = group_col,
                            isotype_col = isotype_col)
  save_table(cm, paste0(prefix, "c_clone_mean_mutation"), run_dir,
             "mean mutation frequency per clone, group and isotype")
  p3 <- clone_mutation_density_plot(cm, group_col, isotype_col)
  # A panel with fewer than two clones stays empty, as in the original
  # figure; ggplot2's notices about such panels are not repeated.
  withCallingHandlers(
    save_figure(p3, paste0(prefix, "c_clone_mutation_density"), run_dir, cfg,
                width = 10, height = 8,
                legend = paste(
                  "Density of the mean mutation frequency of each clone, by",
                  "isotype (rows) and", group_col, "(columns). Panels with",
                  "fewer than two clones are empty.")),
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("fewer than two data points", msg, fixed = TRUE) ||
          grepl("no non-missing arguments to max", msg, fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    })

  attr(prof, "total_mutation_col") <- total_col
  invisible(prof)
}

#' IGHV gene by CDR3 length overview (bulk)
#'
#' One point per clone within each group, by IGHV gene and junction length,
#' coloured by the median of a mutation measure and sized by the number of
#' sequences, one panel per group and isotype. The defaults reproduce the
#' original CDR3 replacement figure (colour scale 0-30%); IgE and points
#' above 250 sequences are not drawn (`bcr$overview`).
#'
#' @inheritParams section_mutation_load
#' @param value_col Mutation measure (default the CDR3 replacement frequency
#'   added by `section_mutation_by_isotype()`).
#' @param group_col Sample or condition column (panel rows).
#' @param group_levels Optional order of the groups; others are left out.
#' @param isotype_col Isotype column (panel columns).
#' @param clone_col,subject_col Clone and subject columns.
#' @param clones Optional subject-qualified clone keys to restrict to.
#' @param colour_title Colour legend title.
#' @param colour_values,colour_limits Positions of the four colours of the
#'   scale, and its limits, in percent.
#' @param size_breaks Breaks of the size legend.
#' @param unit What a row counts, "sequences" (bulk) or "cells".
#' @return The overview table, invisibly.
#' @export
section_vgene_cdr3_overview <- function(db, value_col = "mu_freq_cdr3_r",
                                        group_col = "sample_id",
                                        group_levels = NULL,
                                        isotype_col = "c_gene",
                                        clone_col = "clone_id",
                                        subject_col = "subject_id",
                                        clones = NULL,
                                        title = "Frequency of CDR3 replacement mutations",
                                        colour_title = "Median CDR3 replacement SHM %",
                                        colour_values = c(0, 10, 20, 30),
                                        colour_limits = c(0, 30),
                                        size_breaks = c(2, 20, 50, 100, 180, 290),
                                        unit = "sequences",
                                        name = "vgene_cdr3_overview",
                                        run_dir, cfg, report, prefix = "11") {
  report$h1(title)
  if (!value_col %in% names(db)) {
    stop("Column '", value_col, "' not found. The CDR3 replacement frequency ",
         "is added by section_mutation_by_isotype().", call. = FALSE)
  }
  write_vgene_cdr3_overview(
    db, value_col = value_col, clone_col = clone_col, subject_col = subject_col,
    group_col = group_col, group_levels = group_levels,
    isotype_col = isotype_col, clones = clones, title = title,
    colour_title = colour_title, unit = unit, colour_values = colour_values,
    colour_limits = colour_limits, size_breaks = size_breaks,
    name = paste0(prefix, "_", name), run_dir = run_dir, cfg = cfg,
    report = report)
}

#' The largest clones and their mutation load (bulk)
#'
#' Summary of the `n` largest heavy-chain clones (100 by default,
#' `bcr$overview$top_n_clones`): size, mean and median mutation frequency,
#' V and J genes, groups. Then the IGHV by CDR3 length overview restricted
#' to them, coloured by median VH mutation (scale 0-15%), as in the
#' original analysis.
#'
#' @inheritParams section_vgene_cdr3_overview
#' @param shm_col Mutation frequency column. By default the total frequency
#'   recorded by `section_mutation_by_isotype()` (attribute
#'   `total_mutation_col`), else `mu_freq`.
#' @param n Number of clones; the configured number by default.
#' @return The summary of the largest clones, invisibly.
#' @export
section_top_clones <- function(db, shm_col = NULL, n = NULL,
                               group_col = "sample_id", group_levels = NULL,
                               isotype_col = "c_gene", clone_col = "clone_id",
                               subject_col = "subject_id", unit = "sequences",
                               run_dir, cfg, report, prefix = "12") {
  s <- bcr_advanced_settings(cfg)
  n <- as.integer(n %||% s$top_n)
  shm_col <- shm_col %||% attr(db, "total_mutation_col") %||% "mu_freq"
  report$h1(sprintf("The %d largest clones", n))

  top <- top_clones_summary(db, n = n, shm_col = shm_col, clone_col = clone_col,
                            subject_col = subject_col, group_col = group_col)
  report$p(sprintf(paste0(
    "Heavy-chain clones ranked by number of %s; mutation frequency from %s."),
    unit, shm_col))
  report$p("")
  report$tbl(as.data.frame(utils::head(top, 20)))
  if (nrow(top) > 20) report$p(sprintf("(first 20 of %d shown)", nrow(top)))
  save_table(top, paste0(prefix, "a_top_clones"), run_dir,
             sprintf("the %d largest clones with mutation load and genes", n))

  report$h2(sprintf("IGHV by CDR3 length, %d largest clones", n))
  write_vgene_cdr3_overview(
    db, value_col = shm_col, clone_col = clone_col, subject_col = subject_col,
    group_col = group_col, group_levels = group_levels,
    isotype_col = isotype_col,
    clones = subject_clone_key(top, clone_col, subject_col),
    title = sprintf("Top %d clones - Frequency of VH mutations", n),
    colour_title = "Median VH SHM %", unit = unit,
    colour_values = c(0, 10, 15), colour_limits = c(0, 15),
    size_breaks = c(2, 20, 50, 100, 180, 260),
    name = paste0(prefix, "b_top_clones_overview"), run_dir = run_dir,
    cfg = cfg, report = report)
  invisible(top)
}
