# ---------------------------------------------------------------------------
# Clone tracking: following clonotypes across samples, time and compartments
#
# Three views, each answering a different question.
#
#   1. DOMINANT CLONOTYPES ACROSS SAMPLES. The union of each sample's top-N
#      clonotypes, followed across an ordered series of samples (conditions
#      or timepoints of one subject). This is the clonalCompare() view of
#      scRepertoire (Borcherding et al. 2020, F1000Research 9:47,
#      PMID 32789006): proportions are taken over each sample's whole
#      repertoire before any selection; the top N is taken per sample, with
#      ties broken by order rather than by keeping extra clonotypes; and a
#      clonotype absent from a sample simply has no lode there.
#
#   2. BASELINE AND EXPANDED CLONOTYPES OVER TIME, within one subject. The
#      clonotypes that dominated the baseline sample, and the clonotypes
#      called significantly expanded against it by test_clonal_abundance(),
#      de novo and pre-existing kept apart, followed across every timepoint
#      the subject has. Frequencies are completed with zero where a
#      clonotype is absent: without the zero, a ribbon breaks exactly where
#      the clonotype's absence is the finding.
#
#   3. CLONOTYPES SHARED BETWEEN COMPARTMENTS of one subject (sorted
#      populations, tissues, culture conditions), optionally highlighting
#      the clonotypes also found in one designated compartment.
#
# The alluvial figures are drawn with ggalluvial (Brunson 2020, Journal of
# Open Source Software 5:2017, doi:10.21105/joss.02017). Ribbon height is a
# frequency or a count in the pipeline's unit (molecules in bulk data, cells
# in single-cell data), and the axis label says which.
# ---------------------------------------------------------------------------

# --- internal helpers -------------------------------------------------------

#' Stop with an informative message when columns are missing
#' @noRd
tracking_check_columns <- function(df, cols, what) {
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols)) {
    stop(what, " needs column(s): ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Display order of a sample-level variable
#'
#' Factor levels when the column is a factor, else the configured order,
#' else first appearance. Values present but absent from the declared order
#' are appended rather than dropped, so nothing disappears from a figure.
#' @noRd
tracking_levels <- function(x, configured = NULL) {
  present <- unique(as.character(x[!is.na(x)]))
  lv <- if (is.factor(x)) levels(x) else if (length(configured))
    as.character(configured) else present
  c(intersect(lv, present), setdiff(present, lv))
}

#' File-name-safe version of a label
#' @noRd
tracking_stem <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
}

#' Tracking parameters from the configuration, with their defaults
#'
#' Read from cfg$repertoire$tracking. The defaults are the values of the
#' original analyses: top 10 and top 20 clonotypes per sample, top 30 in raw
#' counts for a pairwise comparison, the 15 dominant baseline clonotypes,
#' the 15 expanded clonotypes with the largest fold change per subject and
#' type, and clonotypes present in at least two compartments.
#' @noRd
tracking_config <- function(cfg) {
  tr <- cfg$repertoire$tracking %||% list()
  out <- list(
    compare_top        = as.integer(unlist(tr$compare_top %||% c(10L, 20L))),
    compare_proportion = isTRUE(tr$compare_proportion %||% TRUE),
    compare_pairs      = lapply(tr$compare_pairs %||% list(), as.character),
    compare_pairs_top  = as.integer(tr$compare_pairs_top %||% 30L),
    baseline_top       = as.integer(tr$baseline_top %||% 15L),
    expanded_top       = as.integer(tr$expanded_top %||% 15L),
    expanded_rank_by   = as.character(tr$expanded_rank_by %||% "first"),
    compartments       = if (length(tr$compartments))
      as.character(unlist(tr$compartments)) else NULL,
    highlight_compartment = if (length(tr$highlight_compartment))
      as.character(tr$highlight_compartment)[1] else NULL,
    min_compartments   = as.integer(tr$min_compartments %||% 2L)
  )
  ints <- c(out$compare_top, out$compare_pairs_top, out$baseline_top,
            out$expanded_top, out$min_compartments)
  if (!length(out$compare_top) || anyNA(ints) || any(ints < 1L)) {
    stop("repertoire$tracking: compare_top, compare_pairs_top, baseline_top, ",
         "expanded_top and min_compartments must be positive integers.",
         call. = FALSE)
  }
  if (!out$expanded_rank_by %in% c("first", "peak")) {
    stop("repertoire$tracking$expanded_rank_by must be 'first' or 'peak'.",
         call. = FALSE)
  }
  bad_pair <- vapply(out$compare_pairs, length, integer(1)) != 2L
  if (any(bad_pair)) {
    stop("repertoire$tracking$compare_pairs must be a list of pairs of ",
         "sample labels, e.g. [[label_a, label_b]].", call. = FALSE)
  }
  out
}

#' Human-readable label for an expansion category
#' @noRd
tracking_type_label <- function(type) {
  ifelse(type == "de_novo", "de novo", "pre-existing expanded")
}

# --- 1. dominant clonotypes across samples ---------------------------------

#' Dominant clonotypes followed across an ordered series of samples
#'
#' The table behind a clonalCompare-style alluvial (scRepertoire): for every
#' sample, the abundance of each clonotype and its proportion of the sample;
#' then the union of each sample's `top_n` most abundant clonotypes, kept in
#' every selected sample where they occur.
#'
#' Proportions are computed over each sample's whole repertoire, before any
#' sample or clonotype is set aside, so they do not depend on which other
#' samples are compared. Within a sample, ties at the `top_n` boundary are
#' broken by clonotype name (as `slice_max(with_ties = FALSE)` does in
#' scRepertoire), so exactly `top_n` clonotypes come from each sample.
#'
#' @param clono Clonotype table: one row per clonotype per sample, as from
#'   `collapse_clonotypes()`, or one row per cell or sequence with
#'   `abundance_col = NULL`.
#' @param samples Samples to compare, in display order. NULL keeps every
#'   sample, in factor-level order when the sample column is a factor and in
#'   order of first appearance otherwise.
#' @param top_n Number of clonotypes taken from each sample.
#' @param abundance_col Column holding each row's abundance, or NULL to count
#'   rows.
#' @param clone_col,sample_col Names of the clonotype and sample columns.
#' @return Data frame: sample_id (factor, display order), clonotype,
#'   count, proportion, rank (within the sample, 1 = most abundant) and
#'   in_top (whether the clonotype is among that sample's own `top_n`).
#' @export
clone_compare_table <- function(clono, samples = NULL, top_n = 10,
                                abundance_col = "abundance",
                                clone_col = "clonotype",
                                sample_col = "sample_id") {
  need("dplyr", "clone_compare_table()")
  tracking_check_columns(clono, c(sample_col, clone_col, abundance_col),
                         "clone_compare_table()")
  stopifnot(length(top_n) == 1L, top_n >= 1)

  d <- data.frame(
    sample_id = as.character(clono[[sample_col]]),
    clonotype = as.character(clono[[clone_col]]),
    w = if (is.null(abundance_col)) 1 else as.numeric(clono[[abundance_col]]),
    stringsAsFactors = FALSE)
  d <- d[!is.na(d$sample_id) & !is.na(d$clonotype) & nzchar(d$clonotype) &
           !is.na(d$w), , drop = FALSE]
  if (!nrow(d)) stop("clone_compare_table(): no clonotype left.", call. = FALSE)

  order_lv <- if (is.null(samples)) {
    tracking_levels(clono[[sample_col]])
  } else {
    absent <- setdiff(as.character(samples), d$sample_id)
    if (length(absent)) {
      warning("Sample(s) not found and skipped: ",
              paste(absent, collapse = ", "), call. = FALSE)
    }
    intersect(as.character(samples), d$sample_id)
  }
  if (!length(order_lv)) {
    stop("clone_compare_table(): none of the requested samples is present.",
         call. = FALSE)
  }

  agg <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(d, .data$sample_id, .data$clonotype),
    count = sum(.data$w), .groups = "drop"))
  agg <- agg[order(agg$sample_id, -agg$count, agg$clonotype,
                   method = "radix"), , drop = FALSE]
  agg <- dplyr::ungroup(dplyr::mutate(
    dplyr::group_by(agg, .data$sample_id),
    proportion = .data$count / sum(.data$count),
    rank = dplyr::row_number()))
  agg$in_top <- agg$rank <= top_n

  agg <- agg[agg$sample_id %in% order_lv, , drop = FALSE]
  keep <- unique(agg$clonotype[agg$in_top])
  out <- agg[agg$clonotype %in% keep, , drop = FALSE]

  out$sample_id <- factor(out$sample_id, levels = order_lv)
  out$clonotype <- factor(out$clonotype,
                          levels = sort(unique(out$clonotype), method = "radix"))
  out <- out[order(out$sample_id, out$rank), , drop = FALSE]
  as.data.frame(out[, c("sample_id", "clonotype", "count", "proportion",
                        "rank", "in_top")])
}

#' Alluvial plot of dominant clonotypes across samples
#'
#' Strata stacked per sample and flows joining the same clonotype from one
#' sample to the next (`geom_stratum()` + `geom_flow(stat = "alluvium")`),
#' filled along the inferno palette, as scRepertoire's clonalCompare() draws
#' them.
#'
#' @param tab Output of `clone_compare_table()`.
#' @param value "proportion" (default) or "count".
#' @param title Plot title.
#' @param unit Unit of `count`, used in the axis label.
#' @param show_legend Show the clonotype legend. Off by default: clonotype
#'   keys are long, and the identities are in the table written beside the
#'   figure.
#' @return A ggplot object.
#' @export
plot_clone_compare <- function(tab, value = c("proportion", "count"),
                               title = NULL, unit = "cells",
                               show_legend = FALSE) {
  value <- match.arg(value)
  need("ggalluvial", "plot_clone_compare()")
  tab$clonotype <- droplevels(as.factor(tab$clonotype))
  n_clones <- nlevels(tab$clonotype)

  p <- ggplot2::ggplot(tab, ggplot2::aes(x = .data$sample_id,
                                         y = .data[[value]],
                                         fill = .data$clonotype,
                                         group = .data$clonotype,
                                         stratum = .data$clonotype,
                                         alluvium = .data$clonotype)) +
    ggalluvial::geom_stratum() +
    ggalluvial::geom_flow(stat = ggalluvial::StatAlluvium) +
    ggplot2::scale_fill_manual(
      values = grDevices::hcl.colors(n_clones, "inferno", fixup = TRUE)) +
    ggplot2::labs(x = NULL,
                  y = if (value == "proportion") "Proportion" else
                    sprintf("Count (%s)", unit),
                  fill = "Clonotype", title = title) +
    theme_pub(grid = "none") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  if (!show_legend) p <- p + ggplot2::theme(legend.position = "none")
  p
}

# --- 2. baseline and expanded clonotypes over time --------------------------

#' Clonotypes dominant at baseline, followed over each subject's timepoints
#'
#' For each subject, the `top_n` most abundant clonotypes of the baseline
#' sample, with their frequency in every sample of that subject, and zero
#' where a clonotype is absent. This is the trackClonotypes() view of
#' immunarch with the baseline sample as reference.
#'
#' @param clono Output of `collapse_clonotypes()` joined to the sample
#'   metadata: sample_id, subject_id, clonotype, abundance, frequency and the
#'   time column.
#' @param top_n Number of baseline clonotypes followed.
#' @param baseline Baseline level of the time column. NULL, or a level the
#'   subject does not have, takes that subject's earliest timepoint.
#' @param time_col Column ordering a subject's samples.
#' @param time_levels Order of the time column when it is not a factor.
#' @return Data frame: subject_id, sample_id, the time column (factor),
#'   clonotype, baseline_rank, abundance and frequency, one row per followed
#'   clonotype per sample of the subject.
#' @export
track_baseline_clonotypes <- function(clono, top_n = 15, baseline = NULL,
                                      time_col = "timepoint",
                                      time_levels = NULL) {
  tracking_check_columns(clono, c("sample_id", "subject_id", "clonotype",
                                  "abundance", time_col),
                         "track_baseline_clonotypes()")
  lv <- tracking_levels(clono[[time_col]], time_levels)
  if (!"frequency" %in% names(clono)) {
    clono$frequency <- clono$abundance /
      stats::ave(clono$abundance, as.character(clono$sample_id), FUN = sum)
  }

  out <- list()
  for (subj in tracking_levels(clono$subject_id)) {
    sub <- clono[as.character(clono$subject_id) == subj, , drop = FALSE]
    t_chr <- as.character(sub[[time_col]])
    times <- lv[lv %in% t_chr]
    if (!length(times)) next
    base_t <- if (!is.null(baseline) && baseline %in% times) baseline else times[1]

    base <- sub[t_chr == base_t, , drop = FALSE]
    if (length(unique(as.character(base$sample_id))) > 1L) {
      stop("Subject '", subj, "' has more than one sample at the baseline ",
           "timepoint '", base_t, "'. Resolve which is the baseline first.",
           call. = FALSE)
    }
    base <- base[order(-base$abundance, as.character(base$clonotype),
                       method = "radix"), , drop = FALSE]
    top <- utils::head(as.character(base$clonotype), top_n)

    smp <- unique(sub[, c("sample_id", time_col), drop = FALSE])
    smp <- smp[order(match(as.character(smp[[time_col]]), lv)), , drop = FALSE]
    for (i in seq_len(nrow(smp))) {
      s_rows <- sub[as.character(sub$sample_id) == as.character(smp$sample_id[i]),
                    , drop = FALSE]
      idx <- match(top, as.character(s_rows$clonotype))
      ab <- s_rows$abundance[idx]
      fr <- s_rows$frequency[idx]
      ab[is.na(ab)] <- 0
      fr[is.na(fr)] <- 0
      piece <- data.frame(subject_id = subj,
                          sample_id = as.character(smp$sample_id[i]),
                          time = as.character(smp[[time_col]][i]),
                          clonotype = top, baseline_rank = seq_along(top),
                          abundance = ab, frequency = fr,
                          stringsAsFactors = FALSE)
      out[[length(out) + 1L]] <- piece
    }
  }
  if (!length(out)) {
    return(data.frame(subject_id = character(0), sample_id = character(0),
                      time = factor(character(0), levels = lv),
                      clonotype = character(0), baseline_rank = integer(0),
                      abundance = numeric(0), frequency = numeric(0)))
  }
  res <- do.call(rbind, out)
  res$time <- factor(res$time, levels = lv)
  names(res)[names(res) == "time"] <- time_col
  res
}

#' Significantly expanded clonotypes chosen for tracking
#'
#' From `test_clonal_abundance()` calls, the clonotypes called significantly
#' expanded, split into de novo and pre-existing expanded, and within each
#' subject and type the `top_n` with the largest log2 fold change. Ties at
#' the boundary are kept, so a few more than `top_n` can be returned.
#'
#' @param calls Output of `test_clonal_abundance()`.
#' @param top_n Clonotypes kept per subject and type.
#' @param rank_by "first" (default, as in the original analysis): a
#'   clonotype significant at several timepoints is ranked by its fold
#'   change at the earliest of them. "peak": by its largest fold change.
#' @param time_levels Order of the timepoints when `calls$timepoint` is not a
#'   factor.
#' @return Data frame: subject_id, type ("de_novo" or
#'   "pre_existing_expanded"), clonotype, log2_fold_change (the value ranked
#'   on), timepoint (where that value was taken) and
#'   n_timepoints_significant.
#' @export
select_tracked_expanded <- function(calls, top_n = 15,
                                    rank_by = c("first", "peak"),
                                    time_levels = NULL) {
  rank_by <- match.arg(rank_by)
  tracking_check_columns(calls, c("subject_id", "timepoint", "clonotype",
                                  "category", "log2_fold_change",
                                  "sig_expansion"),
                         "select_tracked_expanded()")
  empty <- data.frame(subject_id = character(0), type = character(0),
                      clonotype = character(0), log2_fold_change = numeric(0),
                      timepoint = character(0),
                      n_timepoints_significant = integer(0),
                      stringsAsFactors = FALSE)
  keep <- which(calls$sig_expansion &
                  calls$category %in% c("de_novo", "pre_existing_expanded"))
  if (!length(keep)) return(empty)
  sig <- calls[keep, , drop = FALSE]

  lv <- tracking_levels(sig$timepoint, time_levels)
  sig$subject_id <- as.character(sig$subject_id)
  sig$clonotype  <- as.character(sig$clonotype)
  sig$type       <- ifelse(sig$category == "de_novo", "de_novo",
                           "pre_existing_expanded")
  sig$t_order    <- match(as.character(sig$timepoint), lv)

  out <- list()
  for (key in unique(paste(sig$subject_id, sig$type, sep = "\r"))) {
    g <- sig[paste(sig$subject_id, sig$type, sep = "\r") == key, , drop = FALSE]
    g <- if (rank_by == "first") {
      g[order(g$clonotype, g$t_order, -g$log2_fold_change, method = "radix"), ]
    } else {
      g[order(g$clonotype, -g$log2_fold_change, g$t_order, method = "radix"), ]
    }
    n_sig <- table(g$clonotype)
    one <- g[!duplicated(g$clonotype), , drop = FALSE]
    cutoff <- sort(one$log2_fold_change, decreasing = TRUE)[min(top_n, nrow(one))]
    one <- one[one$log2_fold_change >= cutoff, , drop = FALSE]
    one <- one[order(-one$log2_fold_change, one$clonotype, method = "radix"), ]
    out[[length(out) + 1L]] <- data.frame(
      subject_id = one$subject_id, type = one$type, clonotype = one$clonotype,
      log2_fold_change = one$log2_fold_change,
      timepoint = as.character(one$timepoint),
      n_timepoints_significant = as.integer(n_sig[one$clonotype]),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

#' Frequencies of chosen clonotypes across each subject's samples
#'
#' @param clono Clonotype table: sample_id, subject_id, clonotype, abundance,
#'   frequency and the time column.
#' @param selection Data frame with subject_id and clonotype, and optionally
#'   type, such as the output of `select_tracked_expanded()`.
#' @param time_col Column ordering a subject's samples.
#' @param time_levels Order of the time column when it is not a factor.
#' @return Long data frame: subject_id, type (when given), clonotype,
#'   sample_id, the time column (factor), abundance and frequency, with zero
#'   for every sample of the subject in which the clonotype is absent.
#' @export
track_clonotypes <- function(clono, selection, time_col = "timepoint",
                             time_levels = NULL) {
  tracking_check_columns(clono, c("sample_id", "subject_id", "clonotype",
                                  "abundance", time_col), "track_clonotypes()")
  tracking_check_columns(selection, c("subject_id", "clonotype"),
                         "track_clonotypes()")
  lv <- tracking_levels(clono[[time_col]], time_levels)
  if (!"frequency" %in% names(clono)) {
    clono$frequency <- clono$abundance /
      stats::ave(clono$abundance, as.character(clono$sample_id), FUN = sum)
  }
  has_type <- "type" %in% names(selection)
  selection <- unique(selection[, intersect(c("subject_id", "type", "clonotype"),
                                            names(selection)), drop = FALSE])

  out <- list()
  for (subj in unique(as.character(selection$subject_id))) {
    sel <- selection[as.character(selection$subject_id) == subj, , drop = FALSE]
    sub <- clono[as.character(clono$subject_id) == subj, , drop = FALSE]
    if (!nrow(sub)) next
    smp <- unique(sub[, c("sample_id", time_col), drop = FALSE])
    smp <- smp[order(match(as.character(smp[[time_col]]), lv)), , drop = FALSE]
    for (i in seq_len(nrow(smp))) {
      s_rows <- sub[as.character(sub$sample_id) == as.character(smp$sample_id[i]),
                    , drop = FALSE]
      idx <- match(as.character(sel$clonotype), as.character(s_rows$clonotype))
      ab <- s_rows$abundance[idx]
      fr <- s_rows$frequency[idx]
      ab[is.na(ab)] <- 0
      fr[is.na(fr)] <- 0
      piece <- data.frame(subject_id = subj,
                          clonotype = as.character(sel$clonotype),
                          sample_id = as.character(smp$sample_id[i]),
                          time = as.character(smp[[time_col]][i]),
                          abundance = ab, frequency = fr,
                          stringsAsFactors = FALSE)
      if (has_type) piece$type <- as.character(sel$type)
      out[[length(out) + 1L]] <- piece
    }
  }
  if (!length(out)) return(NULL)
  res <- do.call(rbind, out)
  res$time <- factor(res$time, levels = lv)
  names(res)[names(res) == "time"] <- time_col
  cols <- c("subject_id", if (has_type) "type", "clonotype", "sample_id",
            time_col, "abundance", "frequency")
  res[, cols]
}

#' Alluvial plot of clonotype frequencies over time
#'
#' Ribbons join the same clonotype from one timepoint to the next
#' (`geom_alluvium()` over `geom_stratum()`, unsorted strata), one colour per
#' clonotype along viridis, without a legend.
#'
#' @param track Rows of `track_clonotypes()` or
#'   `track_baseline_clonotypes()` for one subject (and one type).
#' @param time_col Column on the x axis.
#' @param title,subtitle Plot title and subtitle.
#' @param y_label Axis label for the frequency.
#' @return A ggplot object.
#' @export
plot_clone_tracking <- function(track, time_col = "timepoint", title = NULL,
                                subtitle = NULL, y_label = "Frequency") {
  need("ggalluvial", "plot_clone_tracking()")
  track[[time_col]] <- droplevels(as.factor(track[[time_col]]))
  # A timepoint at which none of the clonotypes is detected (the baseline,
  # for de novo clonotypes) cannot be stacked by ggalluvial, whose lodes
  # need at least one non-zero height per axis. Its rows are dropped and the
  # axis kept, empty, so the absence still shows.
  all_zero <- tapply(track$frequency, track[[time_col]], function(v) all(v == 0))
  empty_axes <- names(all_zero)[which(all_zero)]
  track <- track[!as.character(track[[time_col]]) %in% empty_axes, , drop = FALSE]
  ggplot2::ggplot(track, ggplot2::aes(x = .data[[time_col]],
                                      y = .data$frequency,
                                      stratum = .data$clonotype,
                                      alluvium = .data$clonotype,
                                      fill = .data$clonotype)) +
    ggalluvial::geom_alluvium(alpha = 0.7, decreasing = FALSE) +
    ggalluvial::geom_stratum(alpha = 0.9, linewidth = 0.1, decreasing = FALSE) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_fill_viridis_d(option = "viridis", guide = "none") +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = y_label) +
    theme_pub(grid = "none")
}

# --- 3. clonotypes shared between compartments ------------------------------

#' Clonotypes shared between compartments of the same subject
#'
#' Within each subject, the abundance (cells, in single-cell data) of every
#' clonotype in every compartment, restricted to clonotypes present in at
#' least `min_compartments` of the listed compartments. A clonotype is
#' flagged as highlighted when the subject has it in the `highlight`
#' compartment.
#'
#' @param clono Clonotype table: subject_id, the clonotype and compartment
#'   columns, and the abundance column.
#' @param compartment_col Column naming the compartment.
#' @param compartments Compartments to include, in display order. NULL: all,
#'   in factor-level order or order of first appearance.
#' @param highlight Compartment whose clonotypes are highlighted, or NULL.
#' @param min_compartments Minimum number of compartments a clonotype must
#'   occupy to be shown.
#' @param clone_col,abundance_col Names of the clonotype and abundance
#'   columns.
#' @return Data frame: subject_id, compartment (factor), clonotype,
#'   abundance, n_compartments, highlighted.
#' @export
compartment_flows <- function(clono, compartment_col = "condition",
                              compartments = NULL, highlight = NULL,
                              min_compartments = 2,
                              clone_col = "clonotype",
                              abundance_col = "abundance") {
  need("dplyr", "compartment_flows()")
  tracking_check_columns(clono, c("subject_id", compartment_col, clone_col,
                                  abundance_col), "compartment_flows()")
  lv <- compartments %||% tracking_levels(clono[[compartment_col]])

  d <- data.frame(subject_id = as.character(clono$subject_id),
                  compartment = as.character(clono[[compartment_col]]),
                  clonotype = as.character(clono[[clone_col]]),
                  abundance = as.numeric(clono[[abundance_col]]),
                  stringsAsFactors = FALSE)
  d <- d[!is.na(d$clonotype) & !is.na(d$compartment), , drop = FALSE]

  hl <- if (is.null(highlight)) character(0) else
    unique(paste(d$subject_id, d$clonotype, sep = "\r")[d$compartment == highlight])

  d <- d[d$compartment %in% lv, , drop = FALSE]
  agg <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(d, .data$subject_id, .data$compartment, .data$clonotype),
    abundance = sum(.data$abundance), .groups = "drop"))
  agg <- dplyr::ungroup(dplyr::mutate(
    dplyr::group_by(agg, .data$subject_id, .data$clonotype),
    n_compartments = dplyr::n_distinct(.data$compartment)))
  agg <- agg[agg$n_compartments >= min_compartments, , drop = FALSE]
  agg$highlighted <- paste(agg$subject_id, agg$clonotype, sep = "\r") %in% hl
  agg$compartment <- factor(agg$compartment, levels = lv)
  agg <- agg[order(agg$subject_id, agg$compartment, -agg$abundance), , drop = FALSE]
  as.data.frame(agg[, c("subject_id", "compartment", "clonotype", "abundance",
                        "n_compartments", "highlighted")])
}

#' Alluvial plot of clonotypes shared between compartments
#'
#' One stratum per compartment and one flow per shared clonotype, orange for
#' highlighted clonotypes and light blue for the others.
#'
#' @param flows Rows of `compartment_flows()` for one subject.
#' @param title Plot title.
#' @param unit Unit of the abundance, for the axis label.
#' @param highlight_label Legend label of the highlighted clonotypes; NULL
#'   when nothing is highlighted (no legend).
#' @return A ggplot object.
#' @export
plot_compartment_flows <- function(flows, title = NULL, unit = "cells",
                                   highlight_label = NULL) {
  need("ggalluvial", "plot_compartment_flows()")
  flows$compartment <- droplevels(as.factor(flows$compartment))
  if (is.null(highlight_label)) {
    flows$group <- "Shared"
    values <- c(Shared = "#9ecae1")
  } else {
    flows$group <- ifelse(flows$highlighted, highlight_label, "Other")
    values <- stats::setNames(c("#e6550d", "#9ecae1"),
                              c(highlight_label, "Other"))
  }
  flows$group <- factor(flows$group, levels = names(values))
  p <- ggplot2::ggplot(flows, ggplot2::aes(x = .data$compartment,
                                           stratum = .data$compartment,
                                           alluvium = .data$clonotype,
                                           y = .data$abundance,
                                           fill = .data$group)) +
    ggalluvial::geom_flow(alpha = 0.6, width = 0.4) +
    ggalluvial::geom_stratum(width = 0.4, fill = "grey90", colour = "grey50") +
    ggplot2::scale_fill_manual(values = values, name = NULL, drop = FALSE) +
    ggplot2::labs(title = title, x = NULL,
                  y = sprintf("%s (shared clonotypes)",
                              paste0(toupper(substr(unit, 1, 1)),
                                     substring(unit, 2)))) +
    theme_pub()
  if (is.null(highlight_label)) p <- p + ggplot2::theme(legend.position = "none")
  p
}

# --- sections -----------------------------------------------------------------

#' Dominant clonotypes across samples (clonalCompare-style alluvials)
#'
#' For each subject (or other grouping), the union of each sample's top
#' clonotypes followed across the subject's samples in configured order, for
#' every value of `repertoire$tracking$compare_top` (10 and 20 by default),
#' and the top `compare_pairs_top` (30) in raw counts for each configured
#' pair of samples. Runs in both pipelines.
#'
#' @param clono Clonotype table: sample_id, clonotype, abundance, frequency
#'   and the grouping and ordering columns.
#' @param run_dir Run directory from `init_run()`.
#' @param cfg Configuration list.
#' @param report Report writer from `open_report()`.
#' @param within Column defining one alluvial per group (default
#'   "subject_id"); NULL compares all samples together.
#' @param order_col Column ordering and labelling the samples on the x axis.
#'   NULL takes "timepoint" when present, else "condition", else sample_id.
#' @param unit Unit of the abundance ("cells" or "molecules").
#' @param abundance_col Abundance column of `clono`, or NULL to count rows
#'   (for instance one row per sequence, as clonalCompare() counts a
#'   Change-O table keyed on clone_id).
#' @param prefix File name prefix.
#' @return The combined table of tracked clonotypes, invisibly.
#' @export
section_clone_compare <- function(clono, run_dir, cfg, report,
                                  within = "subject_id", order_col = NULL,
                                  unit = "cells", abundance_col = "abundance",
                                  prefix = "10") {
  need(c("dplyr", "ggalluvial"), "section_clone_compare()")
  report$h1("Dominant clonotypes across samples")
  tr <- tracking_config(cfg)
  tracking_check_columns(clono, c("sample_id", "clonotype", abundance_col),
                         "section_clone_compare()")

  if (is.null(order_col)) {
    usable <- function(col) col %in% names(clono) && !all(is.na(clono[[col]]))
    order_col <- if (usable("timepoint")) "timepoint" else
      if (usable("condition")) "condition" else "sample_id"
  }
  configured <- switch(order_col, timepoint = cfg$cohort$timepoints,
                       condition = cfg$cohort$conditions, NULL)

  grouped <- !is.null(within) && within %in% names(clono) &&
    !all(is.na(clono[[within]]))
  group_lv <- if (grouped) tracking_levels(clono[[within]]) else "all samples"

  report$p(paste0(
    "For each sample, its most abundant clonotypes; the union of these is ",
    "followed across the samples", if (grouped) paste0(" of each ",
                                                        sub("_id$", "", within))
    else "", ", in ", order_col, " order. Proportions are of each sample's ",
    "whole repertoire, computed before any selection, so a clonotype's height ",
    "does not depend on which other samples are compared. A clonotype absent ",
    "from a sample has no lode there."))
  report$p("")

  tables <- list(); panels <- list(); share_rows <- list()
  for (g in group_lv) {
    sub <- if (grouped) {
      clono[as.character(clono[[within]]) == g, , drop = FALSE]
    } else clono
    smp <- unique(data.frame(sample_id = as.character(sub$sample_id),
                             label = as.character(sub[[order_col]]),
                             stringsAsFactors = FALSE))
    if (anyNA(smp$label) || anyDuplicated(smp$label) ||
        anyDuplicated(smp$sample_id)) {
      smp <- data.frame(sample_id = unique(as.character(sub$sample_id)),
                        stringsAsFactors = FALSE)
      smp$label <- smp$sample_id
      lab_lv <- tracking_levels(sub$sample_id)
    } else {
      lab_lv <- tracking_levels(sub[[order_col]], configured)
    }
    if (nrow(smp) < 2L) {
      report$p(sprintf("%s: a single sample, nothing to compare.", g))
      next
    }
    sub$.label <- smp$label[match(as.character(sub$sample_id), smp$sample_id)]

    for (n in tr$compare_top) {
      tab <- clone_compare_table(sub, samples = lab_lv, top_n = n,
                                 abundance_col = abundance_col,
                                 sample_col = ".label")
      value <- if (tr$compare_proportion) "proportion" else "count"
      p <- plot_clone_compare(tab, value = value, unit = unit,
                              title = sprintf("%s - top %d", g, n))
      save_figure(p, sprintf("%sa_clone_compare_top%d_%s", prefix, n,
                             tracking_stem(g)), run_dir, cfg,
                  width = 12, height = 7,
                  legend = sprintf(paste0(
                    "Union of the %d most abundant clonotypes of each sample ",
                    "of %s, followed across samples. Height: %s of the ",
                    "sample's repertoire (%s). Flows join the same clonotype ",
                    "between consecutive samples; a clonotype absent from a ",
                    "sample has no lode there."), n, g,
                    if (value == "proportion") "proportion" else "count", unit))
      panels[[length(panels) + 1L]] <- list(group = g, top = n, plot = p)
      tab$group <- g
      tab$top_n <- n
      tables[[length(tables) + 1L]] <- tab
      held <- tapply(tab$proportion, tab$sample_id, sum)
      share_rows[[length(share_rows) + 1L]] <- data.frame(
        group = g, top_n = n, sample = names(held),
        n_clonotypes_followed = length(unique(tab$clonotype)),
        repertoire_held_pct = round(100 * as.numeric(held), 1),
        stringsAsFactors = FALSE)
    }

    for (pr in tr$compare_pairs) {
      if (!all(pr %in% smp$label)) next
      tab <- clone_compare_table(sub, samples = pr, top_n = tr$compare_pairs_top,
                                 abundance_col = abundance_col,
                                 sample_col = ".label")
      p <- plot_clone_compare(tab, value = "count", unit = unit,
                              title = sprintf("%s - %s vs %s, top %d", g, pr[1],
                                              pr[2], tr$compare_pairs_top))
      save_figure(p, sprintf("%sb_clone_compare_%s_vs_%s_%s", prefix,
                             tracking_stem(pr[1]), tracking_stem(pr[2]),
                             tracking_stem(g)), run_dir, cfg,
                  width = 7, height = 5,
                  legend = sprintf(paste0(
                    "The %d most abundant clonotypes of %s and of %s (%s), ",
                    "followed between the two samples in raw counts (%s)."),
                    tr$compare_pairs_top, pr[1], pr[2], g, unit))
      tab$group <- g
      tab$top_n <- tr$compare_pairs_top
      tables[[length(tables) + 1L]] <- tab
    }
  }

  if (!length(tables)) {
    report$p("No group had two or more samples; no comparison drawn.")
    return(invisible(NULL))
  }

  all_tab <- do.call(rbind, lapply(tables, function(t) {
    t$sample_id <- as.character(t$sample_id)
    t$clonotype <- as.character(t$clonotype)
    t
  }))
  names(all_tab)[names(all_tab) == "sample_id"] <- "sample_label"
  save_table(all_tab, paste0(prefix, "_clone_compare"), run_dir,
             "clonotypes followed across samples, with count, proportion and rank")

  shares <- do.call(rbind, share_rows)
  report$h2("Share of each sample held by the followed clonotypes")
  report$tbl(shares)

  # Combined figure: one column per group, one row per top-N value.
  if (length(panels) >= 2L && length(panels) <= 8L &&
      requireNamespace("patchwork", quietly = TRUE)) {
    p_group <- vapply(panels, function(x) as.character(x$group), character(1))
    p_top <- vapply(panels, function(x) as.numeric(x$top), numeric(1))
    n_groups <- length(unique(p_group))
    ord <- order(p_top, match(p_group, group_lv))
    plots <- lapply(panels[ord], `[[`, "plot")
    comp <- patchwork::wrap_plots(plots, ncol = n_groups) +
      patchwork::plot_annotation(
        title = "Clonal dynamics across samples (top clonotypes per sample)")
    n_rows <- ceiling(length(plots) / n_groups)
    save_figure(comp, paste0(prefix, "c_clone_compare_combined"), run_dir, cfg,
                width = 6 * n_groups, height = 4.5 * n_rows,
                legend = paste0("Combined view of the dominant-clonotype ",
                                "alluvials: one column per group, one row per ",
                                "number of clonotypes taken from each sample."))
  } else if (length(panels) >= 2L && !requireNamespace("patchwork", quietly = TRUE)) {
    report$p("(Package patchwork not installed: combined figure not drawn.)")
  }
  invisible(all_tab)
}

#' Clonotype tracking over time: baseline-dominant and expanded clonotypes
#'
#' Bulk pipeline, after differential abundance. Per subject: the
#' `baseline_top` (15) clonotypes dominant at baseline followed over the
#' subject's timepoints; and, when `calls` is given, the `expanded_top` (15)
#' significantly expanded clonotypes with the largest log2 fold change, de
#' novo and pre-existing expanded drawn separately, one panel each, with a
#' combined figure of all panels.
#'
#' @param clono Clonotype table with sample_id, subject_id, timepoint,
#'   clonotype, abundance and frequency.
#' @param calls Output of `test_clonal_abundance()`, or NULL.
#' @param run_dir Run directory from `init_run()`.
#' @param cfg Configuration list.
#' @param report Report writer from `open_report()`.
#' @param locus Locus label for the axis (e.g. "TRB"), or NULL.
#' @param unit Unit of the abundance ("molecules" in bulk data).
#' @param prefix File name prefix.
#' @return Invisibly, a list with the baseline and expanded tracking tables.
#' @export
section_clone_tracking <- function(clono, calls = NULL, run_dir, cfg, report,
                                   locus = NULL, unit = "molecules",
                                   prefix = "11") {
  need(c("dplyr", "ggalluvial"), "section_clone_tracking()")
  report$h1("Clonotype tracking over time")
  if (!"timepoint" %in% names(clono) || all(is.na(clono$timepoint))) {
    report$p("No timepoint in the data; tracking over time skipped.")
    return(invisible(NULL))
  }
  tr <- tracking_config(cfg)
  lv <- tracking_levels(clono$timepoint, cfg$cohort$timepoints)
  ref <- cfg$differential_abundance$reference_timepoint
  def <- cfg$repertoire$clonotype_definition %||% "strict"
  y_lab <- sprintf("%sfrequency (share of %s)",
                   if (length(locus)) paste0(paste(locus, collapse = "/"), " ")
                   else "", unit)

  # --- baseline-dominant clonotypes ---------------------------------------
  base_track <- track_baseline_clonotypes(clono, top_n = tr$baseline_top,
                                          baseline = ref, time_col = "timepoint",
                                          time_levels = lv)
  report$h2(sprintf("The %d clonotypes dominant at baseline", tr$baseline_top))
  report$p(paste0(
    "Each subject's most abundant baseline clonotypes, followed over the ",
    "subject's timepoints, with zero where a clonotype is not detected. Their ",
    "combined share of the repertoire at each timepoint is below."))
  held_rows <- list()
  for (subj in unique(base_track$subject_id)) {
    d <- base_track[base_track$subject_id == subj, , drop = FALSE]
    if (length(unique(d$sample_id)) < 2L) next
    p <- plot_clone_tracking(d, title = sprintf(paste0(
      "Tracking of the top %d baseline clonotypes - %s"), tr$baseline_top, subj),
      subtitle = sprintf("clonotype definition: %s", def), y_label = y_lab)
    save_figure(p, sprintf("%sa_baseline_tracking_%s", prefix,
                           tracking_stem(subj)), run_dir, cfg,
                width = 10, height = 6,
                legend = sprintf(paste0(
                  "Frequency of the %d most abundant baseline clonotypes of ",
                  "%s at each timepoint (share of %s; zero where not ",
                  "detected). Ribbons join the same clonotype over time."),
                  tr$baseline_top, subj, unit))
    held <- tapply(d$frequency, d$timepoint, sum)
    held <- held[!is.na(held)]
    held_rows[[length(held_rows) + 1L]] <- data.frame(
      subject_id = subj, timepoint = names(held),
      top_baseline_share_pct = round(100 * as.numeric(held), 2),
      stringsAsFactors = FALSE)
  }
  if (nrow(base_track)) {
    save_table(base_track, paste0(prefix, "a_baseline_tracking"), run_dir,
               "baseline-dominant clonotypes: frequency at every timepoint")
  }
  if (length(held_rows)) report$tbl(do.call(rbind, held_rows)) else
    report$p("No subject has two or more timepoints.")

  # --- significantly expanded clonotypes ----------------------------------
  report$h2("Significantly expanded clonotypes over time")
  exp_track <- NULL
  if (is.null(calls) || !any(calls$sig_expansion, na.rm = TRUE)) {
    report$p("No clonotype called significantly expanded; nothing to track.")
    return(invisible(list(baseline = base_track, expanded = NULL)))
  }
  sel <- select_tracked_expanded(calls, top_n = tr$expanded_top,
                                 rank_by = tr$expanded_rank_by, time_levels = lv)
  exp_track <- track_clonotypes(clono, sel, time_col = "timepoint",
                                time_levels = lv)
  report$p(sprintf(paste0(
    "Per subject and type (de novo, pre-existing expanded), the %d ",
    "significantly expanded clonotypes with the largest log2 fold change ",
    "(ranked at the %s), followed over every timepoint of the subject with ",
    "zero where absent. Ties at the boundary are kept."), tr$expanded_top,
    if (tr$expanded_rank_by == "first")
      "earliest timepoint at which each is significant" else
        "timepoint of each clonotype's largest fold change"))

  panels <- list()
  for (subj in tracking_levels(sel$subject_id, levels(calls$subject_id))) {
    for (type in c("de_novo", "pre_existing_expanded")) {
      d <- exp_track[exp_track$subject_id == subj & exp_track$type == type, ,
                     drop = FALSE]
      if (!nrow(d)) next
      type_lab <- tracking_type_label(type)
      solo <- plot_clone_tracking(d, title = sprintf(paste0(
        "Tracking of significantly expanded clonotypes (%s) - %s"), type_lab, subj),
        subtitle = sprintf(paste0("clonotype definition: %s; ribbons connect ",
                                  "the same clonotype over time"), def),
        y_label = y_lab)
      save_figure(solo, sprintf("%sb_expanded_tracking_%s_%s", prefix, type,
                                tracking_stem(subj)), run_dir, cfg,
                  width = 8, height = 6,
                  legend = sprintf(paste0(
                    "Frequency over time of the %s clonotypes of %s called ",
                    "significantly expanded against baseline with the largest ",
                    "log2 fold change (share of %s; zero where not detected)."),
                    type_lab, subj, unit))
      panels[[length(panels) + 1L]] <- plot_clone_tracking(
        d, title = sprintf("%s - %s", subj,
                           if (type == "de_novo") "de novo" else "pre-exist."),
        y_label = y_lab) +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 8),
                       axis.text = ggplot2::element_text(size = 6),
                       axis.title.y = ggplot2::element_text(size = 7))
    }
  }
  save_table(exp_track, paste0(prefix, "b_expanded_tracking"), run_dir,
             "tracked expanded clonotypes: frequency at every timepoint")
  counts <- as.data.frame(table(subject_id = sel$subject_id, type = sel$type),
                          stringsAsFactors = FALSE)
  names(counts)[names(counts) == "Freq"] <- "n_clonotypes_tracked"
  report$tbl(counts[counts$n_clonotypes_tracked > 0, , drop = FALSE])

  if (length(panels) && requireNamespace("patchwork", quietly = TRUE)) {
    comp <- patchwork::wrap_plots(panels, ncol = 2) +
      patchwork::plot_annotation(
        tag_levels = "A",
        title = "Tracking of significantly expanded clonotypes over time",
        theme = ggplot2::theme(plot.title = ggplot2::element_text(
          face = "bold", size = 12))) &
      ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 10))
    n_rows <- ceiling(length(panels) / 2)
    save_figure(comp, paste0(prefix, "c_expanded_tracking_combined"), run_dir,
                cfg, width = 12, height = max(4, 3 * n_rows),
                legend = paste0("All expanded-clonotype tracking panels, ",
                                "subject by type; panels without any expanded ",
                                "clonotype are omitted."))
  } else if (length(panels)) {
    report$p("(Package patchwork not installed: combined figure not drawn.)")
  }
  invisible(list(baseline = base_track, expanded = exp_track))
}

#' Clonotypes shared between compartments of the same subject
#'
#' Single-cell pipeline (sorted populations, tissues or culture conditions
#' of one subject). One alluvial per subject of the clonotypes found in at
#' least `repertoire$tracking$min_compartments` (2) compartments, with the
#' clonotypes of `highlight_compartment`, when configured, in orange.
#'
#' @param clono Clonotype table: subject_id, clonotype, abundance and the
#'   compartment column.
#' @param run_dir Run directory from `init_run()`.
#' @param cfg Configuration list.
#' @param report Report writer from `open_report()`.
#' @param compartment_col Column naming the compartment (default
#'   "condition").
#' @param clone_col Clonotype column; for B cells, a clone identifier from
#'   clonal inference shares lineages that differ by somatic mutation.
#' @param unit Unit of the abundance ("cells").
#' @param prefix File name prefix.
#' @return The flow table, invisibly.
#' @export
section_compartment_tracking <- function(clono, run_dir, cfg, report,
                                         compartment_col = "condition",
                                         clone_col = "clonotype",
                                         unit = "cells", prefix = "12") {
  need(c("dplyr", "ggalluvial"), "section_compartment_tracking()")
  report$h1("Clonotypes shared between compartments")
  if (!compartment_col %in% names(clono) || all(is.na(clono[[compartment_col]]))) {
    report$p(sprintf("No '%s' column; compartment tracking skipped.",
                     compartment_col))
    return(invisible(NULL))
  }
  tr <- tracking_config(cfg)
  configured <- if (identical(compartment_col, "condition"))
    cfg$cohort$conditions else NULL
  comps <- tr$compartments %||% tracking_levels(clono[[compartment_col]], configured)
  hl <- tr$highlight_compartment
  if (!is.null(hl) && !hl %in% as.character(clono[[compartment_col]])) {
    report$p(sprintf("Highlight compartment '%s' not found in the data; no highlight.",
                     hl))
    hl <- NULL
  }

  flows <- compartment_flows(clono, compartment_col = compartment_col,
                             compartments = comps, highlight = hl,
                             min_compartments = tr$min_compartments,
                             clone_col = clone_col)
  report$p(sprintf(paste0(
    "Clonotypes found in at least %d of the compartments (%s), per subject. ",
    "Flow height: %s of the clonotype in each compartment.%s"),
    tr$min_compartments, paste(comps, collapse = ", "), unit,
    if (is.null(hl)) "" else sprintf(
      " Orange: clonotypes the subject also has in '%s'.", hl)))
  if (!nrow(flows)) {
    report$p("No clonotype is shared between compartments.")
    return(invisible(flows))
  }

  hl_label <- if (is.null(hl)) NULL else sprintf("Found in %s", hl)
  summary_rows <- list()
  for (subj in unique(flows$subject_id)) {
    d <- flows[flows$subject_id == subj, , drop = FALSE]
    p <- plot_compartment_flows(
      d, unit = unit, highlight_label = hl_label,
      title = sprintf("Clone tracking across compartments - %s", subj))
    save_figure(p, sprintf("%s_compartment_flows_%s", prefix, tracking_stem(subj)),
                run_dir, cfg, width = 10, height = 5,
                legend = sprintf(paste0(
                  "Each ribbon is a clonotype of %s present in at least %d ",
                  "compartments; height in %s.%s"), subj, tr$min_compartments,
                  unit, if (is.null(hl)) "" else
                    sprintf(" Orange: clonotypes also found in %s.", hl)))
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      subject_id = subj,
      n_shared_clonotypes = length(unique(d$clonotype)),
      n_highlighted = length(unique(d$clonotype[d$highlighted])),
      abundance_in_shared = sum(d$abundance), stringsAsFactors = FALSE)
  }
  save_table(flows, paste0(prefix, "_compartment_flows"), run_dir,
             "clonotypes shared between compartments, abundance per compartment")
  report$tbl(do.call(rbind, summary_rows))
  invisible(flows)
}
