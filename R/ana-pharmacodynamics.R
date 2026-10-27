# ---------------------------------------------------------------------------
# Pharmacodynamic view: repertoire focusing and clonotype appearance
#
# Two complementary readouts, each subject compared with ITS OWN baseline:
#
#   (i)  FOCUSING. Does the repertoire concentrate into fewer clonotypes
#        after the intervention? Measured as the change from baseline in
#        Shannon clonality and in the Gini index.
#   (ii) APPEARANCE AND EXPANSION. How much of the repertoire is held by
#        clonotypes that were not there before, and how many of them are
#        there?
#
# The layout follows the pharmacodynamic figure of
#   Weber JS, Carlino MS, Khattak A, et al. Individualised neoantigen therapy mRNA-4157
#   (V940) plus pembrolizumab versus pembrolizumab monotherapy in resected
#   melanoma (KEYNOTE-942): a randomised, phase 2b study. Lancet.
#   2024;403(10427):632-644. doi:10.1016/S0140-6736(23)02268-7. PMID 38246194.
#
# WHY TWO FOCUSING METRICS. Shannon clonality is 1 - Pielou's evenness,
# which is what the commercial repertoire platforms report and therefore
# what the published figures use; because log(p) gives the rare clonotypes
# appreciable weight, it is sensitive to the long tail and so to sequencing
# depth and to the number of singletons. The Gini index is the Lorenz-curve
# inequality of the abundance distribution, which depends more on the
# overall shape and less on the ultra-rare clonotypes. Running both is the
# test of the conclusion: if the two tell the same story, it does not depend
# on the index; if they diverge, the long tail -- and hence depth -- is
# driving the result, which is why section_pharmacodynamics() also reports
# the correlation of each metric with depth.
#
# WHY TWO CLONOTYPE VIEWS. "Novel" counts only clonotypes absent from the
# baseline above the detection thresholds. "All expanded" adds the
# clonotypes already present at baseline whose frequency rose by at least
# the fold-change threshold. The first is the narrower claim, the second
# the fuller description of the response; the composition panel shows how
# much of the response is fresh and how much is recall.
#
# RELATION TO ana-differential.R. These calls are DETECTION-based: a
# clonotype qualifies on thresholds, with no test. test_clonal_abundance()
# instead tests every clonotype with a Fisher exact test and controls the
# FDR, and is the stricter criterion. The two are complementary and are
# deliberately kept apart; neither produces a label asserting that the
# intervention caused the change.
#
# DETECTION THRESHOLDS. A clonotype counts as present at baseline only at
# or above `min_count` molecules. That single definition is shared by both
# views, which is what keeps a clonotype's de novo status the same in both,
# and it has a second effect worth stating: a baseline singleton does not
# enter the baseline, so it is called de novo if it reappears, and no fold
# change is ever computed against a singleton's frequency, which would
# produce an explosive ratio out of one molecule.
#   Shugay M, et al. Nat Methods. 2014;11:653-655 (MIGEC).
#   Egorov ES, et al. J Immunol. 2015;194:6155-6163.
#   Greiff V, et al. Genome Med. 2015;7:49.
# ---------------------------------------------------------------------------

#' Repertoire focusing metrics per sample
#'
#' Both focusing metrics are computed on the same collapsed table, so they
#' are strictly comparable with each other and with every other number in
#' the run. Which of them is plotted is a separate decision.
#'
#' @param clono Output of [collapse_clonotypes()], with a `timepoint` column.
#' @return Data frame: `subject_id`, `timepoint`, `sample_id`,
#'   `total_abundance`, `richness`, `shannon`, `shannon_clonality`, `gini`.
#' @export
focalisation_metrics <- function(clono) {
  need("dplyr", "focalisation_metrics()")
  for (cl in c("sample_id", "subject_id", "abundance")) {
    if (!cl %in% names(clono)) {
      stop("focalisation_metrics() needs column '", cl, "'.", call. = FALSE)
    }
  }
  group_cols <- intersect(c("subject_id", "timepoint", "sample_id"),
                          names(clono))

  dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(clono, dplyr::across(dplyr::all_of(group_cols))),
    total_abundance   = sum(.data$abundance),
    richness          = dplyr::n_distinct(.data$clonotype),
    shannon           = shannon_entropy(.data$abundance),
    # 1 - H/log(richness). NA at a richness of one, where log(richness) is
    # zero: a single-clonotype sample has no evenness to normalise by, and
    # returning a number there would put a division by zero into a figure.
    shannon_clonality = shannon_clonality(.data$abundance),
    gini              = gini_index(.data$abundance),
    .groups = "drop"))
}

#' Change in each focusing metric from the subject's own baseline
#'
#' Comparing each subject with their own baseline neutralises the
#' between-subject differences in the starting repertoire, which are large:
#' receptor repertoires are largely private, so a cohort-level baseline
#' would be a mixture of individual repertoires rather than a reference for
#' any of them.
#'
#' @param focal Output of [focalisation_metrics()].
#' @param reference_timepoint Baseline level of `timepoint`.
#' @param metrics Metric columns to follow.
#' @return Long data frame: `subject_id`, `timepoint`, `sample_id`, `metric`,
#'   `value`, `baseline_value`, `delta`.
#' @export
focalisation_change <- function(focal, reference_timepoint,
                                metrics = c("shannon_clonality", "gini")) {
  need(c("dplyr", "tidyr"), "focalisation_change()")
  if (!"timepoint" %in% names(focal)) {
    stop("focalisation_change() needs a 'timepoint' column.", call. = FALSE)
  }
  metrics <- intersect(metrics, names(focal))
  if (!length(metrics)) {
    stop("None of the requested metrics is present in the input.",
         call. = FALSE)
  }

  long <- tidyr::pivot_longer(
    focal[, c("subject_id", "timepoint", "sample_id", metrics)],
    cols = dplyr::all_of(metrics), names_to = "metric", values_to = "value")

  base <- long[as.character(long$timepoint) == reference_timepoint,
               c("subject_id", "metric", "value"), drop = FALSE]
  names(base)[names(base) == "value"] <- "baseline_value"

  out <- dplyr::left_join(long, base, by = c("subject_id", "metric"))
  out$delta <- out$value - out$baseline_value
  out
}

#' Baseline reference shared by both clonotype views
#'
#' A clonotype counts as present at baseline when it is detected there at or
#' above `min_count`. Defining this once and using it for both views is what
#' makes a clonotype's de novo status identical in both.
#'
#' @param clono Output of [collapse_clonotypes()], with `timepoint`.
#' @param reference_timepoint Baseline level of `timepoint`.
#' @param min_count Minimum abundance for a baseline detection to count.
#' @return Data frame: `subject_id`, `clonotype`, `frequency_baseline`,
#'   `abundance_baseline`.
#' @export
baseline_reference <- function(clono, reference_timepoint, min_count = 2) {
  need("dplyr", "baseline_reference()")
  base <- clono[as.character(clono$timepoint) == reference_timepoint &
                  clono$abundance >= min_count, , drop = FALSE]
  if (!nrow(base)) {
    warning("No clonotype reaches min_count = ", min_count,
            " at the reference timepoint; every post-baseline clonotype ",
            "will be called de novo.", call. = FALSE)
  }
  out <- unique(data.frame(
    subject_id = as.character(base$subject_id),
    clonotype = as.character(base$clonotype),
    frequency_baseline = base$frequency,
    abundance_baseline = base$abundance,
    stringsAsFactors = FALSE))
  out[!duplicated(out[, c("subject_id", "clonotype")]), , drop = FALSE]
}

#' Clonotypes absent from the baseline ("novel" view)
#'
#' There is no de novo column in the data: de novo is the ABSENCE of the
#' clonotype from the subject's baseline. The join below therefore gives a
#' missing baseline to exactly the clonotypes never seen at baseline, and
#' that missingness is what the flag records. `min_count` and `min_freq`
#' only exclude noise; they are not what defines de novo.
#'
#' @param clono Output of [collapse_clonotypes()], with `timepoint`.
#' @param baseline Output of [baseline_reference()].
#' @param reference_timepoint Baseline level of `timepoint`.
#' @param min_count,min_freq Detection thresholds applied post-baseline.
#' @return `clono` with `present_baseline` and `novel_expanded` columns.
#' @export
classify_novel_clonotypes <- function(clono, baseline, reference_timepoint,
                                      min_count = 2, min_freq = 1e-5) {
  need("dplyr", "classify_novel_clonotypes()")
  key <- unique(baseline[, c("subject_id", "clonotype")])
  key$present_baseline <- TRUE

  out <- dplyr::left_join(
    dplyr::mutate(clono, subject_id = as.character(.data$subject_id),
                  clonotype = as.character(.data$clonotype)),
    key, by = c("subject_id", "clonotype"))
  out$present_baseline[is.na(out$present_baseline)] <- FALSE
  out$novel_expanded <- as.character(out$timepoint) != reference_timepoint &
    !out$present_baseline &
    out$abundance >= min_count &
    out$frequency >= min_freq
  out
}

#' Clonotypes expanded or newly detected ("all expanded" view)
#'
#' A clonotype already present at baseline cannot be new; it is called
#' expanded when its frequency has risen by at least `min_fold`. Clonotypes
#' absent from the baseline are expanded by construction.
#'
#' @param clono Output of [collapse_clonotypes()], with `timepoint`.
#' @param baseline Output of [baseline_reference()].
#' @param reference_timepoint Baseline level of `timepoint`.
#' @param min_count,min_freq Detection thresholds applied post-baseline.
#' @param min_fold Fold change in frequency above which a pre-existing
#'   clonotype counts as expanded.
#' @return `clono` with `present_baseline`, `fold_change`, `clonotype_type`
#'   and `is_expanded` columns.
#' @export
classify_expanded_clonotypes <- function(clono, baseline, reference_timepoint,
                                         min_count = 2, min_freq = 1e-5,
                                         min_fold = 2) {
  need("dplyr", "classify_expanded_clonotypes()")
  base <- baseline[, c("subject_id", "clonotype", "frequency_baseline")]

  out <- dplyr::left_join(
    dplyr::mutate(clono, subject_id = as.character(.data$subject_id),
                  clonotype = as.character(.data$clonotype)),
    base, by = c("subject_id", "clonotype"))

  out$present_baseline <- !is.na(out$frequency_baseline)
  # NA for the clonotypes absent at baseline, by design: there is no ratio
  # to form, and substituting a floor here would turn an absence into a
  # finite fold change and mix the two categories.
  out$fold_change <- out$frequency / out$frequency_baseline

  post <- as.character(out$timepoint) != reference_timepoint
  detected <- out$abundance >= min_count & out$frequency >= min_freq
  out$clonotype_type <- ifelse(
    post & !out$present_baseline & detected, "de novo",
    ifelse(post & out$present_baseline & detected &
             !is.na(out$fold_change) & out$fold_change >= min_fold,
           "pre-existing (expanded)", "not expanded"))
  out$is_expanded <- out$clonotype_type %in% c("de novo",
                                               "pre-existing (expanded)")
  out
}

#' Weight and number of the novel clonotypes per sample
#'
#' Two quantities, not one: the summed frequency is the magnitude of the
#' response (how much of the repertoire it occupies) and the count is its
#' breadth (how many clonotypes). A response can be large in one and small
#' in the other, and reporting only one of them hides that.
#'
#' @param novel Output of [classify_novel_clonotypes()].
#' @return Data frame: `subject_id`, `timepoint`, `sample_id`,
#'   `summed_frequency_novel`, `n_novel`.
#' @export
summarise_novel_clonotypes <- function(novel) {
  need("dplyr", "summarise_novel_clonotypes()")
  dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(novel, .data$subject_id, .data$timepoint,
                    .data$sample_id),
    summed_frequency_novel = sum(.data$frequency[.data$novel_expanded],
                                 na.rm = TRUE),
    n_novel = dplyr::n_distinct(.data$clonotype[.data$novel_expanded]),
    .groups = "drop"))
}

#' Weight, number and composition of the expanded clonotypes per sample
#'
#' @param expanded Output of [classify_expanded_clonotypes()].
#' @return List with `summary` (one row per sample) and `by_type` (one row
#'   per sample per clonotype type), for the composition panel.
#' @export
summarise_expanded_clonotypes <- function(expanded) {
  need("dplyr", "summarise_expanded_clonotypes()")
  summary <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(expanded, .data$subject_id, .data$timepoint,
                    .data$sample_id),
    summed_frequency_expanded = sum(.data$frequency[.data$is_expanded],
                                    na.rm = TRUE),
    n_expanded = dplyr::n_distinct(.data$clonotype[.data$is_expanded]),
    n_de_novo = dplyr::n_distinct(
      .data$clonotype[.data$clonotype_type == "de novo"]),
    n_pre_existing = dplyr::n_distinct(
      .data$clonotype[.data$clonotype_type == "pre-existing (expanded)"]),
    .groups = "drop"))

  by_type <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(expanded[expanded$is_expanded, , drop = FALSE],
                    .data$subject_id, .data$timepoint,
                    .data$clonotype_type),
    n = dplyr::n_distinct(.data$clonotype), .groups = "drop"))

  list(summary = summary, by_type = by_type)
}

#' Correlation of the focusing metrics with sequencing depth
#'
#' The control that decides how the focusing panel may be read. Shannon
#' clonality is not independent of depth, because richness there is OBSERVED
#' richness; if depth varies strongly between a subject's samples, part of
#' the change can be technical. The Gini index is less exposed, so comparing
#' the two correlations is the diagnostic. The agreement between the two
#' metrics on the changes themselves is reported alongside: a high
#' correlation there means the conclusion does not depend on the index.
#'
#' @param focal Output of [focalisation_metrics()].
#' @param change Output of [focalisation_change()].
#' @param reference_timepoint Baseline level, excluded from the agreement
#'   correlation because its change is zero by construction.
#' @return One-row data frame of Spearman coefficients.
#' @export
focalisation_depth_correlations <- function(focal, change,
                                            reference_timepoint = NULL) {
  need(c("dplyr", "tidyr"), "focalisation_depth_correlations()")
  rho <- function(a, b) {
    ok <- stats::complete.cases(a, b)
    if (sum(ok) < 3L) return(NA_real_)
    suppressWarnings(stats::cor(a[ok], b[ok], method = "spearman"))
  }

  wide <- tidyr::pivot_wider(
    change[, c("subject_id", "timepoint", "metric", "delta")],
    names_from = "metric", values_from = "delta")
  if (!is.null(reference_timepoint)) {
    wide <- wide[as.character(wide$timepoint) != reference_timepoint, ,
                 drop = FALSE]
  }

  data.frame(
    rho_depth_shannon_clonality =
      if ("shannon_clonality" %in% names(focal))
        rho(focal$total_abundance, focal$shannon_clonality) else NA_real_,
    rho_depth_gini =
      if ("gini" %in% names(focal))
        rho(focal$total_abundance, focal$gini) else NA_real_,
    rho_delta_shannon_vs_gini =
      if (all(c("shannon_clonality", "gini") %in% names(wide)))
        rho(wide$shannon_clonality, wide$gini) else NA_real_,
    n_samples = nrow(focal),
    stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Section
# ---------------------------------------------------------------------------

#' Pharmacodynamic section: focusing and clonotype appearance over time
#'
#' Panel A is the change in focusing from baseline, drawn for each
#' configured metric. Panels B and C are the magnitude and the breadth of
#' the response, once for the novel clonotypes alone and once for all
#' expanded clonotypes, and panel D is the composition of the latter. The
#' depth control table and the correlations are written to the report
#' because panel A cannot responsibly be read without them.
#'
#' @param clono Output of [collapse_clonotypes()], with a `timepoint` column.
#' @param run_dir,cfg,report Run directory, configuration and report writer.
#' @param unit Unit of observation, used in labels.
#' @param prefix File-name prefix.
#' @return List with `focal`, `change`, `novel` and `expanded` summaries,
#'   invisibly.
#' @export
section_pharmacodynamics <- function(clono, run_dir, cfg, report,
                                     unit = "molecules", prefix = "08b") {
  need(c("dplyr", "tidyr"), "section_pharmacodynamics()")
  report$h1("Pharmacodynamics: repertoire focusing and clonotype appearance")

  pd <- cfg$differential_abundance$pharmacodynamics %||% list()
  focal_metric <- pd$focal_metric %||% "both"
  if (!focal_metric %in% c("shannon", "gini", "both")) {
    stop("differential_abundance$pharmacodynamics$focal_metric must be ",
         "'shannon', 'gini' or 'both' (got '", focal_metric, "').",
         call. = FALSE)
  }
  min_count <- as.numeric(pd$min_count %||% 2)
  min_freq  <- as.numeric(pd$min_freq %||% 1e-5)
  # Aligned with the project-wide expansion threshold unless the
  # configuration overrides it: two different fold-change thresholds in one
  # project would call different clonotypes expanded in two figures.
  min_fold  <- as.numeric(pd$min_fold %||%
                            cfg$differential_abundance$fold_change_threshold)
  ref_tp <- cfg$differential_abundance$reference_timepoint

  if (!"timepoint" %in% names(clono) || is.null(ref_tp) ||
      !ref_tp %in% as.character(clono$timepoint)) {
    report$p(sprintf(paste0(
      "No baseline timepoint ('%s') in the data. Every quantity in this ",
      "section is a change from the subject's own baseline, so the section ",
      "is skipped."), ref_tp %||% "(unset)"))
    return(invisible(NULL))
  }

  report$p(sprintf(paste0(
    "Detection thresholds: at least %g %s and a frequency of at least %g; ",
    "a pre-existing clonotype counts as expanded from a %g-fold rise in ",
    "frequency. These are detection criteria, not a test: the Fisher and ",
    "FDR view of the same data is in the differential abundance section, ",
    "and the two are complementary."),
    min_count, unit, min_freq, min_fold))
  report$p("")

  # --- focusing ------------------------------------------------------------
  focal <- focalisation_metrics(clono)
  change <- focalisation_change(focal, reference_timepoint = ref_tp)
  save_table(focal, paste0(prefix, "_focalisation_metrics"), run_dir,
             sprintf("depth, richness, Shannon clonality and Gini per sample (%s)",
                     unit))
  save_table(change, paste0(prefix, "_focalisation_change"), run_dir,
             "change in each focusing metric from the subject's own baseline")

  metric_labels <- c(
    shannon_clonality = "Change in Shannon clonality\nfrom baseline",
    gini = "Change in clonal inequality\n(Gini) from baseline")
  metric_subtitles <- c(
    shannon_clonality = paste("Normalised Shannon clonality (1 - Pielou J);",
                              "sensitive to the long tail and so to depth"),
    gini = paste("Gini index of clonal inequality; less exposed to the",
                 "ultra-rare clonotypes"))
  metrics_used <- switch(focal_metric,
                         "shannon" = "shannon_clonality",
                         "gini"    = "gini",
                         "both"    = c("shannon_clonality", "gini"))

  panel_a <- list()
  for (m in metrics_used) {
    dat <- change[change$metric == m, , drop = FALSE]
    panel_a[[m]] <- ggplot2::ggplot(
      dat, ggplot2::aes(.data$timepoint, .data$delta,
                        group = .data$subject_id,
                        colour = .data$subject_id)) +
      # The dashed zero line is the reference the panel is read against:
      # above it the repertoire has focused relative to that subject's own
      # baseline, below it has diversified.
      ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
      ggplot2::geom_line(linewidth = 0.7, alpha = 0.85) +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(x = NULL, y = metric_labels[[m]],
                    subtitle = metric_subtitles[[m]], colour = "Subject") +
      theme_pub(10)
    save_figure(panel_a[[m]], sprintf("%sa_focalisation_%s", prefix, m),
                run_dir, cfg, width = 8, height = 4,
                legend = sprintf(paste0("Change from each subject's own ",
                                        "baseline in %s. The dashed line is ",
                                        "no change."), m))
  }

  if (identical(focal_metric, "both") &&
      requireNamespace("patchwork", quietly = TRUE)) {
    p_cmp <- patchwork::wrap_plots(panel_a, nrow = 1) +
      patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = "bottom")
    save_figure(p_cmp, paste0(prefix, "a_focalisation_shannon_vs_gini"),
                run_dir, cfg, width = 12, height = 4.5,
                legend = paste("The two focusing metrics side by side. Where",
                               "they agree, the conclusion does not depend",
                               "on the index; where they diverge, the long",
                               "tail is driving it."))
  }

  # --- novel view ----------------------------------------------------------
  baseline <- baseline_reference(clono, reference_timepoint = ref_tp,
                                 min_count = min_count)
  novel <- classify_novel_clonotypes(clono, baseline,
                                     reference_timepoint = ref_tp,
                                     min_count = min_count,
                                     min_freq = min_freq)
  novel_summary <- summarise_novel_clonotypes(novel)
  save_table(novel_summary, paste0(prefix, "_novel_clonotypes"), run_dir,
             "weight and number of clonotypes absent from the baseline")

  report$h2("Clonotypes absent from the baseline")
  report$tbl(as.data.frame(novel_summary), digits = 4)
  report$p("")
  report$p(paste0(
    "De novo here means absent from the subject's own baseline above the ",
    "detection thresholds; no test is involved. The summed frequency is the ",
    "magnitude of the response, the count its breadth, and the two can move ",
    "in opposite directions."))

  pb <- trajectory_panel(novel_summary, "summed_frequency_novel",
                         "Summed frequency of\nnovel clonotypes")
  pc <- trajectory_panel(novel_summary, "n_novel",
                         "Novel clonotypes (n)")
  save_figure(pb, paste0(prefix, "b_novel_frequency"), run_dir, cfg,
              width = 8, height = 4,
              legend = paste("Fraction of the repertoire held by clonotypes",
                             "absent from the subject's own baseline."))
  save_figure(pc, paste0(prefix, "c_novel_count"), run_dir, cfg,
              width = 8, height = 4,
              legend = paste("Number of distinct clonotypes absent from the",
                             "subject's own baseline."))

  # --- all expanded view ---------------------------------------------------
  expanded <- classify_expanded_clonotypes(
    clono, baseline, reference_timepoint = ref_tp, min_count = min_count,
    min_freq = min_freq, min_fold = min_fold)
  exp_summary <- summarise_expanded_clonotypes(expanded)
  save_table(exp_summary$summary, paste0(prefix, "_expanded_clonotypes"),
             run_dir,
             "weight and number of clonotypes newly detected or expanded")
  save_table(exp_summary$by_type,
             paste0(prefix, "_expanded_clonotype_composition"), run_dir,
             "expanded clonotypes split into de novo and pre-existing")

  report$h2("Clonotypes newly detected or expanded")
  report$tbl(as.data.frame(exp_summary$summary), digits = 4)

  pb_all <- trajectory_panel(exp_summary$summary,
                             "summed_frequency_expanded",
                             "Summed frequency of\nall expanded clonotypes")
  pc_all <- trajectory_panel(exp_summary$summary, "n_expanded",
                             "All expanded\nclonotypes (n)")
  save_figure(pb_all, paste0(prefix, "d_expanded_frequency"), run_dir, cfg,
              width = 8, height = 4,
              legend = paste("Fraction of the repertoire held by clonotypes",
                             "newly detected or expanded against the",
                             "subject's own baseline."))
  save_figure(pc_all, paste0(prefix, "e_expanded_count"), run_dir, cfg,
              width = 8, height = 4,
              legend = paste("Number of distinct clonotypes newly detected",
                             "or expanded."))

  if (nrow(exp_summary$by_type)) {
    type_cols <- palette_categorical(2)
    names(type_cols) <- c("de novo", "pre-existing (expanded)")
    p_comp <- ggplot2::ggplot(exp_summary$by_type,
                              ggplot2::aes(.data$timepoint, .data$n,
                                           fill = .data$clonotype_type)) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::facet_wrap(~subject_id, nrow = 1) +
      ggplot2::scale_fill_manual(values = type_cols, name = NULL) +
      ggplot2::labs(x = NULL, y = "Expanded clonotypes (n)",
                    title = "Composition of the expanded repertoire",
                    subtitle = paste("How much of the response is new and",
                                     "how much is pre-existing")) +
      theme_pub(9) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                         hjust = 1, size = 7),
                     legend.position = "bottom")
    save_figure(p_comp, paste0(prefix, "f_expanded_composition"), run_dir,
                cfg, width = 11, height = 4,
                legend = paste("Expanded clonotypes per timepoint split into",
                               "those absent from the baseline and those",
                               "already present and expanded."))
  }

  # --- depth control -------------------------------------------------------
  control <- focal[, intersect(c("subject_id", "timepoint", "sample_id",
                                 "total_abundance", "richness", "shannon",
                                 "shannon_clonality", "gini"),
                               names(focal)), drop = FALSE]
  control <- control[order(control$subject_id, control$timepoint), ,
                     drop = FALSE]
  save_table(control, paste0(prefix, "_depth_control"), run_dir,
             sprintf("depth, richness and both focusing metrics per sample (%s)",
                     unit))

  cors <- focalisation_depth_correlations(focal, change,
                                          reference_timepoint = ref_tp)
  save_table(cors, paste0(prefix, "_focalisation_correlations"), run_dir,
             "Spearman correlations of the focusing metrics with depth and with each other")

  report$h2("Depth control")
  report$tbl(as.data.frame(control), digits = 3)
  report$p("")
  report$p(sprintf(paste0(
    "Spearman correlation with depth: %.2f for Shannon clonality, %.2f for ",
    "Gini. Agreement between the two metrics on the changes: %.2f."),
    cors$rho_depth_shannon_clonality, cors$rho_depth_gini,
    cors$rho_delta_shannon_vs_gini))
  report$p("")
  report$p(paste0(
    "Read this before the focusing panel. Shannon clonality normalises by ",
    "OBSERVED richness, so it is not independent of depth: where a ",
    "subject's depth varies between timepoints, part of the change can be ",
    "technical. A high correlation with depth for Shannon clonality but not ",
    "for Gini means the Shannon trajectory is partly a depth artefact; a ",
    "high agreement between the two metrics means the conclusion does not ",
    "depend on which was chosen."))

  invisible(list(focal = focal, change = change, novel = novel_summary,
                 expanded = exp_summary))
}

#' One trajectory panel: a quantity against time, one line per subject
#'
#' Shared by every panel of the pharmacodynamic figure so that the panels
#' are visually identical and only the quantity differs.
#' @noRd
trajectory_panel <- function(df, value_col, y_label) {
  ggplot2::ggplot(df, ggplot2::aes(.data$timepoint, .data[[value_col]],
                                   group = .data$subject_id,
                                   colour = .data$subject_id)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
    ggplot2::geom_line(linewidth = 0.7, alpha = 0.85) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(x = NULL, y = y_label, colour = "Subject") +
    theme_pub(10)
}
