# ---------------------------------------------------------------------------
# Differential clonal abundance against a within-subject baseline
#
# Each post-baseline sample is compared with the same subject's baseline
# sample. Between-subject comparison is not attempted: repertoires are
# private, and the only meaningful reference for a subject's clonotype is
# that subject's own earlier sample.
#
# TWO AXES, KEPT SEPARATE
#
#   Axis 1, category (descriptive): was the clonotype detected at baseline?
#     - de_novo       not detected at baseline, detected after
#     - pre_existing  detected at baseline, higher after
#     - stable        detected at baseline, unchanged
#     - contracting   detected at baseline, lower after
#
#   Axis 2, significance (statistical): is the change larger than counting
#     noise? Fisher exact test on the 2x2 of clonotype count against the
#     rest of the repertoire, at both timepoints, with FDR control.
#
# These are different questions and the answers do not coincide. A clonotype
# can be de novo and non-significant (one molecule, absent before: the test
# has no power, and absence at baseline is exactly what a detection floor
# produces at random). It can be pre-existing and strongly significant.
#
# Collapsing the two into a single "induced" label is the standard error in
# this analysis: it attributes a causal origin on the strength of a
# descriptive comparison. The functions below return both axes and a
# combined status, and never a label asserting causation.
#
# The categorical scheme follows the logic of reverse-immunology clonotype
# tracking as described for individualised mRNA vaccines in
#   Sahin U, Schmidt M, Derhovanessian E, et al. Nature. 2026;651:1088-1096.
#   doi:10.1038/s41586-025-10004-2. PMID 41708868.
# In that work candidate selection combined bulk tracking with paired
# single-cell TRA/TRB and functional validation, and only two categories
# (de novo, pre-existing) were used. The additional categories here, the
# fold-change threshold, the Fisher/FDR axis and the "significantly
# expanded" label are this implementation's own and should not be
# attributed to that paper.
# ---------------------------------------------------------------------------

#' Categorise and test clonotypes against a within-subject baseline
#'
#' @param clono Output of `collapse_clonotypes()`, with a `timepoint` column.
#' @param reference_timepoint Baseline level of `timepoint`.
#' @param fold_change_threshold Minimum ratio of frequencies to call
#'   expansion or contraction.
#' @param detection_floor Minimum abundance for a clonotype to count as
#'   detected. A clonotype below the floor is "not detected", which is not
#'   the same claim as "absent".
#' @param fdr_threshold,fdr_method FDR control.
#' @return Data frame, one row per clonotype per post-baseline sample.
test_clonal_abundance <- function(clono,
                                  reference_timepoint,
                                  fold_change_threshold = 2,
                                  detection_floor = 1,
                                  fdr_threshold = 0.05,
                                  fdr_method = "BH") {
  need("dplyr", "test_clonal_abundance()")
  for (cl in c("sample_id", "subject_id", "timepoint", "clonotype",
               "abundance", "frequency")) {
    if (!cl %in% names(clono)) {
      stop("test_clonal_abundance() needs column '", cl, "'.", call. = FALSE)
    }
  }
  if (!reference_timepoint %in% clono$timepoint) {
    stop("Reference timepoint '", reference_timepoint,
         "' is not present in the data.", call. = FALSE)
  }

  depths <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(clono, .data$sample_id),
    depth = sum(.data$abundance), .groups = "drop"))

  results <- list()

  for (subj in unique(clono$subject_id)) {
    sub_db <- clono[clono$subject_id == subj, , drop = FALSE]

    base <- sub_db[sub_db$timepoint == reference_timepoint, , drop = FALSE]
    if (!nrow(base)) {
      warning("Subject '", subj, "' has no sample at the reference timepoint; ",
              "skipped. A subject without a baseline cannot contribute to a ",
              "within-subject comparison.", call. = FALSE)
      next
    }
    if (length(unique(base$sample_id)) > 1L) {
      stop("Subject '", subj, "' has more than one sample at the reference ",
           "timepoint. Resolve which is the baseline before testing.",
           call. = FALSE)
    }
    base_depth <- depths$depth[depths$sample_id == base$sample_id[1]]

    post_samples <- setdiff(unique(sub_db$sample_id), base$sample_id[1])

    for (smp in post_samples) {
      post <- sub_db[sub_db$sample_id == smp, , drop = FALSE]
      post_depth <- depths$depth[depths$sample_id == smp]

      all_clones <- union(base$clonotype, post$clonotype)

      b_ab <- base$abundance[match(all_clones, base$clonotype)]
      p_ab <- post$abundance[match(all_clones, post$clonotype)]
      b_ab[is.na(b_ab)] <- 0
      p_ab[is.na(p_ab)] <- 0

      b_freq <- b_ab / base_depth
      p_freq <- p_ab / post_depth

      detected_base <- b_ab >= detection_floor
      detected_post <- p_ab >= detection_floor

      # Fold change on frequencies. Where the clonotype is undetected at
      # baseline the ratio is undefined; the detection floor divided by
      # baseline depth is used as a conservative stand-in, which is the
      # smallest frequency the assay could have reported.
      b_freq_adj <- ifelse(detected_base, b_freq, detection_floor / base_depth)
      fold_change <- p_freq / b_freq_adj

      # Assigned in a fixed order, expansion before contraction, and the
      # two are mutually exclusive by construction for any threshold above
      # 1. At a threshold of exactly 1 every clonotype would satisfy both,
      # so that degenerate setting is rejected in validate_config() rather
      # than silently resolved by assignment order.
      category <- rep(NA_character_, length(all_clones))
      both <- detected_base & detected_post
      category[!detected_base & detected_post] <- "de_novo"
      category[both & fold_change >= fold_change_threshold] <-
        "pre_existing_expanded"
      category[both & is.na(category) &
                 fold_change <= 1 / fold_change_threshold] <- "contracting"
      category[both & is.na(category)] <- "stable"
      category[detected_base & !detected_post] <- "lost"
      # A clonotype present in neither sample above the detection floor can
      # arise whenever detection_floor exceeds the abundance floor used to
      # build the table. It is labelled rather than left as NA, which would
      # otherwise surface as the string "not_significant:NA".
      category[!detected_base & !detected_post] <- "below_detection"

      # Fisher exact test on the 2x2 table of this clonotype's molecules
      # against the rest of the repertoire, baseline versus post.
      p_value <- rep(NA_real_, length(all_clones))
      testable <- (p_ab + b_ab) > 0
      for (i in which(testable)) {
        tab <- matrix(c(p_ab[i], post_depth - p_ab[i],
                        b_ab[i], base_depth - b_ab[i]), nrow = 2)
        p_value[i] <- stats::fisher.test(tab, alternative = "two.sided")$p.value
      }
      fdr <- stats::p.adjust(p_value, method = fdr_method)

      results[[length(results) + 1L]] <- data.frame(
        subject_id = subj,
        sample_id = smp,
        timepoint = post$timepoint[1],
        clonotype = all_clones,
        abundance_baseline = b_ab,
        abundance_post = p_ab,
        frequency_baseline = b_freq,
        frequency_post = p_freq,
        fold_change = fold_change,
        log2_fold_change = log2(fold_change),
        category = category,
        p_value = p_value,
        fdr = fdr,
        sig_expansion = !is.na(fdr) & fdr < fdr_threshold &
          fold_change >= fold_change_threshold,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(results)) {
    stop("No subject had both a baseline and a post-baseline sample.",
         call. = FALSE)
  }

  out <- dplyr::bind_rows(results)

  # The combined status names the two axes side by side, never merged.
  out$status <- ifelse(
    out$sig_expansion,
    paste0("significantly_expanded:", out$category),
    paste0("not_significant:", out$category)
  )
  out
}

#' Sweep the abundance floor and report its effect on the calls
#'
#' The abundance floor is the most consequential free parameter in bulk
#' repertoire analysis, and it acts asymmetrically: de novo calls are defined
#' by absence at baseline, so single-molecule clonotypes are over-represented
#' among them and are removed first when the floor is raised.
#'
#' Run this before fixing `repertoire$min_abundance`, report the chosen value
#' in Methods, and show the sweep as a supplementary figure.
#'
#' @param db AIRR data frame.
#' @param abundance_col Abundance column name.
#' @param floors Numeric vector of floors to compare.
#' @param ... Passed to `test_clonal_abundance()`.
#' @return Data frame, one row per floor per sample.
sweep_abundance_floor <- function(db, abundance_col, floors = c(1, 2, 3, 5),
                                  definition = "strict",
                                  reference_timepoint, ...) {
  need("dplyr", "sweep_abundance_floor()")
  rows <- list()

  for (f in floors) {
    clono <- collapse_clonotypes(db, abundance_col = abundance_col,
                                 definition = definition, min_abundance = f)
    # Carry the timepoint column across the collapse.
    meta <- unique(db[, c("sample_id", "timepoint"), drop = FALSE])
    clono <- dplyr::left_join(clono, meta, by = "sample_id")

    retained <- dplyr::ungroup(dplyr::summarise(
      dplyr::group_by(clono, .data$sample_id, .data$subject_id,
                      .data$timepoint),
      n_clonotypes = dplyr::n(),
      total_abundance = sum(.data$abundance),
      .groups = "drop"))

    calls <- try(
      test_clonal_abundance(clono, reference_timepoint = reference_timepoint,
                            detection_floor = f, ...),
      silent = TRUE)

    if (inherits(calls, "try-error")) {
      called <- data.frame(sample_id = retained$sample_id,
                           n_de_novo = NA_integer_,
                           n_sig_expanded = NA_integer_,
                           stringsAsFactors = FALSE)
    } else {
      called <- dplyr::ungroup(dplyr::summarise(
        dplyr::group_by(calls, .data$sample_id),
        n_de_novo = sum(.data$category == "de_novo", na.rm = TRUE),
        n_sig_expanded = sum(.data$sig_expansion, na.rm = TRUE),
        .groups = "drop"))
    }

    rows[[length(rows) + 1L]] <- dplyr::mutate(
      dplyr::left_join(retained, called, by = "sample_id"),
      min_abundance = f)
  }

  out <- dplyr::bind_rows(rows)

  # Express everything relative to the most permissive floor, which is what
  # makes the cost of each threshold legible.
  ref <- out[out$min_abundance == min(floors),
             c("sample_id", "n_clonotypes", "total_abundance", "n_de_novo")]
  names(ref) <- c("sample_id", "ref_clonotypes", "ref_abundance", "ref_de_novo")

  out <- dplyr::left_join(out, ref, by = "sample_id")
  out$pct_clonotypes_retained <- 100 * out$n_clonotypes / out$ref_clonotypes
  out$pct_abundance_retained  <- 100 * out$total_abundance / out$ref_abundance
  out$pct_de_novo_retained    <- 100 * out$n_de_novo / out$ref_de_novo
  out
}

#' Summarise differential abundance per sample
#'
#' @param calls Output of `test_clonal_abundance()`.
#' @return One row per sample: counts by category and by significance, and
#'   the fraction of the repertoire the significantly expanded clonotypes
#'   occupy (magnitude).
summarise_abundance_calls <- function(calls) {
  need("dplyr", "summarise_abundance_calls()")
  dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(calls, .data$subject_id, .data$sample_id, .data$timepoint),
    n_tested = dplyr::n(),
    n_de_novo = sum(.data$category == "de_novo", na.rm = TRUE),
    n_pre_existing_expanded = sum(.data$category == "pre_existing_expanded",
                                  na.rm = TRUE),
    n_contracting = sum(.data$category == "contracting", na.rm = TRUE),
    n_sig_expanded = sum(.data$sig_expansion, na.rm = TRUE),
    n_sig_expanded_de_novo = sum(.data$sig_expansion &
                                   .data$category == "de_novo", na.rm = TRUE),
    n_sig_expanded_pre_existing = sum(.data$sig_expansion &
                                        .data$category == "pre_existing_expanded",
                                      na.rm = TRUE),
    magnitude_pct = 100 * sum(.data$frequency_post[.data$sig_expansion],
                              na.rm = TRUE),
    .groups = "drop"
  ))
}

#' Persistence of significantly expanded clonotypes across timepoints
#'
#' A response detected at a single timepoint and a response sustained across
#' several are different findings; reporting only the peak count conflates
#' them.
clonotype_persistence <- function(calls) {
  need("dplyr", "clonotype_persistence()")
  sig <- calls[calls$sig_expansion, , drop = FALSE]
  if (!nrow(sig)) {
    return(data.frame(subject_id = character(0), clonotype = character(0),
                      n_timepoints = integer(0), timepoints = character(0),
                      stringsAsFactors = FALSE))
  }
  dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(sig, .data$subject_id, .data$clonotype),
    n_timepoints = dplyr::n_distinct(.data$timepoint),
    timepoints = paste(sort(unique(as.character(.data$timepoint))),
                       collapse = ","),
    peak_log2fc = max(.data$log2_fold_change, na.rm = TRUE),
    .groups = "drop"
  ))
}
