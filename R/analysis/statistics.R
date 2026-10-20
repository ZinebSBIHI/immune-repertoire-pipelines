# ---------------------------------------------------------------------------
# Paired comparisons at small sample size
#
# Immunology cohorts are often three to six subjects. A paired rank test on
# that many pairs has a floor on the p-value it can return, set by the number
# of distinct sign patterns and independent of the effect size:
#
#     n = 3   ->  smallest two-sided p = 0.25
#     n = 4   ->  0.125
#     n = 5   ->  0.0625
#     n = 6   ->  0.03125
#
# Below n = 6, no result can reach p < 0.05 however large and consistent the
# difference. Reporting "n.s." in that regime describes the test, not the
# data, and reporting an uncorrected p that happens to fall under 0.05 after
# testing several metrics describes neither.
#
# What carries evidence at this size is the consistency of the direction
# across subjects: a difference in the same direction in every subject is a
# statement about reproducibility that does not depend on a test the design
# cannot support. The functions here compute the attainable bound, return it
# alongside every p-value, and report direction concordance as a first-class
# result.
# ---------------------------------------------------------------------------

#' Smallest attainable two-sided p-value for a paired signed-rank test
#'
#' The most extreme outcome is all n differences sharing a sign, which has
#' probability 2^-n under the null; two-sided, 2 * 2^-n.
#'
#' @param n Number of pairs (subjects with both values).
#' @return The bound, or NA for n < 1.
attainable_p_paired <- function(n) {
  n <- as.integer(n)
  ifelse(n < 1L, NA_real_, pmin(1, 2^(-n + 1)))
}

#' Is a paired test capable of reaching the chosen alpha?
test_is_powered <- function(n, alpha = 0.05) {
  bound <- attainable_p_paired(n)
  !is.na(bound) & bound <= alpha
}

#' Paired comparison of two groups across subjects
#'
#' @param df Long data frame with one value per subject per group.
#' @param value_col Column of values.
#' @param group_col Column identifying the two groups.
#' @param subject_col Column identifying the subject (the pairing).
#' @param groups Length-2 character vector naming the groups to compare.
#' @param alpha Significance level, used only to report whether the design
#'   can reach it.
#' @return One-row data frame with the paired difference, the test result,
#'   the attainable bound and direction concordance.
paired_comparison <- function(df, value_col, group_col, subject_col,
                              groups, alpha = 0.05) {
  stopifnot(length(groups) == 2L)

  pair <- pair_by_subject(df, value_col, group_col, subject_col, groups)
  a <- pair$a; b <- pair$b
  n <- length(a)

  empty <- data.frame(
    group_a = groups[1], group_b = groups[2], n_pairs = n,
    mean_a = if (n) mean(a) else NA_real_,
    mean_b = if (n) mean(b) else NA_real_,
    mean_difference = NA_real_, median_difference = NA_real_,
    p_value = NA_real_, exact_test = NA, n_effective = NA_integer_,
    attainable_p = NA_real_, powered = FALSE,
    n_same_direction = NA_integer_, direction_concordant = FALSE,
    direction = NA_character_, stringsAsFactors = FALSE)
  if (n < 2L) return(empty)

  d <- b - a
  # A signed-rank test discards zero differences, so the sample size that
  # governs both the p-value and its floor is the number of NON-ZERO
  # differences, not the number of pairs. Using n_pairs for the bound would
  # make the bound too optimistic and could report `powered = TRUE` for a
  # test that cannot reach alpha.
  d_nz <- d[d != 0 & !is.na(d)]
  n_eff <- length(d_nz)
  if (n_eff < 1L) {
    empty$mean_difference <- mean(d); empty$median_difference <- stats::median(d)
    empty$n_effective <- 0L; empty$direction <- "no change"
    return(empty)
  }

  # The exact signed-rank distribution is only available without ties in
  # |d|. With ties, wilcox.test() silently falls back to the normal
  # approximation, whose p-value can fall BELOW the exact test's floor --
  # which would otherwise be printed next to that floor in the same
  # sentence, contradicting it. Which test ran is therefore recorded, and
  # the attainable bound is reported only for the exact test, to which it
  # applies.
  exact_ok <- !any(duplicated(abs(d_nz))) && n_eff < 50L
  p <- suppressWarnings(
    stats::wilcox.test(b, a, paired = TRUE, exact = exact_ok)$p.value)

  n_pos <- sum(d_nz > 0)
  n_neg <- sum(d_nz < 0)

  data.frame(
    group_a = groups[1],
    group_b = groups[2],
    n_pairs = n,
    mean_a = mean(a),
    mean_b = mean(b),
    mean_difference = mean(d),
    median_difference = stats::median(d),
    p_value = p,
    exact_test = exact_ok,
    n_effective = n_eff,
    attainable_p = if (exact_ok) attainable_p_paired(n_eff) else NA_real_,
    powered = exact_ok && test_is_powered(n_eff, alpha),
    n_same_direction = max(n_pos, n_neg),
    direction_concordant = (n_pos == n_eff) | (n_neg == n_eff),
    direction = if (n_pos == n_eff) "increase"
                else if (n_neg == n_eff) "decrease" else "mixed",
    stringsAsFactors = FALSE
  )
}

#' Align two groups' values by subject
#'
#' Built with `match()` rather than `reshape()`. `reshape()` prefixes the
#' value column onto the new names, and stripping that prefix with a greedy
#' regex destroys any group label containing a dot -- `Day.0`, `wk.2`, `V1.0`
#' are ordinary timepoint labels -- which yields a silent all-NA result
#' rather than an error.
#'
#' @return A list of `a`, `b` (paired, complete cases only) and `subjects`.
pair_by_subject <- function(df, value_col, group_col, subject_col, groups) {
  for (cl in c(value_col, group_col, subject_col)) {
    if (!cl %in% names(df)) {
      stop("Column '", cl, "' not found in the input.", call. = FALSE)
    }
  }
  grp <- as.character(df[[group_col]])
  sbj <- as.character(df[[subject_col]])
  val <- df[[value_col]]

  rows_a <- which(grp == groups[1])
  rows_b <- which(grp == groups[2])

  if (anyDuplicated(sbj[rows_a]) || anyDuplicated(sbj[rows_b])) {
    stop("More than one value per subject per group for '", value_col,
         "' in the comparison ", groups[1], " vs ", groups[2],
         ".\nAggregate to one value per subject per group before testing: ",
         "silently taking the first would report a result from an ",
         "arbitrary sample.", call. = FALSE)
  }

  subjects <- intersect(sbj[rows_a], sbj[rows_b])
  a <- val[rows_a][match(subjects, sbj[rows_a])]
  b <- val[rows_b][match(subjects, sbj[rows_b])]

  ok <- !is.na(a) & !is.na(b)
  list(a = a[ok], b = b[ok], subjects = subjects[ok])
}

#' All pairwise paired comparisons across the levels of a group
#'
#' @inheritParams paired_comparison
#' @param fdr_method Correction applied across the comparisons made here.
#'   Several metrics each tested across several arms multiplies quickly; the
#'   correction belongs in the same table as the raw p-values.
paired_comparison_all <- function(df, value_col, group_col, subject_col,
                                  alpha = 0.05, fdr_method = "BH") {
  need("dplyr", "paired_comparison_all()")
  lv <- if (is.factor(df[[group_col]])) levels(df[[group_col]])
        else sort(unique(as.character(df[[group_col]])))
  if (length(lv) < 2L) {
    stop("Need at least two groups to compare.", call. = FALSE)
  }

  combos <- utils::combn(lv, 2L, simplify = FALSE)
  out <- dplyr::bind_rows(lapply(combos, function(g) {
    paired_comparison(df, value_col, group_col, subject_col, g, alpha)
  }))

  out$fdr <- stats::p.adjust(out$p_value, method = fdr_method)
  out$metric <- value_col
  out
}

#' Format a paired result for a manuscript
#'
#' Produces the sentence the result actually supports. Where the design
#' cannot reach alpha, the p-value is reported with its bound and the
#' direction concordance is given as the primary evidence, rather than a
#' bare "not significant".
format_paired_result <- function(res) {
  vapply(seq_len(nrow(res)), function(i) {
    r <- res[i, ]
    if (is.na(r$p_value)) {
      return(sprintf("%s vs %s: insufficient paired observations (n = %d).",
                     r$group_a, r$group_b, r$n_pairs))
    }

    base <- sprintf("%s vs %s (n = %d pairs%s): mean difference %.3f",
                    r$group_a, r$group_b, r$n_pairs,
                    if (!is.na(r$n_effective) && r$n_effective < r$n_pairs)
                      sprintf(", %d non-tied", r$n_effective) else "",
                    r$mean_difference)

    concordance <- if (isTRUE(r$direction_concordant))
      sprintf(", concordant in all %d", r$n_effective)
    else sprintf(", concordant in %d of %d", r$n_same_direction, r$n_effective)

    if (!isTRUE(r$exact_test)) {
      # Without the exact distribution the floor does not apply, so it is
      # not quoted: printing a bound beside a p-value from a different test
      # would invite exactly the comparison it cannot support.
      sprintf(paste0("%s; paired Wilcoxon p = %.3g, normal approximation ",
                     "(tied absolute differences, so the exact test is ",
                     "unavailable and its floor does not apply). ",
                     "Direction: %s%s."),
              base, r$p_value, r$direction, concordance)
    } else if (isTRUE(r$powered)) {
      sprintf("%s; exact paired Wilcoxon p = %.3g%s.",
              base, r$p_value,
              if (!is.null(r$fdr) && !is.na(r$fdr))
                sprintf(" (FDR %.3g)", r$fdr) else "")
    } else {
      sprintf(paste0("%s; exact paired Wilcoxon p = %.3g, against a smallest ",
                     "attainable p of %.3g at this sample size, so the test ",
                     "cannot reach the nominal level. Direction: %s%s."),
              base, r$p_value, r$attainable_p, r$direction, concordance)
    }
  }, character(1))
}

#' Direction concordance across subjects
#'
#' @return One row per comparison: how many subjects move in each direction.
direction_concordance <- function(df, value_col, group_col, subject_col,
                                  groups) {
  pair <- pair_by_subject(df, value_col, group_col, subject_col, groups)
  d <- pair$b - pair$a

  data.frame(
    group_a = groups[1], group_b = groups[2],
    n_pairs = length(d),
    n_increase = sum(d > 0), n_decrease = sum(d < 0), n_tied = sum(d == 0),
    all_same_direction = length(d) > 0 &&
      (all(d > 0) || all(d < 0)),
    stringsAsFactors = FALSE
  )
}
