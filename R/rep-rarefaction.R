# ---------------------------------------------------------------------------
# Rarefaction to a common sequencing depth
#
# Richness, clone counts and clonal concentration all scale with depth.
# Comparing them between samples of unequal depth measures the depth
# difference, so every such comparison is made at a common depth obtained
# by resampling without replacement.
# ---------------------------------------------------------------------------

#' Rarefy every sample to a common depth
#'
#' Multinomial resampling without replacement from each sample's abundance
#' vector, repeated `n_bootstrap` times, at a depth common to the samples
#' being compared.
#'
#' `within` names the grouping inside which the common depth is taken. For a
#' within-subject comparison of arms, the right common depth is the minimum
#' across that subject's samples, not the global minimum: using the global
#' minimum throws away depth from every subject to accommodate the shallowest
#' sample in the cohort.
#'
#' @param clono Output of `collapse_clonotypes()`.
#' @param depth Target depth. Default: the minimum depth within each group.
#' @param within Grouping column for the common depth, or NULL for global.
#' @param n_bootstrap Number of resamples.
#' @param statistics Named list of functions taking an abundance vector.
#' @return Data frame with the mean and quantiles of each statistic.
#' @export
rarefy_to_common_depth <- function(clono,
                                   depth = NULL,
                                   within = "subject_id",
                                   n_bootstrap = 200,
                                   statistics = list(
                                     n_clonotypes = function(a) sum(a > 0),
                                     top1_pct     = function(a) 100 * max(a) / sum(a),
                                     shannon      = shannon_entropy,
                                     inv_simpson  = inverse_simpson
                                   )) {
  need("dplyr", "rarefy_to_common_depth()")

  depths <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(clono, .data$sample_id),
    depth = sum(.data$abundance), .groups = "drop"))

  key <- if (is.null(within)) NULL else within
  meta <- unique(clono[, c("sample_id", intersect(key, names(clono))),
                       drop = FALSE])
  depths <- dplyr::left_join(depths, meta, by = "sample_id")

  if (is.null(depth)) {
    if (is.null(key)) {
      depths$target_depth <- min(depths$depth)
    } else {
      depths <- dplyr::ungroup(dplyr::mutate(
        dplyr::group_by(depths, dplyr::across(dplyr::all_of(key))),
        target_depth = min(.data$depth)))
    }
  } else {
    depths$target_depth <- depth
  }

  results <- vector("list", nrow(depths))

  for (i in seq_len(nrow(depths))) {
    s   <- depths$sample_id[i]
    tgt <- depths$target_depth[i]
    ab  <- clono$abundance[clono$sample_id == s]

    # Molecule counts are integers; a non-integer abundance (a normalised
    # or imputed value) cannot be expanded into a molecule pool, and
    # truncating it would make the guard below compare against a different
    # total than the pool actually holds.
    if (any(ab != floor(ab))) {
      stop("rarefy_to_common_depth() needs integer abundances (molecule or ",
           "cell counts); sample '", s, "' has fractional values. Rarefaction ",
           "resamples discrete observations and has no meaning on ",
           "normalised abundances.", call. = FALSE)
    }

    if (sum(ab) < tgt) {
      # Should not happen when the target is a within-group minimum, but a
      # user-supplied depth can exceed a sample's depth. Up-sampling would
      # invent molecules, so the sample is reported as NA.
      results[[i]] <- data.frame(sample_id = s, target_depth = tgt,
                                 statistic = names(statistics),
                                 mean = NA_real_, q025 = NA_real_,
                                 q975 = NA_real_,
                                 stringsAsFactors = FALSE)
      next
    }

    # Expand to a molecule vector once, then resample indices: for the
    # depths typical of bulk repertoire data this is both faster and more
    # obviously correct than repeated multivariate hypergeometric draws.
    pool <- rep.int(seq_along(ab), times = as.integer(ab))

    draws <- replicate(n_bootstrap, {
      idx <- sample.int(length(pool), size = as.integer(tgt), replace = FALSE)
      a <- tabulate(pool[idx], nbins = length(ab))
      vapply(statistics, function(f) f(a[a > 0]), numeric(1))
    })
    if (is.null(dim(draws))) draws <- matrix(draws, nrow = length(statistics))

    results[[i]] <- data.frame(
      sample_id    = s,
      target_depth = tgt,
      statistic    = names(statistics),
      mean = apply(draws, 1, mean, na.rm = TRUE),
      q025 = apply(draws, 1, stats::quantile, probs = 0.025, na.rm = TRUE),
      q975 = apply(draws, 1, stats::quantile, probs = 0.975, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  out <- dplyr::bind_rows(results)
  dplyr::left_join(out, meta, by = "sample_id")
}
