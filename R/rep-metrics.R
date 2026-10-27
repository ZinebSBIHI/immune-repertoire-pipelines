# ---------------------------------------------------------------------------
# Repertoire structure, diversity and overlap
#
# The indices are written out rather than taken from a package, so that the
# exact definition behind each number in a manuscript is visible and
# citable. All take a vector of clonotype abundances.
#
# SEQUENCING DEPTH. Richness, clone counts and clonal concentration are all
# monotone in depth, so any comparison of them between samples of unequal
# depth measures the depth difference. Such comparisons go through
# rarefy_to_common_depth() in rarefaction.R.
# ---------------------------------------------------------------------------

#' Per-sample repertoire summary
#'
#' @param clono Output of `collapse_clonotypes()`.
#' @return Data frame of depth, richness and clonality metrics per sample.
#' @export
repertoire_summary <- function(clono) {
  need("dplyr", "repertoire_summary()")
  dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(clono, .data$sample_id, .data$subject_id),
    depth         = sum(.data$abundance),
    n_clonotypes  = dplyr::n(),
    n_singletons  = sum(.data$abundance == 1),
    pct_singleton = 100 * sum(.data$abundance == 1) / dplyr::n(),
    top1_pct      = 100 * max(.data$frequency),
    top10_pct     = 100 * sum(sort(.data$frequency, decreasing = TRUE)[
      seq_len(min(10L, dplyr::n()))]),
    shannon       = shannon_entropy(.data$abundance),
    clonality     = shannon_clonality(.data$abundance),
    inv_simpson   = inverse_simpson(.data$abundance),
    gini          = gini_index(.data$abundance),
    chao1         = chao1_richness(.data$abundance),
    d50           = d50_index(.data$abundance),
    .groups = "drop"
  ))
}

# --- Diversity and clonality indices ---------------------------------------
#
# Written out rather than taken from a package so that the exact definition
# behind each number in a manuscript is visible and citable. All take a
# vector of clonotype abundances.

#' Shannon entropy (natural log)
#' @param abundance Abundance vector of one sample; zeros are ignored.
#' @export
shannon_entropy <- function(abundance) {
  p <- abundance[abundance > 0]
  if (!length(p)) return(NA_real_)
  p <- p / sum(p)
  -sum(p * log(p))
}

#' Shannon clonality
#'
#' 1 - normalised Shannon entropy: 0 for a perfectly even repertoire,
#' approaching 1 as the repertoire becomes dominated by one clonotype.
#' Normalising by log(richness) makes it less depth-dependent than raw
#' entropy, but it is not depth-free: compare at common depth.
#' @param abundance Abundance vector of one sample; zeros are ignored.
#' @export
shannon_clonality <- function(abundance) {
  p <- abundance[abundance > 0]
  if (length(p) < 2L) return(NA_real_)
  1 - shannon_entropy(p) / log(length(p))
}

#' Inverse Simpson index
#' @param abundance Abundance vector of one sample; zeros are ignored.
#' @export
inverse_simpson <- function(abundance) {
  p <- abundance[abundance > 0]
  if (!length(p)) return(NA_real_)
  p <- p / sum(p)
  1 / sum(p^2)
}

#' Gini index of abundance inequality
#' @param abundance Abundance vector of one sample; zeros are ignored.
#' @export
gini_index <- function(abundance) {
  x <- sort(abundance[abundance > 0])
  n <- length(x)
  if (n < 2L) return(NA_real_)
  sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x))
}

#' Chao1 richness estimator
#'
#' Chao A. Nonparametric estimation of the number of classes in a
#' population. Scand J Stat. 1984;11:265-270.
#'
#' Chao1 extrapolates unseen classes from the singletons and doubletons, so
#' it is only meaningful on UNFILTERED abundances. If an abundance floor has
#' removed the singletons, f1 is zero by construction and the estimate
#' collapses onto the observed richness, which looks like a result and is
#' an artefact of the filter. The guard returns NA in that case rather than
#' a number that invites interpretation.
#' @param abundance Abundance vector of one sample; zeros are ignored.
#' @export
chao1_richness <- function(abundance) {
  obs <- sum(abundance > 0)
  f1 <- sum(abundance == 1)
  f2 <- sum(abundance == 2)
  if (f1 == 0 && f2 == 0) return(NA_real_)
  if (f2 > 0) obs + f1^2 / (2 * f2) else obs + f1 * (f1 - 1) / 2
}

#' D50: percentage of clonotypes accounting for half the repertoire
#'
#' Reported as a percentage of observed clonotypes, which is the convention
#' in the commercial repertoire platforms where the index originates. The
#' raw count is also strongly depth-dependent; the percentage less so,
#' though neither is depth-free.
#' @param abundance Abundance vector of one sample; zeros are ignored.
#' @export
d50_index <- function(abundance) {
  x <- sort(abundance[abundance > 0], decreasing = TRUE)
  if (!length(x)) return(NA_real_)
  100 * which(cumsum(x) >= sum(x) / 2)[1] / length(x)
}

#' Clonal space homeostasis
#'
#' Fraction of the repertoire held by clonotypes in each size class. Reports
#' the shape of the distribution rather than a single index, which is what
#' makes a shift between a few large clones and many small ones visible.
#'
#' @param clono Output of `collapse_clonotypes()`.
#' @param breaks Frequency boundaries defining the classes.
#' @export
clonal_homeostasis <- function(clono,
                               breaks = c(0, 1e-5, 1e-4, 1e-3, 1e-2, 1)) {
  need("dplyr", "clonal_homeostasis()")
  labels <- sprintf("(%g, %g]", breaks[-length(breaks)], breaks[-1])
  clono$size_class <- cut(clono$frequency, breaks = breaks, labels = labels,
                          include.lowest = TRUE)
  # .drop = FALSE keeps size classes with no clonotypes as explicit zero
  # rows. Dropped, an empty class is simply absent from the figure, which
  # reads as "not measured" rather than "none".
  dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(clono, .data$sample_id, .data$subject_id,
                    .data$size_class, .drop = FALSE),
    n_clonotypes = dplyr::n(),
    repertoire_fraction = sum(.data$frequency),
    .groups = "drop"
  ))
}

#' Pairwise repertoire overlap
#'
#' @param clono Output of `collapse_clonotypes()`.
#' @param method "jaccard" (presence), "morisita" (abundance-weighted) or
#'   "shared" (raw count of shared clonotypes).
#' @return Long data frame of pairwise values.
#'
#' Jaccard counts clonotypes and is therefore strongly depth-dependent;
#' Morisita-Horn weights by abundance and is comparatively robust. Where
#' depths differ, report Morisita-Horn, or Jaccard at rarefied depth.
#' @export
repertoire_overlap <- function(clono, method = c("morisita", "jaccard", "shared")) {
  method <- match.arg(method)

  # sample_id is coerced to character up front. Left as a factor, the names
  # below would index the split list by LEVEL CODE rather than by label: with
  # an unused level that is a subscript-out-of-bounds error, and with a
  # different level ordering it silently returns another pair's overlap.
  sid <- as.character(clono$sample_id)
  samples <- unique(sid)
  n <- length(samples)
  if (n < 2L) stop("Overlap needs at least two samples.", call. = FALSE)

  by_sample <- split(clono[, c("clonotype", "abundance", "frequency")],
                     sid)[samples]

  res <- expand.grid(sample_a = samples, sample_b = samples,
                     stringsAsFactors = FALSE)
  res$sample_a <- as.character(res$sample_a)
  res$sample_b <- as.character(res$sample_b)
  res <- res[match(res$sample_a, samples) < match(res$sample_b, samples), ,
             drop = FALSE]
  res$value <- NA_real_

  for (k in seq_len(nrow(res))) {
    A <- by_sample[[res$sample_a[k]]]
    B <- by_sample[[res$sample_b[k]]]
    shared <- intersect(A$clonotype, B$clonotype)

    res$value[k] <- switch(
      method,
      "shared"  = length(shared),
      "jaccard" = length(shared) /
        length(union(A$clonotype, B$clonotype)),
      "morisita" = {
        if (!length(shared)) 0 else {
          pa <- A$frequency[match(shared, A$clonotype)]
          pb <- B$frequency[match(shared, B$clonotype)]
          2 * sum(pa * pb) / (sum(A$frequency^2) + sum(B$frequency^2))
        }
      }
    )
  }
  res$method <- method
  res
}
