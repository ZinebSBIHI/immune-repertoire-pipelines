# ---------------------------------------------------------------------------
# Convergent recombination and public clonotype groups
#
# A "public" group is a set of CDR3 sequences that share a V gene and are
# similar enough to each other to suggest selection on the same epitope,
# found in more than one subject.
#
# Three decisions make or break this analysis, and all three are explicit
# here because the literature is inconsistent about them:
#
#   1. NORMALISED DISTANCE. An absolute edit distance of 1 is a different
#      biological statement on a 6-residue CDR3 than on a 22-residue one. On
#      short CDR3s an absolute threshold typically yields no edges at all,
#      and the analysis silently reports no convergence. Distance is
#      therefore normalised by sequence length.
#
#   2. GERMLINE PROXIMITY. Unmutated, unswitched sequences sit close to
#      germline and therefore close to each other. A convergence analysis
#      that includes naive sequences measures shared germline, not shared
#      selection. Filtering to mutated or class-switched sequences is the
#      default for BCR.
#
#   3. A NULL. The number of public groups found depends on repertoire
#      depth, the number of subjects and the threshold. Without a null it
#      cannot be read. Permuting subject labels while keeping the sequences
#      fixed gives the number of groups expected when no sequence is shared
#      for biological reasons.
# ---------------------------------------------------------------------------

#' Normalised Levenshtein identity between two sets of sequences
#'
#' @param a,b Character vectors of sequences.
#' @param normalise "shorter" (default), "longer" or "mean": the length used
#'   to convert edit distance into identity.
#' @return Matrix of identities in [0, 1], rows a, columns b.
cdr3_identity_matrix <- function(a, b = a, normalise = c("shorter", "longer", "mean")) {
  normalise <- match.arg(normalise)
  need("stringdist", "cdr3_identity_matrix()")

  d <- stringdist::stringdistmatrix(a, b, method = "lv")
  la <- nchar(a); lb <- nchar(b)

  denom <- switch(
    normalise,
    "shorter" = outer(la, lb, pmin),
    "longer"  = outer(la, lb, pmax),
    "mean"    = outer(la, lb, function(x, y) (x + y) / 2)
  )
  # With normalise = "shorter" the edit distance can exceed the shorter
  # length, so the raw value can go below zero for a very unequal pair. It
  # is floored at 0: a negative identity has no meaning, and leaving it
  # unbounded would break any later rescaling even though thresholding is
  # unaffected.
  # The matrix goes first in pmax(): with the scalar first, the dim
  # attribute is not carried over and the result comes back as a plain
  # vector, which then fails as soon as it is indexed as a matrix.
  identity <- pmax(1 - d / denom, 0)
  identity[denom == 0] <- NA_real_
  dim(identity) <- c(length(a), length(b))
  identity
}

#' Detect public CDR3 groups
#'
#' Within each V gene, sequences are linked when their normalised identity
#' reaches `identity_threshold`; connected components of that graph are the
#' candidate groups. A component is public when it spans at least
#' `min_subjects` subjects.
#'
#' Single linkage is used deliberately: a lineage of related CDR3s forms a
#' chain, and requiring every pair in a group to meet the threshold
#' (complete linkage) would fragment it. The cost is that long chains can
#' link distant sequences, so `mean_identity` is reported per group and
#' should be inspected.
#'
#' @param df Data frame with one row per sequence.
#' @param cdr3_col,v_col,subject_col Column names.
#' @param identity_threshold Minimum normalised identity to link two CDR3s.
#' @param min_subjects Minimum distinct subjects for a group to be public.
#' @param normalise Passed to `cdr3_identity_matrix()`.
#' @return List with `groups` (one row per sequence in a public group) and
#'   `summary` (one row per group).
detect_public_groups <- function(df,
                                 cdr3_col = "junction_aa",
                                 v_col = "v_gene",
                                 subject_col = "subject_id",
                                 identity_threshold = 0.75,
                                 min_subjects = 2L,
                                 normalise = "shorter",
                                 components = NULL) {
  need(c("dplyr", "stringdist"), "detect_public_groups()")
  for (cl in c(cdr3_col, v_col, subject_col)) {
    if (!cl %in% names(df)) {
      stop("Column '", cl, "' not found in the input.", call. = FALSE)
    }
  }

  df <- df[!is.na(df[[cdr3_col]]) & nzchar(df[[cdr3_col]]), , drop = FALSE]
  if (!nrow(df)) {
    return(list(groups = df[0, , drop = FALSE],
                summary = empty_group_summary()))
  }

  # Components depend only on the sequences and their V genes, never on the
  # subject labels, so they can be computed once and reused — which is what
  # makes the permutation null affordable.
  if (is.null(components)) {
    components <- cdr3_components(df, cdr3_col, v_col, identity_threshold,
                                 normalise)
  }
  if (!length(components$component_id)) {
    return(list(groups = df[0, , drop = FALSE],
                summary = empty_group_summary()))
  }

  group_rows <- list()
  gid <- 0L

  for (k in seq_along(components$component_id)) {
    rows <- components$rows[[k]]
    subjects <- unique(df[[subject_col]][rows])
    if (length(subjects) < min_subjects) next

    gid <- gid + 1L
    # The V gene is already a column of `df` (named by `v_col`), so it is
    # deliberately not repeated here: binding it again would produce two
    # identically named columns.
    group_rows[[length(group_rows) + 1L]] <- data.frame(
      public_group = sprintf("PG%04d", gid),
      row_id = rows,
      n_members = components$n_members[k],
      n_subjects = length(subjects),
      mean_identity = components$mean_identity[k],
      min_identity = components$min_identity[k],
      stringsAsFactors = FALSE
    )
  }

  if (!length(group_rows)) {
    return(list(groups = df[0, , drop = FALSE],
                summary = empty_group_summary()))
  }

  gr <- dplyr::bind_rows(group_rows)
  groups <- cbind(df[gr$row_id, , drop = FALSE],
                  gr[, setdiff(names(gr), "row_id"), drop = FALSE])
  rownames(groups) <- NULL

  summary <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(groups, .data$public_group,
                    dplyr::across(dplyr::all_of(v_col))),
    n_sequences = dplyr::n(),
    n_distinct_cdr3 = dplyr::n_distinct(.data[[cdr3_col]]),
    n_subjects = dplyr::n_distinct(.data[[subject_col]]),
    subjects = paste(sort(unique(as.character(.data[[subject_col]]))),
                     collapse = ","),
    mean_identity = dplyr::first(.data$mean_identity),
    min_identity = dplyr::first(.data$min_identity),
    mean_cdr3_length = mean(nchar(.data[[cdr3_col]])),
    .groups = "drop"
  ))
  summary <- summary[order(-summary$n_subjects, -summary$n_sequences), ]

  list(groups = groups, summary = summary)
}

#' Permutation null for the number of public groups
#'
#' Subject labels are shuffled while the sequences themselves are held fixed,
#' so the sequence composition and the V gene structure are preserved and
#' only the assignment of sequences to subjects is randomised. The resulting
#' distribution is the number of public groups expected with no
#' between-subject sharing beyond chance.
#'
#' @inheritParams detect_public_groups
#' @param n_permutations Number of label shuffles.
#' @return List with `observed`, `null` (vector of counts), `p_empirical`
#'   and `fold_enrichment`.
permute_public_groups <- function(df,
                                  cdr3_col = "junction_aa",
                                  v_col = "v_gene",
                                  subject_col = "subject_id",
                                  identity_threshold = 0.75,
                                  min_subjects = 2L,
                                  normalise = "shorter",
                                  n_permutations = 200L) {
  df <- df[!is.na(df[[cdr3_col]]) & nzchar(df[[cdr3_col]]), , drop = FALSE]

  # The permutation holds sequences and V genes fixed and shuffles only the
  # subject labels, so the identity matrices and the connected components
  # are identical in every replicate. Computing them once and re-counting
  # subjects per component gives exactly the same null at a fraction of the
  # cost: the component search is the expensive part, and it was being
  # repeated n_permutations + 1 times for no reason.
  components <- cdr3_components(df, cdr3_col, v_col, identity_threshold,
                               normalise)

  count_public <- function(labels) {
    if (!length(components$component_id)) return(0L)
    sum(vapply(components$rows,
               function(rows) length(unique(labels[rows])) >= min_subjects,
               logical(1)))
  }

  observed <- count_public(df[[subject_col]])

  null_counts <- integer(n_permutations)
  for (i in seq_len(n_permutations)) {
    null_counts[i] <- count_public(sample(df[[subject_col]]))
  }

  # Add-one empirical p-value: with n permutations the smallest attainable
  # value is 1/(n+1), and reporting p = 0 would overstate the resolution of
  # the test.
  p_emp <- (1 + sum(null_counts >= observed)) / (n_permutations + 1)

  list(
    observed = observed,
    null = null_counts,
    null_mean = mean(null_counts),
    p_empirical = p_emp,
    fold_enrichment = if (mean(null_counts) > 0) observed / mean(null_counts) else NA_real_,
    n_permutations = n_permutations
  )
}

#' Distance from each sequence to its nearest neighbour in another subject
#'
#' A per-sequence view of convergence that does not depend on a grouping
#' threshold: for every sequence, the identity to the most similar sequence
#' carrying the same V gene in any other subject.
#'
#' @return The input data frame with `nn_identity` and `nn_subject` columns.
nearest_neighbour_other_subject <- function(df,
                                            cdr3_col = "junction_aa",
                                            v_col = "v_gene",
                                            subject_col = "subject_id",
                                            normalise = "shorter") {
  need("stringdist", "nearest_neighbour_other_subject()")

  df$nn_identity <- NA_real_
  df$nn_subject  <- NA_character_

  by_v <- split(seq_len(nrow(df)), df[[v_col]])
  for (idx in by_v) {
    if (length(idx) < 2L) next
    seqs <- df[[cdr3_col]][idx]
    subj <- df[[subject_col]][idx]

    ident <- cdr3_identity_matrix(seqs, seqs, normalise = normalise)
    # Mask same-subject comparisons, including the diagonal.
    same <- outer(subj, subj, `==`)
    ident[same] <- NA_real_

    best <- apply(ident, 1, function(r) {
      if (all(is.na(r))) return(c(NA_real_, NA_real_))
      j <- which.max(r)
      c(r[j], j)
    })
    df$nn_identity[idx] <- best[1, ]
    hit <- !is.na(best[2, ])
    df$nn_subject[idx[hit]] <- subj[best[2, hit]]
  }
  df
}

#' Filter to sequences informative for convergence
#'
#' Keeps sequences that are mutated above `shm_threshold` or class-switched.
#' Unmutated unswitched sequences are excluded because their similarity
#' reflects shared germline rather than shared selection.
#'
#' @param df Data frame.
#' @param shm_col Column of mutation frequency (fraction, not percent).
#' @param isotype_col Column of isotype calls.
#' @param unswitched Isotypes counted as unswitched.
filter_informative_sequences <- function(df,
                                         shm_col = "mu_freq",
                                         isotype_col = "c_gene",
                                         shm_threshold = 0.02,
                                         unswitched = c("IGHM", "IGHD")) {
  has_shm <- shm_col %in% names(df)
  has_iso <- isotype_col %in% names(df)
  if (!has_shm && !has_iso) {
    warning("Neither '", shm_col, "' nor '", isotype_col, "' present; ",
            "convergence will include germline-proximal sequences and is ",
            "likely to be inflated.", call. = FALSE)
    return(df)
  }

  mutated  <- if (has_shm) !is.na(df[[shm_col]]) & df[[shm_col]] >= shm_threshold
              else rep(FALSE, nrow(df))
  # An absent isotype call is not evidence of switching. Treating it as
  # such would retain exactly the germline-proximal sequences this filter
  # exists to remove.
  switched <- if (has_iso) {
                g <- sub("\\*.*$", "", sub(",.*$", "",
                                           as.character(df[[isotype_col]])))
                !is.na(g) & nzchar(trimws(g)) & !g %in% unswitched
              } else rep(FALSE, nrow(df))

  keep <- mutated | switched
  message(sprintf("Convergence filter: %d of %d sequences retained (%.1f%%) ",
                  sum(keep), nrow(df), 100 * mean(keep)),
          "as mutated or class-switched.")
  df[keep, , drop = FALSE]
}

# --- internals -------------------------------------------------------------

#' Connected components of the CDR3 similarity graph, within each V gene
#'
#' Depends only on the sequences and their V genes, never on subject labels.
#' Separating this from the public-group test is what lets the permutation
#' null reuse one component search instead of repeating it per replicate.
#'
#' @return A list with, per component: `component_id`, `rows` (indices into
#'   `df`), `n_members` (distinct CDR3s), `mean_identity`, `min_identity`.
#'   Components of a single CDR3 are not returned.
cdr3_components <- function(df, cdr3_col, v_col, identity_threshold,
                            normalise = "shorter",
                            max_sequences_per_gene = 20000L) {
  need("stringdist", "cdr3_components()")

  by_v <- split(seq_len(nrow(df)), df[[v_col]])
  out <- list(component_id = character(0), rows = list(),
              n_members = integer(0), mean_identity = numeric(0),
              min_identity = numeric(0))
  cid_global <- 0L

  for (v in names(by_v)) {
    idx <- by_v[[v]]
    seqs <- unique(df[[cdr3_col]][idx])
    if (length(seqs) < 2L) next

    # The identity matrix is dense and quadratic in the number of distinct
    # CDR3s carrying one V gene. Deep bulk repertoires can exceed what fits
    # in memory, so the limit is checked and reported rather than allowed to
    # exhaust the machine mid-run.
    if (length(seqs) > max_sequences_per_gene) {
      warning(sprintf(
        paste0("V gene %s carries %s distinct CDR3s, above the limit of %s; ",
               "it is skipped. A dense identity matrix at that size needs ",
               "roughly %.0f GB. Reduce the input (for example by analysing ",
               "clonotypes rather than sequences, or by raising the ",
               "abundance floor) before interpreting convergence."),
        v, format(length(seqs), big.mark = " "),
        format(max_sequences_per_gene, big.mark = " "),
        3 * 8 * length(seqs)^2 / 1e9), call. = FALSE)
      next
    }

    ident <- cdr3_identity_matrix(seqs, seqs, normalise = normalise)
    adj <- !is.na(ident) & ident >= identity_threshold
    diag(adj) <- FALSE
    if (!any(adj)) next

    comp <- connected_components(adj)

    for (cid in unique(comp)) {
      in_comp <- comp == cid
      if (sum(in_comp) < 2L) next

      members <- seqs[in_comp]
      sub_ident <- ident[in_comp, in_comp, drop = FALSE]
      off_diag <- sub_ident[lower.tri(sub_ident)]

      cid_global <- cid_global + 1L
      out$component_id <- c(out$component_id, sprintf("C%06d", cid_global))
      out$rows[[length(out$rows) + 1L]] <-
        idx[df[[cdr3_col]][idx] %in% members]
      out$n_members     <- c(out$n_members, length(members))
      out$mean_identity <- c(out$mean_identity, mean(off_diag, na.rm = TRUE))
      out$min_identity  <- c(out$min_identity, min(off_diag, na.rm = TRUE))
    }
  }
  out
}

#' Connected components of a logical adjacency matrix
#'
#' Breadth-first search. Avoids a hard dependency on igraph for what is a
#' few lines, and keeps the component definition visible.
connected_components <- function(adj) {
  n <- nrow(adj)
  comp <- integer(n)
  current <- 0L

  for (start in seq_len(n)) {
    if (comp[start] != 0L) next
    current <- current + 1L
    queue <- start
    comp[start] <- current

    while (length(queue)) {
      node <- queue[1L]
      queue <- queue[-1L]
      nbrs <- which(adj[node, ] & comp == 0L)
      if (length(nbrs)) {
        comp[nbrs] <- current
        queue <- c(queue, nbrs)
      }
    }
  }
  comp
}

empty_group_summary <- function() {
  data.frame(public_group = character(0), v_gene = character(0),
             n_sequences = integer(0), n_distinct_cdr3 = integer(0),
             n_subjects = integer(0), subjects = character(0),
             mean_identity = numeric(0), min_identity = numeric(0),
             mean_cdr3_length = numeric(0), stringsAsFactors = FALSE)
}
