# ---------------------------------------------------------------------------
# Repertoire metrics on AIRR-formatted tables
#
# Input throughout this file is a long data frame in AIRR Rearrangement
# format (Vander Heiden et al. 2018, PMID 30323809), one row per
# rearrangement, with at least:
#
#   sample_id, subject_id, v_call, j_call, junction, junction_aa, locus
#   and an abundance column (duplicate_count or umi_count)
#
# Two points of method are enforced here rather than left to the caller,
# because getting either wrong invalidates every comparison downstream.
#
# ABUNDANCE COLUMN. UMI-based protocols and read-based protocols are not
# interchangeable. `duplicate_count` in an AIRR file produced by one tool may
# be reads and by another may be UMI-collapsed molecules; using reads where
# UMIs exist inflates abundance by the PCR duplication factor, which differs
# per sample. The column is therefore named explicitly, never guessed.
#
# SEQUENCING DEPTH. Richness, clone counts and clonal concentration are all
# monotone in depth. Comparing them between samples of unequal depth measures
# the depth difference. Any such comparison goes through
# `rarefy_to_common_depth()`.
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

#' Per-sample repertoire summary
#'
#' @param clono Output of `collapse_clonotypes()`.
#' @return Data frame of depth, richness and clonality metrics per sample.
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
shannon_clonality <- function(abundance) {
  p <- abundance[abundance > 0]
  if (length(p) < 2L) return(NA_real_)
  1 - shannon_entropy(p) / log(length(p))
}

#' Inverse Simpson index
inverse_simpson <- function(abundance) {
  p <- abundance[abundance > 0]
  if (!length(p)) return(NA_real_)
  p <- p / sum(p)
  1 / sum(p^2)
}

#' Gini index of abundance inequality
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
