# ---------------------------------------------------------------------------
# CDR3 structure, physicochemistry, sequence similarity and diversity spectra
#
# This file answers one question in several ways: what do the CDR3s look
# like, and do the ones that expanded look different from the rest?
#
# WHY LENGTH AND CHARGE. A selected clonotype is not a random draw from the
# repertoire. Selection on an epitope acts on the loop that contacts it, so
# a shift in CDR3 length, net charge or hydrophobicity between the expanded
# clonotypes and the rest of the same subject's repertoire is a fingerprint
# of selection rather than of recombination. The comparison is made within
# subject, against that subject's own repertoire, because the distributions
# differ between people.
#
# WHY WEIGHT BY ABUNDANCE. A length or charge distribution over distinct
# clonotypes and the same distribution over molecules are different
# statements: the first describes the recombination machinery, the second
# describes the cells actually present. Every distribution here is weighted
# by abundance and says so in its legend; the unweighted clonotype-level
# view is in section_gene_usage().
#
# WHY NORMALISED EDIT DISTANCE IN THE NETWORK. An absolute Levenshtein
# distance of 1 is a far stronger statement about two 20-residue CDR3s than
# about two 9-residue ones, and on real data an absolute threshold of 1
# typically returns a graph with no edges at all, which reads as "no
# convergence" when it is really "wrong threshold". Distance is therefore
# divided by the length of the LONGER sequence of each pair, which is the
# conservative choice (it cannot flatter a pair of unequal length), and
# identical sequences are excluded so that the edges carry information
# about similarity rather than about duplicate rows.
#
# WHY A THRESHOLD DIAGNOSTIC. The number of edges is extremely sensitive to
# that fraction, and a single figure drawn at one threshold hides it.
# cdr3_network_diagnostic() reports the edge count at several thresholds so
# that the chosen value can be seen as a choice and reported in Methods.
# ---------------------------------------------------------------------------

#' CDR3 length spectratype, weighted by abundance
#'
#' The expected shape is a near-Gaussian bell; a human TCR beta repertoire
#' peaks around 14-15 amino acids. A sharp peak at one length is an oligoclonal
#' expansion large enough to deform the distribution of the whole
#' repertoire, which is itself the finding.
#'
#' Related but not the same as [cdr3_length_spectrum()], which works from
#' collapsed clonotypes and reports counts per sample; this works from the
#' rearrangement table and reports a within-group frequency, so that groups
#' of unequal depth are comparable.
#'
#' @param db AIRR data frame with `junction_aa` and an abundance column.
#' @param abundance_col Name of the abundance column (UMI or read count).
#' @param group_cols Columns defining the groups within which frequency is
#'   computed, e.g. subject and timepoint.
#' @param trim Drop the conserved first and last residue before measuring
#'   length. Off by default, so that the lengths match the untrimmed
#'   convention [cdr3_properties()] reports.
#' @return Data frame: group columns, `cdr3_length`, `abundance`,
#'   `frequency`.
#' @export
cdr3_spectratype <- function(db, abundance_col,
                             group_cols = c("subject_id", "timepoint"),
                             trim = FALSE) {
  need("dplyr", "cdr3_spectratype()")
  group_cols <- intersect(group_cols, names(db))
  if (!length(group_cols)) {
    stop("cdr3_spectratype() needs at least one grouping column present in ",
         "the input.", call. = FALSE)
  }
  if (!abundance_col %in% names(db)) {
    stop("Abundance column '", abundance_col, "' not found.", call. = FALSE)
  }

  x <- db[!is.na(db$junction_aa) & nzchar(db$junction_aa), , drop = FALSE]
  len <- nchar(as.character(x$junction_aa))
  # Trimming removes the conserved cysteine and the terminal F/W, which are
  # present in every sequence and therefore carry no information; it shifts
  # every length by two and is only applied when asked for, so that a
  # spectratype is never silently on a different scale from a published one.
  x$cdr3_length <- if (isTRUE(trim)) pmax(len - 2L, 0L) else len
  x$.abund <- as.numeric(x[[abundance_col]])

  out <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(x, dplyr::across(dplyr::all_of(c(group_cols,
                                                     "cdr3_length")))),
    abundance = sum(.data$.abund, na.rm = TRUE),
    n_rearrangements = dplyr::n(),
    .groups = "drop"))

  dplyr::ungroup(dplyr::mutate(
    dplyr::group_by(out, dplyr::across(dplyr::all_of(group_cols))),
    frequency = .data$abundance / sum(.data$abundance)))
}

#' Physicochemical properties of the CDR3 amino acid sequence
#'
#' Thin wrapper over `alakazam::aminoAcidProperties()`, kept here so that
#' every pipeline computes the properties the same way and with the same
#' column names. Length, net charge and hydrophobicity (GRAVY) are the three
#' used downstream; the function returns everything alakazam computes.
#'
#' @param db AIRR data frame.
#' @param seq_col Column holding the CDR3 amino acid sequence.
#' @param trim Passed to alakazam. Note what it does and does not do: in
#'   alakazam 1.3.0 `trim` is applied only while translating NUCLEOTIDE
#'   input, so with the amino acid input used here it has no effect and the
#'   conserved first and last residues are INCLUDED in every property. That
#'   is the behaviour the original analyses ran with, so it is documented
#'   here rather than silently corrected: trimming the sequences would shift
#'   every length, charge and GRAVY value and make them incomparable with
#'   numbers already reported. Trim beforehand if the trimmed convention is
#'   wanted, and say so in Methods.
#' @param label Prefix for the generated columns, giving `cdr3_aa_length`,
#'   `cdr3_aa_charge`, `cdr3_aa_gravy` and the rest.
#' @return The input data frame with the property columns added.
#' @export
cdr3_properties <- function(db, seq_col = "junction_aa", trim = TRUE,
                            label = "cdr3") {
  need("alakazam", "cdr3_properties()")
  if (!seq_col %in% names(db)) {
    stop("Sequence column '", seq_col, "' not found.", call. = FALSE)
  }
  db <- db[!is.na(db[[seq_col]]) & nzchar(db[[seq_col]]), , drop = FALSE]
  alakazam::aminoAcidProperties(as.data.frame(db), seq = seq_col, nt = FALSE,
                                trim = trim, label = label)
}

#' Abundance-weighted summary of the CDR3 properties per group
#'
#' One row per group with the abundance-weighted mean length, charge and
#' hydrophobicity. Weighted rather than plain means, for the reason given at
#' the top of this file: the unweighted mean describes the clonotypes, the
#' weighted mean describes the molecules.
#'
#' @param props Output of [cdr3_properties()].
#' @param abundance_col Name of the abundance column.
#' @param group_cols Grouping columns.
#' @param label Property column prefix used in [cdr3_properties()].
#' @return Data frame: group columns, `cdr3_length`, `cdr3_charge`,
#'   `cdr3_gravy`, `n_rearrangements`, `abundance`.
#' @export
cdr3_property_summary <- function(props, abundance_col,
                                  group_cols = c("subject_id", "timepoint"),
                                  label = "cdr3") {
  need("dplyr", "cdr3_property_summary()")
  group_cols <- intersect(group_cols, names(props))
  cols <- paste0(label, c("_aa_length", "_aa_charge", "_aa_gravy"))
  missing_cols <- setdiff(cols, names(props))
  if (length(missing_cols)) {
    stop("cdr3_property_summary() needs column(s): ",
         paste(missing_cols, collapse = ", "),
         "\nRun cdr3_properties() first, with the same `label`.",
         call. = FALSE)
  }

  props$.abund <- as.numeric(props[[abundance_col]])
  props$.len    <- props[[cols[1]]]
  props$.charge <- props[[cols[2]]]
  props$.gravy  <- props[[cols[3]]]

  dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(props, dplyr::across(dplyr::all_of(group_cols))),
    n_rearrangements = dplyr::n(),
    abundance    = sum(.data$.abund, na.rm = TRUE),
    cdr3_length  = stats::weighted.mean(.data$.len, w = .data$.abund,
                                        na.rm = TRUE),
    cdr3_charge  = stats::weighted.mean(.data$.charge, w = .data$.abund,
                                        na.rm = TRUE),
    cdr3_gravy   = stats::weighted.mean(.data$.gravy, w = .data$.abund,
                                        na.rm = TRUE),
    .groups = "drop"))
}

#' Convergent recombination within a sample: one CDR3, several junctions
#'
#' The fraction of distinct CDR3 amino acid sequences that are encoded by
#' more than one nucleotide junction in the same sample. Independent
#' recombination events converging on the same amino acid sequence is
#' evidence of selection on that sequence, and unlike the between-subject
#' convergence in `ana-convergence.R` it needs no threshold and no null:
#' it is a count of exact matches.
#'
#' Read with depth in mind. The chance of observing two junctions for one
#' CDR3 rises with sequencing depth, so this index is comparable across
#' timepoints of one subject at similar depth and not across samples of very
#' unequal depth.
#'
#' @param db AIRR data frame with `junction` and `junction_aa`.
#' @param group_cols Columns defining the unit within which convergence is
#'   counted; must identify a sample.
#' @return Data frame: group columns, `n_cdr3_aa`, `n_convergent`,
#'   `pct_convergent`.
#' @export
nucleotide_convergence <- function(db,
                                   group_cols = c("sample_id", "subject_id",
                                                  "timepoint")) {
  need("dplyr", "nucleotide_convergence()")
  for (cl in c("junction", "junction_aa")) {
    if (!cl %in% names(db)) {
      stop("nucleotide_convergence() needs column '", cl, "'.", call. = FALSE)
    }
  }
  group_cols <- intersect(group_cols, names(db))

  per_cdr3 <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(db, dplyr::across(dplyr::all_of(c(group_cols,
                                                      "junction_aa")))),
    n_junctions = dplyr::n_distinct(.data$junction),
    .groups = "drop"))

  dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(per_cdr3, dplyr::across(dplyr::all_of(group_cols))),
    n_cdr3_aa = dplyr::n(),
    n_convergent = sum(.data$n_junctions > 1L),
    pct_convergent = 100 * mean(.data$n_junctions > 1L),
    .groups = "drop"))
}

# --- Diversity spectra -----------------------------------------------------

#' Hill numbers of a clonotype abundance distribution
#'
#' Hill MO. Diversity and evenness: a unifying notation and its
#' consequences. Ecology. 1973;54:427-432.
#'
#' `D_q = (sum p^q)^(1/(1-q))`, with the limiting case `D_1 = exp(H)`. One
#' curve replaces the argument about which diversity index to use: q = 0 is
#' richness, which counts every clonotype equally and is the most
#' depth-dependent; q = 1 is the exponential of Shannon entropy; q = 2 is
#' inverse Simpson; higher q weights the dominant clonotypes more heavily.
#' A high, flat profile is a rich and even repertoire; a profile that drops
#' steeply with q is one where a few clonotypes carry the repertoire.
#'
#' @param abundance Numeric vector of clonotype abundances.
#' @param q Orders at which to evaluate the profile.
#' @return Data frame with `q` and `hill`.
#' @export
hill_numbers <- function(abundance, q = 0:6) {
  p <- abundance[!is.na(abundance) & abundance > 0]
  if (!length(p)) {
    return(data.frame(q = q, hill = NA_real_))
  }
  p <- p / sum(p)
  value <- vapply(q, function(qq) {
    # q = 1 is a removable singularity: the exponent 1/(1-q) is undefined
    # there, and the limit is exp(H). Evaluating the general formula at
    # q = 1 would return Inf or NaN rather than the intended value.
    if (isTRUE(all.equal(qq, 1))) exp(-sum(p * log(p)))
    else sum(p^qq)^(1 / (1 - qq))
  }, numeric(1))
  data.frame(q = q, hill = value)
}

#' Hill profile for every sample
#'
#' @param clono Output of [collapse_clonotypes()].
#' @param q Orders at which to evaluate the profile.
#' @return Data frame: `sample_id`, `subject_id`, `q`, `hill`.
#' @export
hill_spectrum <- function(clono, q = 0:6) {
  need("dplyr", "hill_spectrum()")
  sid <- as.character(clono$sample_id)
  rows <- lapply(unique(sid), function(s) {
    sub <- clono[sid == s, , drop = FALSE]
    out <- hill_numbers(sub$abundance, q = q)
    out$sample_id  <- s
    out$subject_id <- as.character(sub$subject_id[1])
    out
  })
  dplyr::bind_rows(rows)
}

#' Expected clonotype richness at a reduced sequencing depth
#'
#' The exact expectation of the number of clonotypes observed in a random
#' subsample of `n` molecules drawn without replacement:
#'
#'   `E[S_n] = sum_i ( 1 - C(N - N_i, n) / C(N, n) )`
#'
#' Hurlbert SH. The nonconcept of species diversity. Ecology. 1971;52:577-586.
#'
#' Computed analytically rather than by resampling: the expectation is the
#' quantity a rarefaction curve is meant to show, and an exact value removes
#' both the Monte Carlo noise and the dependence on a resampling seed. The
#' binomial coefficients are evaluated on the log scale, because the
#' coefficients themselves overflow at the depths bulk repertoires reach.
#'
#' @param abundance Numeric vector of clonotype abundances.
#' @param depths Depths at which to evaluate the curve. By default, a grid
#'   from a small depth up to the observed depth.
#' @param n_points Number of depths in the default grid.
#' @return Data frame with `depth` and `expected_richness`.
#' @export
rarefaction_curve <- function(abundance, depths = NULL, n_points = 20L) {
  x <- abundance[!is.na(abundance) & abundance > 0]
  total <- sum(x)
  if (!length(x) || total < 1) {
    return(data.frame(depth = numeric(0), expected_richness = numeric(0)))
  }
  if (is.null(depths)) {
    depths <- unique(round(seq(max(1, floor(total / n_points)), total,
                               length.out = n_points)))
  }
  depths <- sort(unique(depths[depths >= 1 & depths <= total]))

  expected <- vapply(depths, function(n) {
    # P(clonotype i unseen in a sample of n) = C(N - N_i, n) / C(N, n),
    # which is zero whenever the clonotype is too abundant to be missed.
    log_miss <- lchoose(total - x, n) - lchoose(total, n)
    miss <- exp(log_miss)
    miss[(total - x) < n] <- 0
    sum(1 - miss)
  }, numeric(1))

  data.frame(depth = depths, expected_richness = expected)
}

# --- CDR3 similarity network -----------------------------------------------

#' Pairwise CDR3 similarity network
#'
#' Nodes are distinct CDR3 amino acid sequences; two nodes are linked when
#' their Levenshtein distance, divided by the length of the longer of the
#' two, is at most `max_distance_fraction`, and when they are not identical.
#' See the note at the top of this file for why the distance is normalised
#' and why identical sequences are excluded.
#'
#' @param sequences Character vector of CDR3 amino acid sequences; duplicates
#'   are collapsed.
#' @param abundance Abundance per sequence, used only to choose which
#'   sequences to keep when the cap below applies.
#' @param focus Logical vector marking sequences of interest (for example the
#'   significantly expanded clonotypes). Every focus sequence is kept
#'   whatever the cap.
#' @param max_distance_fraction Distance threshold as a fraction of the
#'   longer sequence.
#' @param min_length Sequences shorter than this are dropped as
#'   uninformative.
#' @param max_sequences Cap on the number of nodes. The identity matrix is
#'   dense and quadratic, so the cap is what keeps the figure computable on
#'   a deep repertoire; focus sequences are kept first, then the most
#'   abundant of the rest.
#' @param layout_seed Seed for the force-directed layout, so a rerun
#'   reproduces the same picture.
#' @return List with `nodes` (`sequence`, `abundance`, `focus`, `x`, `y`,
#'   `degree`), `edges` (`x_from`, `y_from`, `x_to`, `y_to`), `n_nodes`,
#'   `n_edges` and `n_candidates`. `n_nodes` is zero when no pair met the
#'   threshold, which is a result and not an error.
#' @export
cdr3_similarity_network <- function(sequences, abundance = NULL, focus = NULL,
                                    max_distance_fraction = 0.20,
                                    min_length = 8L, max_sequences = 400L,
                                    layout_seed = 1234L) {
  need(c("stringdist", "igraph"), "cdr3_similarity_network()")

  seqs <- as.character(sequences)
  if (is.null(abundance)) abundance <- rep(1, length(seqs))
  if (is.null(focus))     focus <- rep(FALSE, length(seqs))
  abundance <- as.numeric(abundance)
  focus <- as.logical(focus)
  focus[is.na(focus)] <- FALSE

  keep <- !is.na(seqs) & nchar(seqs) >= min_length
  seqs <- seqs[keep]; abundance <- abundance[keep]; focus <- focus[keep]
  if (!length(seqs)) return(empty_cdr3_network())

  # Collapse duplicates: a node is a sequence, so the same CDR3 appearing in
  # several rearrangements must not become several nodes, which would create
  # zero-distance pairs and a cluster of artefacts.
  ord <- order(seqs)
  seqs <- seqs[ord]; abundance <- abundance[ord]; focus <- focus[ord]
  agg <- stats::aggregate(list(abundance = abundance, focus = as.integer(focus)),
                          by = list(sequence = seqs), FUN = sum)
  agg$focus <- agg$focus > 0L
  n_candidates <- nrow(agg)

  if (nrow(agg) > max_sequences) {
    foc <- agg[agg$focus, , drop = FALSE]
    oth <- agg[!agg$focus, , drop = FALSE]
    room <- max(0L, max_sequences - nrow(foc))
    oth <- oth[order(-oth$abundance), , drop = FALSE]
    agg <- rbind(foc, utils::head(oth, room))
  }
  if (nrow(agg) < 3L) return(empty_cdr3_network(n_candidates))

  s <- agg$sequence
  d <- stringdist::stringdistmatrix(s, s, method = "lv")
  len <- nchar(s)
  denom <- outer(len, len, pmax)
  frac <- d / denom
  adj <- frac <= max_distance_fraction & d > 0
  adj[is.na(adj)] <- FALSE
  diag(adj) <- FALSE
  if (!any(adj)) return(empty_cdr3_network(n_candidates))

  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected",
                                           diag = FALSE)
  igraph::V(g)$name  <- s
  igraph::V(g)$focus <- agg$focus
  g <- igraph::delete_vertices(g, igraph::degree(g) == 0)
  if (igraph::gorder(g) == 0L) return(empty_cdr3_network(n_candidates))

  # The layout is stochastic, so it is seeded: two runs of the same analysis
  # should produce the same figure, not two pictures of the same graph.
  set.seed(layout_seed)
  lay <- igraph::layout_with_fr(g)
  nodes <- data.frame(sequence = igraph::V(g)$name,
                      focus = igraph::V(g)$focus,
                      degree = as.integer(igraph::degree(g)),
                      x = lay[, 1], y = lay[, 2],
                      stringsAsFactors = FALSE)
  nodes$abundance <- agg$abundance[match(nodes$sequence, agg$sequence)]

  e <- igraph::as_data_frame(g, what = "edges")
  edges <- data.frame(
    x_from = nodes$x[match(e$from, nodes$sequence)],
    y_from = nodes$y[match(e$from, nodes$sequence)],
    x_to   = nodes$x[match(e$to,   nodes$sequence)],
    y_to   = nodes$y[match(e$to,   nodes$sequence)],
    stringsAsFactors = FALSE)

  list(nodes = nodes, edges = edges, n_nodes = nrow(nodes),
       n_edges = nrow(edges), n_candidates = n_candidates)
}

#' Save a figure whose statistic consumes the weight aesthetic
#'
#' Every distribution in this file is weighted by abundance, for the reason
#' given at the top. ggplot2 3.4 does not list `weight` among the aesthetics
#' `stat_ydensity()` and `stat_density()` declare, so each weighted violin
#' and density emits "the following aesthetics were dropped during
#' statistical transformation: weight" even though the weights are applied:
#' a weighted and an unweighted build of the same plot differ. The kernel
#' density additionally reports that its bandwidth was selected without the
#' weights, which is true and is stated in the report text instead, where a
#' reader will see it once rather than once per panel. Those two are muffled
#' here; nothing else is.
#' @noRd
save_weighted_figure <- function(...) {
  quiet <- c("aesthetics were dropped", "Selecting bandwidth")
  withCallingHandlers(
    save_figure(...),
    warning = function(w) {
      msg <- conditionMessage(w)
      if (any(vapply(quiet, grepl, logical(1), x = msg, fixed = TRUE))) {
        invokeRestart("muffleWarning")
      }
    })
}

#' An empty network result
#' @noRd
empty_cdr3_network <- function(n_candidates = 0L) {
  list(nodes = data.frame(sequence = character(0), focus = logical(0),
                          degree = integer(0), x = numeric(0), y = numeric(0),
                          abundance = numeric(0), stringsAsFactors = FALSE),
       edges = data.frame(x_from = numeric(0), y_from = numeric(0),
                          x_to = numeric(0), y_to = numeric(0)),
       n_nodes = 0L, n_edges = 0L, n_candidates = n_candidates)
}

#' Edges obtained at several similarity thresholds
#'
#' The edge count of a CDR3 similarity graph is extremely sensitive to the
#' threshold, and a graph drawn at one value gives no indication of that.
#' This reports the count at each candidate threshold, including the
#' absolute distance of 1 that the literature often uses, so the chosen
#' value can be reported as a choice.
#'
#' @param sequences Character vector of CDR3 amino acid sequences.
#' @param fractions Normalised thresholds to evaluate.
#' @param absolute_distances Absolute Levenshtein thresholds to evaluate.
#' @param min_length,max_sequences As in [cdr3_similarity_network()].
#' @return Data frame: `criterion`, `threshold`, `n_sequences`, `n_edges`,
#'   `n_connected_nodes`.
#' @export
cdr3_network_diagnostic <- function(sequences, fractions = c(0.10, 0.15, 0.20),
                                    absolute_distances = 1L,
                                    min_length = 8L, max_sequences = 400L) {
  need("stringdist", "cdr3_network_diagnostic()")
  s <- unique(as.character(sequences))
  s <- s[!is.na(s) & nchar(s) >= min_length]
  if (length(s) > max_sequences) s <- utils::head(s, max_sequences)
  if (length(s) < 2L) {
    return(data.frame(criterion = character(0), threshold = numeric(0),
                      n_sequences = integer(0), n_edges = integer(0),
                      n_connected_nodes = integer(0),
                      stringsAsFactors = FALSE))
  }

  d <- stringdist::stringdistmatrix(s, s, method = "lv")
  len <- nchar(s)
  frac <- d / outer(len, len, pmax)

  count <- function(adj) {
    adj[is.na(adj)] <- FALSE
    diag(adj) <- FALSE
    c(n_edges = sum(adj) / 2, n_connected = sum(rowSums(adj) > 0))
  }

  rows <- list()
  for (f in fractions) {
    k <- count(frac <= f & d > 0)
    rows[[length(rows) + 1L]] <- data.frame(
      criterion = "fraction of longer length", threshold = f,
      n_sequences = length(s), n_edges = as.integer(k[["n_edges"]]),
      n_connected_nodes = as.integer(k[["n_connected"]]),
      stringsAsFactors = FALSE)
  }
  for (a in absolute_distances) {
    k <- count(d <= a & d > 0)
    rows[[length(rows) + 1L]] <- data.frame(
      criterion = "absolute edit distance", threshold = a,
      n_sequences = length(s), n_edges = as.integer(k[["n_edges"]]),
      n_connected_nodes = as.integer(k[["n_connected"]]),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

# --- Sequence logos --------------------------------------------------------

#' Sequence logos of the CDR3 motifs of shared groups
#'
#' One logo per group, over the distinct CDR3 amino acid sequences of that
#' group. Conserved positions are the convergent motif; a logo in which
#' every position is equally mixed is a group held together by the
#' similarity threshold rather than by a motif.
#'
#' ALIGNMENT. Logos compare positions, so the sequences must be on a common
#' coordinate system. `align = "msa"` aligns them with DECIPHER, which is
#' the right answer when that package is installed. `align = "by_length"`
#' instead draws one logo per CDR3 length within the group, so that every
#' logo is over sequences of equal length and no gap character is invented;
#' it shows less at once but cannot misalign. `align = "auto"` uses the
#' multiple alignment when DECIPHER is available and falls back to the
#' per-length view otherwise, and the chosen route is returned so the caller
#' can state it.
#'
#' @param df Data frame with one row per sequence.
#' @param group_col Column identifying the group.
#' @param cdr3_col Column holding the CDR3 amino acid sequence.
#' @param label_col Optional column whose first value labels each panel
#'   (typically the V gene).
#' @param min_sequences Minimum distinct sequences for a group to get a logo.
#' @param max_groups Maximum number of groups to draw.
#' @param align One of "auto", "msa", "by_length".
#' @return List with `plot` (a ggplot, or NULL when nothing qualified),
#'   `alignment` (the route actually used), `groups` (the groups drawn).
#' @export
cdr3_motif_logos <- function(df, group_col = "public_group",
                             cdr3_col = "junction_aa", label_col = NULL,
                             min_sequences = 4L, max_groups = 9L,
                             align = c("auto", "msa", "by_length")) {
  align <- match.arg(align)
  need("ggseqlogo", "cdr3_motif_logos()")
  if (!nrow(df) || !group_col %in% names(df)) {
    return(list(plot = NULL, alignment = NA_character_,
                groups = character(0)))
  }

  have_msa <- requireNamespace("DECIPHER", quietly = TRUE) &&
    requireNamespace("Biostrings", quietly = TRUE)
  route <- switch(align,
                  "auto" = if (have_msa) "msa" else "by_length",
                  "msa" = { need(c("DECIPHER", "Biostrings"),
                                 "cdr3_motif_logos(align = \"msa\")")
                            "msa" },
                  "by_length" = "by_length")

  groups <- unique(as.character(df[[group_col]]))
  sets <- list()
  for (g in groups) {
    rows <- df[as.character(df[[group_col]]) == g, , drop = FALSE]
    seqs <- unique(as.character(rows[[cdr3_col]]))
    seqs <- seqs[!is.na(seqs) & nzchar(seqs)]
    if (length(seqs) < min_sequences) next

    tag <- if (!is.null(label_col) && label_col %in% names(rows)) {
      sprintf("%s\n%s (n = %d)", g, as.character(rows[[label_col]][1]),
              length(seqs))
    } else {
      sprintf("%s (n = %d)", g, length(seqs))
    }

    if (route == "msa") {
      aligned <- as.character(
        DECIPHER::AlignSeqs(Biostrings::AAStringSet(seqs), verbose = FALSE))
      sets[[tag]] <- unname(aligned)
    } else {
      # One logo per length: positions are then comparable by construction,
      # and a group spanning several lengths shows as several panels rather
      # than as one misaligned panel.
      by_len <- split(seqs, nchar(seqs))
      for (L in names(by_len)) {
        if (length(by_len[[L]]) < min_sequences) next
        sets[[sprintf("%s | %s aa (n = %d)", tag, L,
                      length(by_len[[L]]))]] <- by_len[[L]]
      }
    }
    if (length(sets) >= max_groups) break
  }

  if (!length(sets)) {
    return(list(plot = NULL, alignment = route, groups = character(0)))
  }
  p <- ggseqlogo::ggseqlogo(sets, method = "prob",
                            ncol = min(3L, length(sets))) +
    ggplot2::labs(title = "Convergent CDR3 motifs") +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  list(plot = p, alignment = route, groups = names(sets))
}

# ---------------------------------------------------------------------------
# Sections
# ---------------------------------------------------------------------------

#' Diversity spectrum: Hill profile, rarefaction curves and D50
#'
#' Three views of the same distribution, chosen because each fails
#' differently. The Hill profile shows whether a conclusion depends on the
#' choice of diversity index. The rarefaction curves show whether richness
#' has saturated: a curve still climbing at the observed depth means the
#' sample is undersampled and its richness is a lower bound. D50 reports
#' how few clonotypes carry half the repertoire, which is the quantity
#' clinical repertoire reports usually quote.
#'
#' @param clono Output of [collapse_clonotypes()].
#' @param run_dir,cfg,report Run directory, configuration and report writer.
#' @param unit Unit of observation, used in labels.
#' @param prefix File-name prefix.
#' @return The Hill spectrum, invisibly.
#' @export
section_diversity_spectrum <- function(clono, run_dir, cfg, report,
                                       unit = "molecules", prefix = "04b") {
  need("dplyr", "section_diversity_spectrum()")
  report$h1("Diversity spectrum")

  q <- as.numeric(cfg$repertoire$hill_q %||% 0:6)
  hill <- hill_spectrum(clono, q = q)
  save_table(hill, paste0(prefix, "_hill_spectrum"), run_dir,
             "Hill number at each order q, per sample")

  report$p(paste0(
    "The Hill profile evaluates diversity at several orders q: q = 0 is ",
    "richness, q = 1 the exponential of Shannon entropy, q = 2 inverse ",
    "Simpson, and higher orders weight the dominant clonotypes more. A ",
    "high flat profile is a rich and even repertoire; a steep drop with q ",
    "means a few clonotypes carry the repertoire. Reporting the profile ",
    "rather than one index removes the question of which index was chosen."))
  report$p("")

  p_hill <- ggplot2::ggplot(hill, ggplot2::aes(.data$q, .data$hill,
                                               group = .data$sample_id,
                                               colour = .data$subject_id)) +
    ggplot2::geom_line(linewidth = 0.6, alpha = 0.85) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "Order q", y = "Hill number (log scale)",
                  colour = "Subject", title = "Hill diversity profile") +
    theme_pub()
  save_figure(p_hill, paste0(prefix, "a_hill_spectrum"), run_dir, cfg,
              width = 6.5, height = 4,
              legend = paste("Hill numbers from q = 0 (richness) upwards,",
                             "one line per sample. A profile that falls",
                             "steeply with q is dominated by a few",
                             "clonotypes."))

  # Rarefaction curves, one per sample, at the exact expected richness.
  sid <- as.character(clono$sample_id)
  n_points <- as.integer(cfg$repertoire$rarefaction$curve_points %||% 20L)
  curves <- dplyr::bind_rows(lapply(unique(sid), function(s) {
    sub <- clono[sid == s, , drop = FALSE]
    out <- rarefaction_curve(sub$abundance, n_points = n_points)
    if (!nrow(out)) return(NULL)
    out$sample_id  <- s
    out$subject_id <- as.character(sub$subject_id[1])
    out
  }))
  save_table(curves, paste0(prefix, "_rarefaction_curves"), run_dir,
             sprintf("expected clonotype richness against depth in %s", unit))

  p_rar <- ggplot2::ggplot(curves,
                           ggplot2::aes(.data$depth,
                                        .data$expected_richness,
                                        group = .data$sample_id,
                                        colour = .data$subject_id)) +
    ggplot2::geom_line(linewidth = 0.6, alpha = 0.9) +
    ggplot2::labs(x = sprintf("Depth (%s)", unit),
                  y = "Expected clonotypes",
                  colour = "Subject", title = "Rarefaction curves",
                  subtitle = paste("Exact expectation without replacement;",
                                   "a curve still rising at its end point",
                                   "is undersampled")) +
    theme_pub()
  save_figure(p_rar, paste0(prefix, "b_rarefaction_curves"), run_dir, cfg,
              width = 6.5, height = 4,
              legend = paste("Expected number of clonotypes against",
                             "sequencing depth, computed exactly rather than",
                             "by resampling. Richness is a lower bound for",
                             "any sample whose curve has not flattened."))

  # D50: the headline number, kept beside the curves it has to be read with.
  d50 <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(clono, .data$sample_id, .data$subject_id),
    depth = sum(.data$abundance),
    n_clonotypes = dplyr::n(),
    d50_pct = d50_index(.data$abundance),
    .groups = "drop"))
  d50$d50_clonotypes <- round(d50$d50_pct / 100 * d50$n_clonotypes)
  report$h2("D50")
  report$tbl(as.data.frame(d50), digits = 2)
  save_table(d50, paste0(prefix, "_d50"), run_dir,
             "clonotypes accounting for half the repertoire, per sample")
  report$p("")
  report$p(paste0(
    "D50 is reported as a percentage of observed clonotypes and as a count. ",
    "A low D50 means few clonotypes carry half the repertoire, which is ",
    "high clonality. Both forms depend on depth, so read them against the ",
    "rarefaction curves above."))
  invisible(hill)
}

#' CDR3 structure: spectratype, physicochemistry and convergent coding
#'
#' @param db AIRR data frame, already restricted to the locus of interest.
#' @param abundance_col Abundance column name.
#' @param locus Locus label, used in titles.
#' @param run_dir,cfg,report Run directory, configuration and report writer.
#' @param group_cols Grouping columns for the spectratype and the summary.
#' @param unit Unit of observation, used in labels.
#' @param prefix File-name prefix.
#' @return The property summary, invisibly.
#' @export
section_cdr3_structure <- function(db, abundance_col, locus, run_dir, cfg,
                                   report,
                                   group_cols = c("subject_id", "timepoint"),
                                   unit = "molecules", prefix = "05c") {
  need("dplyr", "section_cdr3_structure()")
  report$h1("CDR3 structure and physicochemistry")

  db_locus <- filter_locus_productive(db, locus)
  group_cols <- intersect(group_cols, names(db_locus))

  # --- spectratype ---------------------------------------------------------
  spectra <- cdr3_spectratype(db_locus, abundance_col = abundance_col,
                              group_cols = group_cols)
  save_table(spectra, paste0(prefix, "_cdr3_spectratype"), run_dir,
             sprintf("CDR3 length distribution weighted by %s", unit))

  p_spec <- ggplot2::ggplot(spectra, ggplot2::aes(.data$cdr3_length,
                                                  .data$frequency)) +
    ggplot2::geom_col(fill = "grey60", width = 0.8) +
    ggplot2::labs(x = "CDR3 length (amino acids)",
                  y = sprintf("Frequency (%s)", unit),
                  title = sprintf("%s CDR3 spectratype",
                                  paste(locus, collapse = "/")),
                  subtitle = paste("A near-Gaussian bell is expected; a",
                                   "sharp peak at one length is an expansion",
                                   "large enough to deform the whole",
                                   "distribution")) +
    theme_pub()
  if (length(group_cols) >= 2L) {
    p_spec <- p_spec + ggplot2::facet_grid(
      stats::as.formula(sprintf("%s ~ %s", group_cols[1], group_cols[2])))
  } else if (length(group_cols) == 1L) {
    p_spec <- p_spec + ggplot2::facet_wrap(
      stats::as.formula(paste("~", group_cols[1])))
  }
  save_figure(p_spec, paste0(prefix, "a_cdr3_spectratype"), run_dir, cfg,
              width = 11, height = 7,
              legend = sprintf(paste0("CDR3 amino acid length distribution, ",
                                      "weighted by %s, within each group."),
                               unit))

  # --- physicochemistry ----------------------------------------------------
  props <- try(cdr3_properties(db_locus), silent = TRUE)
  if (inherits(props, "try-error")) {
    report$p("")
    report$p(paste0(
      "CDR3 physicochemical properties were not computed: ",
      conditionMessage(attr(props, "condition"))))
  } else {
    summ <- cdr3_property_summary(props, abundance_col = abundance_col,
                                  group_cols = group_cols)
    report$h2("Abundance-weighted CDR3 properties")
    report$tbl(as.data.frame(summ), digits = 2)
    save_table(summ, paste0(prefix, "_cdr3_properties_summary"), run_dir,
               sprintf("mean CDR3 length, charge and GRAVY weighted by %s",
                       unit))
    report$p("")
    report$p(paste0(
      "Net charge and hydrophobicity (GRAVY) of the CDR3, weighted by ",
      unit, ". Selected clonotypes are not a random draw from the ",
      "repertoire, so a shift in these properties over time, within a ",
      "subject, is a fingerprint of selection rather than of recombination. ",
      "Properties are computed over the whole junction, conserved first and ",
      "last residues included."))

    # The violins are drawn from a frame with fixed column names, so the
    # mapping below is literal: an aes() built from variable column names
    # leaves ggplot to infer the grouping, and a mis-inferred group silently
    # merges the panels instead of failing.
    abund <- as.numeric(props[[abundance_col]])
    x_col <- if (length(group_cols)) group_cols[length(group_cols)] else NULL
    for (prop in c("cdr3_aa_charge", "cdr3_aa_gravy")) {
      if (!prop %in% names(props)) next
      dat <- data.frame(
        value = props[[prop]],
        weight = abund,
        x = if (is.null(x_col)) "all" else as.character(props[[x_col]]),
        facet = if (length(group_cols) >= 2L)
          as.character(props[[group_cols[1]]]) else "all",
        stringsAsFactors = FALSE)

      p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$x, y = .data$value,
                                             weight = .data$weight,
                                             group = .data$x,
                                             fill = .data$x)) +
        ggplot2::geom_violin(alpha = 0.7, linewidth = 0.3) +
        ggplot2::scale_fill_manual(
          values = palette_sequential(length(unique(dat$x)), direction = -1),
          guide = "none") +
        ggplot2::labs(x = NULL,
                      y = if (prop == "cdr3_aa_charge") "Net charge"
                          else "GRAVY (hydrophobicity)",
                      title = sprintf("CDR3 %s",
                                      if (prop == "cdr3_aa_charge")
                                        "net charge" else "hydrophobicity"),
                      subtitle = sprintf("Weighted by abundance in %s",
                                         unit)) +
        theme_pub() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                           hjust = 1,
                                                           size = 7))
      if (length(group_cols) >= 2L) {
        p <- p + ggplot2::facet_wrap(~facet, nrow = 1)
      }
      save_weighted_figure(p, paste0(prefix, "b_", prop), run_dir, cfg,
                           width = 12, height = 4,
                           legend = sprintf(paste0("Distribution of CDR3 %s, ",
                                                   "each rearrangement ",
                                                   "weighted by its abundance ",
                                                   "in %s."), prop, unit))
    }
  }

  # --- convergent coding ---------------------------------------------------
  conv <- nucleotide_convergence(db_locus)
  if (nrow(conv)) {
    report$h2("Convergent coding within sample")
    report$tbl(as.data.frame(conv), digits = 2)
    save_table(conv, paste0(prefix, "_nucleotide_convergence"), run_dir,
               "percentage of CDR3 amino acid sequences encoded by more than one junction")
    report$p("")
    report$p(paste0(
      "Percentage of distinct CDR3 amino acid sequences encoded by more ",
      "than one nucleotide junction in the same sample. Independent ",
      "recombinations converging on one amino acid sequence is evidence of ",
      "selection on that sequence. The index rises with depth, so compare ",
      "samples of similar depth."))

    if (all(c("timepoint", "subject_id") %in% names(conv))) {
      p_conv <- ggplot2::ggplot(conv,
                                ggplot2::aes(.data$timepoint,
                                             .data$pct_convergent,
                                             group = .data$subject_id,
                                             colour = .data$subject_id)) +
        ggplot2::geom_line(linewidth = 0.7, alpha = 0.85) +
        ggplot2::geom_point(size = 2) +
        ggplot2::labs(x = NULL, y = "Convergent CDR3 (%)", colour = "Subject",
                      title = "Convergent coding of the same CDR3",
                      subtitle = "% of CDR3 amino acid sequences encoded by more than one junction") +
        theme_pub()
      save_figure(p_conv, paste0(prefix, "c_nucleotide_convergence"), run_dir,
                  cfg, width = 7, height = 4,
                  legend = paste("Percentage of CDR3 amino acid sequences",
                                 "encoded by more than one nucleotide",
                                 "junction, one line per subject."))
    }
  }
  invisible(conv)
}

#' CDR3 features of the expanded clonotypes against the rest
#'
#' Length, charge and sequence similarity of a marked set of clonotypes
#' compared with the remainder of the same subject's repertoire. The
#' comparison is within subject for the reason given at the top of this
#' file, and the marked set is supplied by the caller rather than recomputed
#' here, so that the figure and the differential-abundance table can never
#' disagree about which clonotypes are expanded.
#'
#' @param db AIRR data frame, restricted to the locus of interest, with a
#'   clonotype key column.
#' @param abundance_col Abundance column name.
#' @param focus_clonotypes Clonotype keys of the marked set (for example the
#'   significantly expanded clonotypes).
#' @param clonotype_col Column in `db` holding the clonotype key.
#' @param run_dir,cfg,report Run directory, configuration and report writer.
#' @param unit Unit of observation, used in labels.
#' @param prefix File-name prefix.
#' @return The per-subject network summary, invisibly.
#' @export
section_cdr3_features_expanded <- function(db, abundance_col,
                                           focus_clonotypes, clonotype_col,
                                           run_dir, cfg, report,
                                           unit = "molecules",
                                           prefix = "07c") {
  need("dplyr", "section_cdr3_features_expanded()")
  report$h1("CDR3 features of the expanded clonotypes")

  if (!length(focus_clonotypes)) {
    report$p(paste0(
      "No clonotype was marked as expanded, so there is nothing to compare ",
      "with the rest of the repertoire. This is a result, not a failure: a ",
      "response can be polyclonal with no clonotype passing both the ",
      "fold-change and the FDR criterion."))
    return(invisible(NULL))
  }
  if (!clonotype_col %in% names(db)) {
    stop("section_cdr3_features_expanded() needs the clonotype key column '",
         clonotype_col, "' in the input.", call. = FALSE)
  }

  db <- as.data.frame(db)
  db$.abund <- as.numeric(db[[abundance_col]])
  db$.focus <- as.character(db[[clonotype_col]]) %in%
    as.character(focus_clonotypes)
  db$group <- ifelse(db$.focus, "Expanded", "Other")
  db$cdr3_length <- nchar(as.character(db$junction_aa))

  group_cols <- intersect("subject_id", names(db))
  cols <- palette_categorical(2)
  names(cols) <- c("Expanded", "Other")

  report$p(sprintf(paste0(
    "%s clonotypes are compared with the rest of the same subject's ",
    "repertoire. Distributions are weighted by %s, so they describe the ",
    "molecules present rather than the distinct clonotypes; the kernel ",
    "bandwidth of the density is chosen from the unweighted sample, which ",
    "affects only how smooth the curve is and not where its mass sits."),
    format(length(unique(focus_clonotypes)), big.mark = " "), unit))
  report$p("")

  # --- length --------------------------------------------------------------
  p_len <- ggplot2::ggplot(db, ggplot2::aes(.data$cdr3_length,
                                            weight = .data$.abund,
                                            colour = .data$group)) +
    ggplot2::geom_density(linewidth = 0.7) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::labs(x = "CDR3 length (amino acids)", y = "Density",
                  title = "CDR3 length: expanded clonotypes against the rest") +
    theme_pub()
  if (length(group_cols)) {
    p_len <- p_len + ggplot2::facet_wrap(~subject_id, nrow = 1)
  }
  save_weighted_figure(p_len, paste0(prefix, "a_cdr3_length_expanded"),
                       run_dir, cfg, width = 12, height = 4,
                       legend = sprintf(paste0("CDR3 length density, weighted ",
                                               "by %s, for the expanded ",
                                               "clonotypes and for the rest ",
                                               "of each subject's ",
                                               "repertoire."), unit))

  # --- charge --------------------------------------------------------------
  props <- try(cdr3_properties(db), silent = TRUE)
  if (!inherits(props, "try-error") && "cdr3_aa_charge" %in% names(props)) {
    p_chg <- ggplot2::ggplot(props, ggplot2::aes(.data$group,
                                                 .data$cdr3_aa_charge,
                                                 weight = .data$.abund,
                                                 fill = .data$group)) +
      ggplot2::geom_violin(alpha = 0.7, linewidth = 0.3) +
      ggplot2::scale_fill_manual(values = cols, guide = "none") +
      ggplot2::labs(x = NULL, y = "Net charge",
                    title = "CDR3 net charge: expanded clonotypes against the rest") +
      theme_pub()
    if (length(group_cols)) {
      p_chg <- p_chg + ggplot2::facet_wrap(~subject_id, nrow = 1)
    }
    save_weighted_figure(p_chg, paste0(prefix, "b_cdr3_charge_expanded"),
                         run_dir, cfg, width = 12, height = 4,
                         legend = sprintf(paste0("CDR3 net charge, weighted ",
                                                 "by %s, for the expanded ",
                                                 "clonotypes and for the rest ",
                                                 "of each subject's ",
                                                 "repertoire."), unit))
  } else {
    report$p(paste0(
      "CDR3 charge was not computed; alakazam is required for the ",
      "physicochemical properties."))
  }

  # --- similarity network --------------------------------------------------
  net_cfg <- cfg$repertoire$cdr3_network %||% list()
  frac      <- as.numeric(net_cfg$max_distance_fraction %||% 0.20)
  min_len   <- as.integer(net_cfg$min_length %||% 8L)
  max_seq   <- as.integer(net_cfg$max_sequences %||% 400L)
  diag_frac <- as.numeric(net_cfg$diagnostic_fractions %||%
                            c(0.10, 0.15, 0.20))

  subjects <- if (length(group_cols)) unique(as.character(db$subject_id)) else
    "all subjects"
  panels <- list(); summaries <- list(); diagnostics <- list()

  for (s in subjects) {
    sub <- if (length(group_cols))
      db[as.character(db$subject_id) == s, , drop = FALSE] else db
    agg <- dplyr::ungroup(dplyr::summarise(
      dplyr::group_by(sub, .data$junction_aa),
      abundance = sum(.data$.abund, na.rm = TRUE),
      focus = any(.data$.focus), .groups = "drop"))

    net <- cdr3_similarity_network(
      agg$junction_aa, abundance = agg$abundance, focus = agg$focus,
      max_distance_fraction = frac, min_length = min_len,
      max_sequences = max_seq, layout_seed = cfg$project$seed %||% 1234L)

    diagnostics[[length(diagnostics) + 1L]] <- data.frame(
      subject_id = s,
      cdr3_network_diagnostic(agg$junction_aa, fractions = diag_frac,
                              min_length = min_len, max_sequences = max_seq),
      stringsAsFactors = FALSE)

    summaries[[length(summaries) + 1L]] <- data.frame(
      subject_id = s, n_candidate_cdr3 = net$n_candidates,
      n_nodes = net$n_nodes, n_edges = net$n_edges,
      stringsAsFactors = FALSE)

    panels[[s]] <- if (net$n_nodes) {
      net$nodes$status <- ifelse(net$nodes$focus, "Expanded", "Other")
      ggplot2::ggplot() +
        ggplot2::geom_segment(
          data = net$edges,
          ggplot2::aes(x = .data$x_from, y = .data$y_from,
                       xend = .data$x_to, yend = .data$y_to),
          colour = "grey80", linewidth = 0.2) +
        ggplot2::geom_point(data = net$nodes,
                            ggplot2::aes(.data$x, .data$y,
                                         colour = .data$status), size = 1) +
        ggplot2::scale_colour_manual(values = cols, name = NULL) +
        ggplot2::labs(title = sprintf("%s | n = %d (distance <= %s of length)",
                                      s, net$n_nodes,
                                      scales::percent(frac))) +
        ggplot2::theme_void(base_size = 8) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                          size = 8),
                       legend.position = "none")
    } else {
      # No pair met the threshold. The panel says so rather than being
      # omitted: an absent panel reads as a failed step, where "no edges at
      # this threshold" is a finding about the repertoire.
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0,
                          label = sprintf("%s\n(no edge at %s)", s,
                                          scales::percent(frac)),
                          size = 2.5) +
        ggplot2::theme_void()
    }
  }

  net_summary <- dplyr::bind_rows(summaries)
  diag_table  <- dplyr::bind_rows(diagnostics)
  save_table(net_summary, paste0(prefix, "_cdr3_network_summary"), run_dir,
             "nodes and edges of the CDR3 similarity network, per subject")
  save_table(diag_table, paste0(prefix, "_cdr3_network_thresholds"), run_dir,
             "edges obtained at each similarity threshold")

  report$h2("CDR3 similarity network")
  report$tbl(as.data.frame(net_summary))
  report$p("")
  report$p(sprintf(paste0(
    "Edges link CDR3s whose Levenshtein distance is at most %s of the ",
    "longer sequence, excluding identical sequences. Normalising by length ",
    "is what makes one threshold comparable across short and long CDR3s: an ",
    "absolute distance of 1 is a much stronger statement about two long ",
    "CDR3s than about two short ones, and on real repertoires it frequently ",
    "yields no edges at all."), scales::percent(frac)))
  report$p("")
  report$h2("Edges at each candidate threshold")
  report$tbl(as.data.frame(diag_table))
  report$p("")
  report$p(paste0(
    "The edge count is highly sensitive to the threshold, so the table ",
    "above is the evidence for the value used. Report it in Methods; a ",
    "network figure at a single threshold does not let a reader judge it."))

  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_net <- patchwork::wrap_plots(panels, nrow = 1)
  } else {
    p_net <- panels[[1]]
    report$p("")
    report$p(paste0(
      "patchwork is not installed, so only the first subject's network was ",
      "written as a figure; the per-subject counts are in the table."))
  }
  save_figure(p_net, paste0(prefix, "c_cdr3_network"), run_dir, cfg,
              width = max(6, 3 * length(panels)), height = 3.4,
              legend = paste("CDR3 similarity network per subject. Nodes are",
                             "distinct CDR3 amino acid sequences, coloured by",
                             "whether the clonotype was called expanded;",
                             "edges link sequences within the configured",
                             "normalised edit distance. Isolated nodes are",
                             "not drawn."))
  invisible(net_summary)
}

#' Sequence logos of the public CDR3 groups
#'
#' @param public_groups The `groups` element of [detect_public_groups()].
#' @param run_dir,cfg,report Run directory, configuration and report writer.
#' @param group_col,cdr3_col,label_col Column names.
#' @param prefix File-name prefix.
#' @return The logo result, invisibly.
#' @export
section_cdr3_logos <- function(public_groups, run_dir, cfg, report,
                               group_col = "public_group",
                               cdr3_col = "junction_aa",
                               label_col = "v_gene", prefix = "07b") {
  report$h1("CDR3 motifs of the public groups")

  if (is.null(public_groups) || !nrow(public_groups)) {
    report$p("No public group was detected, so no motif logo is drawn.")
    return(invisible(NULL))
  }
  if (!requireNamespace("ggseqlogo", quietly = TRUE)) {
    report$p(paste0(
      "ggseqlogo is not installed, so the motif logos were skipped. The ",
      "public group membership table carries the sequences."))
    return(invisible(NULL))
  }

  logo_cfg <- cfg$convergence$logos %||% list()
  res <- cdr3_motif_logos(
    as.data.frame(public_groups), group_col = group_col, cdr3_col = cdr3_col,
    label_col = label_col,
    min_sequences = as.integer(logo_cfg$min_sequences %||% 4L),
    max_groups = as.integer(logo_cfg$max_groups %||% 9L),
    align = logo_cfg$align %||% "auto")

  if (is.null(res$plot)) {
    report$p(sprintf(paste0(
      "No public group had at least %d distinct CDR3 sequences, which is ",
      "the minimum for a logo to be readable."),
      as.integer(logo_cfg$min_sequences %||% 4L)))
    return(invisible(res))
  }

  save_figure(res$plot, paste0(prefix, "_cdr3_logos"), run_dir, cfg,
              width = 11, height = 7,
              legend = paste("CDR3 amino acid logos of the public groups.",
                             "Conserved positions are the motif shared",
                             "between subjects."))
  report$p(sprintf("Panels drawn: %d.", length(res$groups)))
  report$p("")
  if (identical(res$alignment, "msa")) {
    report$p(paste0(
      "Sequences within a group were aligned with DECIPHER before the logo ",
      "was drawn, so positions are comparable across sequences of different ",
      "length."))
  } else {
    report$p(paste0(
      "DECIPHER is not installed, so no multiple alignment was computed. ",
      "One logo is drawn per CDR3 length within each group instead: every ",
      "panel is then over sequences of equal length, which is comparable by ",
      "construction and invents no gap, at the cost of splitting a group ",
      "across panels. Install DECIPHER for the aligned view."))
  }
  invisible(res)
}
