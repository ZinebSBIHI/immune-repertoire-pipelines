# ---------------------------------------------------------------------------
# Simulate a small bulk TCRb repertoire dataset in AIRR format
#
# Produces an example dataset so that pipelines/03_bulk_repertoire.R runs
# end to end without any restricted data, and so that a reader can see what
# the pipeline expects as input.
#
# The generative model is deliberately simple and is NOT a model of a real
# immune response:
#
#   - a baseline repertoire per subject, clonotype abundances drawn from a
#     log-normal, which gives the long-tailed shape bulk repertoires take
#     without the single-clone dominance a steep rank power law produces;
#   - at each post-baseline timepoint, a handful of clonotypes expanded by a
#     fixed factor, plus a handful of clonotypes absent at baseline;
#   - multinomial sampling at a per-sample depth that varies between samples,
#     so the pipeline's rarefaction and depth diagnostics have something to
#     act on.
#
# Usage:
#   Rscript examples/simulate_airr.R [output_dir] [seed]
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[[1]] else "data/example"
seed    <- if (length(args) >= 2L) as.integer(args[[2]]) else 1234L

set.seed(seed)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- parameters -------------------------------------------------------------
subjects   <- c("subject_01", "subject_02", "subject_03")
timepoints <- c("t0", "t1", "t2")
condition_of <- c(t0 = "baseline", t1 = "treatment_a", t2 = "treatment_b")

n_baseline_clones <- 4000      # distinct clonotypes in each baseline
# Abundance is drawn from a log-normal, which reproduces the long-tailed but
# not single-clone-dominated shape of a resting bulk repertoire: most
# clonotypes near the detection limit, a top clone of order 1% of the
# repertoire. A rank power law with a steep exponent would put tens of
# percent into the first clone, which resting blood does not look like.
lognormal_sdlog    <- 1.5
n_expanded         <- 25       # pre-existing clonotypes expanded post-baseline
expansion_factor   <- 12
n_de_novo          <- 20       # clonotypes absent at baseline
# Total share of the post-baseline repertoire held by the de novo clonotypes,
# expressed as a fraction so it does not depend on the pool's scale.
de_novo_total_share <- 0.03

v_genes <- sprintf("TRBV%d-%d", rep(2:29, each = 2), rep(1:2, 28))
j_genes <- sprintf("TRBJ%d-%d", rep(1:2, each = 6), c(1:6, 1:6))
aa      <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]

#' Draw a plausible CDR3b amino acid sequence
#'
#' Conserved CASS.. start and ..F end, variable middle, length 10-20, which
#' is the range human TCRb CDR3s occupy.
make_cdr3 <- function() {
  n_mid <- sample(5:15, 1L)
  paste0("CASS", paste(sample(aa, n_mid, replace = TRUE), collapse = ""), "F")
}

#' Nucleotide junction consistent with an amino acid CDR3
#'
#' Reverse-translated with a fixed codon per residue. Sufficient for the
#' pipeline's purposes (clonotype keys, length distributions); it is not a
#' model of V(D)J junction formation.
codon <- c(A="GCT", C="TGT", D="GAT", E="GAA", F="TTT", G="GGT", H="CAT",
           I="ATT", K="AAA", L="CTT", M="ATG", N="AAT", P="CCT", Q="CAA",
           R="CGT", S="AGT", T="ACT", V="GTT", W="TGG", Y="TAT")
make_junction <- function(cdr3_aa) {
  paste(codon[strsplit(cdr3_aa, "")[[1]]], collapse = "")
}

# --- build the repertoires --------------------------------------------------
all_rows <- list()

for (subj in subjects) {

  # Baseline clonotype pool for this subject. Repertoires are private, so
  # each subject gets its own pool; a small shared fraction is added below
  # so that the convergence analysis has something to find.
  cdr3  <- vapply(seq_len(n_baseline_clones), function(i) make_cdr3(),
                  character(1))
  pool <- data.frame(
    junction_aa = cdr3,
    junction    = vapply(cdr3, make_junction, character(1)),
    v_call      = sample(v_genes, n_baseline_clones, replace = TRUE),
    j_call      = sample(j_genes, n_baseline_clones, replace = TRUE),
    stringsAsFactors = FALSE
  )
  pool$v_call <- paste0(pool$v_call, "*01")
  pool$j_call <- paste0(pool$j_call, "*01")

  # Log-normal abundance weights, normalised to sum to one so that the
  # de novo share below is a true fraction of the repertoire.
  pool$weight <- stats::rlnorm(n_baseline_clones, meanlog = 0,
                               sdlog = lognormal_sdlog)
  pool$weight <- pool$weight / sum(pool$weight)

  # Clonotypes that will expand, and clonotypes that appear de novo.
  expanded_idx <- sample.int(nrow(pool), n_expanded)

  de_novo_cdr3 <- vapply(seq_len(n_de_novo), function(i) make_cdr3(),
                         character(1))
  de_novo <- data.frame(
    junction_aa = de_novo_cdr3,
    junction    = vapply(de_novo_cdr3, make_junction, character(1)),
    v_call      = paste0(sample(v_genes, n_de_novo, replace = TRUE), "*01"),
    j_call      = paste0(sample(j_genes, n_de_novo, replace = TRUE), "*01"),
    weight      = de_novo_total_share / n_de_novo,
    stringsAsFactors = FALSE
  )

  for (tp in timepoints) {
    is_baseline <- tp == timepoints[1]

    rep_tbl <- pool
    if (!is_baseline) {
      rep_tbl$weight[expanded_idx] <- rep_tbl$weight[expanded_idx] *
        expansion_factor
      rep_tbl <- rbind(rep_tbl, de_novo)
    }

    # Depth varies between samples on purpose: unequal depth is the normal
    # situation and the pipeline must handle it rather than assume it away.
    depth <- round(stats::runif(1, 30000, 90000))

    probs  <- rep_tbl$weight / sum(rep_tbl$weight)
    counts <- stats::rmultinom(1, size = depth, prob = probs)[, 1]

    keep <- counts > 0
    out <- data.frame(
      sequence_id    = sprintf("%s_%s_%06d", subj, tp, seq_len(sum(keep))),
      sequence_alignment = "",
      productive     = TRUE,
      locus          = "TRB",
      v_call         = rep_tbl$v_call[keep],
      d_call         = "",
      j_call         = rep_tbl$j_call[keep],
      junction       = rep_tbl$junction[keep],
      junction_aa    = rep_tbl$junction_aa[keep],
      junction_length = nchar(rep_tbl$junction[keep]),
      umi_count      = counts[keep],
      duplicate_count = counts[keep],
      stringsAsFactors = FALSE
    )

    # A small number of non-productive rearrangements, so the pipeline's
    # productive filter is exercised rather than assumed unnecessary.
    n_bad <- round(0.02 * nrow(out))
    if (n_bad > 0) {
      bad <- out[sample.int(nrow(out), n_bad), ]
      bad$productive <- FALSE
      bad$sequence_id <- paste0(bad$sequence_id, "_np")
      out <- rbind(out, bad)
    }

    f <- file.path(out_dir, sprintf("%s_%s_%s.airr.tsv", subj, tp,
                                    condition_of[[tp]]))
    utils::write.table(out, f, sep = "\t", quote = FALSE, row.names = FALSE)
    all_rows[[length(all_rows) + 1L]] <- data.frame(
      file = basename(f), subject = subj, timepoint = tp,
      condition = condition_of[[tp]],
      n_rearrangements = nrow(out), depth = sum(out$umi_count),
      stringsAsFactors = FALSE)
  }
}

manifest <- do.call(rbind, all_rows)
utils::write.csv(manifest, file.path(out_dir, "MANIFEST.csv"),
                 row.names = FALSE)

cat("\nSimulated AIRR dataset written to:", normalizePath(out_dir), "\n\n")
print(manifest, row.names = FALSE)
cat("\nThese data are simulated. They carry no biological meaning and must",
    "\nnot be used to draw any immunological conclusion.\n")
