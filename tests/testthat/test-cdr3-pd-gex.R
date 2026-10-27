# ---------------------------------------------------------------------------
# Tests for rep-cdr3.R, ana-pharmacodynamics.R and gex-markers.R
#
# The repertoire fixtures are built here rather than read from a file, so
# that every expectation below is checked against a repertoire whose answer
# is known by construction: which clonotypes are present at baseline, which
# appear afterwards, which pair of CDR3s differs by one substitution.
#
# The transcriptome fixture is SeuratObject::pbmc_small, which ships with
# the package and is small enough to test on. It is 230 genes, so it also
# exercises the case a full dataset never reaches: an object too small for
# AddModuleScore's default number of control genes.
# ---------------------------------------------------------------------------

# --- fixtures --------------------------------------------------------------

#' A deterministic toy repertoire: 3 subjects x 3 timepoints
#'
#' Baseline holds clonotypes c001..c040, with c040 a singleton so that the
#' detection floor has something to act on. At each later timepoint:
#'   - c001 expands, by a subject-specific factor so that the subjects are
#'     not identical and a rank correlation over them is defined;
#'   - new_a and new_b appear, absent from the baseline;
#'   - c040 rises above the floor, which is the case that tests the
#'     documented consequence of the shared floor: a baseline singleton is
#'     not in the baseline, so it counts as de novo when it reappears.
toy_clono <- function() {
  subjects <- c("s1", "s2", "s3")
  timepoints <- c("t0", "t1", "t2")
  factors <- c(s1 = 10, s2 = 8, s3 = 6)
  rows <- list()

  for (s in subjects) {
    for (tp in timepoints) {
      clonotypes <- sprintf("c%03d", 1:40)
      abundance <- c(50, rep(10, 19), rep(4, 19), 1)
      if (tp != "t0") {
        abundance[1] <- 50 * factors[[s]]              # c001 expands
        abundance[40] <- 20                            # the singleton returns
        clonotypes <- c(clonotypes, "new_a", "new_b")  # absent at baseline
        abundance   <- c(abundance, 30, 3)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        sample_id = paste(s, tp, sep = "_"), subject_id = s, timepoint = tp,
        clonotype = clonotypes, abundance = abundance,
        frequency = abundance / sum(abundance), stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

#' A toy AIRR rearrangement table with CDR3s whose relationships are known
#'
#' CASSLAPGATNEQFF and CASSLAPGATNEQFY differ by one substitution;
#' CAWWWWWWWWWWWWF is unrelated to both; CASSQDR is below the length floor.
#' The first CDR3 is encoded by two different junctions, so the convergent
#' coding index has a known value.
toy_airr <- function() {
  data.frame(
    sample_id = rep(c("s1_t0", "s1_t1"), each = 4),
    subject_id = "s1",
    timepoint = rep(c("t0", "t1"), each = 4),
    productive = TRUE,
    locus = "TRB",
    v_call = rep(c("TRBV20-1*01", "TRBV20-1*02", "TRBV5-1*01",
                   "TRBV5-1*01"), 2),
    j_call = "TRBJ2-1*01",
    junction_aa = rep(c("CASSLAPGATNEQFF", "CASSLAPGATNEQFF",
                        "CASSLAPGATNEQFY", "CAWWWWWWWWWWWWF"), 2),
    junction = rep(c("TGTGCTAGTAGTTTA", "TGTGCTAGCAGCTTA", "TGTGCTAGTAGTTAT",
                     "TGTGCTTGGTGGTTT"), 2),
    umi_count = c(100, 50, 10, 5, 200, 50, 10, 5),
    stringsAsFactors = FALSE)
}

#' A configuration and a fresh run directory for the section tests
toy_run <- function(extra = list()) {
  data_dir <- withr::local_tempdir(.local_envir = parent.frame())
  out_dir  <- withr::local_tempdir(.local_envir = parent.frame())
  base <- list(
    project = list(name = "test", seed = 1L),
    paths = list(data_dir = data_dir, output_dir = out_dir),
    cohort = list(subjects = c("s1", "s2", "s3"),
                  # t9 is declared but never present in the data, which is
                  # what the "no baseline in the data" test needs.
                  timepoints = c("t0", "t1", "t2", "t9"),
                  conditions = c("g1", "g2")),
    repertoire = list(locus = "TRB", abundance_column = "umi_count"),
    differential_abundance = list(reference_timepoint = "t0",
                                  fold_change_threshold = 2),
    output = list(formats = "png", dpi = 72))
  cfg <- validate_config(utils::modifyList(base, extra))
  run_dir <- init_run(cfg, "test")
  list(cfg = cfg, run_dir = run_dir,
       report = open_report(file.path(run_dir, "REPORT.txt"), "test"))
}

run_files <- function(run_dir) list.files(run_dir, recursive = TRUE)

# ===========================================================================
# Diversity spectra
# ===========================================================================

test_that("Hill numbers reproduce the indices they generalise", {
  set.seed(1)
  ab <- c(rpois(200, 5) + 1, 500, 300)

  h <- hill_numbers(ab, q = c(0, 1, 2))
  expect_equal(h$hill[h$q == 0], sum(ab > 0))
  expect_equal(h$hill[h$q == 1], exp(shannon_entropy(ab)))
  expect_equal(h$hill[h$q == 2], inverse_simpson(ab))
})

test_that("the Hill profile decreases with q and survives degenerate input", {
  set.seed(2)
  h <- hill_numbers(c(rpois(100, 3) + 1, 400), q = 0:6)
  expect_true(all(diff(h$hill) <= 1e-8))

  expect_true(is.na(hill_numbers(numeric(0), q = 1)$hill))
  # One clonotype: every order says "one effective clonotype".
  expect_equal(hill_numbers(10, q = c(0, 1, 2))$hill, c(1, 1, 1))
})

test_that("the rarefaction curve is exact at full depth and monotone", {
  set.seed(3)
  ab <- c(rpois(300, 4) + 1, 200)
  curve <- rarefaction_curve(ab, n_points = 8L)

  expect_true(nrow(curve) > 1L)
  expect_true(all(diff(curve$expected_richness) > 0))
  expect_true(all(curve$expected_richness <= sum(ab > 0) + 1e-8))
  # At the observed depth every clonotype is certain to be seen, so the
  # expectation must be the observed richness exactly, not approximately.
  expect_equal(curve$expected_richness[nrow(curve)], sum(ab > 0))
  expect_equal(curve$depth[nrow(curve)], sum(ab))
})

test_that("rarefaction handles an evenly spread repertoire and empty input", {
  # Ten clonotypes of one molecule each: a subsample of n molecules sees
  # exactly n clonotypes, with no randomness left.
  curve <- rarefaction_curve(rep(1, 10), depths = c(1, 5, 10))
  expect_equal(curve$expected_richness, c(1, 5, 10))
  expect_equal(nrow(rarefaction_curve(numeric(0))), 0L)
})

test_that("the diversity spectrum section writes its tables and figures", {
  r <- toy_run(list(repertoire = list(locus = "TRB",
                                      abundance_column = "umi_count",
                                      hill_q = 0:3,
                                      rarefaction = list(enabled = TRUE,
                                                         curve_points = 5))))
  out <- section_diversity_spectrum(toy_clono(), r$run_dir, r$cfg, r$report)
  r$report$close()

  expect_s3_class(out, "data.frame")
  files <- run_files(r$run_dir)
  expect_true(any(grepl("hill_spectrum[.]csv$", files)))
  expect_true(any(grepl("rarefaction_curves[.]csv$", files)))
  expect_true(any(grepl("d50[.]csv$", files)))
  expect_true(any(grepl("hill_spectrum[.]png$", files)))
  expect_true(any(grepl("rarefaction_curves[.]png$", files)))
})

# ===========================================================================
# CDR3 descriptors
# ===========================================================================

test_that("the spectratype is a within-group frequency", {
  db <- toy_airr()
  sp <- cdr3_spectratype(db, "umi_count", group_cols = c("subject_id",
                                                         "timepoint"))
  totals <- tapply(sp$frequency, sp$timepoint, sum)
  expect_true(all(abs(totals - 1) < 1e-10))
  expect_equal(sum(sp$abundance), sum(db$umi_count))

  # Trimming the conserved residues shifts every length by exactly two.
  trimmed <- cdr3_spectratype(db, "umi_count", trim = TRUE)
  expect_equal(sort(unique(trimmed$cdr3_length)),
               sort(unique(sp$cdr3_length)) - 2L)
})

test_that("convergent coding counts CDR3s with more than one junction", {
  conv <- nucleotide_convergence(toy_airr())
  # Three distinct CDR3s per sample, one of which has two junctions.
  expect_equal(unique(conv$n_cdr3_aa), 3L)
  expect_equal(unique(conv$n_convergent), 1L)
  expect_equal(unique(round(conv$pct_convergent, 6)), round(100 / 3, 6))

  # One junction per CDR3 is zero convergence, not a missing value.
  single <- toy_airr()
  single$junction <- single$junction_aa
  expect_equal(unique(nucleotide_convergence(single)$pct_convergent), 0)
})

test_that("CDR3 properties are summarised weighted by abundance", {
  skip_if_not_installed("alakazam")
  db <- toy_airr()
  props <- cdr3_properties(db)
  expect_true(all(c("cdr3_aa_length", "cdr3_aa_charge", "cdr3_aa_gravy") %in%
                    names(props)))
  # The properties cover the whole junction: alakazam applies `trim` only
  # when translating nucleotide input, which the documentation of
  # cdr3_properties() states and this pins down, since trimming would shift
  # every value reported from these columns.
  expect_equal(unique(props$cdr3_aa_length),
               unique(nchar(db$junction_aa)))
  expect_equal(cdr3_properties(db, trim = FALSE)$cdr3_aa_charge,
               props$cdr3_aa_charge)

  summ <- cdr3_property_summary(props, "umi_count",
                                group_cols = c("subject_id", "timepoint"))
  expect_equal(nrow(summ), 2L)
  one <- props[props$timepoint == "t0", ]
  expect_equal(summ$cdr3_charge[summ$timepoint == "t0"],
               stats::weighted.mean(one$cdr3_aa_charge, one$umi_count))
  expect_equal(summ$abundance[summ$timepoint == "t0"], sum(one$umi_count))
})

test_that("cdr3_property_summary refuses to guess the property columns", {
  expect_error(cdr3_property_summary(toy_airr(), "umi_count"),
               "cdr3_properties")
})

# ===========================================================================
# CDR3 similarity network
# ===========================================================================

test_that("the network links near neighbours and ignores identical sequences", {
  skip_if_not_installed("igraph")
  seqs <- c("CASSLAPGATNEQFF", "CASSLAPGATNEQFF", "CASSLAPGATNEQFY",
            "CAWWWWWWWWWWWWF")
  net <- cdr3_similarity_network(seqs, abundance = c(10, 10, 5, 1),
                                 focus = c(TRUE, TRUE, FALSE, FALSE))

  # Duplicates collapse to one node, so the candidate count is 3, not 4.
  expect_equal(net$n_candidates, 3L)
  # One edge only: the single-substitution pair. The identical pair is
  # excluded by construction and the unrelated sequence is far away.
  expect_equal(net$n_edges, 1L)
  expect_equal(net$n_nodes, 2L)
  expect_setequal(net$nodes$sequence, c("CASSLAPGATNEQFF",
                                        "CASSLAPGATNEQFY"))
  expect_true(any(net$nodes$focus))
  # Abundance is summed over the collapsed duplicates.
  expect_equal(net$nodes$abundance[net$nodes$sequence == "CASSLAPGATNEQFF"],
               20)
})

test_that("an edgeless network is a result rather than an error", {
  skip_if_not_installed("igraph")
  far <- c("CASSLAPGATNEQFF", "CAWWWWWWWWWWWWF", "CFFFFFFFFFFFFFF")
  net <- cdr3_similarity_network(far)
  expect_equal(net$n_nodes, 0L)
  expect_equal(net$n_edges, 0L)
  expect_s3_class(net$nodes, "data.frame")
})

test_that("the network honours the length floor and the node cap", {
  skip_if_not_installed("igraph")
  short <- c("CASSF", "CASSY", "CASSLAPGATNEQFF", "CASSLAPGATNEQFY")
  net <- cdr3_similarity_network(short, min_length = 8L)
  expect_false(any(c("CASSF", "CASSY") %in% net$nodes$sequence))

  set.seed(4)
  many <- vapply(1:60, function(i) paste0("CASS",
                                          paste(sample(LETTERS[1:20], 10,
                                                       replace = TRUE),
                                                collapse = ""), "F"),
                 character(1))
  capped <- cdr3_similarity_network(c(many, "CASSLAPGATNEQFF",
                                      "CASSLAPGATNEQFY"),
                                    focus = c(rep(FALSE, 60), TRUE, TRUE),
                                    max_sequences = 5L)
  # The two focus sequences are kept whatever the cap, so their edge
  # survives a cap far below the number of candidates.
  expect_equal(capped$n_candidates, 62L)
  expect_true(all(c("CASSLAPGATNEQFF", "CASSLAPGATNEQFY") %in%
                    capped$nodes$sequence))
})

test_that("the layout is reproducible for a given seed", {
  skip_if_not_installed("igraph")
  seqs <- c("CASSLAPGATNEQFF", "CASSLAPGATNEQFY", "CASSLAPGATNEQFW",
            "CASSLAPGATNEQFL")
  a <- cdr3_similarity_network(seqs, layout_seed = 7L)
  b <- cdr3_similarity_network(seqs, layout_seed = 7L)
  expect_equal(a$nodes$x, b$nodes$x)
  expect_equal(a$nodes$y, b$nodes$y)
})

test_that("the threshold diagnostic reports more edges at looser thresholds", {
  set.seed(5)
  seqs <- unique(c("CASSLAPGATNEQFF", "CASSLAPGATNEQFY", "CASSLAPGATNEQFW",
                   "CASSLAPGQTNEQFF",
                   vapply(1:20, function(i)
                     paste0("CASS", paste(sample(LETTERS[1:20], 10,
                                                 replace = TRUE),
                                          collapse = ""), "F"),
                     character(1))))
  diag <- cdr3_network_diagnostic(seqs, fractions = c(0.05, 0.10, 0.20),
                                  absolute_distances = 1L)
  frac <- diag[diag$criterion == "fraction of longer length", ]
  frac <- frac[order(frac$threshold), ]
  expect_true(all(diff(frac$n_edges) >= 0))
  expect_true(any(diag$criterion == "absolute edit distance"))
  expect_true(all(diag$n_edges >= 0))
})

# ===========================================================================
# Sequence logos
# ===========================================================================

test_that("logos fall back to one panel per CDR3 length without DECIPHER", {
  skip_if_not_installed("ggseqlogo")
  groups <- data.frame(
    public_group = "PG0001",
    v_gene = "TRBV20-1",
    junction_aa = c("CASSLAPGATNEQFF", "CASSLAPGATNEQFY", "CASSLAPGATNEQFW",
                    "CASSLAPGATNEQFL", "CASSLAPGATNEQF", "CASSLAPGATNEQY"),
    stringsAsFactors = FALSE)

  res <- cdr3_motif_logos(groups, min_sequences = 2L, align = "by_length")
  expect_identical(res$alignment, "by_length")
  # Two lengths are present, so the group is split into two panels, each
  # over sequences of equal length.
  expect_equal(length(res$groups), 2L)
  expect_s3_class(res$plot, "ggplot")

  # A group below the minimum yields no plot and says so rather than erroring.
  thin <- cdr3_motif_logos(groups[1, , drop = FALSE], min_sequences = 4L,
                           align = "by_length")
  expect_null(thin$plot)
})

test_that("the logo section tolerates an empty set of public groups", {
  r <- toy_run()
  expect_null(section_cdr3_logos(NULL, r$run_dir, r$cfg, r$report))
  expect_null(section_cdr3_logos(data.frame(), r$run_dir, r$cfg, r$report))
  r$report$close()
})

# ===========================================================================
# Pharmacodynamics
# ===========================================================================

test_that("focusing metrics agree with the indices they are built from", {
  clono <- toy_clono()
  focal <- focalisation_metrics(clono)
  expect_equal(nrow(focal), 9L)

  one <- clono[clono$sample_id == "s1_t1", ]
  row <- focal[focal$sample_id == "s1_t1", ]
  expect_equal(row$total_abundance, sum(one$abundance))
  expect_equal(row$richness, nrow(one))
  expect_equal(row$shannon, shannon_entropy(one$abundance))
  expect_equal(row$shannon_clonality, shannon_clonality(one$abundance))
  expect_equal(row$gini, gini_index(one$abundance))

  # Both metrics must rise when one clonotype takes over, which is the
  # direction the whole section is read in.
  expect_gt(row$shannon_clonality,
            focal$shannon_clonality[focal$sample_id == "s1_t0"])
  expect_gt(row$gini, focal$gini[focal$sample_id == "s1_t0"])
})

test_that("the change from baseline is zero at baseline and signed after it", {
  focal <- focalisation_metrics(toy_clono())
  change <- focalisation_change(focal, reference_timepoint = "t0")

  expect_setequal(unique(change$metric), c("shannon_clonality", "gini"))
  expect_equal(nrow(change), 18L)
  expect_true(all(change$delta[change$timepoint == "t0"] == 0))
  expect_true(all(change$delta[change$timepoint != "t0"] > 0))
  # Each subject is compared with their own baseline, never another's.
  s1 <- change[change$subject_id == "s1" &
                 change$metric == "gini", ]
  expect_equal(unique(s1$baseline_value),
               focal$gini[focal$sample_id == "s1_t0"])
})

test_that("the baseline reference applies the detection floor", {
  clono <- toy_clono()
  loose <- baseline_reference(clono, "t0", min_count = 1)
  strict <- baseline_reference(clono, "t0", min_count = 2)

  expect_equal(nrow(loose), 120L)             # 40 clonotypes x 3 subjects
  expect_equal(nrow(strict), 117L)            # the singleton c040 drops out
  expect_true("c040" %in% loose$clonotype)
  expect_false("c040" %in% strict$clonotype)
  expect_setequal(names(strict), c("subject_id", "clonotype",
                                   "frequency_baseline",
                                   "abundance_baseline"))
})

test_that("novel clonotypes are those absent from the subject's baseline", {
  clono <- toy_clono()
  base <- baseline_reference(clono, "t0", min_count = 2)
  novel <- classify_novel_clonotypes(clono, base, "t0", min_count = 2,
                                     min_freq = 1e-5)

  # The baseline itself is never novel, whatever its contents.
  expect_false(any(novel$novel_expanded[novel$timepoint == "t0"]))
  called <- unique(novel$clonotype[novel$novel_expanded])
  # new_a and new_b are above both thresholds. c040 is a singleton at
  # baseline, so it is not in the baseline at min_count = 2 and its return
  # above the floor counts as de novo: the documented consequence of using
  # one floor for both timepoints.
  expect_true(all(c("new_a", "new_b") %in% called))
  expect_true("c040" %in% called)
  # c001 was present at baseline, so however much it expands it is never
  # novel. That is the line between this view and the expansion view.
  expect_false("c001" %in% called)
  # At a floor of 1 the singleton IS in the baseline, so it stops being
  # called de novo; nothing else changes.
  loose <- classify_novel_clonotypes(
    clono, baseline_reference(clono, "t0", min_count = 1), "t0",
    min_count = 1, min_freq = 1e-5)
  expect_false("c040" %in% loose$clonotype[loose$novel_expanded])
  expect_true(all(c("new_a", "new_b") %in%
                    loose$clonotype[loose$novel_expanded]))

  # A floor above a clonotype's post-baseline abundance excludes it.
  tight <- classify_novel_clonotypes(clono, base, "t0", min_count = 5,
                                     min_freq = 1e-5)
  expect_false("new_b" %in% tight$clonotype[tight$novel_expanded])
})

test_that("expanded clonotypes separate de novo from pre-existing", {
  clono <- toy_clono()
  base <- baseline_reference(clono, "t0", min_count = 2)
  exp <- classify_expanded_clonotypes(clono, base, "t0", min_count = 2,
                                      min_freq = 1e-5, min_fold = 2)

  expect_setequal(unique(exp$clonotype_type),
                  c("de novo", "pre-existing (expanded)", "not expanded"))
  # No fold change is formed for a clonotype with no baseline frequency.
  expect_true(all(is.na(exp$fold_change[!exp$present_baseline])))
  # c001 is expanded tenfold in abundance and rises in frequency, so it is
  # the pre-existing expanded clonotype; it is never called de novo.
  c001 <- exp[exp$clonotype == "c001" & exp$timepoint == "t1", ]
  expect_true(all(c001$clonotype_type == "pre-existing (expanded)"))
  expect_true(all(c001$fold_change >= 2))
  expect_true(all(exp$clonotype_type[exp$clonotype %in% c("new_a") &
                                       exp$timepoint != "t0"] == "de novo"))
  # The baseline is never expanded.
  expect_false(any(exp$is_expanded[exp$timepoint == "t0"]))

  # A fold-change threshold above the observed rise removes the
  # pre-existing calls without touching the de novo ones.
  strict <- classify_expanded_clonotypes(clono, base, "t0", min_fold = 1e6)
  expect_false(any(strict$clonotype_type == "pre-existing (expanded)"))
  expect_true(any(strict$clonotype_type == "de novo"))
})

test_that("magnitude and breadth are summarised separately", {
  clono <- toy_clono()
  base <- baseline_reference(clono, "t0", min_count = 2)
  novel <- summarise_novel_clonotypes(
    classify_novel_clonotypes(clono, base, "t0"))
  expect_equal(nrow(novel), 9L)
  expect_true(all(novel$summed_frequency_novel[novel$timepoint == "t0"] == 0))
  expect_true(all(novel$n_novel[novel$timepoint == "t0"] == 0L))
  expect_true(all(novel$summed_frequency_novel[novel$timepoint != "t0"] > 0))

  ex <- summarise_expanded_clonotypes(
    classify_expanded_clonotypes(clono, base, "t0", min_fold = 2))
  expect_named(ex, c("summary", "by_type"))
  # The two components add up to the total, which is what makes the
  # composition panel readable.
  expect_equal(ex$summary$n_expanded,
               ex$summary$n_de_novo + ex$summary$n_pre_existing)
  expect_setequal(unique(ex$by_type$clonotype_type),
                  c("de novo", "pre-existing (expanded)"))
})

test_that("the depth control reports both correlations and their agreement", {
  focal <- focalisation_metrics(toy_clono())
  change <- focalisation_change(focal, "t0")
  cors <- focalisation_depth_correlations(focal, change,
                                          reference_timepoint = "t0")

  expect_equal(nrow(cors), 1L)
  expect_true(all(c("rho_depth_shannon_clonality", "rho_depth_gini",
                    "rho_delta_shannon_vs_gini") %in% names(cors)))
  expect_equal(cors$n_samples, 9L)
  expect_false(is.na(cors$rho_delta_shannon_vs_gini))
  expect_true(abs(cors$rho_delta_shannon_vs_gini) <= 1)
  # Both metrics rise with the same expansion in this fixture, so they agree
  # perfectly; a disagreement here would be the signal the report warns
  # about.
  expect_equal(cors$rho_delta_shannon_vs_gini, 1)

  # Too few complete pairs for a rank correlation returns NA, not an error.
  tiny <- focal[1:2, ]
  expect_true(is.na(focalisation_depth_correlations(
    tiny, change[change$sample_id %in% tiny$sample_id, ])$rho_depth_gini))
})

test_that("the pharmacodynamic section writes every panel and the control", {
  r <- toy_run(list(differential_abundance = list(
    reference_timepoint = "t0", fold_change_threshold = 2,
    pharmacodynamics = list(focal_metric = "both", min_count = 2,
                            min_freq = 1e-5))))
  out <- section_pharmacodynamics(toy_clono(), r$run_dir, r$cfg, r$report)
  r$report$close()

  expect_named(out, c("focal", "change", "novel", "expanded"))
  files <- run_files(r$run_dir)
  for (stem in c("focalisation_metrics[.]csv", "focalisation_change[.]csv",
                 "novel_clonotypes[.]csv", "expanded_clonotypes[.]csv",
                 "expanded_clonotype_composition[.]csv",
                 "depth_control[.]csv", "focalisation_correlations[.]csv",
                 "focalisation_shannon_clonality[.]png",
                 "focalisation_gini[.]png", "novel_frequency[.]png",
                 "novel_count[.]png", "expanded_frequency[.]png",
                 "expanded_count[.]png", "expanded_composition[.]png")) {
    expect_true(any(grepl(stem, files)), info = stem)
  }

  # The depth caveat is in the report, because panel A cannot be read
  # without it.
  txt <- paste(readLines(file.path(r$run_dir, "REPORT.txt"), warn = FALSE),
               collapse = " ")
  expect_true(grepl("OBSERVED richness", txt))
  expect_true(grepl("Spearman correlation with depth", txt))
})

test_that("the pharmacodynamic section skips cleanly with no baseline", {
  r <- toy_run(list(differential_abundance = list(
    reference_timepoint = "t9", fold_change_threshold = 2)))
  clono <- toy_clono()
  expect_null(section_pharmacodynamics(clono, r$run_dir, r$cfg, r$report))
  r$report$close()
  txt <- paste(readLines(file.path(r$run_dir, "REPORT.txt"), warn = FALSE),
               collapse = " ")
  expect_true(grepl("No baseline timepoint", txt))
})

test_that("an unknown focusing metric is rejected before anything is written", {
  r <- toy_run(list(differential_abundance = list(
    reference_timepoint = "t0", fold_change_threshold = 2,
    pharmacodynamics = list(focal_metric = "entropy"))))
  expect_error(section_pharmacodynamics(toy_clono(), r$run_dir, r$cfg,
                                        r$report),
               "focal_metric")
  r$report$close()
})

# ===========================================================================
# CDR3 features of the expanded clonotypes
# ===========================================================================

test_that("the expanded-CDR3 section says so when nothing was called", {
  r <- toy_run()
  db <- toy_airr()
  db$clonotype <- clonotype_key(db, "strict")
  expect_null(section_cdr3_features_expanded(db, "umi_count",
                                             character(0), "clonotype",
                                             r$run_dir, r$cfg, r$report))
  r$report$close()
  txt <- paste(readLines(file.path(r$run_dir, "REPORT.txt"), warn = FALSE),
               collapse = " ")
  expect_true(grepl("This is a result, not a failure", txt))
})

test_that("the expanded-CDR3 section writes the network and its diagnostic", {
  skip_if_not_installed("igraph")
  r <- toy_run()
  db <- toy_airr()
  db$clonotype <- clonotype_key(db, "strict")
  out <- section_cdr3_features_expanded(db, "umi_count",
                                        unique(db$clonotype)[1], "clonotype",
                                        r$run_dir, r$cfg, r$report)
  r$report$close()

  expect_s3_class(out, "data.frame")
  expect_true(all(c("n_candidate_cdr3", "n_nodes", "n_edges") %in% names(out)))
  files <- run_files(r$run_dir)
  expect_true(any(grepl("cdr3_network_summary[.]csv$", files)))
  expect_true(any(grepl("cdr3_network_thresholds[.]csv$", files)))
  expect_true(any(grepl("cdr3_length_expanded[.]png$", files)))
  expect_true(any(grepl("cdr3_network[.]png$", files)))
})

test_that("the CDR3 structure section writes its three views", {
  skip_if_not_installed("alakazam")
  r <- toy_run()
  section_cdr3_structure(toy_airr(), "umi_count", "TRB", r$run_dir, r$cfg,
                         r$report)
  r$report$close()
  files <- run_files(r$run_dir)
  expect_true(any(grepl("cdr3_spectratype[.]csv$", files)))
  expect_true(any(grepl("cdr3_properties_summary[.]csv$", files)))
  expect_true(any(grepl("nucleotide_convergence[.]csv$", files)))
  expect_true(any(grepl("cdr3_aa_charge[.]png$", files)))
})

# ===========================================================================
# Transcriptome: panels, markers, contrasts and scores
# ===========================================================================

#' pbmc_small with the metadata the sections expect
toy_seurat <- function() {
  obj <- SeuratObject::pbmc_small
  obj$seurat_clusters <- obj$RNA_snn_res.0.8
  obj$condition  <- obj$groups
  obj$sample_id  <- obj$orig.ident
  obj$subject_id <- obj$groups
  obj
}

test_that("panel coverage names the genes the object does not carry", {
  skip_if_not_installed("Seurat")
  obj <- toy_seurat()
  cov <- panel_coverage(obj, list(bcell = c("MS4A1", "CD79A", "NOT_A_GENE"),
                                  absent = c("NOPE1", "NOPE2")))

  expect_equal(cov$coverage$n_present, c(2L, 0L))
  expect_equal(cov$coverage$missing[1], "NOT_A_GENE")
  expect_equal(cov$coverage$pct_present[2], 0)
  # A panel with no gene present is dropped from what will be plotted,
  # rather than being passed to Seurat, which would abort.
  expect_named(cov$kept, "bcell")
})

test_that("cluster markers come back tidy and rankable", {
  skip_if_not_installed("Seurat")
  obj <- toy_seurat()
  mk <- suppressWarnings(cluster_markers(obj, min_pct = 0.25,
                                         logfc_threshold = 0.25))

  expect_true(nrow(mk) > 0L)
  expect_true(all(c("gene", "cluster", "avg_log2FC", "p_val_adj") %in%
                    names(mk)))
  # only.pos: every marker is enriched in its own cluster.
  expect_true(all(mk$avg_log2FC > 0))

  top <- top_markers(mk, n = 3L)
  expect_true(all(table(top$cluster) <= 3L))
  strong <- top_markers(mk, n = 3L, min_log2fc = 1)
  expect_true(all(strong$avg_log2FC > 1))
  by_p <- top_markers(mk, n = 3L, by = "p_val_adj")
  expect_true(nrow(by_p) <= nrow(mk))

  expect_error(cluster_markers(obj, group_by = "not_a_column"),
               "not_a_column")
})

test_that("a contrast is labelled, ordered and checked", {
  skip_if_not_installed("Seurat")
  obj <- toy_seurat()
  res <- suppressWarnings(compare_conditions(obj, "g1", "g2",
                                             group_by = "condition"))

  expect_true(nrow(res) > 0L)
  expect_equal(unique(res$comparison), "g1_vs_g2")
  expect_equal(unique(res$group1), "g1")
  expect_true(all(diff(res$p_val_adj) >= 0))
  expect_true("gene" %in% names(res))

  # A group that is not in the data is an error naming the levels that are,
  # rather than an empty result.
  expect_error(compare_conditions(obj, "g1", "g9", group_by = "condition"),
               "g9")
  expect_error(compare_conditions(obj, "g1", "g2", group_by = "nope"), "nope")
})

test_that("a set of contrasts skips the impossible ones and keeps the rest", {
  skip_if_not_installed("Seurat")
  obj <- toy_seurat()
  res <- suppressWarnings(compare_condition_pairs(
    obj, list(list(group1 = "g1", group2 = "g2"),
              list(group1 = "g1", group2 = "absent")),
    group_by = "condition"))

  expect_equal(unique(res$comparison), "g1_vs_g2")
  expect_warning(compare_condition_pairs(obj, list(list(group1 = "g1")),
                                         group_by = "condition"),
                 "does not name two groups")
})

test_that("module scores are named after their panel and report coverage", {
  skip_if_not_installed("Seurat")
  obj <- toy_seurat()
  res <- module_scores(obj, list(bcell = c("MS4A1", "CD79A", "CD79B"),
                                 tcell = c("CD3D", "CD3E", "IL32"),
                                 thin = c("MS4A1", "CD79A")))

  # Named columns, not the positional sig_1/sig_2 AddModuleScore produces:
  # a dropped panel would otherwise shift every later panel onto the wrong
  # column.
  expect_setequal(res$scored, c("bcell", "tcell"))
  expect_true(all(c("bcell_score", "tcell_score") %in%
                    names(res$object@meta.data)))
  expect_false(any(grepl("_score1$", names(res$object@meta.data))))
  expect_equal(res$coverage$score_column[res$coverage$panel == "thin"],
               NA_character_)
  expect_false(res$coverage$scored[res$coverage$panel == "thin"])
  expect_true(is.numeric(res$object$bcell_score))

  # The B-cell score must be higher in the cluster expressing the B-cell
  # genes than in the other one, or the score is not measuring the module.
  means <- tapply(res$object$bcell_score, res$object$seurat_clusters, mean)
  expect_gt(max(means), min(means))
})

test_that("the dot plot and the heatmap return plots", {
  skip_if_not_installed("Seurat")
  obj <- toy_seurat()
  p <- suppressWarnings(marker_dotplot(
    obj, list(`B cell` = c("MS4A1", "CD79A"), other = c("CD3D", "IL32")),
    group_by = "seurat_clusters", title = "panel"))
  expect_s3_class(p, "ggplot")

  hm <- suppressWarnings(suppressMessages(marker_heatmap(
    obj, c("MS4A1", "CD79A", "CD3D"), group_by = "seurat_clusters",
    downsample = 20L)))
  expect_s3_class(hm, "ggplot")
  # No requested gene is present: NULL rather than an error from Seurat.
  expect_null(marker_heatmap(obj, c("NOPE1", "NOPE2")))
})

test_that("the transcriptome sections write their tables and figures", {
  skip_if_not_installed("Seurat")
  obj <- toy_seurat()
  r <- toy_run(list(gex = list(
    marker_panels = list(bcell = list(`B cell` = c("MS4A1", "CD79A"),
                                      naive = c("TCL1A", "NOT_A_GENE"))),
    markers = list(min_pct = 0.25, logfc_threshold = 0.25, top_n = 5,
                   heatmap_downsample = 20),
    comparisons = list(group_variable = "condition",
                       contrasts = list(list(group1 = "g1", group2 = "g2")),
                       top_n = 5, heatmap_downsample = 20),
    module_scores = list(bcell = c("MS4A1", "CD79A", "CD79B")))))

  suppressWarnings(suppressMessages({
    section_marker_panels(obj, r$run_dir, r$cfg, r$report)
    section_cluster_markers(obj, r$run_dir, r$cfg, r$report)
    section_condition_contrasts(obj, r$run_dir, r$cfg, r$report)
    scores <- section_module_scores(obj, r$run_dir, r$cfg, r$report)
  }))
  r$report$close()

  expect_true("bcell_score" %in% names(scores$object@meta.data))
  files <- run_files(r$run_dir)
  for (stem in c("panel_coverage[.]csv", "cluster_markers[.]csv",
                 "cluster_markers_top[.]csv", "condition_contrasts[.]csv",
                 "condition_contrasts_top[.]csv",
                 "module_score_coverage[.]csv",
                 "module_scores_by_cluster[.]csv",
                 "bcell_by_seurat_clusters[.]png",
                 "bcell_by_condition[.]png")) {
    expect_true(any(grepl(stem, files)), info = stem)
  }

  # The absent gene is recorded rather than silently dropped.
  cov <- utils::read.csv(file.path(r$run_dir, "tables",
                                   "06b_panel_coverage.csv"))
  expect_true(any(grepl("NOT_A_GENE", cov$missing)))
})

test_that("the transcriptome sections say so when nothing is configured", {
  skip_if_not_installed("Seurat")
  obj <- toy_seurat()
  r <- toy_run()
  expect_null(section_marker_panels(obj, r$run_dir, r$cfg, r$report))
  expect_null(section_condition_contrasts(obj, r$run_dir, r$cfg, r$report))
  r$report$close()

  txt <- paste(readLines(file.path(r$run_dir, "REPORT.txt"), warn = FALSE),
               collapse = " ")
  expect_true(grepl("No marker panel in the configuration", txt))
  expect_true(grepl("No contrast in the configuration", txt))
})
