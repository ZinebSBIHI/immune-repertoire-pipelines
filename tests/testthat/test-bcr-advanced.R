# Advanced B-cell receptor methods: mutation profiles, R/S, BASELINe,
# switched/unswitched load, bulk overviews and lineage trees.
#
# Hand-built tables check the arithmetic exactly; the Immcantation example
# data (alakazam::ExampleDb, dowser::ExampleAirr) check the steps that call
# shazam and dowser.

# A small subset of alakazam::ExampleDb: the first clones of 3 to 10
# sequences in each sample, so that every clone has several members and the
# shazam steps stay fast.
small_db <- function(n_clones = 10) {
  e <- as.data.frame(alakazam::ExampleDb)
  keep <- unlist(lapply(split(e, e$sample_id), function(s) {
    size <- table(s$clone_id)
    ok <- names(size)[size >= 3 & size <= 10]
    s$sequence_id[s$clone_id %in% utils::head(sort(as.numeric(ok)), n_clones)]
  }))
  e <- e[e$sequence_id %in% keep, ]
  e$subject_id <- "S1"
  e$compartment <- ifelse(e$sample_id == "+7d", "B", "A")
  rownames(e) <- NULL
  e
}

# A run directory, a configuration and a report, as init_run() would give.
adv_run <- function(bcr = list()) {
  run_dir <- tempfile("bcr_adv_")
  dir.create(file.path(run_dir, "figures"), recursive = TRUE)
  dir.create(file.path(run_dir, "tables"))
  list(run_dir = run_dir,
       cfg = list(output = list(formats = "png", dpi = 40), bcr = bcr),
       report = open_report(file.path(run_dir, "REPORT.txt"), "test"))
}

out_file <- function(run, kind, name) {
  file.path(run$run_dir, kind, name)
}

# --- mutation profiles -------------------------------------------------------

test_that("the mutation profile adds sub-region and total columns and keeps mu_freq", {
  skip_if_not_installed("shazam")
  skip_if_not_installed("alakazam")
  db <- small_db(4)
  db$mu_freq <- 0.123
  prof <- observed_mutation_profile(db, region = "IMGT_V")

  counts <- c("mu_count_cdr_r", "mu_count_cdr_s", "mu_count_fwr_r",
              "mu_count_fwr_s")
  expect_true(all(c(counts, "mu_count_imgt_v", "mu_freq_cdr_r",
                    "mu_freq_fwr_s", "mu_freq_imgt_v") %in% names(prof)))
  # Every mutation of the region falls in exactly one sub-region and is
  # either R or S: the identity the original analysis checked.
  expect_equal(prof$mu_count_imgt_v, unname(rowSums(prof[, counts])))
  expect_true(all(prof$mu_freq_imgt_v >= 0 & prof$mu_freq_imgt_v <= 1))
  expect_identical(prof$sequence_id, db$sequence_id)
  expect_true(all(prof$mu_freq == 0.123))
  expect_identical(attr(prof, "total_mutation_col"), "mu_freq_imgt_v")
  expect_identical(attr(prof, "germline_col"), "germline_alignment_d_mask")

  # Recomputing on a profiled table replaces the columns without warnings.
  expect_no_warning(again <- observed_mutation_profile(prof, region = "IMGT_V",
                                                       counts_only = TRUE))
  expect_equal(again$mu_count_cdr_r, prof$mu_count_cdr_r)
  expect_null(attr(again, "total_mutation_col"))
})

test_that("the V(D)J total does not depend on how the region is partitioned", {
  skip_if_not_installed("shazam")
  skip_if_not_installed("alakazam")
  db <- small_db(3)
  by_regions <- observed_mutation_profile(db, region = "IMGT_VDJ_BY_REGIONS")
  whole <- observed_mutation_profile(db, region = "IMGT_VDJ")
  expect_true(all(c("mu_freq_cdr3_r", "mu_count_fwr4_s") %in% names(by_regions)))
  expect_equal(by_regions$mu_freq_imgt_vdj_by_regions, whole$mu_freq_imgt_vdj)
  expect_equal(by_regions$mu_count_imgt_vdj_by_regions, whole$mu_count_imgt_vdj)

  # V(D)J boundaries need the junction length: rows without one are removed
  # for those definitions only.
  db$junction_length[1] <- NA
  expect_equal(nrow(observed_mutation_profile(db, region = "IMGT_VDJ",
                                              counts_only = TRUE)),
               nrow(db) - 1L)
  expect_equal(nrow(observed_mutation_profile(db, region = "IMGT_V",
                                              counts_only = TRUE)),
               nrow(db))
})

test_that("unknown region definitions and missing germlines are reported", {
  db <- data.frame(sequence_alignment = "ACGT", locus = "IGH")
  expect_error(selection_germline_col(db), "No germline column")
  skip_if_not_installed("shazam")
  expect_error(shazam_region("IMGT_NOPE"), "Unknown shazam region")
})

# --- R/S -----------------------------------------------------------------------

test_that("R/S sums mutations over the group before dividing, flooring S at 1", {
  d <- data.frame(subject_id = "S1", compartment = c("A", "A", "B"),
                  mu_count_cdr_r = c(2, 4, 3), mu_count_cdr_s = c(1, 1, 0),
                  mu_count_fwr_r = c(1, 1, 0), mu_count_fwr_s = c(2, 2, 0))
  rs <- replacement_silent_ratio(d)
  a <- rs[rs$compartment == "A", ]
  b <- rs[rs$compartment == "B", ]
  expect_equal(a$n, 2L)
  expect_equal(a$rs_cdr, 6 / 2)
  expect_equal(a$rs_fwr, 2 / 4)
  expect_equal(b$rs_cdr, 3)      # silent count 0, floored at 1
  expect_equal(b$rs_fwr, 0)
  expect_error(replacement_silent_ratio(d[, -3]), "mu_count_cdr_r")
})

# --- switched and unswitched ---------------------------------------------------

test_that("switched and unswitched load follows the explicit switched list", {
  d <- data.frame(
    locus = "IGH",
    subject_id  = c("S1", "S1", "S1", "S2", "S2", "S2", "S2", "S1"),
    wave        = c("w1", "w1", "w1", "w2", "w2", "w2", "w2", "w1"),
    compartment = "M",
    c_gene      = c("IGHG1", "IGHM", "IGHA2", "IGHG", "IGHD", "", NA, "IGHG1"),
    mu_freq     = c(0.10, 0.00, 0.08, 0.05, 0.01, 0.20, 0.30, 0.06),
    stimulated  = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE)
  res <- shm_by_switch_status(d, subject_group_col = "wave")

  expect_equal(res$excluded[["no_isotype"]], 2)
  expect_equal(res$excluded[["stimulated"]], 1)

  bs <- as.data.frame(res$by_subject)
  s1_sw <- bs[bs$subject_id == "S1" & bs$isotype_class == "class-switched", ]
  expect_equal(s1_sw$cells, 2L)
  expect_equal(s1_sw$shm, 9)                         # mean of 10% and 8%
  s2_un <- bs[bs$subject_id == "S2" & bs$isotype_class == "unswitched", ]
  expect_equal(s2_un$cells, 2L)                      # IGHG without subclass counts as unswitched
  expect_equal(s2_un$shm, 3)

  sh <- as.data.frame(res$unswitched_share)
  expect_equal(sh$pct_unswitched[sh$subject_id == "S1"], 100 / 3)
  expect_equal(sh$shm_all[sh$subject_id == "S1"], 6)
  expect_equal(res$outside_lists$isotype, "IGHG")
  expect_equal(res$outside_lists$n_cells, 1L)

  # Without subject groups, subjects are summarised together.
  pooled <- as.data.frame(shm_by_switch_status(d)$summary)
  un <- pooled[pooled$isotype_class == "unswitched", ]
  expect_equal(un$n_subjects, 2L)
  expect_equal(c(un$shm_mean, un$shm_min, un$shm_max), c(1.5, 0, 3))
})

# --- bulk aggregations ---------------------------------------------------------

bulk_toy <- function() {
  data.frame(
    locus = "IGH", subject_id = "S1",
    clone_id = c(1, 1, 1, 2, 2, 3),
    sample_id = c("x", "x", "y", "x", "x", "y"),
    c_gene = c("IGHG1", "IGHG1", "IGHG1", "IGHM", "IGHM", "IGHA1"),
    v_call = c("IGHV1-2*01", "IGHV1-2*02", "IGHV1-2*01", "IGHV3-23*01",
               "IGHV3-23*01", "IGHV4-34*01"),
    j_call = "IGHJ4*02",
    junction_aa = c("CARDYW", "CARDYW", "CARDYW", "CAKW", "CAKW", "CARGGW"),
    mu_freq = c(0.1, 0.2, 0.3, 0, 0.02, 0.05),
    mu_freq_cdr3_r = c(0.01, 0.03, 0.02, 0, 0, 0.1),
    stringsAsFactors = FALSE)
}

test_that("per-clone means, largest clones and the overview aggregate as in the bulk analysis", {
  d <- bulk_toy()

  cm <- as.data.frame(clone_mean_mutation(d))
  expect_equal(nrow(cm), 4L)
  expect_equal(cm$mean_mut_freq[cm$clone_id == 1 & cm$sample_id == "x"], 0.15)
  expect_equal(cm$n_sequences[cm$clone_id == 1 & cm$sample_id == "y"], 1L)

  top <- as.data.frame(top_clones_summary(d, n = 2))
  expect_equal(top$clone_id, c(1, 2))
  expect_equal(top$n_sequences, c(3L, 2L))
  expect_equal(top$mean_mu_freq[1], 0.2)
  expect_equal(top$v_genes[1], "IGHV1-2")
  expect_equal(top$present_in[1], "x; y")

  ov <- as.data.frame(vgene_cdr3_overview(d))
  expect_equal(nrow(ov), 4L)
  row <- ov[ov$clone_id == 1 & ov$sample_id == "x", ]
  expect_equal(row$median_mut, 2)                    # median(1%, 3%)
  expect_equal(row$seq_count, 2L)
  expect_equal(row$junction_aa_length, 6L)
  only1 <- vgene_cdr3_overview(d, clones = "S1_1")
  expect_true(all(only1$clone_id == 1))
})

# --- BASELINe --------------------------------------------------------------------

test_that("BASELINe estimates Sigma per group and records skipped groups", {
  skip_if_not_installed("shazam")
  skip_if_not_installed("alakazam")
  db <- small_db(6)
  db$compartment[1:2] <- "C"
  res <- suppressMessages(baseline_selection(db, group_col = "compartment"))
  expect_s4_class(res$baseline, "Baseline")
  expect_setequal(unique(as.character(res$stats$compartment)), c("A", "B"))
  expect_true(all(c("baseline_sigma", "baseline_ci_lower",
                    "baseline_ci_upper") %in% names(res$stats)))
  st <- res$groups
  expect_match(st$status[st$compartment == "C"], "^skipped")
  expect_equal(st$status[st$compartment != "C"], c("computed", "computed"))
})

# --- lineage -------------------------------------------------------------------

test_that("lineage clones must reach every set and are kept apart by subject", {
  d <- data.frame(
    locus = "IGH",
    subject_id  = c(rep("S1", 7), rep("S2", 4)),
    clone_id    = c(1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1),
    compartment = c("mem", "mem", "pb", "mem", "mem", "mem", "mem",
                    "mem", "pb", "pb", "pb"),
    condition   = c("c1", "c1", "c1", "c1", "c1", "c1", "c1",
                    "c1", "c2", "c2", "c2"),
    stringsAsFactors = FALSE)
  sets <- list(list(compartment = "mem"), list(compartment = "pb"))

  sel <- select_lineage_clones(d, sets = sets, min_cells = 3)
  expect_equal(sel$clone_key, c("S2_1", "S1_1"))     # largest first
  expect_equal(sel$n_sequences, c(4L, 3L))
  expect_true(all(sel$in_set_1 & sel$in_set_2))

  expect_equal(select_lineage_clones(d, sets, min_cells = 4)$clone_key, "S2_1")
  # A set may constrain several columns at once.
  two_cols <- list(list(compartment = "mem"),
                   list(compartment = "pb", condition = "c2"))
  expect_equal(select_lineage_clones(d, two_cols)$clone_key, "S2_1")
  # Without sets, size alone decides.
  expect_equal(nrow(select_lineage_clones(d, NULL, min_cells = 4)), 2L)
  expect_error(select_lineage_clones(d, list(list("mem"))), "named list")
})

test_that("clone labels take the most frequent V, J and CDR3", {
  d <- data.frame(
    locus = "IGH", subject_id = "S1", clone_id = c(5, 5, 5, 6, 6),
    v_call = c("IGHV1-2*01", "IGHV1-2*02", "IGHV1-69*01", "IGHV3-7*01",
               "IGHV3-23*01"),
    j_call = c("IGHJ4*02", "IGHJ6*02", "IGHJ6*03", "IGHJ4*02", "IGHJ4*02"),
    junction_aa = c("CARW", "CARW", "CAKW", "CTTW", "CSSW"),
    stringsAsFactors = FALSE)
  lab <- clone_dominant_labels(d)
  expect_equal(lab$label[lab$clone_key == "S1_5"], "IGHV1-2 | IGHJ6 | CARW")
  # Ties go to the first value in alphabetical order.
  expect_equal(lab$label[lab$clone_key == "S1_6"], "IGHV3-23 | IGHJ4 | CSSW")
})

example_airr <- function() {
  a <- as.data.frame(dowser::ExampleAirr)
  a$locus <- "IGH"
  a$subject_id <- "S1"
  # The largest clone is left out to keep tree building quick.
  a[a$clone_id != 3128, ]
}

test_that("one maximum-likelihood tree is built per selected clone", {
  skip_if_not_installed("dowser")
  skip_if_not_installed("phangorn")
  skip_if_not_installed("alakazam")
  a <- example_airr()
  sel <- select_lineage_clones(a, sets = list(list(location = "7d"),
                                              list(location = "14d")))
  expect_gt(nrow(sel), 1L)
  keys <- utils::head(sel$clone_key, 2)

  trees <- suppressWarnings(
    build_lineage_trees(a, keys, traits = c("location", "subject_id")))
  expect_true(nrow(trees) >= 1L)
  expect_true(all(trees$clone_id %in% keys))
  expect_true(all(vapply(trees$trees, inherits, logical(1), what = "phylo")))

  bulk <- suppressWarnings(
    build_lineage_trees(a, keys, mode = "bulk",
                        traits = c("subject_id", "c_call")))
  expect_true(all(bulk$clone_id %in% keys))
})

# --- sections ------------------------------------------------------------------

test_that("the single-cell sections write their tables and figures", {
  skip_if_not_installed("shazam")
  skip_if_not_installed("alakazam")
  run <- adv_run()
  db <- observed_mutation_profile(small_db(6), region = "IMGT_V")
  db$mu_freq <- db$mu_freq_imgt_v
  db$c_gene <- strip_allele(db$c_call)

  section_mutation_load(db, x_col = "compartment", run_dir = run$run_dir,
                        cfg = run$cfg, report = run$report)
  section_mutation_load(db, x_col = "c_gene", facet_rows = "subject_id",
                        facet_cols = "compartment",
                        fill_values = c(IGHA = "grey70", IGHD = "grey40",
                                        IGHG = "#cb181d", IGHM = "grey90"),
                        name = "mutation_load_isotype",
                        run_dir = run$run_dir, cfg = run$cfg,
                        report = run$report)
  rs <- section_rs_ratio(db, group_col = "compartment", run_dir = run$run_dir,
                         cfg = run$cfg, report = run$report)
  sw <- section_shm_switched(db, run_dir = run$run_dir, cfg = run$cfg,
                             report = run$report)
  run$report$close()

  for (f in c("10_mutation_load", "10_mutation_load_isotype", "11_rs_ratio")) {
    expect_true(file.exists(out_file(run, "tables", paste0(f, ".csv"))))
    expect_true(file.exists(out_file(run, "figures", paste0(f, ".png"))))
  }
  expect_true(file.exists(out_file(run, "tables",
                                   "13_shm_switched_summary.csv")))
  expect_setequal(as.character(rs$compartment), c("A", "B"))
  # ExampleDb isotypes carry no subclass, so none is in the switched list.
  expect_true(all(sw$by_subject$isotype_class == "unswitched"))
  expect_match(paste(readLines(file.path(run$run_dir, "REPORT.txt")),
                     collapse = "\n"), "counted as unswitched")
})

test_that("the BASELINe section writes Sigma and its figure", {
  skip_if_not_installed("shazam")
  skip_if_not_installed("alakazam")
  run <- adv_run(bcr = list(selection = list(baseline_min_clones = 3)))
  res <- suppressMessages(section_baseline_selection(
    small_db(5), group_col = "compartment", label = "test view",
    run_dir = run$run_dir, cfg = run$cfg, report = run$report))
  run$report$close()
  expect_true(file.exists(out_file(run, "tables",
                                   "12_baseline_selection_sigma.csv")))
  expect_true(file.exists(out_file(run, "figures",
                                   "12_baseline_selection.png")))
  expect_equal(res$groups$status, c("computed", "computed"))
})

test_that("the bulk sections reproduce the isotype, overview and largest-clone steps", {
  skip_if_not_installed("shazam")
  skip_if_not_installed("alakazam")
  # ExampleDb isotypes have no subclass; the order is set in the
  # configuration, as a cohort would.
  run <- adv_run(bcr = list(
    mutation_by_isotype = list(isotype_levels = c("IGHD", "IGHM", "IGHG",
                                                  "IGHA")),
    overview = list(top_n_clones = 5)))
  db <- small_db(6)
  prof <- section_mutation_by_isotype(db, group_col = "sample_id",
                                      group_levels = c("-1h", "+7d"),
                                      run_dir = run$run_dir, cfg = run$cfg,
                                      report = run$report)
  expect_identical(attr(prof, "total_mutation_col"),
                   "mu_freq_imgt_vdj_by_regions")
  ov <- section_vgene_cdr3_overview(prof, run_dir = run$run_dir,
                                    cfg = run$cfg, report = run$report)
  top <- section_top_clones(prof, run_dir = run$run_dir, cfg = run$cfg,
                            report = run$report)
  run$report$close()

  for (f in c("10a_vh_mutation_by_isotype", "10b_cdr3_replacement_by_isotype",
              "10c_clone_mutation_density", "11_vgene_cdr3_overview",
              "12b_top_clones_overview")) {
    expect_true(file.exists(out_file(run, "figures", paste0(f, ".png"))),
                info = f)
  }
  for (f in c("10a_mutation_by_isotype", "10c_clone_mean_mutation",
              "11_vgene_cdr3_overview", "12a_top_clones")) {
    expect_true(file.exists(out_file(run, "tables", paste0(f, ".csv"))),
                info = f)
  }
  expect_equal(nrow(top), 5L)
  expect_true(all(ov$median_mut >= 0))
})

test_that("the lineage section saves trees and states why figures are missing", {
  skip_if_not_installed("dowser")
  skip_if_not_installed("phangorn")
  skip_if_not_installed("alakazam")
  run <- adv_run(bcr = list(lineage = list(max_clones = 2)))
  trees <- suppressWarnings(section_lineage_trees(
    example_airr(), sets = list(list(location = "7d"), list(location = "14d")),
    tip_col = "location", run_dir = run$run_dir, cfg = run$cfg,
    report = run$report))
  run$report$close()

  expect_true(file.exists(file.path(run$run_dir, "14_lineage_trees.rds")))
  expect_true(file.exists(out_file(run, "tables", "14_lineage_trees.csv")))
  expect_lte(nrow(trees), 2L)
  report <- paste(readLines(file.path(run$run_dir, "REPORT.txt")),
                  collapse = "\n")
  if (length(lineage_plot_missing())) {
    expect_match(report, "Tree figures not drawn")
  } else {
    expect_true(file.exists(out_file(run, "figures",
                                     "14_lineage_trees_p01.png")))
  }

  # No qualifying clone: reported, nothing built.
  run2 <- adv_run()
  none <- section_lineage_trees(
    example_airr(), sets = list(list(location = "nowhere")),
    run_dir = run2$run_dir, cfg = run2$cfg, report = run2$report)
  run2$report$close()
  expect_null(none)
})
