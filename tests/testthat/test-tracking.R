# Tests for clone tracking (rep-tracking.R) and V-J chord diagrams
# (rep-circos.R).
#
# Every dataset is simulated here. The bulk simulator follows the generative
# model of examples/simulate_airr.R (log-normal baseline, a set of
# clonotypes expanded after baseline, a set absent at baseline, multinomial
# sampling) at a size that keeps the tests fast. The data carry no
# biological meaning.

# --- simulated data -----------------------------------------------------------

trk_cdr3 <- function(n) {
  aa <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
  vapply(seq_len(n), function(i) {
    paste0("CASS", paste(sample(aa, sample(6:11, 1L), replace = TRUE),
                         collapse = ""), "F")
  }, character(1))
}

trk_sim_bulk <- function(seed = 1, n_clones = 250,
                         subjects = c("s1", "s2"),
                         timepoints = c("t0", "t1", "t2")) {
  set.seed(seed)
  v_genes <- sprintf("TRBV%d-%d", rep(2:13, each = 2), rep(1:2, 12))
  j_genes <- sprintf("TRBJ%d-%d", rep(1:2, each = 6), c(1:6, 1:6))
  rows <- list()
  for (s in subjects) {
    pool <- data.frame(junction_aa = trk_cdr3(n_clones),
                       v_call = paste0(sample(v_genes, n_clones, TRUE), "*01"),
                       j_call = paste0(sample(j_genes, n_clones, TRUE), "*01"),
                       stringsAsFactors = FALSE)
    pool$w <- stats::rlnorm(n_clones, 0, 1.5)
    expanded <- sample.int(n_clones, 10)
    dn <- data.frame(junction_aa = trk_cdr3(8),
                     v_call = paste0(sample(v_genes, 8, TRUE), "*01"),
                     j_call = paste0(sample(j_genes, 8, TRUE), "*01"),
                     stringsAsFactors = FALSE)
    dn$w <- 0.05 * sum(pool$w) / 8
    for (tp in timepoints) {
      tab <- pool
      if (tp != timepoints[1]) {
        tab$w[expanded] <- tab$w[expanded] * 12
        tab <- rbind(tab, dn)
      }
      counts <- stats::rmultinom(1, 6000, tab$w / sum(tab$w))[, 1]
      keep <- counts > 0
      rows[[length(rows) + 1L]] <- data.frame(
        sample_id = paste(s, tp, sep = "_"), subject_id = s, timepoint = tp,
        locus = "TRB", productive = TRUE,
        v_call = tab$v_call[keep], j_call = tab$j_call[keep],
        junction = tab$junction_aa[keep], junction_aa = tab$junction_aa[keep],
        umi_count = counts[keep], stringsAsFactors = FALSE)
    }
  }
  db <- do.call(rbind, rows)
  db$subject_id <- factor(db$subject_id, levels = subjects)
  db$timepoint <- factor(db$timepoint, levels = timepoints)
  db
}

# One subject, three compartments; some clonotypes span compartments; some
# cells carry two TRB contigs; a few TRA and non-productive contigs.
trk_sim_sc <- function(seed = 2) {
  set.seed(seed)
  v_genes <- sprintf("TRBV%d-%d", rep(2:13, each = 2), rep(1:2, 12))
  j_genes <- sprintf("TRBJ%d-%d", rep(1:2, each = 6), c(1:6, 1:6))
  n_cl <- 60
  clones <- data.frame(junction_aa = trk_cdr3(n_cl),
                       v_call = sample(v_genes, n_cl, TRUE),
                       j_call = sample(j_genes, n_cl, TRUE),
                       stringsAsFactors = FALSE)
  comps <- c("comp_a", "comp_b", "comp_c")
  rows <- list()
  for (k in seq_along(comps)) {
    w <- stats::rlnorm(n_cl, 0, 1.2)
    w[sample.int(n_cl, 30)] <- 0
    w[1:5] <- 3
    idx <- sample.int(n_cl, 150, replace = TRUE, prob = w / sum(w))
    cells <- sprintf("%s_cell%03d", comps[k], seq_along(idx))
    rows[[k]] <- data.frame(sample_id = paste0("p1_", comps[k]),
                            subject_id = "p1", condition = comps[k],
                            cell_id = cells, locus = "TRB", productive = TRUE,
                            v_call = clones$v_call[idx], j_call = clones$j_call[idx],
                            junction = clones$junction_aa[idx],
                            junction_aa = clones$junction_aa[idx],
                            umi_count = 3L, stringsAsFactors = FALSE)
  }
  db <- do.call(rbind, rows)
  extra <- db[1:10, ]
  extra$locus <- "TRA"
  extra$v_call <- "TRAV1-2"
  extra$j_call <- "TRAJ33"
  np <- db[11:15, ]
  np$productive <- FALSE
  rbind(db, extra, np)
}

trk_sc_clono <- function(db, definition = "strict") {
  cells <- filter_locus_productive(db, "TRB")
  cells$clonotype <- clonotype_key(cells, definition)
  u <- unique(cells[, c("sample_id", "clonotype", "cell_id")])
  cl <- stats::aggregate(list(abundance = rep(1, nrow(u))),
                         by = list(sample_id = u$sample_id,
                                   clonotype = u$clonotype), FUN = sum)
  cl$frequency <- cl$abundance / stats::ave(cl$abundance, cl$sample_id, FUN = sum)
  meta <- unique(db[, c("sample_id", "subject_id", "condition")])
  merge(cl, meta, by = "sample_id")
}

trk_cfg <- function(tracking = list(), circos = list(),
                    formats = c("pdf", "png")) {
  validate_config(list(
    project = list(name = "test", seed = 1L),
    paths = list(data_dir = tempdir(), output_dir = tempdir()),
    cohort = list(subjects = c("s1", "s2", "p1"),
                  timepoints = c("t0", "t1", "t2"),
                  conditions = c("comp_a", "comp_b", "comp_c")),
    repertoire = list(locus = "TRB", clonotype_definition = "strict",
                      abundance_column = "umi_count",
                      tracking = tracking, circos = circos),
    differential_abundance = list(reference_timepoint = "t0"),
    output = list(formats = formats, dpi = 72)))
}

trk_run_dir <- function() {
  d <- tempfile("trk_run_")
  dir.create(file.path(d, "figures"), recursive = TRUE)
  dir.create(file.path(d, "tables"))
  d
}

# --- clonalCompare-style table ---------------------------------------------

test_that("clone_compare_table takes the union of per-sample top clonotypes", {
  toy <- data.frame(sample_id = c("A", "A", "A", "B", "B", "B"),
                    clonotype = c("a1", "a2", "x", "x", "b1", "a1"),
                    abundance = c(50, 30, 20, 60, 30, 10))
  tab <- clone_compare_table(toy, top_n = 1)
  expect_setequal(as.character(unique(tab$clonotype)), c("a1", "x"))
  expect_equal(nrow(tab), 4L)
  # proportions are of the whole sample, computed before selection
  expect_equal(tab$proportion[tab$sample_id == "A" & tab$clonotype == "x"], 0.2)
  expect_equal(tab$proportion[tab$sample_id == "B" & tab$clonotype == "a1"], 0.1)
  expect_true(tab$in_top[tab$sample_id == "A" & tab$clonotype == "a1"])
  expect_false(tab$in_top[tab$sample_id == "A" & tab$clonotype == "x"])
  expect_equal(tab$rank[tab$sample_id == "A" & tab$clonotype == "x"], 3L)
})

test_that("clone_compare_table orders samples and breaks ties without extras", {
  toy <- data.frame(sample_id = c("A", "A", "B", "B", "B"),
                    clonotype = c("c2", "c1", "q", "c1", "c2"),
                    abundance = c(10, 10, 5, 3, 1))
  tab <- clone_compare_table(toy, samples = c("B", "A"), top_n = 1)
  expect_equal(levels(tab$sample_id), c("B", "A"))
  expect_equal(sum(tab$in_top[tab$sample_id == "A"]), 1L)
  expect_true(tab$in_top[tab$sample_id == "A" & tab$clonotype == "c1"])
  # the requested samples restrict the comparison without changing proportions
  only_a <- clone_compare_table(toy, samples = "A", top_n = 2)
  expect_equal(sum(only_a$proportion), 1)
  expect_warning(clone_compare_table(toy, samples = c("A", "Z"), top_n = 1),
                 "not found")
})

test_that("clone_compare_table counts rows when no abundance column is given", {
  cells <- data.frame(sample_id = rep(c("A", "B"), c(5, 4)),
                      clonotype = c("k", "k", "k", "m", "n", "k", "m", "m", "m"))
  tab <- clone_compare_table(cells, top_n = 1, abundance_col = NULL)
  expect_equal(tab$count[tab$sample_id == "A" & tab$clonotype == "k"], 3)
  expect_equal(tab$proportion[tab$sample_id == "B" & tab$clonotype == "m"], 0.75)
})

test_that("plot_clone_compare and plot_clone_tracking build", {
  skip_if_not_installed("ggalluvial")
  toy <- data.frame(sample_id = c("A", "A", "B", "B", "C"),
                    clonotype = c("a", "b", "a", "c", "a"),
                    abundance = c(5, 3, 2, 6, 4))
  p <- plot_clone_compare(clone_compare_table(toy, top_n = 2))
  expect_s3_class(p, "ggplot")
  expect_error(ggplot2::ggplot_build(p), NA)
  p2 <- plot_clone_compare(clone_compare_table(toy, top_n = 2), value = "count")
  expect_match(p2$labels$y, "Count")

  # a timepoint where no tracked clonotype is detected (the baseline, for de
  # novo clonotypes) stays on the axis, empty, without breaking the stacking
  tr <- data.frame(timepoint = factor(c("t0", "t0", "t1", "t1", "t2", "t2"),
                                      levels = c("t0", "t1", "t2")),
                   clonotype = c("a", "b", "a", "b", "a", "b"),
                   frequency = c(0, 0, 0.2, 0.1, 0.05, 0))
  p3 <- plot_clone_tracking(tr)
  expect_warning(b3 <- ggplot2::ggplot_build(p3), NA)
  expect_equal(b3$layout$panel_params[[1]]$x$get_limits(), c("t0", "t1", "t2"))
})

# --- tracking over time ---------------------------------------------------

test_that("track_baseline_clonotypes follows baseline top N with zeros", {
  clono <- data.frame(
    sample_id = c("s_t0", "s_t0", "s_t0", "s_t1", "s_t2", "s_t2"),
    subject_id = "s",
    timepoint = factor(c("t0", "t0", "t0", "t1", "t2", "t2"),
                       levels = c("t0", "t1", "t2")),
    clonotype = c("c1", "c2", "c3", "c1", "c2", "c9"),
    abundance = c(30, 20, 10, 5, 8, 2))
  tr <- track_baseline_clonotypes(clono, top_n = 2)
  expect_equal(nrow(tr), 6L)
  expect_setequal(unique(tr$clonotype), c("c1", "c2"))
  expect_equal(tr$frequency[tr$timepoint == "t1" & tr$clonotype == "c2"], 0)
  expect_equal(tr$frequency[tr$timepoint == "t0" & tr$clonotype == "c1"], 0.5)
  expect_equal(levels(tr$timepoint), c("t0", "t1", "t2"))
  # a declared baseline the subject has is used instead of the first timepoint
  tr2 <- track_baseline_clonotypes(clono, top_n = 1, baseline = "t2")
  expect_equal(unique(tr2$clonotype), "c2")
})

test_that("select_tracked_expanded ranks within subject and type", {
  calls <- data.frame(
    subject_id = "s", timepoint = factor(c("t1", "t2", "t1", "t1", "t2"),
                                         levels = c("t0", "t1", "t2")),
    clonotype = c("k", "k", "m", "d", "e"),
    category = c("pre_existing_expanded", "pre_existing_expanded",
                 "pre_existing_expanded", "de_novo", "de_novo"),
    log2_fold_change = c(2, 5, 4, 3, 3),
    sig_expansion = TRUE, stringsAsFactors = FALSE)
  first <- select_tracked_expanded(calls, top_n = 1, rank_by = "first")
  expect_equal(first$clonotype[first$type == "pre_existing_expanded"], "m")
  peak <- select_tracked_expanded(calls, top_n = 1, rank_by = "peak")
  expect_equal(peak$clonotype[peak$type == "pre_existing_expanded"], "k")
  expect_equal(peak$n_timepoints_significant[peak$clonotype == "k"], 2L)
  # ties at the boundary are kept
  expect_setequal(first$clonotype[first$type == "de_novo"], c("d", "e"))
  none <- select_tracked_expanded(transform(calls, sig_expansion = FALSE))
  expect_equal(nrow(none), 0L)
})

test_that("tracking works on simulated bulk data with real calls", {
  db <- trk_sim_bulk()
  clono <- collapse_clonotypes(db, "umi_count", definition = "strict")
  meta <- unique(db[, c("sample_id", "timepoint")])
  clono <- merge(clono, meta, by = "sample_id")
  calls <- test_clonal_abundance(clono, reference_timepoint = "t0")
  sel <- select_tracked_expanded(calls, top_n = 5)
  expect_true(all(sel$type %in% c("de_novo", "pre_existing_expanded")))
  expect_true(all(table(sel$subject_id, sel$type) >= 1))
  tr <- track_clonotypes(clono, sel)
  # every tracked clonotype has one row per sample of its subject
  per <- table(paste(tr$subject_id, tr$type, tr$clonotype))
  expect_true(all(per == 3L))
  de_novo_t0 <- tr$frequency[tr$type == "de_novo" & tr$timepoint == "t0"]
  expect_true(all(de_novo_t0 == 0))
})

test_that("compartment_flows keeps shared clonotypes and flags the highlighted", {
  clono <- data.frame(subject_id = "p", condition = c("a", "b", "a", "c", "c", "b"),
                      clonotype = c("x", "x", "y", "y", "z", "w"),
                      abundance = c(3, 1, 2, 2, 5, 1))
  fl <- compartment_flows(clono, highlight = "c")
  expect_setequal(unique(fl$clonotype), c("x", "y"))
  expect_true(all(fl$highlighted[fl$clonotype == "y"]))
  expect_false(any(fl$highlighted[fl$clonotype == "x"]))
  expect_equal(levels(fl$compartment), c("a", "b", "c"))
  fl3 <- compartment_flows(clono, min_compartments = 3)
  expect_equal(nrow(fl3), 0L)
})

# --- V-J chord diagrams -------------------------------------------------------

test_that("vj_circos_data keeps the locus, the top V genes and their J genes", {
  db <- data.frame(
    sample_id = "A",
    locus = c(rep("TRB", 9), "TRA", "TRB"),
    v_call = c("TRBV20-1*01", "TRBV20-1", "TRBV20-1", "TRBV5-1", "TRBV5-1",
               "TRBV9", "TRBV9", "TRBV2", "TRBV2", "TRAV1-2", "TRBV20-1"),
    j_call = c("TRBJ2-7", "TRBJ2-7", "TRBJ1-1", "TRBJ2-7", "TRBJ2-1",
               "TRBJ1-1", "TRBJ2-7", "TRBJ1-1", "TRBJ2-7", "TRAJ33", "None"),
    stringsAsFactors = FALSE)
  cd <- vj_circos_data(db, sample = "A", top_v = 2, locus = "TRB")
  expect_equal(cd$total, 9L)
  expect_equal(rownames(cd$matrix)[1], "V_V20-1")
  expect_equal(nrow(cd$matrix), 2L)
  expect_equal(colnames(cd$matrix)[1], "J_J2-7")
  expect_equal(sum(cd$matrix), 5)
  expect_true(all(diff(rowSums(cd$matrix)) <= 0))
  expect_null(vj_circos_data(db, sample = "A", locus = "IGH"))
  igh <- data.frame(sample_id = "B", v_call = "IGHV3-23*01", j_call = "IGHJ4*02")
  expect_equal(rownames(vj_circos_data(igh)$matrix), "V_V3-23")
})

test_that("vj_colour_map keeps the fixed palettes and extends them", {
  cm <- vj_colour_map(paste0("V_V", 1:3), paste0("J_J", 1:2))
  expect_equal(unname(cm$v_colors[1]), "#1f77b4")
  expect_equal(unname(cm$j_colors[1]), "#8dd3c7")
  big <- vj_colour_map(paste0("V_V", 1:25), paste0("J_J", 1:17))
  expect_length(big$v_colors, 25L)
  expect_false(anyNA(big$v_colors))
  expect_length(unique(big$j_colors), 17L)
  ov <- vj_colour_map(c("V_V20-1", "V_V9"), "J_J1",
                      override_v = c(`TRBV20-1` = "#123456", V9 = "#abcdef"))
  expect_equal(unname(ov$v_colors[c("V_V20-1", "V_V9")]), c("#123456", "#abcdef"))
})

test_that("save_base_figure writes each format and records the legend", {
  skip_if_not_installed("circlize")
  run_dir <- trk_run_dir()
  cfg <- trk_cfg()
  db <- trk_sim_sc()
  cd <- vj_circos_data(db, sample = "p1_comp_a", locus = "TRB")
  cm <- vj_colour_map(sort(cd$top_v), sort(cd$all_j))
  paths <- save_base_figure(function(fmt) draw_vj_circos(cd, cm, "title", "sub"),
                            "chord", run_dir, cfg, width = 8, height = 8.5,
                            legend = "a legend")
  expect_length(paths, 2L)
  expect_true(all(file.size(paths) > 1000))
  expect_match(paste(readLines(file.path(run_dir, "FIGURE_LEGENDS.txt")),
                     collapse = "\n"), "chord\na legend")
  expect_warning(bad <- save_base_figure(function() stop("boom"), "bad", run_dir,
                                         trk_cfg(formats = "pdf")), "boom")
  expect_length(bad, 0L)
  expect_false(file.exists(file.path(run_dir, "figures", "bad.pdf")))
})

test_that("draw_vj_chord keeps the top genes with ties", {
  skip_if_not_installed("circlize")
  pairs <- data.frame(v_gene = c("TRBV1", "TRBV2", "TRBV3", "TRBV4"),
                      j_gene = c("TRBJ1", "TRBJ1", "TRBJ2", "TRBJ3"),
                      n = c(5, 3, 3, 1))
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  m <- draw_vj_chord(pairs, c(TRBV1 = "red", TRBV2 = "blue", TRBV3 = "green",
                              TRBV4 = "grey85"),
                     c(TRBJ1 = "black", TRBJ2 = "orange", TRBJ3 = "pink"),
                     title = "test", top_v = 2, top_j = 3)
  empty <- draw_vj_chord(pairs[0, ], c(), c(), title = "empty")
  grDevices::dev.off()
  expect_equal(rownames(m), c("TRBV1", "TRBV2", "TRBV3"))
  expect_null(empty)
})

# --- sections, end to end -------------------------------------------------

test_that("bulk sections run end to end on simulated data", {
  skip_if_not_installed("ggalluvial")
  skip_if_not_installed("circlize")
  run_dir <- trk_run_dir()
  # PDF only here, to keep the test fast; the single-cell test covers PNG.
  cfg <- trk_cfg(tracking = list(compare_pairs = list(c("t1", "t2"))),
                 formats = "pdf")
  report <- open_report(file.path(run_dir, "REPORT.txt"), "test")
  on.exit(report$close(), add = TRUE)

  db <- trk_sim_bulk()
  clono <- collapse_clonotypes(db, "umi_count", definition = "strict")
  clono <- merge(clono[, setdiff(names(clono), "subject_id")],
                 unique(db[, c("sample_id", "subject_id", "timepoint")]),
                 by = "sample_id")
  calls <- test_clonal_abundance(clono, reference_timepoint = "t0")

  cmp <- section_clone_compare(clono, run_dir, cfg, report, unit = "molecules")
  expect_true(all(c(10, 20, 30) %in% cmp$top_n))
  trk <- section_clone_tracking(clono, calls, run_dir, cfg, report, locus = "TRB")
  expect_true(nrow(trk$expanded) > 0)
  vj <- section_expanded_vj_circos(db, calls, "umi_count", "TRB", run_dir, cfg,
                                   report)
  expect_setequal(unique(vj$panel), c("timepoint", "baseline", "subject"))

  figs <- list.files(file.path(run_dir, "figures"))
  for (f in c("10a_clone_compare_top10_s1.pdf", "10a_clone_compare_top20_s2.pdf",
              "10b_clone_compare_t1_vs_t2_s1.pdf",
              "11a_baseline_tracking_s1.pdf",
              "11b_expanded_tracking_de_novo_s1.pdf",
              "11b_expanded_tracking_pre_existing_expanded_s2.pdf",
              "13a_vj_chord_t1.pdf", "13b_vj_chord_baseline.pdf",
              "13c_vj_chord_s2.pdf", "13d_vj_chord_combined.pdf")) {
    expect_true(f %in% figs, info = f)
  }
  tabs <- list.files(file.path(run_dir, "tables"))
  expect_true(all(c("10_clone_compare.csv", "11a_baseline_tracking.csv",
                    "11b_expanded_tracking.csv", "13_vj_chord_colours.csv",
                    "13_vj_pairs_expanded.csv") %in% tabs))
  if (requireNamespace("patchwork", quietly = TRUE)) {
    expect_true("10c_clone_compare_combined.pdf" %in% figs)
    expect_true("11c_expanded_tracking_combined.pdf" %in% figs)
  }
  legends <- readLines(file.path(run_dir, "FIGURE_LEGENDS.txt"))
  expect_true("13b_vj_chord_baseline" %in% legends)
})

test_that("single-cell sections run end to end on simulated data", {
  skip_if_not_installed("ggalluvial")
  skip_if_not_installed("circlize")
  run_dir <- trk_run_dir()
  cfg <- trk_cfg(tracking = list(highlight_compartment = "comp_c"),
                 circos = list(v_gene_colours = list(`TRBV2-1` = "#000000")))
  report <- open_report(file.path(run_dir, "REPORT.txt"), "test")
  on.exit(report$close(), add = TRUE)

  db <- trk_sim_sc()
  clono <- trk_sc_clono(db)
  meta <- unique(db[, c("sample_id", "subject_id", "condition")])

  cmp <- section_clone_compare(clono, run_dir, cfg, report, unit = "cells")
  expect_equal(sort(unique(cmp$sample_label)), c("comp_a", "comp_b", "comp_c"))
  fl <- section_compartment_tracking(clono, run_dir, cfg, report)
  expect_true(any(fl$highlighted))
  res <- section_vj_circos(db, "TRB", run_dir, cfg, report, sample_meta = meta)
  expect_length(res$circos, 3L)
  expect_true(all(res$usage$count[res$usage$in_top] > 0))
  # the colour mapping is shared: one colour per gene across samples
  expect_false(any(duplicated(names(res$colour_map$v_colors))))

  figs <- list.files(file.path(run_dir, "figures"))
  for (f in c("10a_clone_compare_top10_p1.pdf", "12_compartment_flows_p1.png",
              "14a_vj_circos_p1_comp_a.pdf", "14a_vj_circos_p1_comp_c.png",
              "14b_top_v_usage.pdf", "14c_top_v_heatmap.png")) {
    expect_true(f %in% figs, info = f)
  }
  tabs <- list.files(file.path(run_dir, "tables"))
  expect_true(all(c("12_compartment_flows.csv", "14_vj_colour_mapping.csv",
                    "14_top_v_genes.csv") %in% tabs))
})

test_that("tracking and circos settings are validated", {
  expect_error(tracking_config(trk_cfg(tracking = list(compare_top = 0))),
               "positive integers")
  expect_error(tracking_config(trk_cfg(tracking = list(expanded_rank_by = "max"))),
               "first")
  expect_error(tracking_config(trk_cfg(tracking = list(compare_pairs = list("a")))),
               "pairs")
  expect_error(circos_config(trk_cfg(circos = list(top_v = -1))),
               "positive integers")
  d <- tracking_config(trk_cfg())
  expect_equal(d$compare_top, c(10L, 20L))
  expect_equal(d$baseline_top, 15L)
  expect_equal(circos_config(trk_cfg())$top_v, 10L)
})
