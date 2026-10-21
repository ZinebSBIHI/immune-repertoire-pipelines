# Tests for the shared core: configuration, run directories, clonotypes,
# repertoire indices, rarefaction, differential abundance, paired statistics,
# convergence and single-cell BCR resolution.
#
# Each test pins down one analytical decision argued in docs/design_notes.md,
# so that a later change that silently reverses it fails here.

core_cfg <- function(extra = list()) {
  data_dir <- withr::local_tempdir(.local_envir = parent.frame())
  out_dir <- withr::local_tempdir(.local_envir = parent.frame())
  base <- list(
    project = list(name = "core", seed = 1L),
    paths = list(data_dir = data_dir, output_dir = out_dir),
    cohort = list(subjects = c("s1", "s2"), timepoints = c("t0", "t1"),
                  conditions = c("g1", "g2")),
    output = list(formats = "png", dpi = 40))
  utils::modifyList(base, extra)
}

# --- configuration ---------------------------------------------------------------

test_that("validate_config fills defaults and rejects inconsistencies", {
  cfg <- validate_config(core_cfg())
  expect_identical(cfg$repertoire$clonotype_definition, "strict")
  expect_identical(cfg$differential_abundance$fold_change_threshold, 2)
  expect_error(validate_config(core_cfg(list(typo_block = list()))),
               "Unknown configuration block")
  expect_silent(validate_config(core_cfg(list(specificity = list()))))
  # At a fold-change threshold of 1 expansion and contraction overlap.
  expect_error(validate_config(core_cfg(list(
    differential_abundance = list(fold_change_threshold = 1)))),
    "greater than 1")
  # A group confined to one subject is not convergence.
  expect_error(validate_config(core_cfg(list(
    convergence = list(min_subjects = 1)))), "at least 2")
  # Batch correction must name the variable it corrects over.
  expect_error(validate_config(core_cfg(list(
    gex = list(integration = list(method = "harmony"))))), "batch_variable")
  expect_error(validate_config(core_cfg(list(
    output = list(formats = c("png", "jpg"))))), "unsupported")
  expect_error(validate_config(core_cfg(list(
    cohort = list(subjects = c("s1", "s1"))))), "duplicates")
})

test_that("load_config resolves relative paths against the repository root", {
  root <- withr::local_tempdir()
  dir.create(file.path(root, "config"))
  dir.create(file.path(root, "data"))
  writeLines(c("project: {name: demo}",
               "paths: {data_dir: data, output_dir: results}",
               "cohort: {subjects: [s1]}"),
             file.path(root, "config", "config.yml"))
  cfg <- load_config(file.path(root, "config", "config.yml"))
  expect_identical(normalizePath(cfg$paths$data_dir),
                   normalizePath(file.path(root, "data")))
  expect_true(grepl("results$", cfg$paths$output_dir))
})

test_that("unexpected labels are errors, not NA", {
  expect_error(as_configured_factor(c("a", "c"), c("a", "b"), "condition"),
               "not declared")
  expect_identical(levels(as_configured_factor("b", c("a", "b"))),
                   c("a", "b"))
  expect_identical(display_labels(c("a", "b"), list(a = "Arm A")),
                   c("Arm A", "b"))
})

test_that("file names are parsed by capture group", {
  m <- parse_filename_metadata(c("x/P01_t0_g1.tsv", "x/P02_t1_g2.tsv"),
                               list(subject = "^(P[0-9]+)_",
                                    timepoint = "_(t[0-9]+)_"))
  expect_identical(m$subject, c("P01", "P02"))
  expect_identical(m$timepoint, c("t0", "t1"))
  expect_error(parse_filename_metadata("x/bad.tsv",
                                       list(subject = "^(P[0-9]+)_")),
               "Could not parse")
})

# --- run directories -------------------------------------------------------------

test_that("pruning removes only empty run directories of our own", {
  out <- withr::local_tempdir()
  ours_empty <- file.path(out, "run_a")
  ours_full <- file.path(out, "run_b")
  user_dir <- file.path(out, "manuscript_panels")
  for (d in c(ours_empty, ours_full, user_dir)) dir.create(d)
  marker <- immunerepertoire:::RUN_MARKER
  writeLines(marker, file.path(ours_empty, "RUN_INFO.txt"))
  writeLines(marker, file.path(ours_full, "RUN_INFO.txt"))
  writeLines("x", file.path(ours_full, "result.csv"))
  # A folder of the user's, empty apart from a file with the same name.
  writeLines("my notes", file.path(user_dir, "RUN_INFO.txt"))
  removed <- suppressMessages(immunerepertoire:::prune_empty_runs(out))
  expect_equal(removed, 1L)
  expect_false(dir.exists(ours_empty))
  expect_true(dir.exists(ours_full))
  expect_true(dir.exists(user_dir))
})

test_that("tables and reports are written with their manifest", {
  run_dir <- withr::local_tempdir()
  dir.create(file.path(run_dir, "tables"))
  save_table(data.frame(a = 1:2), "t1", run_dir, "two rows")
  man <- readLines(file.path(run_dir, "TABLE_MANIFEST.txt"))
  expect_true(grepl("t1.csv", man) && grepl("two rows", man))
  rep <- open_report(file.path(run_dir, "REPORT.txt"), "title")
  rep$h1("section"); rep$p("café"); rep$tbl(data.frame(x = pi))
  rep$close()
  txt <- readLines(file.path(run_dir, "REPORT.txt"), encoding = "UTF-8")
  expect_true("SECTION" %in% txt)
  expect_true(any(grepl("3.142", txt, fixed = TRUE)))
})

# --- clonotypes ------------------------------------------------------------------

core_db <- function() {
  data.frame(
    sample_id = c("a", "a", "a", "a", "b"),
    subject_id = "s1",
    locus = c("TRB", "TRB", "TRB", "TRA", "TRB"),
    productive = c(TRUE, TRUE, FALSE, TRUE, TRUE),
    v_call = c("TRBV1*01", "TRBV1*02", "TRBV2*01", "TRAV1*01", "TRBV1*01"),
    j_call = "TRBJ1-1*01",
    junction = c("TGT", "TGT", "TGC", "TGA", "TGT"),
    junction_aa = c("C", "C", "C", "X", "C"),
    umi_count = c(3, 2, 7, 4, 1),
    stringsAsFactors = FALSE)
}

test_that("clonotype keys strip alleles and follow the definition", {
  db <- core_db()
  expect_identical(clonotype_key(db[1:2, ], "strict"),
                   rep("TRBV1|TRBJ1-1|TGT", 2))
  expect_identical(clonotype_key(db[1:2, ], "strict", strip_allele = FALSE),
                   c("TRBV1*01|TRBJ1-1*01|TGT", "TRBV1*02|TRBJ1-1*01|TGT"))
  expect_identical(clonotype_key(db[1, ], "aa"), "C")
  expect_error(clonotype_key(db[, c("v_call", "j_call")], "nt"), "junction")
})

test_that("collapse_clonotypes filters before computing frequencies", {
  cl <- collapse_clonotypes(core_db(), "umi_count", locus = "TRB",
                            min_abundance = 2)
  # Sample a: the two TRBV1 rows merge (5); the non-productive row and the
  # TRA row are dropped. Sample b falls below the floor.
  expect_equal(nrow(cl), 1L)
  expect_equal(cl$abundance, 5)
  expect_equal(cl$frequency, 1)
  expect_warning(collapse_clonotypes(core_db()[, -3], "umi_count",
                                     locus = "TRB"), "NO locus filtering")
  expect_error(collapse_clonotypes(core_db(), "reads"), "not interchangeable")
})

# --- indices ---------------------------------------------------------------------

test_that("diversity indices take their documented values", {
  expect_equal(shannon_entropy(c(1, 1, 1, 1)), log(4))
  expect_equal(shannon_clonality(c(1, 1, 1, 1)), 0)
  expect_true(is.na(shannon_clonality(5)))
  expect_equal(inverse_simpson(c(1, 1)), 2)
  expect_equal(gini_index(c(1, 1, 1)), 0)
  expect_equal(gini_index(c(0, 0, 10, 0, 1e-9)), gini_index(c(10, 1e-9)))
  # Chao1 needs singletons or doubletons; without them it is undefined.
  expect_true(is.na(chao1_richness(c(3, 4, 5))))
  expect_equal(chao1_richness(c(1, 1, 2, 5)), 4 + 2^2 / (2 * 1))
  # D50 as a percentage of the observed clonotypes.
  expect_equal(d50_index(c(50, 25, 25)), 100 / 3)
})

test_that("rarefaction samples to the minimum depth within a subject", {
  clono <- data.frame(sample_id = rep(c("a", "b", "c"), c(3, 3, 2)),
                      subject_id = rep(c("s1", "s1", "s2"), c(3, 3, 2)),
                      clonotype = c("x", "y", "z", "x", "y", "w", "x", "y"),
                      abundance = c(10, 5, 5, 40, 30, 30, 3, 3))
  set.seed(1)
  r <- rarefy_to_common_depth(clono, within = "subject_id", n_bootstrap = 5)
  depth <- unique(r[, c("sample_id", "target_depth")])
  expect_equal(depth$target_depth[match(c("a", "b", "c"), depth$sample_id)],
               c(20, 20, 6))
  clono$abundance[1] <- 10.5
  expect_error(rarefy_to_common_depth(clono, n_bootstrap = 2), "integer")
})

# --- differential abundance -------------------------------------------------------

test_that("categories and significance are separate axes", {
  clono <- data.frame(
    sample_id = rep(c("b", "p"), c(3, 3)), subject_id = "s1",
    timepoint = rep(c("t0", "t1"), c(3, 3)),
    clonotype = c("A", "B", "C", "A", "B", "D"),
    abundance = c(500, 500, 9000, 2000, 10, 1))
  clono$frequency <- clono$abundance / ave(clono$abundance, clono$sample_id,
                                           FUN = sum)
  res <- test_clonal_abundance(clono, reference_timepoint = "t0")
  cat_of <- function(k) res$category[res$clonotype == k]
  expect_identical(cat_of("A"), "pre_existing_expanded")
  expect_identical(cat_of("B"), "contracting")
  expect_identical(cat_of("C"), "lost")
  expect_identical(cat_of("D"), "de_novo")
  # One molecule absent at baseline: de novo but not significant.
  expect_false(res$sig_expansion[res$clonotype == "D"])
  expect_true(res$sig_expansion[res$clonotype == "A"])
  expect_identical(res$status[res$clonotype == "D"], "not_significant:de_novo")
  expect_error(test_clonal_abundance(clono, reference_timepoint = "t9"),
               "not present")
})

# --- paired statistics ------------------------------------------------------------

test_that("the attainable p-value bound governs small paired designs", {
  expect_equal(attainable_p_paired(3:6), c(0.25, 0.125, 0.0625, 0.03125))
  expect_false(test_is_powered(5))
  expect_true(test_is_powered(6))
  df <- data.frame(subject = rep(paste0("s", 1:4), 2),
                   group = rep(c("g1", "g2"), each = 4),
                   value = c(1, 2, 3, 4, 2, 3.5, 5, 8))
  r <- paired_comparison(df, "value", "group", "subject", c("g1", "g2"))
  expect_true(r$exact_test)
  expect_equal(r$attainable_p, 0.125)
  expect_equal(r$p_value, 0.125)
  expect_false(r$powered)
  expect_true(r$direction_concordant)
  expect_identical(r$direction, "increase")
  # Zero differences reduce the effective n, and the bound follows.
  df$value[5] <- 1
  r0 <- paired_comparison(df, "value", "group", "subject", c("g1", "g2"))
  expect_equal(r0$n_effective, 3L)
  # Group labels containing dots pair correctly.
  df$group <- rep(c("Day.0", "Day.7"), each = 4)
  r1 <- paired_comparison(df, "value", "group", "subject", c("Day.0", "Day.7"))
  expect_equal(r1$n_pairs, 4L)
})

# --- convergence -----------------------------------------------------------------

test_that("identity is normalised and floored at zero", {
  m <- cdr3_identity_matrix(c("CASSF", "CASSLF"), c("CASSF", "CQ"))
  expect_equal(dim(m), c(2L, 2L))
  expect_equal(m[1, 1], 1)
  expect_equal(m[2, 1], 1 - 1 / 5)
  expect_true(all(m >= 0))
})

test_that("public groups span subjects and are read against a null", {
  df <- data.frame(
    subject_id = c("s1", "s2", "s3", "s1", "s1"),
    v_gene = c("V1", "V1", "V1", "V2", "V2"),
    junction_aa = c("CASSLGQAYEQYF", "CASSLGQAYEQFF", "CASSLGRAYEQYF",
                    "CAWSVGGYTF", "CAWSVGGYTL"),
    stringsAsFactors = FALSE)
  pub <- detect_public_groups(df, identity_threshold = 0.8, min_subjects = 2)
  expect_equal(nrow(pub$summary), 1L)
  expect_equal(pub$summary$n_subjects, 3L)
  # The V2 pair is within one subject: related, not convergent.
  expect_false(any(pub$groups$v_gene == "V2"))
  set.seed(1)
  perm <- permute_public_groups(df, identity_threshold = 0.8,
                                min_subjects = 2, n_permutations = 20)
  expect_equal(perm$observed, 1L)
  expect_length(perm$null, 20L)
})

# --- single-cell BCR ---------------------------------------------------------------

test_that("an empty isotype call is not switched", {
  expect_identical(is_switched(c("IGHG1", "IGHM", "IGHD", "", NA)),
                   c(TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_identical(expansion_threshold(c(TRUE, FALSE), 2L, 4L), c(4L, 2L))
})

test_that("doublets are removed and switched chains preferred", {
  db <- data.frame(
    cell_id = c("c1", "c1", "c2", "c2", "c3", "c3", "c1"),
    locus = c("IGH", "IGH", "IGH", "IGH", "IGH", "IGK", "IGK"),
    c_call = c("IGHM", "IGHG1", "IGHM", "IGHA1", "IGHD", "IGKC", "IGKC"),
    junction = c("J1", "J2", "J3", "J4", "J5", "K1", "K2"),
    umi_count = c(20, 6, 30, 25, 8, 9, 5),
    productive = TRUE, stringsAsFactors = FALSE)
  res <- resolve_one_heavy_per_cell(db, doublet_min_umi = 10)
  heavy <- res$db[res$db$locus == "IGH", ]
  # c2 has two well-supported heavy chains: a doublet, removed.
  expect_true("c2" %in% res$doublets)
  expect_false("c2" %in% heavy$cell_id)
  # c1: the switched chain is preferred over the more abundant IGHM.
  expect_identical(heavy$c_call[heavy$cell_id == "c1"], "IGHG1")
  expect_equal(anyDuplicated(heavy$cell_id), 0L)
})
