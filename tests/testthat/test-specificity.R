# Tests for the specificity annotation (ana-specificity-db.R,
# ana-specificity.R, sc-specificity.R) and the receptor embedding
# (sc-embedding.R).
#
# Reference files are written here with the column layouts of the published
# exports (VDJdb release and web search export, McPAS-TCR, an antibody table
# laid out as the CoV-AbDab CSV); their rows are invented. Repertoires are
# simulated.

# --- helpers -------------------------------------------------------------------

spc_run <- function(extra = list()) {
  data_dir <- withr::local_tempdir(.local_envir = parent.frame())
  out_dir <- withr::local_tempdir(.local_envir = parent.frame())
  ref_dir <- withr::local_tempdir(.local_envir = parent.frame())
  base <- list(
    project = list(name = "test", seed = 1L),
    paths = list(data_dir = data_dir, output_dir = out_dir,
                 reference_dir = ref_dir),
    cohort = list(subjects = c("s1", "s2"), timepoints = c("t0", "t1"),
                  conditions = c("g1", "g2")),
    repertoire = list(locus = "TRB", abundance_column = "umi_count"),
    differential_abundance = list(reference_timepoint = "t0",
                                  fold_change_threshold = 2),
    output = list(formats = "png", dpi = 40))
  cfg <- validate_config(utils::modifyList(base, extra))
  run_dir <- init_run(cfg, "test")
  list(cfg = cfg, run_dir = run_dir, ref_dir = ref_dir,
       report = open_report(file.path(run_dir, "REPORT.txt"), "test"))
}

spc_write_vdjdb_release <- function(path, rows) {
  cols <- c("gene", "cdr3", "species", "antigen.epitope", "antigen.gene",
            "antigen.species", "complex.id", "v.segm", "j.segm", "mhc.a",
            "mhc.b", "mhc.class", "reference.id", "vdjdb.score")
  d <- as.data.frame(matrix("", nrow(rows), length(cols),
                            dimnames = list(NULL, cols)),
                     stringsAsFactors = FALSE)
  for (nm in names(rows)) d[[nm]] <- rows[[nm]]
  utils::write.table(d, path, sep = "\t", quote = FALSE, row.names = FALSE)
  path
}

spc_vdjdb_rows <- function() {
  data.frame(
    gene = c("TRB", "TRB", "TRB", "TRA", "TRB", "TRB"),
    cdr3 = c("CASSLAPGATNEKLFF", "CASSPDRGGYEQYF", "CASSPDRGGYEQYF",
             "CAVRDSNYQLIW", "CASSIRSSYEQYF", "CASRGQGNYGYTF"),
    species = c("HomoSapiens", "HomoSapiens", "HomoSapiens", "HomoSapiens",
                "HomoSapiens", "MusMusculus"),
    antigen.epitope = c("NLVPMVATV", "GILGFVFTL", "TTDPSFLGRY", "GLCTLVAML",
                        "KLGGALQAK", "SSLENFRAYV"),
    antigen.gene = c("pp65", "M", "ORF1ab", "BMLF1", "IE1", "PA"),
    antigen.species = c("CMV", "InfluenzaA", "SARS-CoV-2", "EBV", "CMV",
                        "InfluenzaA"),
    v.segm = c("TRBV7-6*01", "TRBV19*01", "TRBV19*01", "TRAV12-1*01",
               "TRBV19*01", "TRBV2*01"),
    mhc.a = c("HLA-A*02", "HLA-A*02", "HLA-A*02", "HLA-A*02", "HLA-A*03",
              "H-2Db"),
    stringsAsFactors = FALSE)
}

# --- readers -------------------------------------------------------------------

test_that("read_vdjdb reads the release layout and filters host and chain", {
  f <- spc_write_vdjdb_release(tempfile(fileext = ".tsv"), spc_vdjdb_rows())
  ref <- read_vdjdb(f, chain = "TRB", host_species = "HomoSapiens")
  expect_identical(names(ref), c("cdr3", "chain", "v_gene", "j_gene",
                                 "epitope", "epitope_gene", "epitope_species",
                                 "mhc", "tissue", "cell_type", "source"))
  expect_equal(nrow(ref), 4L)
  expect_true(all(ref$chain == "TRB"))
  expect_true(all(ref$source == "VDJdb"))
  # The host column is never taken as the antigen species.
  expect_false(any(ref$epitope_species == "HomoSapiens"))
  expect_setequal(ref$epitope_species, c("CMV", "InfluenzaA", "SARS-CoV-2"))
  expect_equal(nrow(read_vdjdb(f, chain = NULL, host_species = NULL)), 6L)
})

test_that("read_vdjdb reads the web export, JSON columns included", {
  f <- tempfile(fileext = ".tsv")
  header <- c("complex.id", "Gene", "CDR3", "V", "J", "Species", "MHC A",
              "MHC B", "MHC class", "Epitope", "Epitope gene",
              "Epitope species", "Reference", "Method", "Meta", "CDR3fix",
              "Score")
  meta <- '{"cell.subset": "CD8", "clone.id": "", "donor.MHC": "HLA-A*02"}'
  row1 <- c("1", "TRB", "CASSLAPGATNEKLFF", "TRBV7-6*01", "TRBJ1-4*01",
            "HomoSapiens", "HLA-A*02:01", "B2M", "MHCI", "NLVPMVATV", "pp65",
            "CMV", "PMID:1", '{"frequency": "1/11"}', meta,
            '{"cdr3": "CASSLAPGATNEKLFF", "good": true}', "2")
  row2 <- replace(row1, c(2, 3, 12), c("TRA", "CAVRDSNYQLIW", "EBV"))
  writeLines(c(paste(header, collapse = "\t"), paste(row1, collapse = "\t"),
               paste(row2, collapse = "\t")), f)
  ref <- read_vdjdb(f, chain = "TRB")
  expect_equal(nrow(ref), 1L)
  expect_identical(ref$cdr3, "CASSLAPGATNEKLFF")
  expect_identical(ref$epitope_species, "CMV")
  expect_identical(ref$epitope_gene, "pp65")
  expect_identical(ref$mhc, "HLA-A*02:01")
})

test_that("read_mcpas reads both chains and survives Windows-1252 bytes", {
  f <- tempfile(fileext = ".csv")
  header <- paste0('"', c("CDR3.alpha.aa", "CDR3.beta.aa", "Species",
                          "Category", "Pathology", "Antigen.protein",
                          "Epitope.peptide", "MHC", "Tissue", "T.Cell.Type",
                          "TRAV", "TRAJ", "TRBV", "TRBJ"), '"',
                   collapse = ",")
  l1 <- paste0('"CAVSDLEPNSSASKIIF","CASSIRSSYEQYF","Human","Pathogens",',
               '"Cytomegalovirus (CMV)","IE1","KLGGALQAK","HLA-A*03",',
               '"PBMC","CD8","TRAV8-6","TRAJ3","TRBV19","TRBJ2-7"')
  l2 <- paste0('"","CASSPDRGGYEQYF","Mouse","Pathogens","Influenza",',
               '"M","GILGFVFTL","HLA-A*02","PBMC","CD8","","","TRBV19",',
               '"TRBJ2-7"')
  con <- file(f, "wb")
  writeLines(c(header, l1, l2), con, sep = "\n")
  # A dash encoded in Windows-1252 (0x96), as found in the published file.
  writeBin(c(charToRaw('"","CASRZZZ","Human","Pathogens","EBV","X","Y","Z",'),
             charToRaw('"PBMC","CD8","","","TRBV6'), as.raw(0x96),
             charToRaw('2","TRBJ1-1"\n')), con)
  close(con)
  trb <- read_mcpas(f, chain = "TRB", host_species = "Human")
  expect_equal(nrow(trb), 2L)
  expect_true(all(trb$chain == "TRB"))
  expect_identical(trb$epitope_species[1], "Cytomegalovirus (CMV)")
  expect_true(all(validUTF8(trb$v_gene)))
  tra <- read_mcpas(f, chain = "TRA", host_species = "Human")
  expect_identical(tra$cdr3, "CAVSDLEPNSSASKIIF")
})

test_that("read_antibody_reference finds the columns as the original did", {
  f <- tempfile(fileext = ".csv")
  d <- data.frame(Name = c("ab1", "ab2", "ab3", "ab4"),
                  `Ab or Nb` = "Ab",
                  `Heavy V Gene` = c("IGHV3-53 (Human)", "IGHV1-69 (Human)",
                                     "IGHV3-30-3 (Human)", "IGHV3-53 (Human)"),
                  CDRH3 = c("AREGEVEGYYYYGMDV", "ARDLLLGFDY", "ND", ""),
                  check.names = FALSE, stringsAsFactors = FALSE)
  utils::write.csv(d, f, row.names = FALSE)
  ref <- read_antibody_reference(f)
  expect_identical(ref$ab_name, c("ab1", "ab2"))
  expect_identical(ref$ab_vgene, c("IGHV3-53", "IGHV1-69"))
  expect_identical(immunerepertoire:::ighv_gene("IGHV3-30-3*01"), "IGHV3-30")
})

# --- matching ------------------------------------------------------------------

test_that("collapse_reference keeps one row per CDR3", {
  f <- spc_write_vdjdb_release(tempfile(fileext = ".tsv"), spc_vdjdb_rows())
  lk <- collapse_reference(read_vdjdb(f))
  expect_equal(anyDuplicated(lk$cdr3), 0L)
  two <- lk[lk$cdr3 == "CASSPDRGGYEQYF", ]
  expect_identical(two$specificities, "InfluenzaA | SARS-CoV-2")
  expect_equal(two$n_records, 2L)
})

test_that("specificity_status applies the categories in order", {
  s <- c("CMV | EBV", "EBV", "InfluenzaA", "SARS-CoV-2", "HIV-1", NA,
         "CMV;SARS-CoV-2")
  st <- specificity_status(s)
  expect_identical(as.character(st),
                   c("Bystander (CMV)", "Bystander (EBV)",
                     "Bystander (Influenza)", "Bystander (SARS-CoV-2)",
                     "Known specificity (other)", "Unannotated",
                     "Bystander (CMV)"))
  tg <- specificity_status(s, target_species = "SARS",
                           target_label = "SARS-CoV-2")
  expect_identical(as.character(tg)[c(4, 7)], rep("Target (SARS-CoV-2)", 2))
  expect_identical(levels(tg)[1], "Target (SARS-CoV-2)")
})

test_that("a target species leaves the bystander list", {
  cfg <- list(paths = list(reference_dir = "ref"),
              specificity = list(target_species = "SARS"))
  sp <- immunerepertoire:::specificity_config(cfg)
  expect_false("SARS-CoV-2" %in% names(sp$bystanders))
  expect_identical(sp$target_label, "SARS")
  expect_null(immunerepertoire:::specificity_config(
    list(paths = list(reference_dir = "ref")))$vdjdb_file)
  sp2 <- immunerepertoire:::specificity_config(
    list(paths = list(reference_dir = "ref"),
         specificity = list(vdjdb_file = "vdjdb/x.tsv")))
  expect_identical(sp2$vdjdb_file, file.path("ref", "vdjdb/x.tsv"))
  expect_error(immunerepertoire:::specificity_config(
    list(specificity = list(single_cell_reference = "other"))), "one of")
})

test_that("clonotype_cdr3_map is one-to-one for junction-based keys", {
  db <- data.frame(v_call = c("TRBV1*01", "TRBV1*01", "TRBV2*01"),
                   j_call = "TRBJ1-1*01",
                   junction = c("TGTGCC", "TGTGCA", "TGTGCC"),
                   junction_aa = c("CA", "CA", "CA"),
                   stringsAsFactors = FALSE)
  m <- clonotype_cdr3_map(db, "strict")
  expect_equal(nrow(m), 3L)
  expect_true(all(m$junction_aa == "CA"))
  g <- clonotype_cdr3_map(rbind(db, transform(db[1, ], junction_aa = "CB")),
                          "gene")
  expect_true(is.na(g$junction_aa[g$clonotype == "TRBV1|TRBJ1-1"]))
})

test_that("match_antibody_reference requires the V gene and the identity", {
  ref <- data.frame(ab_name = c("a", "b", "c"),
                    ab_vgene = c("IGHV3-53", "IGHV3-53", "IGHV1-69"),
                    ab_cdr3 = c("AREGEVEGYYYYGMDV", "ARDLLLGFDY",
                                "AREGEVEGYYYYGMDV"),
                    stringsAsFactors = FALSE)
  q <- data.frame(public_group = c("PG1", "PG2", "PG3"),
                  v_gene = c("IGHV3-53", "IGHV1-69", "IGHV4-34"),
                  subject_id = "s1",
                  junction_aa = c("CAREGEVEGYYYYGMDVW", "CARDLLLGFDYW",
                                  "CAREGEVEGYYYYGMDVW"),
                  stringsAsFactors = FALSE)
  h <- match_antibody_reference(q, ref, min_identity = 0.80)
  # PG1: same V gene, identical CDRH3 once the anchors are removed.
  # PG2: the identical CDRH3 is under another V gene, so no hit.
  # PG3: no reference antibody with that V gene.
  expect_identical(h$public_group, "PG1")
  expect_identical(h$ab_hit, "a")
  expect_equal(h$identity, 1)
  near <- q[1, ]
  near$junction_aa <- "CAREGEVEGYYYYAMDVW"
  expect_equal(match_antibody_reference(near, ref)$identity, 1 - 1 / 16)
})

# --- single-cell annotation ------------------------------------------------------

test_that("annotate_cells_by_cdr3 propagates only from exact matches", {
  lookup <- data.frame(cdr3 = c("CASSAAAF", "CASSQQQF", "CASSZZZZZF"),
                       specificities = c("CMV", "EBV", "HIV"),
                       epitopes = NA, epitope_genes = NA, sources = "x",
                       tissue = NA, cell_type = NA, n_records = 1L,
                       stringsAsFactors = FALSE)
  cdr3 <- c("CASSAAAF",   # exact (anchor 1)
            "CASSAAAG",   # 1 edit from anchor 1
            "CASSAQQF",   # 1 edit from anchor 2, 2 edits from anchor 1
            "CASSQQQF",   # exact (anchor 2)
            "CASSZZZZYF", # 1 edit from a reference CDR3 no cell carries
            NA, "NA")
  a0 <- annotate_cells_by_cdr3(cdr3, lookup, 0)
  expect_identical(a0$match, c("exact", NA, NA, "exact", NA, NA, NA))
  a1 <- annotate_cells_by_cdr3(cdr3, lookup, 1)
  expect_identical(a1$Epitope.species,
                   c("CMV", "CMV", "EBV", "EBV", NA, NA, NA))
  a2 <- annotate_cells_by_cdr3(cdr3, lookup, 2)
  # Within reach of both anchors: the first anchor in cell order wins, even
  # when the second is nearer -- the rule of Trex::annotateDB().
  expect_identical(a2$Epitope.species[3], "CMV")
  expect_identical(a2$match[3], "distance")
  expect_identical(a2$matched_cdr3[3], "CASSAAAF")
  # A reference CDR3 absent from the data is never an anchor.
  expect_true(is.na(a2$Epitope.species[5]))
})

test_that("annotate_cells_by_cdr3 reproduces Trex::annotateDB", {
  skip_if_not_installed("Trex")
  skip_if_not_installed("Seurat")
  env <- new.env()
  utils::data("trex_example", package = "Trex", envir = env)
  obj <- env$trex_example
  lookup <- collapse_reference(trex_reference("TRB"))
  trb <- vapply(strsplit(as.character(obj$CTaa), "_", fixed = TRUE),
                function(z) if (length(z) >= 2) z[2] else NA_character_,
                character(1))
  # annotateDB() reads Trex.database from the search path.
  withr::local_package("Trex")
  for (k in c(0, 1, 2)) {
    ref <- suppressWarnings(Trex::annotateDB(obj, chains = "TRB",
                                             edit.distance = k))
    ours <- annotate_cells_by_cdr3(trb, lookup, max_distance = k)
    expect_identical(ours$Epitope.species,
                     as.character(ref$TRB_Epitope.species),
                     info = paste("edit distance", k))
  }
})

test_that("cell_receptor_summary writes the scRepertoire layout", {
  db <- data.frame(
    cell_id = c("c1", "c1", "c1", "c2", "c3"),
    sample_id = "s1", subject_id = "p1",
    locus = c("TRB", "TRA", "TRB", "TRB", "TRA"),
    junction_aa = c("CASSA", "CAVB", "CASSC", "CASSD", "CAVE"),
    v_call = c("TRBV9*01", "TRAV14/DV4*01", "TRBV2*01", "TRBV5-1", "TRAV1"),
    d_call = c(NA, NA, "TRBD1", NA, NA),
    j_call = c("TRBJ2-7*01", "TRAJ33", "TRBJ1-1", "TRBJ2-1", "TRAJ1"),
    c_call = c("TRBC2", "TRAC", "TRBC1", NA, "TRAC"),
    umi_count = c(3, 5, 9, 2, 4), productive = TRUE,
    stringsAsFactors = FALSE)
  s <- cell_receptor_summary(db)
  expect_identical(s$cell_id, c("c1", "c2", "c3"))
  # Two TRB chains in c1, the more abundant first.
  expect_identical(s$TRB_cdr3[1], "CASSC;CASSA")
  expect_identical(s$CTaa, c("CAVB_CASSC;CASSA", "NA_CASSD", "CAVE_NA"))
  expect_identical(s$CTgene[1],
                   "TRAV14/DV4.TRAJ33.TRAC_TRBV2.TRBD1.TRBJ1-1.TRBC1;TRBV9.None.TRBJ2-7.TRBC2")
})

test_that("clone size categories follow the original breaks", {
  x <- clone_size_category(c(1, 2, 3, 5, 6, 10, 11, 54))
  expect_identical(as.character(x),
                   c("Single", "Small", "Medium", "Medium", "Large", "Large",
                     "Hyperexpanded", "Hyperexpanded"))
  expect_identical(attr(x, "labels")[["Medium"]], "Medium (2 < X <= 5)")
  expect_identical(attr(x, "labels")[["Hyperexpanded"]],
                   "Hyperexpanded (10 < X <= 54)")
  g <- cell_specificity_groups(x, c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE,
                                    TRUE, FALSE), "CMV")
  expect_identical(as.character(g$clone_expansion_group[c(1, 3, 7, 8)]),
                   c("Other", "Medium clones", "Hyperexpanded CMV",
                     "Hyperexpanded clones"))
  expect_identical(as.character(g$expanded_target_group[c(2, 3, 5)]),
                   c("Other", "Expanded + CMV", "Expanded clone"))
})

test_that("V usage and top CDR3 of a set of cells", {
  v <- cell_v_usage(c("TRBV9;TRBV2", "TRBV9", NA, "NA"), top_n = 1)
  expect_identical(v$v_gene, "TRBV9")
  expect_equal(v$freq, 2 / 3)
  t <- top_cell_cdr3(c("A", "A", "B", "C", NA), top_n = 2)
  expect_identical(t$cdr3, c("A", "B"))
  expect_equal(t$freq_within_top, c(2, 1) / 3)
})

# --- sections ------------------------------------------------------------------

spc_sim_bulk <- function(seed = 2) {
  set.seed(seed)
  aa <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
  cdr3 <- function(n) vapply(seq_len(n), function(i) paste0(
    "CASS", paste(sample(aa, 7, replace = TRUE), collapse = ""), "F"),
    character(1))
  rows <- list()
  for (s in c("s1", "s2")) {
    pool <- cdr3(150)
    w <- stats::rlnorm(150, 0, 1.2)
    expd <- 1:8
    for (tp in c("t0", "t1")) {
      ww <- w
      if (tp == "t1") ww[expd] <- ww[expd] * 20
      n <- stats::rmultinom(1, 5000, ww / sum(ww))[, 1]
      k <- n > 0
      rows[[length(rows) + 1L]] <- data.frame(
        sample_id = paste(s, tp, sep = "_"), subject_id = s, timepoint = tp,
        locus = "TRB", productive = TRUE, v_call = "TRBV19*01",
        j_call = "TRBJ2-7*01", junction = pool[k], junction_aa = pool[k],
        umi_count = n[k], stringsAsFactors = FALSE)
    }
  }
  db <- do.call(rbind, rows)
  db$subject_id <- factor(db$subject_id, levels = c("s1", "s2"))
  db$timepoint <- factor(db$timepoint, levels = c("t0", "t1"))
  db
}

test_that("section_specificity runs end to end on a bulk repertoire", {
  run <- spc_run(list(specificity = list(vdjdb_file = "vdjdb.tsv")))
  on.exit(run$report$close(), add = TRUE)
  db <- spc_sim_bulk()
  exp_cdr3 <- unique(db$junction_aa[db$subject_id == "s1"])[1:3]
  rows <- spc_vdjdb_rows()[c(1, 2, 3), ]
  rows$cdr3 <- exp_cdr3
  spc_write_vdjdb_release(file.path(run$ref_dir, "vdjdb.tsv"), rows)

  clono <- collapse_clonotypes(db, "umi_count", definition = "strict")
  clono <- merge(clono[, setdiff(names(clono), "subject_id")],
                 unique(db[, c("sample_id", "subject_id", "timepoint")]),
                 by = "sample_id")
  calls <- test_clonal_abundance(clono, reference_timepoint = "t0")
  res <- section_specificity(db, "umi_count", "TRB", run$run_dir, run$cfg,
                             run$report, calls = calls)
  expect_true(all(c("per_sample", "composition", "status") %in% names(res)))
  expect_true(any(res$per_sample$n_cdr3_annotated > 0))
  expect_setequal(unique(res$composition$epitope_species),
                  c("CMV", "InfluenzaA", "SARS-CoV-2"))
  fr <- stats::aggregate(fraction ~ subject_id + timepoint,
                         res$composition, sum)
  expect_true(all(abs(fr$fraction - 1) < 1e-9))
  files <- list.files(run$run_dir, recursive = TRUE)
  expect_true(any(grepl("15_species_composition.png", files)))
  expect_true(any(grepl("15_volcano_annotation.png", files)))
  expect_true(any(grepl("15_annotation_per_sample.csv", files)))
  if (nrow(res$status)) {
    expect_true(all(res$status$sig_expansion))
    expect_s3_class(res$status$specificity_status, "factor")
  }
  # Not applicable to an immunoglobulin locus.
  expect_null(section_specificity(db, "umi_count", "IGH", run$run_dir,
                                  run$cfg, run$report))
})

test_that("section_specificity is skipped without a reference file", {
  run <- spc_run()
  on.exit(run$report$close(), add = TRUE)
  expect_null(section_specificity(spc_sim_bulk(), "umi_count", "TRB",
                                  run$run_dir, run$cfg, run$report))
})

spc_sim_sc <- function(seed = 3) {
  set.seed(seed)
  n <- 120
  cells <- sprintf("s%d_cell%03d", rep(1:2, each = n / 2), seq_len(n))
  trb <- c(rep("CASSLAPGATNEKLFF", 12), rep("CASSLAPGATNEKLYF", 4),
           rep("CASSPDRGGYEQYF", 6), sprintf("CASSX%03dF", seq_len(n - 22)))
  tra <- c(rep("CAVRDSNYQLIW", 5), sprintf("CAVX%03dW", seq_len(n - 5)))
  rows <- rbind(
    data.frame(cell_id = cells, locus = "TRB", junction_aa = trb,
               v_call = "TRBV7-6*01", j_call = "TRBJ1-4*01",
               umi_count = 5, stringsAsFactors = FALSE),
    data.frame(cell_id = cells, locus = "TRA", junction_aa = tra,
               v_call = "TRAV12-1*01", j_call = "TRAJ33*01",
               umi_count = 3, stringsAsFactors = FALSE))
  rows$junction <- rows$junction_aa
  rows$productive <- TRUE
  rows$sample_id <- sub("_cell.*$", "", rows$cell_id)
  rows$subject_id <- rows$sample_id
  rows$condition <- ifelse(rows$sample_id == "s1", "g1", "g2")
  rows
}

test_that("section_cell_specificity annotates, flags and groups cells", {
  run <- spc_run(list(
    cohort = list(subjects = c("s1", "s2"), conditions = c("g1", "g2")),
    specificity = list(vdjdb_file = "vdjdb.tsv", mcpas_file = NULL,
                       single_cell_reference = "vdjdb",
                       target_species = "CMV")))
  on.exit(run$report$close(), add = TRUE)
  spc_write_vdjdb_release(file.path(run$ref_dir, "vdjdb.tsv"),
                          spc_vdjdb_rows())
  db <- spc_sim_sc()
  cells <- section_cell_specificity(db, "umi_count", run$run_dir, run$cfg,
                                    run$report)
  expect_equal(nrow(cells), 120L)
  # 12 + 6 exact matches (two reference CDR3s); the 4 one-edit neighbours
  # of the first join at d = 1 and d = 2.
  expect_equal(sum(cells$TRB_match_d1 == "exact", na.rm = TRUE), 18L)
  expect_equal(sum(cells$TRB_match_d1 == "distance", na.rm = TRUE), 4L)
  expect_equal(sum(cells$TRB_target_d1), 16L)
  expect_equal(sum(cells$TRB_target_d2), 16L)
  expect_true(all(cells$TRA_target_d2 == FALSE))
  expect_equal(sum(cells$target_any), 16L)
  # The 12 identical receptors form one Hyperexpanded clone in sample s1.
  expect_true(all(cells$clone_size[cells$TRB_cdr3 == "CASSLAPGATNEKLFF"] ==
                    "Hyperexpanded"))
  expect_true("Hyperexpanded CMV" %in% cells$clone_expansion_group)
  bars <- utils::read.csv(file.path(run$run_dir, "tables",
                                    "15_expanded_target_by_condition.csv"))
  expect_identical(bars$condition, c("g1", "g2"))
  expect_equal(sum(bars$count), 16L)
  files <- list.files(run$run_dir, recursive = TRUE)
  expect_true(any(grepl("15_target_TRBV_usage.png", files)))
  expect_true(any(grepl("15_expanded_target_percentage.png", files)))
})

test_that("section_cell_specificity draws embedding figures", {
  skip_if_not_installed("Seurat")
  run <- spc_run(list(
    cohort = list(subjects = c("s1", "s2"), conditions = c("g1", "g2")),
    specificity = list(vdjdb_file = "vdjdb.tsv",
                       single_cell_reference = "vdjdb",
                       target_species = "CMV")))
  on.exit(run$report$close(), add = TRUE)
  spc_write_vdjdb_release(file.path(run$ref_dir, "vdjdb.tsv"),
                          spc_vdjdb_rows())
  db <- spc_sim_sc()
  db <- db[db$cell_id %in% unique(db$cell_id)[1:80], ]
  obj <- SeuratObject::pbmc_small
  obj <- SeuratObject::RenameCells(obj, new.names = unique(db$cell_id))
  um <- matrix(stats::rnorm(160), 80, 2,
               dimnames = list(colnames(obj), c("umap_1", "umap_2")))
  obj[["umap"]] <- SeuratObject::CreateDimReducObject(um, key = "umap_",
                                                       assay = "RNA")
  section_cell_specificity(db, "umi_count", run$run_dir, run$cfg, run$report,
                           object = obj)
  files <- list.files(run$run_dir, recursive = TRUE)
  expect_true(any(grepl("15_umap_expansion_groups.png", files)))
  expect_true(any(grepl("15_umap_expanded_target.png", files)))
  expect_true(file.exists(file.path(run$run_dir,
                                    "seurat_with_specificity.rds")))
})

test_that("section_antibody_lookup reports the hits", {
  run <- spc_run(list(specificity = list(antibody_reference_file = "ab.csv")))
  on.exit(run$report$close(), add = TRUE)
  utils::write.csv(data.frame(Name = c("ab1", "ab2"),
                              `Heavy V Gene` = c("IGHV3-53 (Human)",
                                                 "IGHV1-69 (Human)"),
                              CDRH3 = c("AREGEVEGYYYYGMDV", "ARDLLLGFDY"),
                              check.names = FALSE),
                   file.path(run$ref_dir, "ab.csv"), row.names = FALSE)
  groups <- data.frame(subject_id = c("s1", "s2"), v_gene = "IGHV3-53",
                       junction_aa = c("CAREGEVEGYYYYGMDVW",
                                       "CAREGEVEGYYFYGMDVW"),
                       public_group = "PG0001", stringsAsFactors = FALSE)
  hits <- section_antibody_lookup(groups, run$run_dir, run$cfg, run$report)
  expect_equal(nrow(hits), 2L)
  expect_true(file.exists(file.path(run$run_dir, "tables",
                                    "07c_antibody_hits_public.csv")))
})

# --- embedding -------------------------------------------------------------------

test_that("multimodal_neighbours builds the joint graph", {
  skip_if_not_installed("Seurat")
  obj <- SeuratObject::pbmc_small
  set.seed(1)
  m <- matrix(stats::rnorm(ncol(obj) * 10), ncol(obj), 10,
              dimnames = list(colnames(obj), paste0("Trex_", 1:10)))
  obj[["Trex"]] <- SeuratObject::CreateDimReducObject(m, key = "Trex_",
                                                       assay = "RNA")
  out <- suppressWarnings(multimodal_neighbours(
    obj, rna_reduction = "pca", tcr_reduction = "Trex", rna_dims = 10,
    tcr_dims = 10, k_nn = 10, knn_range = 50))
  expect_true(all(c("wnn.umap", "rna.umap") %in% names(out@reductions)))
  expect_true("wsnn" %in% names(out@graphs))
  expect_true("RNA.weight" %in% colnames(out[[]]))
  expect_equal(length(out$wsnn_clusters), ncol(obj))
})

test_that("section_tcr_embedding runs with Trex", {
  skip_if_not_installed("Trex")
  skip_if_not_installed("Seurat")
  skip_on_cran()
  run <- spc_run()
  on.exit(run$report$close(), add = TRUE)
  env <- new.env()
  utils::data("trex_example", package = "Trex", envir = env)
  obj <- suppressWarnings(section_tcr_embedding(env$trex_example, NULL,
                                                run$run_dir, run$cfg,
                                                run$report))
  expect_true(all(c("Trex", "Trex.umap", "wnn.umap", "rna.umap") %in%
                    names(obj@reductions)))
  expect_false(anyNA(SeuratObject::Embeddings(obj, "Trex")))
  expect_true(file.exists(file.path(run$run_dir, "seurat_trex_wnn.rds")))
})
