#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Pipeline 02 - Single-cell repertoire  (TCR or BCR)
#
# One axis: this is the pipeline for SINGLE-CELL immune receptor sequencing,
# where the unit of observation is the cell. The receptor is chosen in the
# configuration through repertoire$locus:
#
#   TRB (or TRA) -> single-cell TCR
#   IGH          -> single-cell BCR, which additionally runs the
#                   immunoglobulin-only sections: one heavy chain per cell,
#                   clonal inference, germlines, somatic hypermutation,
#                   isotypes and dominant-clone profiling.
#
# The shared repertoire sections (structure, diversity, gene usage, overlap,
# convergence) are the same ones pipeline 01 runs on bulk data, so the two
# pipelines cannot drift apart in what they offer. The immunoglobulin
# sections switch on automatically for an Ig locus and are skipped for a TCR
# locus, which does not hypermutate.
#
# Input   AIRR Rearrangement TSV per library (Cell Ranger V(D)J -> IgBLAST
#         -> MakeDb), or Cell Ranger filtered_contig_annotations.csv.
# Output  A timestamped run directory with figures, tables, a numeric
#         report, the resolved configuration and sessionInfo().
#
# Run:
#   cp config/config.example.yml config/config.yml   # set repertoire$locus
#   Rscript pipelines/02_single_cell_repertoire.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(tidyr)
})

ROOT <- local({ source("R/setup.R", local = TRUE); find_repo_root() })
source(file.path(ROOT, "R/setup.R")); source_modules(ROOT)

cfg       <- load_config(file.path(ROOT, "config/config.yml"))
locus     <- cfg$repertoire$locus
is_ig     <- grepl("^IG", locus[1])
unit      <- "cells"

# Immunoglobulin work needs the Immcantation stack; TCR work does not, so
# the dependency is only required when the locus calls for it.
if (is_ig) {
  need(c("alakazam", "shazam", "scoper", "dowser"),
       "Single-cell BCR (Ig locus). The Immcantation stack is on CRAN")
}

run_dir <- init_run(cfg, "02_single_cell_repertoire")
use_pub_theme()

report <- open_report(file.path(run_dir, "REPORT.txt"),
                      sprintf("Single-cell %s repertoire - %s",
                              paste(locus, collapse = "/"), cfg$project$name))

# --- read inputs, normalise to AIRR ----------------------------------------
# Two input conventions are accepted and mapped onto AIRR column names here,
# so every section downstream is agnostic to which one produced the data.
csv  <- list.files(cfg$paths$data_dir, recursive = TRUE, full.names = TRUE,
                   pattern = "filtered_contig_annotations[.]csv$")
tsv  <- list.files(cfg$paths$data_dir, recursive = TRUE, full.names = TRUE,
                   pattern = "[.](airr[.])?tsv$")

if (length(csv)) {
  sample_names <- basename(dirname(csv))
  meta <- parse_filename_metadata(sample_names, cfg$cohort$filename_pattern,
                                  required = "subject")
  db <- bind_rows(lapply(seq_along(csv), function(i) {
    x <- utils::read.csv(csv[i], stringsAsFactors = FALSE)
    data.frame(
      sample_id   = sample_names[i],
      subject_id  = meta$subject[i],
      condition   = meta$condition[i] %||% NA_character_,
      cell_id     = paste(sample_names[i], x$barcode, sep = "_"),
      locus       = x$chain,
      v_call      = x$v_gene, j_call = x$j_gene, c_call = x$c_gene %||% NA,
      junction    = x$cdr3_nt, junction_aa = x$cdr3,
      umi_count   = x$umis, productive = x$productive,
      stringsAsFactors = FALSE)
  }))
  abundance <- "umi_count"
} else if (length(tsv)) {
  meta <- parse_filename_metadata(tsv, cfg$cohort$filename_pattern,
                                  required = "subject")
  db <- bind_rows(lapply(seq_len(nrow(meta)), function(i) {
    x <- utils::read.delim(meta$file[i], stringsAsFactors = FALSE,
                           check.names = FALSE)
    # The per-cell barcode must be present in the source. It is checked
    # here, before a synthetic cell_id is built, because a bulk file pointed
    # at this pipeline by mistake has no barcode and would otherwise have
    # every molecule counted as a cell -- a silently wrong report (100%
    # singletons, one cell per clonotype). Bulk data belongs in pipeline 01.
    barcode_col <- intersect(c("cell_id", "barcode"), names(x))
    if (!length(barcode_col)) {
      stop("File has no per-cell barcode column (cell_id or barcode): ",
           meta$basename[i],
           "\nThis is the single-cell pipeline; bulk repertoire data goes ",
           "to pipelines/01_bulk_repertoire.R.", call. = FALSE)
    }
    x$sample_id  <- sub("[.](airr[.])?tsv$", "", meta$basename[i])
    x$subject_id <- meta$subject[i]
    x$condition  <- meta$condition[i] %||% NA_character_
    x$cell_id <- paste(x$sample_id, x[[barcode_col[1]]], sep = "_")
    x
  }))
  abundance <- cfg$repertoire$abundance_column
  if (!abundance %in% names(db)) {
    stop("repertoire$abundance_column is '", abundance, "' but no such ",
         "column is present. Available: ", paste(names(db), collapse = ", "),
         call. = FALSE)
  }
} else {
  stop("No input found under ", cfg$paths$data_dir,
       "\nExpected filtered_contig_annotations.csv per library (TCR) or ",
       "AIRR .tsv per library (BCR).", call. = FALSE)
}

db$subject_id <- as_configured_factor(db$subject_id, cfg$cohort$subjects, "subject")
sample_meta <- unique(db[, c("sample_id", "subject_id", "condition")])

# ===========================================================================
# Immunoglobulin-only preprocessing: one heavy chain per cell, clones,
# germlines, SHM. Produces the mutation and isotype columns that let the
# shared convergence section apply its germline-proximity filter.
# ===========================================================================
shm_col <- NULL; isotype_col <- NULL

if (is_ig) {
  report$h1("B-cell preprocessing (Ig locus)")

  res <- resolve_one_heavy_per_cell(
    db, umi_col = abundance,
    prefer_switched = cfg$bcr$cell_resolution$prefer_switched_isotype,
    aggregate_umi = cfg$bcr$cell_resolution$aggregate_umi,
    doublet_min_umi = cfg$bcr$cell_resolution$doublet_min_umi)
  cells <- res$db
  save_table(res$filter_log, "00_filter_log", run_dir,
             "rows and cells retained at each filtering step")
  report$tbl(res$filter_log)

  heavy <- cells[cells$locus == "IGH", , drop = FALSE]

  if (identical(cfg$bcr$clone_threshold_mode, "estimate")) {
    dn <- shazam::distToNearest(heavy, sequenceColumn = "junction",
                                vCallColumn = "v_call", jCallColumn = "j_call",
                                model = "ham", normalize = "len",
                                fields = cfg$bcr$clone_fields)
    threshold <- shazam::findThreshold(dn$dist_nearest, method = "gmm",
                                       model = "gamma-norm",
                                       cutoff = "optimal")@threshold
    report$p(sprintf("Clonal threshold estimated from the data: %.4f", threshold))
  } else {
    threshold <- cfg$bcr$clone_threshold_fixed
    report$p(sprintf("Clonal threshold fixed by configuration: %.4f", threshold))
  }

  clones <- scoper::hierarchicalClones(
    heavy, threshold = threshold, method = cfg$bcr$clone_method,
    linkage = cfg$bcr$clone_linkage, normalize = "len",
    fields = cfg$bcr$clone_fields, summarize_clones = FALSE)

  imgt_dir <- file.path(cfg$paths$reference_dir, "imgt", "human", "vdj")
  if (!dir.exists(imgt_dir)) {
    stop("IMGT germline references not found at ", imgt_dir,
         "\nSee data/README.md: fetch with Immcantation's fetch_imgtdb.sh.",
         call. = FALSE)
  }
  references <- dowser::readIMGT(dir = imgt_dir)
  clones <- dowser::createGermlines(clones, references, nproc = 1)

  light <- cells[cells$locus %in% c("IGK", "IGL"), , drop = FALSE]
  subclones <- if (nrow(light)) {
    dowser::resolveLightChains(bind_rows(clones, light), cell = "cell_id",
                               locus = "locus")
  } else clones

  shm_heavy <- subclones[subclones$locus == "IGH", , drop = FALSE]
  shm_light <- subclones[subclones$locus != "IGH", , drop = FALSE]
  for (region in cfg$bcr$shm$regions) {
    shm_heavy <- shazam::observedMutations(
      shm_heavy, sequenceColumn = "sequence_alignment",
      germlineColumn = "germline_alignment_d_mask",
      regionDefinition = get(region, envir = asNamespace("shazam")),
      frequency = TRUE, combine = TRUE, nproc = 1)
    names(shm_heavy)[names(shm_heavy) == "mu_freq"] <-
      paste0("mu_freq_", tolower(region))
  }
  shm <- bind_rows(shm_heavy, shm_light)
  shm$mu_freq <- shm[[paste0("mu_freq_", tolower(cfg$bcr$shm$regions[1]))]]
  shm$c_gene  <- strip_allele(shm$c_call)

  # Isotypes, reported with the unswitched compartment combined.
  report$h1("Isotypes")
  report$p(paste0(
    "One isotype is retained per cell, so an IgM-versus-IgD split would ",
    "describe which transcript dominated, not two populations: naive and ",
    "marginal-zone B cells co-transcribe IGHM and IGHD. The combined ",
    "unswitched fraction is reported instead."))
  iso <- summarise_isotypes(shm, group_cols = c("subject_id", "sample_id"),
                            combine_unswitched = cfg$bcr$report_unswitched_combined)
  save_table(iso, "0X_isotype_composition", run_dir, "isotype composition per sample")
  sw <- switched_fraction(shm, group_cols = c("subject_id", "sample_id"))
  report$tbl(as.data.frame(sw))
  save_table(sw, "0X_switched_fraction", run_dir, "class-switched fraction per sample")

  # The shared repertoire sections run on the heavy chain, one row per cell,
  # carrying the mutation and isotype columns.
  db  <- shm[shm$locus == "IGH", , drop = FALSE]
  shm_col     <- "mu_freq"
  isotype_col <- "c_gene"
  saveRDS(shm, file.path(run_dir, "bcr_cells_with_shm.rds"))
}

# ===========================================================================
# Shared repertoire sections. In single-cell data abundance is the cell
# count per clonotype, so the collapse keys on cells, and every section is
# told its unit is cells.
# ===========================================================================
section_composition(db, locus, run_dir, cfg, report, prefix = "01")

# One row per cell per clonotype -> clonotype table in cells.
cell_rows <- filter_locus_productive(db, locus)
cell_rows$clonotype <- clonotype_key(cell_rows, cfg$repertoire$clonotype_definition)
cells_per_clone <- cell_rows %>%
  distinct(sample_id, clonotype, cell_id) %>%
  count(sample_id, clonotype, name = "abundance")
clono <- cells_per_clone %>%
  group_by(sample_id) %>%
  mutate(frequency = abundance / sum(abundance)) %>%
  ungroup() %>%
  left_join(sample_meta, by = "sample_id")

section_clonotype_definitions(db, abundance, locus, run_dir, cfg, report,
                              prefix = "02")
rep_summary <- section_structure(clono, sample_meta, run_dir, cfg, report,
                                 unit = unit, prefix = "03")
# Cross-sectional single-cell panels usually have one sample per subject, so
# the common depth is taken across the cohort rather than within subject.
section_diversity(clono, sample_meta, run_dir, cfg, report,
                  within = NULL, unit = unit, prefix = "04")
section_gene_usage(db, abundance, locus, sample_meta, run_dir, cfg, report,
                   unit = unit, prefix = "05")
section_overlap(clono, run_dir, cfg, report, within = NULL, prefix = "06")
section_convergence(db, locus, run_dir, cfg, report,
                    shm_col = shm_col, isotype_col = isotype_col, prefix = "07")
section_paired_comparisons(rep_summary, run_dir, cfg, report,
                           group_col = "condition", prefix = "08")

saveRDS(clono, file.path(run_dir, "clonotypes.rds"))
report$h1("Notes")
report$p(paste0(
  "Clonotype table saved as clonotypes.rds", if (is_ig)
    "; the resolved cell table with clones, germlines and SHM as bcr_cells_with_shm.rds"
  else "", ". Parameters are in resolved_config.yml and package versions in ",
  "sessionInfo.txt."))
report$close()
finalise_run(cfg, run_dir)
