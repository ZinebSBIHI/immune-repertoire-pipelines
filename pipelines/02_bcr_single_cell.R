#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Pipeline 02 - Single-cell BCR repertoire
#
# Input   AIRR Rearrangement TSV per library, as produced upstream by
#         Cell Ranger V(D)J -> IgBLAST (AssignGenes.py) -> MakeDb.py.
#         Those three steps run outside R and are not reproduced here.
# Output  A timestamped run directory with the resolved cell table, clonal
#         assignments, germlines, SHM, isotypes, convergence and a numeric
#         report.
#
# Run:
#   cp config/config.example.yml config/config.yml   # then edit
#   Rscript pipelines/02_bcr_single_cell.R
#
# Steps
#   1  Read libraries, attach metadata, trace the filtering
#   2  Resolve to one heavy chain per cell; drop doublets
#   3  Clonal assignment (threshold estimated or fixed), within subject
#   4  Germline reconstruction (IMGT, D segment masked)
#   5  Light-chain subclones
#   6  Somatic hypermutation by IMGT region
#   7  Isotypes, with the unswitched compartment reported combined
#   8  Clonal expansion, with a threshold matched to the material
#   9  Dominant-clone profiling (reported, not applied)
#  10  Convergence between subjects, against a permutation null
#  11  Numeric report
#
# WHAT THIS PIPELINE WILL NOT DO FOR YOU
# Clonal assignment groups sequences within a subject; a clone shared
# between subjects is a contradiction in terms and would indicate a sample
# mix-up or a barcode collision. Between-subject similarity is a separate
# question, handled in step 10 as convergence.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

ROOT <- local({
  candidates <- getwd()
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) candidates <- c(dirname(dirname(normalizePath(f[1]))), candidates)
  for (start in candidates) {
    d <- normalizePath(start, mustWork = FALSE)
    for (i in 1:6) {
      if (dir.exists(file.path(d, "R")) &&
          file.exists(file.path(d, "config", "config.example.yml"))) return(d)
      if (identical(dirname(d), d)) break
      d <- dirname(d)
    }
  }
  stop("Could not locate the repository root. Run from the repository ",
       "directory:\n  Rscript pipelines/02_bcr_single_cell.R", call. = FALSE)
})
for (.f in list.files(file.path(ROOT, "R"), pattern = "[.]R$",
                      full.names = TRUE)) source(.f)

# dowser is required, not optional: readIMGT(), createGermlines() and
# resolveLightChains() all live there, not in alakazam/shazam/scoper.
need(c("alakazam", "shazam", "scoper", "dowser"),
     "Pipeline 02 (single-cell BCR). The Immcantation stack is on CRAN")

cfg     <- load_config(file.path(ROOT, "config/config.yml"))
run_dir <- init_run(cfg, "02_bcr_single_cell")
use_pub_theme()

report <- open_report(file.path(run_dir, "REPORT.txt"),
                      sprintf("Single-cell BCR repertoire - %s",
                              cfg$project$name))

# ===========================================================================
# 1  Input and metadata
# ===========================================================================
report$h1("1. Input")

files <- list.files(cfg$paths$data_dir, pattern = "[.](airr[.])?tsv$",
                    full.names = TRUE, recursive = TRUE)
if (!length(files)) {
  stop("No AIRR .tsv files found under ", cfg$paths$data_dir,
       "\nExpected the output of MakeDb.py, one file per library.",
       call. = FALSE)
}
meta <- parse_filename_metadata(files, cfg$cohort$filename_pattern,
                                required = "subject")

db <- bind_rows(lapply(seq_len(nrow(meta)), function(i) {
  x <- utils::read.delim(meta$file[i], stringsAsFactors = FALSE,
                         check.names = FALSE)
  x$sample_id  <- sub("[.](airr[.])?tsv$", "", meta$basename[i])
  x$subject_id <- meta$subject[i]
  x$condition  <- meta$condition[i] %||% NA_character_
  # A cell barcode is unique only within its library. Prefixing with the
  # library is what keeps two cells from different libraries from being
  # merged into one.
  x$cell_id <- paste(x$sample_id, x$cell_id %||% x$barcode, sep = "_")
  x
}))

db$subject_id <- as_configured_factor(db$subject_id, cfg$cohort$subjects,
                                      "subject")
umi_col <- cfg$repertoire$abundance_column
if (!umi_col %in% names(db)) {
  stop("repertoire$abundance_column is '", umi_col, "' but no such column ",
       "is present. Available: ", paste(names(db), collapse = ", "),
       "\nName it explicitly in the configuration: UMI counts and read ",
       "counts are not interchangeable.", call. = FALSE)
}

report$p(sprintf("Libraries:      %d", nrow(meta)))
report$p(sprintf("Rearrangements: %s", format(nrow(db), big.mark = " ")))
report$p(sprintf("UMI column:     %s", umi_col))

# ===========================================================================
# 2  One heavy chain per cell
# ===========================================================================
report$h1("2. Cell resolution")

res <- resolve_one_heavy_per_cell(
  db, umi_col = umi_col,
  prefer_switched = cfg$bcr$cell_resolution$prefer_switched_isotype,
  aggregate_umi = cfg$bcr$cell_resolution$aggregate_umi,
  doublet_min_umi = cfg$bcr$cell_resolution$doublet_min_umi)

cells <- res$db
save_table(res$filter_log, "02_filter_log", run_dir,
           "rows and cells retained at each filtering step")
report$tbl(res$filter_log)
report$p("")
report$p(sprintf(
  paste0("A barcode carrying two or more heavy chains each at or above %g ",
         "UMI is treated as a doublet and removed (%d cells), rather than ",
         "resolved to whichever chain had more support. Where one chain ",
         "dominates, the cell is kept with that chain; class-switched ",
         "isotypes are preferred over IGHM/IGHD when choosing."),
  cfg$bcr$cell_resolution$doublet_min_umi, length(res$doublets)))

# ===========================================================================
# 3  Clonal assignment
# ===========================================================================
report$h1("3. Clonal assignment")

heavy <- cells[cells$locus == "IGH", , drop = FALSE]

if (identical(cfg$bcr$clone_threshold_mode, "estimate")) {
  dist_nearest <- shazam::distToNearest(
    heavy, sequenceColumn = "junction", vCallColumn = "v_call",
    jCallColumn = "j_call", model = "ham", normalize = "len",
    fields = cfg$bcr$clone_fields)

  thr_obj <- shazam::findThreshold(dist_nearest$dist_nearest, method = "gmm",
                                   model = "gamma-norm", cutoff = "optimal")
  threshold <- thr_obj@threshold
  report$p(sprintf("Threshold estimated from the data: %.4f", threshold))

  p_thr <- ggplot(dist_nearest[!is.na(dist_nearest$dist_nearest), ],
                  aes(dist_nearest)) +
    geom_histogram(bins = 60, fill = "grey75", colour = "white",
                   linewidth = 0.2) +
    geom_vline(xintercept = threshold, colour = "#D55E00", linewidth = 0.7) +
    facet_wrap(~subject_id, scales = "free_y") +
    labs(x = "Normalised Hamming distance to nearest neighbour",
         y = "Sequences", title = "Clonal assignment threshold",
         subtitle = sprintf("estimated threshold %.4f", threshold)) +
    theme_pub()
  save_figure(p_thr, "03_clone_threshold", run_dir, cfg, width = 8,
              height = 3.4,
              legend = paste("Distribution of the normalised distance from",
                             "each heavy chain to its nearest neighbour,",
                             "per subject, with the estimated clonal",
                             "assignment threshold."))
} else {
  threshold <- cfg$bcr$clone_threshold_fixed
  report$p(sprintf("Threshold fixed by configuration: %.4f", threshold))
  report$p("A fixed threshold makes runs comparable across batches at the",
           "\ncost of not adapting to each dataset. Either choice is",
           "\ndefensible; both must be stated in Methods.")
}

clones <- scoper::hierarchicalClones(
  heavy, threshold = threshold, method = cfg$bcr$clone_method,
  linkage = cfg$bcr$clone_linkage, normalize = "len",
  fields = cfg$bcr$clone_fields, summarize_clones = FALSE)

clone_counts <- clones %>%
  group_by(subject_id) %>%
  summarise(n_cells = n(), n_clones = n_distinct(clone_id), .groups = "drop")
report$h2("Clones per subject")
report$tbl(as.data.frame(clone_counts))
save_table(clone_counts, "03_clones_per_subject", run_dir,
           "cells and clones per subject")

# ===========================================================================
# 4  Germline reconstruction
# ===========================================================================
report$h1("4. Germlines")

imgt_dir <- file.path(cfg$paths$reference_dir, "imgt", "human", "vdj")
if (!dir.exists(imgt_dir)) {
  stop("IMGT germline references not found at ", imgt_dir,
       "\nSee data/README.md: fetch them with Immcantation's ",
       "fetch_imgtdb.sh, then point paths$reference_dir at the result.",
       call. = FALSE)
}
references <- dowser::readIMGT(dir = imgt_dir)

# The D segment is masked. D assignment is unreliable at the junction, and
# an incorrect D call would be counted as mutation downstream, inflating SHM
# in exactly the region where it is hardest to detect.
clones <- dowser::createGermlines(clones, references, nproc = 1)
report$p("Germlines reconstructed with the D segment masked: D assignment is",
         "\nunreliable at the junction, and a wrong call would be scored as",
         "\nmutation and inflate SHM.")

# ===========================================================================
# 5  Light-chain subclones
# ===========================================================================
report$h1("5. Light chains")

light <- cells[cells$locus %in% c("IGK", "IGL"), , drop = FALSE]
if (nrow(light)) {
  subclones <- dowser::resolveLightChains(
    bind_rows(clones, light), cell = "cell_id", locus = "locus")
  report$p(sprintf("Light chains available for %d of %d cells (%.1f%%).",
                   n_distinct(light$cell_id), n_distinct(clones$cell_id),
                   100 * n_distinct(light$cell_id) /
                     max(1L, n_distinct(clones$cell_id))))
  report$p("")
  report$p("A heavy-chain clone whose cells carry no light chain is more",
           "\noften an assembly artefact than an expansion; the light-chain",
           "\npairing rate per clone is reported in step 9.")
} else {
  subclones <- clones
  report$p("No light chains in the input; subclone resolution skipped.")
}

# ===========================================================================
# 6  Somatic hypermutation
# ===========================================================================
report$h1("6. Somatic hypermutation")

# Mutation counting runs on the heavy chains only. createGermlines() was
# applied to the heavy chains, so every light-chain row has a missing
# germline; calling observedMutations() over the whole table would score
# mutations against NA. The light chains are re-attached afterwards.
shm_heavy <- subclones[subclones$locus == "IGH", , drop = FALSE]
shm_light <- subclones[subclones$locus != "IGH", , drop = FALSE]

shm <- shm_heavy
for (region in cfg$bcr$shm$regions) {
  shm <- shazam::observedMutations(
    shm, sequenceColumn = "sequence_alignment",
    germlineColumn = "germline_alignment_d_mask",
    regionDefinition = get(region, envir = asNamespace("shazam")),
    frequency = TRUE, combine = TRUE, nproc = 1)
  names(shm)[names(shm) == "mu_freq"] <- paste0("mu_freq_", tolower(region))
}
shm <- bind_rows(shm, shm_light)

shm_by_sample <- shm %>%
  filter(locus == "IGH") %>%
  group_by(subject_id, sample_id) %>%
  summarise(n_cells = n(),
            median_shm_v = stats::median(.data[[paste0("mu_freq_",
              tolower(cfg$bcr$shm$regions[1]))]], na.rm = TRUE),
            .groups = "drop")
report$tbl(as.data.frame(shm_by_sample))
save_table(shm_by_sample, "06_shm_by_sample", run_dir,
           "median V-region mutation frequency per sample")

# ===========================================================================
# 7  Isotypes
# ===========================================================================
report$h1("7. Isotypes")
report$p(
  "One isotype is retained per cell, so an IgM versus IgD split would",
  "\ndescribe which transcript dominated rather than two populations: naive",
  "\nand marginal-zone B cells co-transcribe IGHM and IGHD. The combined",
  "\nunswitched fraction is reported instead.")
report$p("")

iso <- summarise_isotypes(
  shm, group_cols = c("subject_id", "sample_id"),
  combine_unswitched = cfg$bcr$report_unswitched_combined)
save_table(iso, "07_isotype_composition", run_dir,
           "isotype composition per sample")

sw <- switched_fraction(shm, group_cols = c("subject_id", "sample_id"))
report$tbl(as.data.frame(sw))
save_table(sw, "07_switched_fraction", run_dir,
           "class-switched fraction per sample")

p_iso <- ggplot(iso, aes(sample_id, fraction, fill = isotype_class)) +
  geom_col(width = 0.75) +
  facet_wrap(~subject_id, scales = "free_x") +
  scale_fill_manual(values = palette_categorical(
    length(unique(iso$isotype_class)))) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(x = NULL, y = "Cells", fill = "Isotype",
       title = "Isotype composition") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_figure(p_iso, "07_isotype_composition", run_dir, cfg, width = 9,
            height = 4,
            legend = paste("Isotype composition per sample. IGHM and IGHD",
                           "are reported as a single unswitched class."))

# ===========================================================================
# 8  Clonal expansion
# ===========================================================================
report$h1("8. Clonal expansion")
report$p(sprintf(
  paste0("Expansion thresholds differ by material: %d cells for resting and ",
         "%d for stimulated. A multi-day stimulation culture produces ",
         "identical daughter cells, so the cut-off that marks genuine ",
         "expansion in resting material is met trivially after ",
         "stimulation."),
  cfg$bcr$shm$expansion_min_cells_resting,
  cfg$bcr$shm$expansion_min_cells_stimulated))
report$p("")

# `stimulated` must be declared by the cohort; it is not inferred from a
# label, because guessing it would silently apply the wrong threshold.
if (!"stimulated" %in% names(shm)) {
  shm$stimulated <- FALSE
  report$p("No 'stimulated' column found; the resting threshold is applied",
           "\nthroughout. Add the column upstream if the cohort contains",
           "\nstimulated material.")
}

clone_sizes <- shm %>%
  filter(locus == "IGH") %>%
  group_by(subject_id, sample_id, stimulated, clone_id) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  mutate(threshold = expansion_threshold(
           stimulated,
           cfg$bcr$shm$expansion_min_cells_resting,
           cfg$bcr$shm$expansion_min_cells_stimulated),
         expanded = n_cells >= threshold)

expansion <- clone_sizes %>%
  group_by(subject_id, sample_id) %>%
  summarise(n_clones = n(),
            n_expanded = sum(expanded),
            pct_cells_in_expanded = 100 * sum(n_cells[expanded]) / sum(n_cells),
            largest_clone_pct = 100 * max(n_cells) / sum(n_cells),
            .groups = "drop")
report$tbl(as.data.frame(expansion))
save_table(expansion, "08_clonal_expansion", run_dir,
           "expanded clones and the share of cells they hold")

# ===========================================================================
# 9  Dominant clones
# ===========================================================================
report$h1("9. Dominant clones")
report$p(
  "Clones above 1% of a subject's cells are profiled because any",
  "\nrepertoire-level statistic computed on that compartment follows them:",
  "\ngene usage, CDR3 length and mutation load all do. Two criteria",
  "\nseparate real expansion from artefact: UMI support at the level of the",
  "\nsubject's other cells, and a paired light chain. Intraclonal diversity",
  "\nis deliberately not a criterion, since a genuine expansion produces",
  "\nidentical daughter cells.",
  "\n",
  "\nNothing is excluded here. The verdict column records which criteria",
  "\neach clone meets; exclusion is a judgement to state, and any statistic",
  "\nit changes should be reported both ways.")
report$p("")

dominant <- profile_dominant_clones(shm, clone_col = "clone_id",
                                    subject_col = "subject_id",
                                    umi_col = umi_col, min_fraction = 0.01)
if (nrow(dominant)) {
  report$tbl(as.data.frame(dominant[, c("subject_id", "clone_id", "n_cells",
                                        "fraction_of_subject", "v_gene",
                                        "isotype", "umi_ratio",
                                        "pct_with_light_chain", "verdict")]))
  save_table(dominant, "09_dominant_clones", run_dir,
             "clones above 1% of a subject, with supporting criteria")
} else {
  report$p("No clone exceeds 1% of any subject's cells.")
}

# ===========================================================================
# 10  Convergence between subjects
# ===========================================================================
report$h1("10. Convergence between subjects")

conv_input <- shm %>%
  filter(locus == "IGH", !is.na(junction_aa), nchar(junction_aa) > 0) %>%
  mutate(v_gene = sub("\\*.*$", "", sub(",.*$", "", v_call)),
         mu_freq = .data[[paste0("mu_freq_", tolower(cfg$bcr$shm$regions[1]))]],
         c_gene = sub("\\*.*$", "", sub(",.*$", "", c_call))) %>%
  as.data.frame()

if (isTRUE(cfg$convergence$require_mutated_or_switched)) {
  conv_input <- filter_informative_sequences(
    conv_input, shm_col = "mu_freq", isotype_col = "c_gene",
    shm_threshold = cfg$convergence$shm_threshold)
  report$p(sprintf(
    paste0("Restricted to sequences mutated above %.0f%% or class-switched: ",
           "%s sequences. Unmutated unswitched sequences sit close to ",
           "germline and therefore close to each other, so including them ",
           "measures shared germline rather than shared selection."),
    100 * cfg$convergence$shm_threshold,
    format(nrow(conv_input), big.mark = " ")))
}

conv_input <- distinct(conv_input, subject_id, v_gene, junction_aa)

if (nrow(conv_input) > 1L && n_distinct(conv_input$subject_id) >= 2L) {
  pub <- detect_public_groups(
    as.data.frame(conv_input), cdr3_col = "junction_aa", v_col = "v_gene",
    subject_col = "subject_id",
    identity_threshold = cfg$convergence$cdr3_identity,
    min_subjects = cfg$convergence$min_subjects,
    normalise = cfg$convergence$distance_normalisation)

  perm <- permute_public_groups(
    as.data.frame(conv_input), cdr3_col = "junction_aa", v_col = "v_gene",
    subject_col = "subject_id",
    identity_threshold = cfg$convergence$cdr3_identity,
    min_subjects = cfg$convergence$min_subjects,
    normalise = cfg$convergence$distance_normalisation,
    n_permutations = cfg$convergence$n_permutations)

  report$h2("Public groups")
  report$p(sprintf("observed            %d", perm$observed))
  report$p(sprintf("expected by chance  %.1f", perm$null_mean))
  report$p(sprintf("fold enrichment     %.2f", perm$fold_enrichment))
  report$p(sprintf("empirical p         %.4f  (%d permutations)",
                   perm$p_empirical, perm$n_permutations))

  if (nrow(pub$summary)) {
    save_table(pub$summary, "10_public_groups", run_dir,
               "CDR3 groups shared between subjects")
    save_table(pub$groups, "10_public_group_members", run_dir,
               "sequences belonging to each public group")
  }
} else {
  report$p("Too few subjects or sequences for a convergence analysis.")
}

# ===========================================================================
# 11  Save and close
# ===========================================================================
saveRDS(shm, file.path(run_dir, "bcr_cells_with_shm.rds"))
report$h1("11. Notes")
report$p("The resolved cell table with clones, germlines and SHM is saved as",
         "\nbcr_cells_with_shm.rds. Parameters are in resolved_config.yml and",
         "\npackage versions in sessionInfo.txt.")
report$close()

finalise_run(cfg, run_dir)
