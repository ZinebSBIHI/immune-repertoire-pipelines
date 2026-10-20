#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Pipeline 01 - Bulk repertoire  (TCR or BCR)
#
# One axis: this is the pipeline for BULK repertoire sequencing, where the
# unit of observation is the molecule (UMI). The receptor is chosen in the
# configuration: set repertoire$locus to TRB for bulk TCR or IGH for bulk
# BCR. The analysis sections are identical either way; only the locus label
# and the input differ.
#
# Input   AIRR Rearrangement TSV per sample, from MiXCR or IgBLAST + MakeDb.
# Output  A timestamped run directory with figures, tables, a numeric
#         report, the resolved configuration and sessionInfo().
#
# Run:
#   Rscript examples/simulate_airr.R          # once, to create example data
#   cp config/config.example.yml config/config.yml
#   Rscript pipelines/01_bulk_repertoire.R
#
# Sections (shared with pipeline 02, so the two offer the same analyses):
#   1 composition QC          5 CDR3 length and gene usage
#   2 clonotype definitions   6 overlap with each subject's baseline
#   3 repertoire structure    7 convergence between subjects
#   4 diversity (rarefied)    8 paired comparisons across conditions
#                             + differential clonal abundance (bulk-specific)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(tidyr)
})

ROOT <- local({ source("R/setup.R", local = TRUE); find_repo_root() })
source(file.path(ROOT, "R/setup.R")); source_modules(ROOT)

cfg       <- load_config(file.path(ROOT, "config/config.yml"))
run_dir   <- init_run(cfg, "01_bulk_repertoire")
use_pub_theme()

locus     <- cfg$repertoire$locus
abundance <- cfg$repertoire$abundance_column
unit      <- "molecules"

report <- open_report(file.path(run_dir, "REPORT.txt"),
                      sprintf("Bulk %s repertoire - %s",
                              paste(locus, collapse = "/"), cfg$project$name))

# --- read inputs -----------------------------------------------------------
files <- list.files(cfg$paths$data_dir, pattern = "[.](airr[.])?tsv$",
                    full.names = TRUE, recursive = TRUE)
if (!length(files)) {
  stop("No AIRR .tsv files under ", cfg$paths$data_dir,
       "\nRun  Rscript examples/simulate_airr.R  to create an example dataset.",
       call. = FALSE)
}
meta <- parse_filename_metadata(files, cfg$cohort$filename_pattern,
                                required = c("subject", "timepoint"))
db <- bind_rows(lapply(seq_len(nrow(meta)), function(i) {
  x <- utils::read.delim(meta$file[i], stringsAsFactors = FALSE,
                         check.names = FALSE)
  x$sample_id  <- sub("[.](airr[.])?tsv$", "", meta$basename[i])
  x$subject_id <- meta$subject[i]
  x$timepoint  <- meta$timepoint[i]
  x$condition  <- meta$condition[i] %||% NA_character_
  x
}))
db$subject_id <- as_configured_factor(db$subject_id, cfg$cohort$subjects, "subject")
db$timepoint  <- as_configured_factor(db$timepoint, cfg$cohort$timepoints, "timepoint")
if (!all(is.na(db$condition))) {
  db$condition <- as_configured_factor(db$condition, cfg$cohort$conditions, "condition")
}
if (!abundance %in% names(db)) {
  stop("repertoire$abundance_column is '", abundance, "' but no such column ",
       "is present. Available: ", paste(names(db), collapse = ", "),
       call. = FALSE)
}
sample_meta <- unique(db[, c("sample_id", "subject_id", "timepoint", "condition")])

# --- shared sections -------------------------------------------------------
section_composition(db, locus, run_dir, cfg, report, prefix = "01")
section_clonotype_definitions(db, abundance, locus, run_dir, cfg, report,
                              prefix = "02")

clono <- collapse_clonotypes(db, abundance_col = abundance,
                             definition = cfg$repertoire$clonotype_definition,
                             min_abundance = cfg$repertoire$min_abundance,
                             locus = locus)
clono <- left_join(select(clono, -subject_id), sample_meta, by = "sample_id")

rep_summary <- section_structure(clono, sample_meta, run_dir, cfg, report,
                                 unit = unit, prefix = "03")
section_diversity(clono, sample_meta, run_dir, cfg, report,
                  within = "subject_id", unit = unit, prefix = "04")
section_gene_usage(db, abundance, locus, sample_meta, run_dir, cfg, report,
                   unit = unit, prefix = "05")
section_overlap(clono, run_dir, cfg, report, within = "subject_id",
                prefix = "06")

# --- differential clonal abundance (bulk-specific) -------------------------
# Needs a within-subject baseline timepoint, which bulk longitudinal designs
# have and a cross-sectional single-cell panel usually does not; that is why
# this section lives in the bulk pipeline rather than in the shared set.
report$h1("Differential clonal abundance")
report$p(paste0(
  "Two axes, reported separately. Category is descriptive, relative to the ",
  "subject's own baseline (de novo, pre-existing expanded, stable, ",
  "contracting, lost). Significance is a Fisher exact test on the molecule ",
  "counts with FDR control. A clonotype can be de novo and non-significant, ",
  "or pre-existing and significantly expanded; no single label asserting ",
  "induction is produced."))
report$p("")

ref_tp <- cfg$differential_abundance$reference_timepoint
if (!is.null(ref_tp) && ref_tp %in% as.character(clono$timepoint)) {
  calls <- test_clonal_abundance(
    clono, reference_timepoint = ref_tp,
    fold_change_threshold = cfg$differential_abundance$fold_change_threshold,
    detection_floor = cfg$differential_abundance$detection_floor,
    fdr_threshold = cfg$differential_abundance$fdr_threshold,
    fdr_method = cfg$differential_abundance$fdr_method)
  save_table(calls[calls$sig_expansion | calls$category == "de_novo", ],
             "07_differential_abundance_calls", run_dir,
             "clonotypes significantly expanded or called de novo")
  call_summary <- summarise_abundance_calls(calls)
  report$h2("Calls per sample")
  report$tbl(as.data.frame(call_summary))
  save_table(call_summary, "07_differential_abundance_summary", run_dir,
             "counts by category and significance, per sample")

  volc <- calls[is.finite(calls$log2_fold_change) & !is.na(calls$fdr), ]
  volc$neglog_fdr <- pmin(-log10(volc$fdr), 50)
  volc <- thin_points(volc, keep = "sig_expansion", max_background = 15000)
  p_v <- ggplot(volc, aes(log2_fold_change, neglog_fdr)) +
    geom_point(aes(colour = sig_expansion), size = 0.8, alpha = 0.6) +
    geom_hline(yintercept = -log10(cfg$differential_abundance$fdr_threshold),
               linetype = "dashed", linewidth = 0.3) +
    geom_vline(xintercept = log2(cfg$differential_abundance$fold_change_threshold),
               linetype = "dashed", linewidth = 0.3) +
    facet_grid(subject_id ~ timepoint) +
    scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#D55E00"),
                        labels = c("not significant", "significant")) +
    labs(x = expression(log[2]~"fold change vs baseline"),
         y = expression(-log[10]~"FDR (capped at 50)"), colour = NULL,
         title = "Differential clonal abundance") +
    theme_pub(grid = "both")
  save_figure(p_v, "07_volcano", run_dir, cfg, width = 9, height = 6,
              legend = paste("Fold change against each subject's own",
                             "baseline versus FDR; dashed lines mark the",
                             "configured thresholds."))
} else {
  report$p(sprintf(paste0("No baseline timepoint ('%s') in the data, so ",
                          "differential abundance is not computed."),
                   ref_tp %||% "(unset)"))
}

section_convergence(db, locus, run_dir, cfg, report,
                    shm_col = NULL, isotype_col = NULL, prefix = "08")
section_paired_comparisons(rep_summary, run_dir, cfg, report,
                           group_col = "condition", prefix = "09")

report$h1("Notes")
report$p(paste0(
  "Parameters are in resolved_config.yml, package versions in ",
  "sessionInfo.txt, figure legends in FIGURE_LEGENDS.txt and the table ",
  "inventory in TABLE_MANIFEST.txt."))
report$close()
finalise_run(cfg, run_dir)
