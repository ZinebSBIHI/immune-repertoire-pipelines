#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Pipeline 01 - Single-cell TCR repertoire
#
# Input   filtered_contig_annotations.csv per library (Cell Ranger VDJ),
#         and optionally a Seurat object to attach the repertoire to.
#         Cell Ranger runs outside R and is not reproduced here.
# Output  A timestamped run directory with clonotype tables, diversity,
#         gene usage, overlap, the repertoire-annotated Seurat object if one
#         was supplied, and a numeric report.
#
# Run:
#   cp config/config.example.yml config/config.yml   # then edit
#   Rscript pipelines/01_tcr_single_cell.R
#
# Steps
#   1  Read contigs, attach metadata, build clonotypes
#   2  Clonotype definition sensitivity: the four standard definitions
#   3  Clonal structure: abundance, size classes, homeostasis
#   4  Diversity, at native and at common depth
#   5  CDR3 length and V/J gene usage
#   6  Overlap between samples
#   7  Attach the repertoire to a Seurat object, if provided
#   8  Numeric report
#
# A NOTE ON CLONOTYPE DEFINITION
# "Clonotype" is used in the literature for at least four different keys,
# and the count of clonotypes in a sample can differ several-fold between
# them. Step 2 reports all four so that the choice is visible rather than
# implicit, and the rest of the pipeline uses the configured one.
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
       "directory:\n  Rscript pipelines/01_tcr_single_cell.R", call. = FALSE)
})
for (.f in list.files(file.path(ROOT, "R"), pattern = "[.]R$",
                      full.names = TRUE)) source(.f)

cfg     <- load_config(file.path(ROOT, "config/config.yml"))
run_dir <- init_run(cfg, "01_tcr_single_cell")
use_pub_theme()

report <- open_report(file.path(run_dir, "REPORT.txt"),
                      sprintf("Single-cell TCR repertoire - %s",
                              cfg$project$name))

# ===========================================================================
# 1  Contigs and clonotypes
# ===========================================================================
report$h1("1. Input and clonotypes")

files <- list.files(cfg$paths$data_dir, recursive = TRUE, full.names = TRUE,
                    pattern = "filtered_contig_annotations[.]csv$")
if (!length(files)) {
  stop("No filtered_contig_annotations.csv found under ", cfg$paths$data_dir,
       "\nExpected one per Cell Ranger VDJ library, in per-library ",
       "subdirectories.", call. = FALSE)
}
# The library name is taken from the containing directory, which is how
# Cell Ranger lays its output out.
sample_names <- basename(dirname(files))
report$p(sprintf("Libraries: %d", length(files)))

contigs <- bind_rows(lapply(seq_along(files), function(i) {
  x <- utils::read.csv(files[i], stringsAsFactors = FALSE)
  x$sample_id <- sample_names[i]
  x
}))

meta <- parse_filename_metadata(
  sample_names, cfg$cohort$filename_pattern, required = "subject")
contigs <- left_join(
  contigs,
  data.frame(sample_id = sample_names, subject_id = meta$subject,
             condition = meta$condition %||% NA_character_,
             stringsAsFactors = FALSE),
  by = "sample_id")

contigs$subject_id <- as_configured_factor(contigs$subject_id,
                                           cfg$cohort$subjects, "subject")

# Cell Ranger columns are mapped onto AIRR names once, here, so that the
# repertoire functions in R/ are agnostic to which assembler produced the
# input.
airr <- contigs %>%
  filter(.data$productive %in% c("TRUE", "True", TRUE),
         !is.na(.data$cdr3), nchar(.data$cdr3) > 0) %>%
  transmute(
    sample_id, subject_id, condition,
    cell_id     = paste(sample_id, barcode, sep = "_"),
    locus       = .data$chain,
    v_call      = .data$v_gene,
    j_call      = .data$j_gene,
    junction    = .data$cdr3_nt,
    junction_aa = .data$cdr3,
    umi_count   = .data$umis,
    productive  = TRUE
  )

report$p(sprintf("Productive contigs: %s", format(nrow(airr), big.mark = " ")))
report$p(sprintf("Cells with a TCR:   %s",
                 format(n_distinct(airr$cell_id), big.mark = " ")))

chain_comp <- airr %>%
  count(sample_id, locus, name = "n_contigs") %>%
  pivot_wider(names_from = locus, values_from = n_contigs, values_fill = 0)
report$h2("Chain composition per library")
report$tbl(as.data.frame(chain_comp))
save_table(chain_comp, "01_chain_composition", run_dir,
           "productive contigs per chain per library")

# ===========================================================================
# 2  Clonotype definition sensitivity
# ===========================================================================
report$h1("2. Clonotype definition")
report$p(
  "The same sample yields a different clonotype count under each of the",
  "\nfour definitions in common use. The configured definition is used",
  "\neverywhere below; this table is what makes that choice legible to a",
  "\nreader, and is worth including as a supplementary table.")
report$p("")

primary <- cfg$repertoire$locus
airr_primary <- airr[airr$locus == primary, , drop = FALSE]
if (!nrow(airr_primary)) {
  stop("No contigs for locus '", primary, "'. Set repertoire$locus to one ",
       "of: ", paste(sort(unique(airr$locus)), collapse = ", "),
       call. = FALSE)
}

defs <- c("aa", "nt", "gene", "strict")
def_counts <- bind_rows(lapply(defs, function(d) {
  cl <- collapse_clonotypes(airr_primary, abundance_col = "umi_count",
                            definition = d, min_abundance = 1,
                            locus = primary)
  cl %>% group_by(sample_id) %>%
    summarise(n_clonotypes = n(), .groups = "drop") %>%
    mutate(definition = d)
}))
def_table <- pivot_wider(def_counts, names_from = definition,
                         values_from = n_clonotypes)
report$tbl(as.data.frame(def_table))
save_table(def_table, "02_clonotype_definition_sensitivity", run_dir,
           "clonotype count under each of the four definitions")

# ===========================================================================
# 3  Clonal structure
# ===========================================================================
report$h1("3. Clonal structure")

clono <- collapse_clonotypes(
  airr_primary, abundance_col = "umi_count",
  definition = cfg$repertoire$clonotype_definition,
  min_abundance = cfg$repertoire$min_abundance, locus = primary)

sample_meta <- distinct(airr, sample_id, subject_id, condition)
clono <- left_join(select(clono, -subject_id), sample_meta, by = "sample_id")

# In single-cell data the natural unit is the cell, not the molecule, so
# clone size is reported in cells as well as in UMIs. A clone of one cell
# with many UMIs is not an expansion.
cells_per_clone <- airr_primary %>%
  mutate(clonotype = clonotype_key(., cfg$repertoire$clonotype_definition)) %>%
  distinct(sample_id, clonotype, cell_id) %>%
  count(sample_id, clonotype, name = "n_cells")
clono <- left_join(clono, cells_per_clone, by = c("sample_id", "clonotype"))

# NOTE on the naming here: the summary column must NOT be called `n_cells`.
# Inside one summarise(), a new column shadows the input column of the same
# name for every expression that follows it, so `n_cells = sum(n_cells)`
# would turn every later reference into that single total — silently
# yielding 0% singletons, a largest clone equal to the whole sample, and a
# Shannon entropy of 0, in a table that is also printed into the report.
rep_summary <- clono %>%
  group_by(sample_id, subject_id, condition) %>%
  summarise(total_cells = sum(n_cells, na.rm = TRUE),
            n_clonotypes = n(),
            pct_singleton_clones = 100 * sum(n_cells == 1, na.rm = TRUE) / n(),
            largest_clone_cells = max(n_cells, na.rm = TRUE),
            largest_clone_pct = 100 * max(n_cells, na.rm = TRUE) /
              sum(n_cells, na.rm = TRUE),
            shannon = shannon_entropy(n_cells),
            inv_simpson = inverse_simpson(n_cells),
            .groups = "drop")
report$tbl(as.data.frame(rep_summary))
save_table(rep_summary, "03_repertoire_summary", run_dir,
           "clonal structure per library, in cells")

homeo <- clono %>%
  mutate(size_class = cut(n_cells, breaks = c(0, 1, 2, 5, 10, Inf),
                          labels = c("1", "2", "3-5", "6-10", ">10"))) %>%
  group_by(sample_id, subject_id, condition, size_class) %>%
  summarise(n_clonotypes = n(), n_cells = sum(n_cells), .groups = "drop") %>%
  group_by(sample_id) %>%
  mutate(cell_fraction = n_cells / sum(n_cells)) %>%
  ungroup()
save_table(homeo, "03_clonal_size_classes", run_dir,
           "cells by clone size class")

p_homeo <- ggplot(homeo, aes(sample_id, cell_fraction, fill = size_class)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = palette_sequential(5, direction = -1)) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(x = NULL, y = "Cells", fill = "Cells per clone",
       title = "Clonal size distribution") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_figure(p_homeo, "03_clonal_size_classes", run_dir, cfg, width = 7,
            height = 4,
            legend = paste("Fraction of cells in clones of each size, per",
                           "library."))

# ===========================================================================
# 4  Diversity
# ===========================================================================
report$h1("4. Diversity")
report$p(
  "Cell number per library varies with loading and with sort yield, and",
  "\nrichness and clonality both scale with it. Values at native depth are",
  "\nreported for completeness; the comparison is made on rarefied values.")
report$p("")

if (isTRUE(cfg$repertoire$rarefaction$enabled)) {
  cell_clono <- clono %>%
    group_by(sample_id) %>%
    # Grouped, so frequency is within sample. Ungrouped it would be a share
    # of the whole cohort, which is not what any downstream metric expects.
    mutate(abundance = n_cells, frequency = n_cells / sum(n_cells)) %>%
    ungroup() %>%
    select(sample_id, subject_id, clonotype, abundance, frequency)

  # within = "subject_id" keeps subject_id on the returned frame, which the
  # join and the figure below both need. With within = NULL the function
  # carries only sample_id, and the join fails.
  rar <- rarefy_to_common_depth(
    cell_clono, within = "subject_id",
    n_bootstrap = cfg$repertoire$rarefaction$n_bootstrap)
  rar <- left_join(rar, select(sample_meta, sample_id, condition),
                   by = "sample_id")
  save_table(rar, "04_rarefied_metrics", run_dir,
             "metrics at common cell number, with 95% intervals")
  report$tbl(as.data.frame(unique(rar[, c("subject_id", "target_depth")])))
  report$tbl(as.data.frame(rar[, c("sample_id", "statistic", "mean",
                                   "q025", "q975")]))

  p_rar <- ggplot(rar, aes(sample_id, mean)) +
    geom_errorbar(aes(ymin = q025, ymax = q975), width = 0.15,
                  linewidth = 0.3, colour = "grey50") +
    geom_point(aes(colour = subject_id), size = 2.2) +
    facet_wrap(~statistic, scales = "free_y") +
    labs(x = NULL, y = "Value at common cell number", colour = "Subject",
         title = "Repertoire metrics at common cell number") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_figure(p_rar, "04_rarefied_metrics", run_dir, cfg, width = 9,
              height = 5,
              legend = paste("Repertoire metrics after resampling every",
                             "library to a common number of cells. Points",
                             "are the resampling mean, bars the 2.5th and",
                             "97.5th percentiles."))
}

# ===========================================================================
# 5  CDR3 length and gene usage
# ===========================================================================
report$h1("5. CDR3 length and gene usage")

len <- airr_primary %>%
  distinct(sample_id, subject_id, cell_id, junction_aa) %>%
  mutate(cdr3_length = nchar(junction_aa))
p_len <- ggplot(len, aes(cdr3_length)) +
  geom_histogram(binwidth = 1, fill = "grey70", colour = "white",
                 linewidth = 0.2) +
  facet_wrap(~sample_id, scales = "free_y") +
  labs(x = "CDR3 length (amino acids)", y = "Cells",
       title = sprintf("%s CDR3 length", primary)) +
  theme_pub()
save_figure(p_len, "05a_cdr3_length", run_dir, cfg, width = 9, height = 5,
            legend = "CDR3 length distribution per library, one cell per count.")

usage <- bind_rows(
  gene_usage(airr_primary, "v_call", level = "clonotype",
             definition = cfg$repertoire$clonotype_definition),
  gene_usage(airr_primary, "v_call", level = "molecule",
             abundance_col = "umi_count",
             definition = cfg$repertoire$clonotype_definition)
) %>% left_join(sample_meta, by = "sample_id")
save_table(usage, "05_v_gene_usage", run_dir,
           "V gene usage at clonotype and molecule level")

top_v <- usage %>% filter(level == "clonotype") %>%
  group_by(gene) %>% summarise(m = mean(frequency), .groups = "drop") %>%
  slice_max(m, n = 20) %>% pull(gene)

p_usage <- usage %>% filter(gene %in% top_v) %>%
  ggplot(aes(reorder(gene, frequency), frequency, colour = level)) +
  stat_summary(fun = mean, geom = "point", size = 1.8,
               position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2,
               linewidth = 0.3, position = position_dodge(width = 0.5)) +
  coord_flip() +
  scale_colour_manual(values = palette_categorical(2)) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(x = NULL, y = "Usage", colour = "Level",
       title = "V gene usage, 20 most used genes",
       subtitle = paste("Clonotype level separates repertoire-wide bias from",
                        "clonal expansion")) +
  theme_pub(grid = "x")
save_figure(p_usage, "05b_v_gene_usage", run_dir, cfg, width = 6.5,
            height = 6,
            legend = paste("Mean V gene usage across libraries at clonotype",
                           "and molecule level, with standard error. Where",
                           "the two disagree, an expanded clone is the",
                           "explanation."))

# ===========================================================================
# 6  Overlap
# ===========================================================================
report$h1("6. Overlap between libraries")

if (n_distinct(clono$sample_id) >= 2L) {
  ov <- bind_rows(
    repertoire_overlap(clono, method = "morisita"),
    repertoire_overlap(clono, method = "jaccard"))
  report$tbl(as.data.frame(ov))
  save_table(ov, "06_overlap", run_dir,
             "pairwise overlap, Morisita-Horn and Jaccard")

  p_ov <- ggplot(ov, aes(sample_a, sample_b, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", value)), size = 2.8,
              colour = "grey15") +
    facet_wrap(~method) +
    scale_fill_gradientn(colours = palette_sequential(64)) +
    labs(x = NULL, y = NULL, fill = "Overlap",
         title = "Repertoire overlap between libraries") +
    theme_pub(grid = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_figure(p_ov, "06_overlap", run_dir, cfg, width = 9, height = 4.5,
              legend = paste("Pairwise repertoire overlap. Morisita-Horn is",
                             "abundance-weighted; Jaccard counts clonotypes",
                             "and is the more depth-dependent of the two."))
}

# ===========================================================================
# 7  Attach to a Seurat object, if one was supplied
# ===========================================================================
report$h1("7. Transcriptome integration")

seurat_path <- file.path(cfg$paths$data_dir, "seurat_annotated.rds")
if (file.exists(seurat_path) &&
    requireNamespace("Seurat", quietly = TRUE)) {

  obj <- readRDS(seurat_path)

  # Barcode conventions differ between Cell Ranger and a merged Seurat
  # object, and a silent mismatch attaches the repertoire to nothing while
  # appearing to succeed. The overlap is therefore checked and reported
  # before anything is written.
  clone_by_cell <- airr_primary %>%
    mutate(clonotype = clonotype_key(., cfg$repertoire$clonotype_definition)) %>%
    group_by(cell_id) %>%
    slice_max(umi_count, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    left_join(cells_per_clone, by = c("sample_id", "clonotype")) %>%
    transmute(cell_id, clonotype, clone_cells = n_cells)

  matched <- sum(colnames(obj) %in% clone_by_cell$cell_id)
  report$p(sprintf("Seurat cells: %d; matched to a TCR: %d (%.1f%%)",
                   ncol(obj), matched, 100 * matched / ncol(obj)))

  if (matched == 0L) {
    report$p("")
    report$p("No barcode matched. The repertoire was NOT attached. Cell",
             "\nnames in the Seurat object and the contig barcodes follow",
             "\ndifferent conventions; reconcile them before integrating.")
  } else {
    idx <- match(colnames(obj), clone_by_cell$cell_id)
    obj$clonotype   <- clone_by_cell$clonotype[idx]
    obj$clone_cells <- clone_by_cell$clone_cells[idx]
    obj$clone_size  <- cut(obj$clone_cells, breaks = c(0, 1, 2, 5, 10, Inf),
                           labels = c("1", "2", "3-5", "6-10", ">10"))
    saveRDS(obj, file.path(run_dir, "seurat_tcr_integrated.rds"))
    report$p("Repertoire attached; object saved as",
             "seurat_tcr_integrated.rds.")

    if ("umap" %in% names(obj@reductions)) {
      p_umap <- Seurat::DimPlot(obj, group.by = "clone_size",
                                order = TRUE) +
        scale_colour_manual(values = palette_sequential(5, direction = -1),
                            na.value = "grey90") +
        labs(title = "Clone size on the transcriptome embedding") +
        theme_pub(grid = "none")
      save_figure(p_umap, "07_umap_clone_size", run_dir, cfg,
                  width = 6, height = 5,
                  legend = paste("Cells on the transcriptome UMAP coloured",
                                 "by the number of cells sharing their",
                                 "clonotype. Grey: no TCR recovered."))
    }
  }
} else {
  report$p(sprintf("No Seurat object at %s; integration skipped.",
                   seurat_path))
  report$p("Pipeline 04 produces a suitable object. Place it, or a symlink",
           "\nto it, at that path to attach the repertoire to it here.")
}

# ===========================================================================
# 8  Close
# ===========================================================================
saveRDS(clono, file.path(run_dir, "clonotypes.rds"))
report$h1("8. Notes")
report$p("Clonotype table saved as clonotypes.rds. Parameters are in",
         "\nresolved_config.yml and package versions in sessionInfo.txt.")
report$close()

finalise_run(cfg, run_dir)
