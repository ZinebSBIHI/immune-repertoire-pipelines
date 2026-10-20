#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Pipeline 03 - Single-cell transcriptome (Seurat v5)
#
# Input   Cell Ranger filtered feature-barcode matrices (.h5), one per
#         library. Cell Ranger count runs outside R.
# Output  A timestamped run directory with per-sample QC, the clustered and
#         scored object, cluster composition diagnostics, markers, and a
#         numeric report.
#
# Run:
#   cp config/config.example.yml config/config.yml   # then edit
#   Rscript pipelines/03_single_cell_gex.R
#
# Steps
#   1  Read matrices, build one object per library
#   2  Per-sample QC with per-sample thresholds
#   3  Normalise, merge, variable features, PCA
#   4  Dimensionality and resolution sweep, then the configured choice
#   5  Cluster composition by subject: the confounding check
#   6  Markers per cluster
#   7  Module scores
#   8  Numeric report
#
# TWO DECISIONS THIS PIPELINE REFUSES TO MAKE SILENTLY
#
# QC THRESHOLDS ARE PER SAMPLE. A resting library and a stimulated one do
# not share an upper nFeature bound: proliferating
# cells are larger and carry more transcripts, so a single global threshold
# preferentially discards the activated cells, which are usually the cells
# of interest. Thresholds come from gex$qc, with one block per sample that
# needs its own.
#
# BATCH CORRECTION IS OFF BY DEFAULT. Integration makes samples look alike,
# which is indistinguishable from removing the effect under study when the
# batch and the condition coincide. It is enabled explicitly in the
# configuration, and whether it was used is recorded in the report.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

ROOT <- local({ source("R/setup.R", local = TRUE); find_repo_root() })
source(file.path(ROOT, "R/setup.R")); source_modules(ROOT)

need(c("Seurat", "SeuratObject"), "Pipeline 03 (single-cell transcriptome)")
suppressPackageStartupMessages(library(Seurat))

cfg     <- load_config(file.path(ROOT, "config/config.yml"))
run_dir <- init_run(cfg, "03_single_cell_gex")
use_pub_theme()

report <- open_report(file.path(run_dir, "REPORT.txt"),
                      sprintf("Single-cell transcriptome - %s",
                              cfg$project$name))

# ===========================================================================
# 1  Read matrices
# ===========================================================================
report$h1("1. Input")

h5 <- list.files(cfg$paths$data_dir, pattern = "[.]h5$", recursive = TRUE,
                 full.names = TRUE)
if (!length(h5)) {
  stop("No .h5 matrices found under ", cfg$paths$data_dir,
       "\nExpected Cell Ranger filtered_feature_bc_matrix.h5, one per ",
       "library.", call. = FALSE)
}
sample_names <- sub("[.]h5$", "", basename(h5))
# Where Cell Ranger's default file name is kept, the library name is the
# containing directory instead.
generic <- sample_names %in% c("filtered_feature_bc_matrix",
                               "raw_feature_bc_matrix", "sample_filtered_feature_bc_matrix")
sample_names[generic] <- basename(dirname(h5[generic]))

if (anyDuplicated(sample_names)) {
  dup <- unique(sample_names[duplicated(sample_names)])
  stop("Two or more matrices resolve to the same library name: ",
       paste(dup, collapse = ", "),
       "\nLibrary names must be unique; they key the QC thresholds and the ",
       "cell-name prefixes. Rename the directories or the files.",
       call. = FALSE)
}
meta <- parse_filename_metadata(sample_names, cfg$cohort$filename_pattern,
                                required = "subject")
report$p(sprintf("Libraries: %d", length(h5)))

objects <- setNames(vector("list", length(h5)), sample_names)
for (i in seq_along(h5)) {
  counts <- Seurat::Read10X_h5(h5[i])
  if (is.list(counts)) counts <- counts[["Gene Expression"]]
  o <- Seurat::CreateSeuratObject(counts, project = sample_names[i],
                                  min.cells = 3, min.features = 200)
  o$sample_id  <- sample_names[i]
  o$subject_id <- meta$subject[i]
  o$condition  <- meta$condition[i] %||% NA_character_
  objects[[i]] <- o
}

# ===========================================================================
# 2  Per-sample QC
# ===========================================================================
report$h1("2. Quality control")
report$p(
  "Thresholds are applied per sample. A single global upper bound on",
  "\nnFeature would preferentially discard proliferating cells, which carry",
  "\nmore transcripts; where those are the cells under study, a global",
  "\nthreshold removes the signal and looks like clean QC while doing it.")
report$p("")

qc_rows <- list()
for (nm in names(objects)) {
  o <- objects[[nm]]
  o[["percent_mt"]] <- Seurat::PercentageFeatureSet(
    o, pattern = cfg$gex$mito_pattern)

  thr <- cfg$gex$qc[[nm]] %||% cfg$gex$qc$default
  n_before <- ncol(o)

  keep <- o$nFeature_RNA > thr$min_features &
          o$nFeature_RNA < thr$max_features &
          o$nCount_RNA   < thr$max_counts &
          o$percent_mt   < thr$max_percent_mt
  o <- o[, keep]

  qc_rows[[nm]] <- data.frame(
    sample_id = nm,
    subject_id = meta$subject[match(nm, sample_names)],
    thresholds_from = if (is.null(cfg$gex$qc[[nm]])) "default" else "sample-specific",
    min_features = thr$min_features, max_features = thr$max_features,
    max_counts = thr$max_counts, max_percent_mt = thr$max_percent_mt,
    cells_before = n_before, cells_after = ncol(o),
    pct_retained = 100 * ncol(o) / n_before,
    stringsAsFactors = FALSE)

  objects[[nm]] <- o
}
qc <- bind_rows(qc_rows)
report$tbl(qc)
save_table(qc, "02_qc_per_sample", run_dir,
           "QC thresholds applied and cells retained per sample")

default_used <- qc$sample_id[qc$thresholds_from == "default"]
if (length(default_used)) {
  report$p("")
  report$p(sprintf(
    paste0("Default thresholds were used for %d of %d libraries (%s). ",
           "Inspect their distributions and add a gex$qc block for any ",
           "library whose metrics differ from the default assumption."),
    length(default_used), nrow(qc),
    paste(utils::head(default_used, 6), collapse = ", ")))
}

# ===========================================================================
# 3  Normalise, merge, PCA
# ===========================================================================
report$h1("3. Normalisation and dimensionality reduction")

# Normalisation is per sample, before the merge. Normalising after merging
# would let one library's depth set the scale factor for all of them.
objects <- lapply(objects, function(o) {
  o <- Seurat::NormalizeData(
    o, normalization.method = cfg$gex$normalisation$method,
    scale.factor = cfg$gex$normalisation$scale_factor, verbose = FALSE)
  Seurat::FindVariableFeatures(
    o, selection.method = "vst",
    nfeatures = cfg$gex$normalisation$n_variable_features, verbose = FALSE)
})

obj <- if (length(objects) == 1L) objects[[1]] else
  merge(objects[[1]], y = objects[-1], add.cell.ids = names(objects))
obj <- SeuratObject::JoinLayers(obj)

obj <- Seurat::FindVariableFeatures(
  obj, selection.method = "vst",
  nfeatures = cfg$gex$normalisation$n_variable_features, verbose = FALSE)
obj <- Seurat::ScaleData(obj, verbose = FALSE)
# Enough components for the widest sweep value, so the sweep never asks
# for dimensions the reduction does not have.
npcs_max <- max(50L, cfg$gex$clustering$n_pcs,
                cfg$gex$clustering$sweep_n_pcs)
obj <- Seurat::RunPCA(obj, npcs = npcs_max, verbose = FALSE)

reduction <- "pca"
if (identical(cfg$gex$integration$method, "harmony")) {
  need("harmony", "Batch correction with Harmony")
  obj <- harmony::RunHarmony(obj, group.by.vars = cfg$gex$integration$batch_variable,
                             verbose = FALSE)
  reduction <- "harmony"
  report$p(sprintf(
    paste0("Batch correction: Harmony on '%s'. Integration makes samples ",
           "resemble each other, so where batch and condition coincide it ",
           "is indistinguishable from removing the effect under study. ",
           "State in Methods why correction was applied and on which ",
           "variable."),
    cfg$gex$integration$batch_variable))
} else {
  report$p("Batch correction: none. The merge is uncorrected, so any",
           "\nstructure that follows subject or library is visible in the",
           "\nembedding rather than removed from it. Section 5 quantifies it.")
}

p_elbow <- Seurat::ElbowPlot(obj, ndims = npcs_max, reduction = reduction) +
  theme_pub()
save_figure(p_elbow, "03_elbow", run_dir, cfg, width = 5, height = 3.5,
            legend = "Standard deviation of each principal component.")

# ===========================================================================
# 4  Dimensionality and resolution sweep
# ===========================================================================
report$h1("4. Clustering")
report$p(
  "The sweep records how many clusters each combination yields, so that the",
  "\nconfigured choice can be seen as a choice. Cluster count is not a",
  "\nproperty of the data alone.")
report$p("")

sweep_rows <- list()
for (npc in cfg$gex$clustering$sweep_n_pcs) {
  tmp <- Seurat::FindNeighbors(obj, reduction = reduction, dims = 1:npc,
                               verbose = FALSE)
  for (res in cfg$gex$clustering$sweep_resolution) {
    tmp <- Seurat::FindClusters(tmp, resolution = res, verbose = FALSE)
    sweep_rows[[length(sweep_rows) + 1L]] <- data.frame(
      n_pcs = npc, resolution = res,
      n_clusters = length(unique(Seurat::Idents(tmp))),
      stringsAsFactors = FALSE)
  }
}
clust_sweep <- bind_rows(sweep_rows)
report$tbl(clust_sweep)
save_table(clust_sweep, "04_clustering_sweep", run_dir,
           "clusters obtained per dimensionality and resolution")

npc <- cfg$gex$clustering$n_pcs
res <- cfg$gex$clustering$resolution
obj <- Seurat::FindNeighbors(obj, reduction = reduction, dims = 1:npc,
                             verbose = FALSE)
obj <- Seurat::FindClusters(obj, resolution = res, verbose = FALSE)
obj <- Seurat::RunUMAP(obj, reduction = reduction, dims = 1:npc,
                       verbose = FALSE)
report$p("")
report$p(sprintf("Retained: %d components, resolution %.2f, %d clusters.",
                 npc, res, length(unique(Seurat::Idents(obj)))))

for (grp in c("seurat_clusters", "sample_id", "subject_id", "condition")) {
  if (!grp %in% names(obj@meta.data)) next
  if (all(is.na(obj@meta.data[[grp]]))) next
  n_lv <- length(unique(stats::na.omit(obj@meta.data[[grp]])))
  p <- Seurat::DimPlot(obj, group.by = grp, label = grp == "seurat_clusters") +
    scale_colour_manual(values = palette_categorical(n_lv),
                        na.value = "grey90") +
    labs(title = grp) +
    theme_pub(grid = "none")
  save_figure(p, paste0("04_umap_", grp), run_dir, cfg, width = 6, height = 5,
              legend = sprintf("Cells on the UMAP embedding, coloured by %s.",
                               grp))
}

# ===========================================================================
# 5  Cluster composition by subject
# ===========================================================================
report$h1("5. Cluster composition")
report$p(
  "A cluster made almost entirely of one subject's cells cannot support a",
  "\ncomparison between conditions: a difference in its size between arms",
  "\nis a difference between subjects. This check is run before any",
  "\ncluster-proportion comparison, and the flagged clusters are named.")
report$p("")

comp <- obj@meta.data %>%
  count(seurat_clusters, subject_id, name = "n_cells") %>%
  group_by(seurat_clusters) %>%
  mutate(fraction = n_cells / sum(n_cells),
         cluster_cells = sum(n_cells)) %>%
  ungroup()
save_table(comp, "05_cluster_composition", run_dir,
           "cells per cluster per subject")

dominated <- comp %>%
  group_by(seurat_clusters) %>%
  summarise(cluster_cells = first(cluster_cells),
            max_fraction = max(fraction),
            dominant_subject = subject_id[which.max(fraction)],
            .groups = "drop") %>%
  mutate(flagged = max_fraction > cfg$gex$max_single_subject_fraction)
report$tbl(as.data.frame(dominated))
save_table(dominated, "05_cluster_dominance", run_dir,
           "largest single-subject fraction per cluster")

flagged <- dominated$seurat_clusters[dominated$flagged]
if (length(flagged)) {
  report$p("")
  report$p(sprintf(
    paste0("Clusters %s exceed %.0f%% from a single subject. No proportion ",
           "based on these clusters can be compared between conditions. ",
           "Either exclude them from such comparisons and say so, or treat ",
           "the subject effect explicitly."),
    paste(flagged, collapse = ", "),
    100 * cfg$gex$max_single_subject_fraction))
} else {
  report$p("")
  report$p(sprintf("No cluster exceeds %.0f%% from a single subject.",
                   100 * cfg$gex$max_single_subject_fraction))
}

p_comp <- ggplot(comp, aes(seurat_clusters, fraction, fill = subject_id)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = cfg$gex$max_single_subject_fraction,
             linetype = "dashed", linewidth = 0.3) +
  scale_fill_manual(values = palette_categorical(
    length(unique(stats::na.omit(comp$subject_id))))) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(x = "Cluster", y = "Cells", fill = "Subject",
       title = "Subject composition of each cluster",
       subtitle = "Dashed line: dominance threshold from the configuration") +
  theme_pub()
save_figure(p_comp, "05_cluster_composition", run_dir, cfg, width = 8,
            height = 4,
            legend = paste("Fraction of each cluster contributed by each",
                           "subject. Clusters above the dashed line cannot",
                           "support between-condition comparison of their",
                           "size."))

# ===========================================================================
# 6  Markers
# ===========================================================================
report$h1("6. Cluster markers")

markers <- Seurat::FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25,
                                  logfc.threshold = 0.25, verbose = FALSE)
save_table(markers, "06_cluster_markers", run_dir,
           "positive markers per cluster (Wilcoxon)")

top <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10) %>%
  ungroup()
save_table(top, "06_cluster_markers_top10", run_dir,
           "ten strongest markers per cluster")
report$p(sprintf("%d markers across %d clusters; top ten per cluster saved.",
                 nrow(markers), n_distinct(markers$cluster)))

# ===========================================================================
# 7  Module scores
# ===========================================================================
report$h1("7. Module scores")

obj <- Seurat::CellCycleScoring(
  obj, s.features = Seurat::cc.genes.updated.2019$s.genes,
  g2m.features = Seurat::cc.genes.updated.2019$g2m.genes,
  set.ident = FALSE)

phase <- obj@meta.data %>%
  count(sample_id, Phase, name = "n_cells") %>%
  group_by(sample_id) %>% mutate(fraction = n_cells / sum(n_cells)) %>%
  ungroup()
report$tbl(as.data.frame(phase))
save_table(phase, "07_cell_cycle_phase", run_dir,
           "cell cycle phase composition per sample")

# Gene signatures are supplied by the project, as a two-column CSV
# (signature, gene). Hard-coding signatures in a pipeline hides which
# publication each came from.
sig_file <- file.path(cfg$paths$reference_dir, "signatures.csv")
if (file.exists(sig_file)) {
  sigs <- utils::read.csv(sig_file, stringsAsFactors = FALSE)
  present <- rownames(obj)
  sig_list <- split(sigs$gene, sigs$signature)

  kept <- lapply(sig_list, function(g) intersect(g, present))
  scorable <- vapply(kept, length, integer(1)) >= 3L

  # AddModuleScore names its output columns positionally (sig_1, sig_2, ...)
  # and ignores names(features), so the mapping from column to signature is
  # the order of the list it was GIVEN. The coverage table therefore records
  # that column name for the signatures actually scored, rather than being
  # built from the unfiltered list, where a dropped signature would shift
  # every row after it onto the wrong column.
  coverage <- data.frame(
    signature = names(sig_list),
    n_genes = vapply(sig_list, length, integer(1)),
    n_present = vapply(kept, length, integer(1)),
    scored = scorable,
    score_column = NA_character_,
    stringsAsFactors = FALSE)
  coverage$pct_present <- 100 * coverage$n_present / coverage$n_genes
  coverage$score_column[scorable] <- sprintf("sig_%d", seq_len(sum(scorable)))

  report$tbl(coverage)
  save_table(coverage, "07_signature_coverage", run_dir,
             "genes of each signature present, and its score column")
  report$p("")
  report$p("A signature scored on a small fraction of its genes is not that",
           "\nsignature. Coverage is reported so that a low-coverage score",
           "\nis not read as a negative result. Signatures with fewer than",
           "\nthree genes present are not scored at all; score_column maps",
           "\neach scored signature onto its column in the object.")

  if (any(scorable)) {
    obj <- Seurat::AddModuleScore(obj, features = kept[scorable],
                                  name = "sig_", seed = cfg$project$seed)
  }
} else {
  report$p(sprintf("No signature file at %s; module scoring limited to the",
                   sig_file))
  report$p("cell cycle. See data/README.md for the expected format.")
}

# ===========================================================================
# 8  Save and close
# ===========================================================================
saveRDS(obj, file.path(run_dir, "seurat_annotated.rds"))
report$h1("8. Notes")
report$p("The clustered and scored object is saved as",
         "seurat_annotated.rds;",
         "\nit carries the clusters and scores the repertoire can be joined",
         "\nto by cell barcode. Parameters are in",
         "\nresolved_config.yml and package versions in sessionInfo.txt.")
report$close()

finalise_run(cfg, run_dir)
