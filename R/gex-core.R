# ---------------------------------------------------------------------------
# Single-cell transcriptome: quality control, reduction, clustering,
# composition and cell cycle (Seurat v5)
#
# TWO DECISIONS THESE SECTIONS REFUSE TO MAKE SILENTLY
#
# QC THRESHOLDS ARE PER SAMPLE. A resting library and a stimulated one do
# not share an upper nFeature bound: proliferating cells are larger and carry
# more transcripts, so a single global threshold preferentially discards the
# activated cells, which are usually the cells of interest. Thresholds come
# from gex$qc, with one block per sample that needs its own.
#
# BATCH CORRECTION IS OFF BY DEFAULT. Integration makes samples look alike,
# which is indistinguishable from removing the effect under study when the
# batch and the condition coincide. It is enabled explicitly in the
# configuration, and whether it was used is recorded in the report.
#
#   Seurat v5: Hao Y et al. Nat Biotechnol. 2024;42(2):293-304.
#     doi:10.1038/s41587-023-01767-y. PMID 37231261.
#   Harmony: Korsunsky I et al. Nat Methods. 2019;16(12):1289-1296.
#     doi:10.1038/s41592-019-0619-0. PMID 31740819.
# ---------------------------------------------------------------------------

#' Per-sample quality control of expression libraries
#'
#' @param objects Named list of Seurat objects, one per library.
#' @param meta Library metadata from `read_gex_inputs()`.
#' @inheritParams section_arguments
#' @return The filtered objects, invisibly; the QC table is attached as the
#'   attribute `qc`.
#' @export
section_gex_qc <- function(objects, meta, run_dir, cfg, report,
                           prefix = "02") {
  need(c("Seurat", "dplyr"), "section_gex_qc()")
  report$h1("Quality control")
  report$p(paste0(
    "Thresholds are applied per sample. A single global upper bound on ",
    "nFeature would preferentially discard proliferating cells, which carry ",
    "more transcripts; where those are the cells under study, a global ",
    "threshold removes the signal and looks like clean QC while doing it."))
  report$p("")

  rows <- list()
  for (nm in names(objects)) {
    o <- objects[[nm]]
    o[["percent_mt"]] <- Seurat::PercentageFeatureSet(
      o, pattern = cfg$gex$mito_pattern)
    thr <- cfg$gex$qc[[nm]] %||% cfg$gex$qc$default
    n_before <- ncol(o)
    keep <- o$nFeature_RNA > thr$min_features &
      o$nFeature_RNA < thr$max_features &
      o$nCount_RNA < thr$max_counts &
      o$percent_mt < thr$max_percent_mt
    o <- o[, keep]
    rows[[nm]] <- data.frame(
      sample_id = nm,
      subject_id = as.character(meta$subject_id[match(nm, meta$sample_id)]),
      thresholds_from = if (is.null(cfg$gex$qc[[nm]])) "default" else
        "sample-specific",
      min_features = thr$min_features, max_features = thr$max_features,
      max_counts = thr$max_counts, max_percent_mt = thr$max_percent_mt,
      cells_before = n_before, cells_after = ncol(o),
      pct_retained = 100 * ncol(o) / n_before,
      stringsAsFactors = FALSE)
    objects[[nm]] <- o
  }
  qc <- dplyr::bind_rows(rows)
  report$tbl(qc)
  save_table(qc, paste0(prefix, "_qc_per_sample"), run_dir,
             "QC thresholds applied and cells retained per sample")
  default_used <- qc$sample_id[qc$thresholds_from == "default"]
  if (length(default_used)) {
    report$p("")
    report$p(sprintf(paste0(
      "Default thresholds were used for %d of %d libraries (%s). Inspect ",
      "their distributions and add a gex$qc block for any library whose ",
      "metrics differ from the default assumption."),
      length(default_used), nrow(qc),
      paste(utils::head(default_used, 6), collapse = ", ")))
  }
  attr(objects, "qc") <- qc
  invisible(objects)
}

#' Normalisation, merge, variable features and PCA
#'
#' Normalisation is per sample, before the merge: normalising after merging
#' would let one library's depth set the scale for all of them.
#'
#' @param objects Named list of filtered Seurat objects.
#' @inheritParams section_arguments
#' @return The merged object, with the reduction to cluster on in the
#'   attribute `reduction` ("pca" or "harmony").
#' @export
section_gex_reduction <- function(objects, run_dir, cfg, report,
                                  prefix = "03") {
  need(c("Seurat", "SeuratObject"), "section_gex_reduction()")
  report$h1("Normalisation and dimensionality reduction")
  nm <- cfg$gex$normalisation
  objects <- lapply(objects, function(o) {
    o <- Seurat::NormalizeData(o, normalization.method = nm$method,
                               scale.factor = nm$scale_factor, verbose = FALSE)
    Seurat::FindVariableFeatures(o, selection.method = "vst",
                                 nfeatures = nm$n_variable_features,
                                 verbose = FALSE)
  })
  obj <- if (length(objects) == 1L) objects[[1]] else
    merge(objects[[1]], y = objects[-1], add.cell.ids = names(objects))
  obj <- SeuratObject::JoinLayers(obj)
  obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst",
                                      nfeatures = nm$n_variable_features,
                                      verbose = FALSE)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  # Enough components for the widest sweep value, so the sweep never asks
  # for dimensions the reduction does not have.
  npcs <- max(50L, cfg$gex$clustering$n_pcs, cfg$gex$clustering$sweep_n_pcs)
  npcs <- min(npcs, ncol(obj) - 1L, length(SeuratObject::VariableFeatures(obj)))
  obj <- Seurat::RunPCA(obj, npcs = npcs, verbose = FALSE)

  reduction <- "pca"
  if (identical(cfg$gex$integration$method, "harmony")) {
    need("harmony", "Batch correction with Harmony")
    obj <- harmony::RunHarmony(
      obj, group.by.vars = cfg$gex$integration$batch_variable,
      verbose = FALSE)
    reduction <- "harmony"
    report$p(sprintf(paste0(
      "Batch correction: Harmony on '%s'. Integration makes samples resemble ",
      "each other, so where batch and condition coincide it is ",
      "indistinguishable from removing the effect under study. State in ",
      "Methods why correction was applied and on which variable."),
      cfg$gex$integration$batch_variable))
  } else {
    report$p(paste0(
      "Batch correction: none. The merge is uncorrected, so any structure ",
      "that follows subject or library is visible in the embedding rather ",
      "than removed from it; the composition section quantifies it."))
  }
  p <- Seurat::ElbowPlot(obj, ndims = npcs, reduction = reduction) + theme_pub()
  save_figure(p, paste0(prefix, "_elbow"), run_dir, cfg, width = 5,
              height = 3.5,
              legend = "Standard deviation of each principal component.")
  attr(obj, "reduction") <- reduction
  obj
}

#' Clustering with a recorded sweep of dimensionality and resolution
#'
#' @param obj Merged Seurat object from `section_gex_reduction()`.
#' @param reduction Reduction to cluster on.
#' @inheritParams section_arguments
#' @return The clustered object with a UMAP.
#' @export
section_gex_clustering <- function(obj, reduction = "pca", run_dir, cfg,
                                   report, prefix = "04") {
  need(c("Seurat", "dplyr", "ggplot2"), "section_gex_clustering()")
  report$h1("Clustering")
  report$p(paste0(
    "The sweep records how many clusters each combination yields, so that ",
    "the configured choice can be seen as a choice. Cluster count is not a ",
    "property of the data alone."))
  report$p("")
  cl <- cfg$gex$clustering
  avail <- ncol(SeuratObject::Embeddings(obj, reduction = reduction))
  rows <- list()
  for (npc in cl$sweep_n_pcs[cl$sweep_n_pcs <= avail]) {
    tmp <- Seurat::FindNeighbors(obj, reduction = reduction,
                                 dims = seq_len(npc), verbose = FALSE)
    for (res in cl$sweep_resolution) {
      tmp <- Seurat::FindClusters(tmp, resolution = res, verbose = FALSE,
                                  random.seed = cfg$project$seed)
      rows[[length(rows) + 1L]] <- data.frame(
        n_pcs = npc, resolution = res,
        n_clusters = length(unique(Seurat::Idents(tmp))))
    }
  }
  sweep <- dplyr::bind_rows(rows)
  report$tbl(sweep)
  save_table(sweep, paste0(prefix, "_clustering_sweep"), run_dir,
             "clusters obtained per dimensionality and resolution")

  npc <- min(cl$n_pcs, avail)
  obj <- Seurat::FindNeighbors(obj, reduction = reduction, dims = seq_len(npc),
                               verbose = FALSE)
  obj <- Seurat::FindClusters(obj, resolution = cl$resolution, verbose = FALSE,
                              random.seed = cfg$project$seed)
  obj <- Seurat::RunUMAP(obj, reduction = reduction, dims = seq_len(npc),
                         seed.use = cfg$project$seed, verbose = FALSE)
  report$p("")
  report$p(sprintf("Retained: %d components, resolution %.2f, %d clusters.",
                   npc, cl$resolution, length(unique(Seurat::Idents(obj)))))

  for (grp in c("seurat_clusters", "sample_id", "subject_id", "condition")) {
    if (!grp %in% names(obj@meta.data)) next
    if (all(is.na(obj@meta.data[[grp]]))) next
    n_lv <- length(unique(stats::na.omit(obj@meta.data[[grp]])))
    cols <- if (n_lv <= 8) palette_categorical(n_lv) else
      grDevices::hcl.colors(n_lv, "Dark 3")
    p <- Seurat::DimPlot(obj, group.by = grp,
                         label = grp == "seurat_clusters") +
      ggplot2::scale_colour_manual(values = cols, na.value = "grey90") +
      ggplot2::labs(title = grp) + theme_pub(grid = "none")
    save_figure(p, paste0(prefix, "_umap_", grp), run_dir, cfg, width = 6,
                height = 5,
                legend = sprintf("Cells on the UMAP embedding, coloured by %s.",
                                 grp))
  }
  obj
}

#' Subject composition of each cluster: the confounding check
#'
#' A cluster made almost entirely of one subject's cells cannot support a
#' comparison between conditions: a difference in its size between arms is
#' a difference between subjects.
#'
#' @param obj Clustered Seurat object.
#' @inheritParams section_arguments
#' @return The dominance table, invisibly.
#' @export
section_cluster_composition <- function(obj, run_dir, cfg, report,
                                        prefix = "05") {
  need(c("dplyr", "ggplot2"), "section_cluster_composition()")
  report$h1("Cluster composition")
  report$p(paste0(
    "A cluster made almost entirely of one subject's cells cannot support a ",
    "comparison between conditions: a difference in its size between arms ",
    "is a difference between subjects. This check is run before any ",
    "cluster-proportion comparison, and the flagged clusters are named."))
  report$p("")
  md <- obj@meta.data
  comp <- dplyr::count(md, .data$seurat_clusters, .data$subject_id,
                       name = "n_cells")
  comp <- dplyr::ungroup(dplyr::mutate(
    dplyr::group_by(comp, .data$seurat_clusters),
    fraction = .data$n_cells / sum(.data$n_cells),
    cluster_cells = sum(.data$n_cells)))
  save_table(comp, paste0(prefix, "_cluster_composition"), run_dir,
             "cells per cluster per subject")
  thr <- cfg$gex$max_single_subject_fraction
  dominated <- dplyr::summarise(
    dplyr::group_by(comp, .data$seurat_clusters),
    cluster_cells = dplyr::first(.data$cluster_cells),
    max_fraction = max(.data$fraction),
    dominant_subject = as.character(.data$subject_id)[which.max(.data$fraction)],
    .groups = "drop")
  dominated$flagged <- dominated$max_fraction > thr
  report$tbl(as.data.frame(dominated))
  save_table(dominated, paste0(prefix, "_cluster_dominance"), run_dir,
             "largest single-subject fraction per cluster")
  flagged <- dominated$seurat_clusters[dominated$flagged]
  report$p("")
  if (length(flagged)) {
    report$p(sprintf(paste0(
      "Clusters %s exceed %.0f%% from a single subject. No proportion based ",
      "on these clusters can be compared between conditions. Either exclude ",
      "them from such comparisons and say so, or treat the subject effect ",
      "explicitly."), paste(flagged, collapse = ", "), 100 * thr))
  } else {
    report$p(sprintf("No cluster exceeds %.0f%% from a single subject.",
                     100 * thr))
  }
  n_sub <- length(unique(stats::na.omit(comp$subject_id)))
  p <- ggplot2::ggplot(comp, ggplot2::aes(.data$seurat_clusters,
                                          .data$fraction,
                                          fill = .data$subject_id)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::geom_hline(yintercept = thr, linetype = "dashed",
                        linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = if (n_sub <= 8)
      palette_categorical(n_sub) else grDevices::hcl.colors(n_sub, "Dark 3")) +
    ggplot2::scale_y_continuous(labels = scales::label_percent()) +
    ggplot2::labs(x = "Cluster", y = "Cells", fill = "Subject",
                  title = "Subject composition of each cluster",
                  subtitle = "Dashed line: dominance threshold") +
    theme_pub()
  save_figure(p, paste0(prefix, "_cluster_composition"), run_dir, cfg,
              width = 8, height = 4,
              legend = paste("Fraction of each cluster contributed by each",
                             "subject. Clusters above the dashed line cannot",
                             "support between-condition comparison of their",
                             "size."))
  invisible(dominated)
}

#' Cell-cycle phase of every cell
#'
#' Scored with the S and G2/M gene lists distributed with Seurat (Tirosh I
#' et al. Science. 2016;352(6282):189-196. doi:10.1126/science.aad0501.
#' PMID 27124452).
#'
#' @param obj Clustered Seurat object.
#' @inheritParams section_arguments
#' @return The object with `S.Score`, `G2M.Score` and `Phase`.
#' @export
section_cell_cycle <- function(obj, run_dir, cfg, report, prefix = "07") {
  need(c("Seurat", "dplyr"), "section_cell_cycle()")
  report$h1("Cell-cycle phase")
  genes <- Seurat::cc.genes.updated.2019
  obj <- Seurat::CellCycleScoring(obj, s.features = genes$s.genes,
                                  g2m.features = genes$g2m.genes,
                                  set.ident = FALSE)
  phase <- dplyr::count(obj@meta.data, .data$sample_id, .data$Phase,
                        name = "n_cells")
  phase <- dplyr::ungroup(dplyr::mutate(
    dplyr::group_by(phase, .data$sample_id),
    fraction = .data$n_cells / sum(.data$n_cells)))
  report$tbl(as.data.frame(phase))
  save_table(phase, paste0(prefix, "_cell_cycle_phase"), run_dir,
             "cell cycle phase composition per sample")
  obj
}

#' Gene signatures from the configuration and the signature file
#'
#' Panels named in `gex$module_scores` come first; a two-column CSV
#' (`signature`, `gene`) at `gex$signature_file`, or at
#' `<reference_dir>/signatures.csv`, adds more. Keeping signatures in a file
#' rather than in code records which publication each came from.
#' @noRd
gex_signature_panels <- function(cfg) {
  panels <- cfg$gex$module_scores %||% list(
    activation = c("IL2RA", "TNFRSF4", "TNFRSF18", "CTLA4", "CD38"),
    exhaustion = c("PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "TOX",
                   "TOX2"))
  f <- cfg$gex$signature_file %||%
    file.path(cfg$paths$reference_dir %||% "reference", "signatures.csv")
  if (!is.null(f) && file.exists(f)) {
    sig <- utils::read.csv(f, stringsAsFactors = FALSE)
    if (!all(c("signature", "gene") %in% names(sig))) {
      stop("The signature file needs columns 'signature' and 'gene': ", f,
           call. = FALSE)
    }
    extra <- split(sig$gene, sig$signature)
    panels <- c(panels, extra[setdiff(names(extra), names(panels))])
  }
  lapply(panels, as.character)
}
