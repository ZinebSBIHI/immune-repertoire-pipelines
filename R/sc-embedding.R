# ---------------------------------------------------------------------------
# Receptor-sequence embedding and joint RNA + TCR neighbourhoods
#
# Trex turns each cell's CDR3 into a numeric vector (here with its
# "geometric" method, which rotates BLOSUM62-derived amino-acid vectors and
# needs no deep-learning backend), stored as a reduction beside the PCA of
# the transcriptome. Three uses follow the original analysis:
#
#   * a UMAP of the receptor embedding alone;
#   * weighted nearest neighbours (WNN) combining the RNA PCA and the
#     receptor embedding, clustered on the joint graph;
#   * a CoNGA-style reduction to one representative cell per clonotype
#     ("dist": the cell closest to the clonotype's centre in PCA space), so
#     that expanded clonotypes do not dominate the receptor neighbourhoods.
#
#   Trex: Borcherding N et al. Nat Immunol. 2024;25(9):1731-1741.
#     doi:10.1038/s41590-024-01888-9. PMID 39164479.
#   WNN: Hao Y et al. Cell. 2021;184(13):3573-3587.e29.
#     doi:10.1016/j.cell.2021.04.048. PMID 34062119.
#   CoNGA: Schattgen SA et al. Nat Biotechnol. 2022;40(1):54-63.
#     doi:10.1038/s41587-021-00989-2. PMID 34426704.
# ---------------------------------------------------------------------------

#' Embedding parameters from the configuration, with their defaults
#'
#' The defaults are the values of the original analysis: TRB, geometric
#' method with theta = pi, 20 dimensions for every UMAP and for both
#' modalities of the WNN graph, WNN clusters at resolution 0.6 with the SLM
#' algorithm (algorithm = 3), CoNGA "dist" representatives clustered with
#' Seurat's default neighbour dimensions (1:10) and resolution (0.8).
#' @noRd
embedding_config <- function(cfg) {
  e <- cfg$repertoire$embedding %||% list()
  list(
    chain = as.character(e$chain %||% "TRB"),
    method = as.character(e$method %||% "geometric"),
    theta = as.numeric(e$theta %||% pi),
    umap_dims = as.integer(e$umap_dims %||% 20L),
    rna_reduction = as.character(e$rna_reduction %||% "pca"),
    wnn_rna_dims = as.integer(e$wnn_rna_dims %||% 20L),
    wnn_tcr_dims = as.integer(e$wnn_tcr_dims %||% 20L),
    wnn_resolution = as.numeric(e$wnn_resolution %||% 0.6),
    cluster_algorithm = as.integer(e$cluster_algorithm %||% 3L),
    k_nn = as.integer(e$k_nn %||% 20L),
    knn_range = as.integer(e$knn_range %||% 200L),
    conga = isTRUE(e$conga %||% TRUE),
    conga_method = as.character(e$conga_method %||% "dist"),
    conga_neighbour_dims = as.integer(e$conga_neighbour_dims %||% 10L),
    conga_resolution = as.numeric(e$conga_resolution %||% 0.8),
    seed = as.integer(cfg$project$seed %||% 1234L)
  )
}

#' Add scRepertoire-style receptor columns to a Seurat object
#'
#' @param object Seurat object.
#' @param cells Output of `cell_receptor_summary()`.
#' @return The object with `CTaa`, `CTgene`, `TRA_cdr3` and `TRB_cdr3` in
#'   its metadata (NA for unmatched cells); the number of matched cells is
#'   in the attribute `matched`.
#' @export
add_receptor_metadata <- function(object, cells) {
  need(c("Seurat", "SeuratObject"), "add_receptor_metadata()")
  idx <- match(colnames(object), cells$cell_id)
  cols <- intersect(c("CTaa", "CTgene", "TRA_cdr3", "TRB_cdr3"),
                    names(cells))
  md <- cells[idx, cols, drop = FALSE]
  rownames(md) <- colnames(object)
  object <- Seurat::AddMetaData(object, md)
  attr(object, "matched") <- sum(!is.na(idx))
  object
}

#' Receptor embedding with Trex, cells with an incomplete embedding removed
#'
#' Trex returns NA coordinates for cells without a usable chain; they are
#' removed before any neighbour search, which would otherwise fail or treat
#' missing values as coordinates.
#'
#' @param object Seurat object with `CTaa` and `CTgene` metadata.
#' @param chain,method,theta Passed to `Trex::runTrex()`.
#' @param reduction_name Name of the reduction created.
#' @param umap_dims Dimensions used for the receptor UMAP ("Trex.umap").
#' @param seed Random seed.
#' @return List: `object` (cells with a complete embedding), `n_removed`.
#' @export
run_trex_embedding <- function(object, chain = "TRB", method = "geometric",
                               theta = pi, reduction_name = "Trex",
                               umap_dims = 20L, seed = 1234L) {
  need(c("Seurat", "SeuratObject", "Trex"), "run_trex_embedding()")
  set.seed(seed)
  obj <- Trex::runTrex(object, chains = chain, method = method,
                       theta = theta, reduction.name = reduction_name)
  emb <- SeuratObject::Embeddings(obj, reduction = reduction_name)
  bad <- rownames(emb)[apply(emb, 1, function(x) any(is.na(x)))]
  # Cells absent from the embedding altogether are removed too.
  bad <- union(bad, setdiff(colnames(obj), rownames(emb)))
  if (length(bad)) obj <- subset(obj, cells = setdiff(colnames(obj), bad))
  emb <- SeuratObject::Embeddings(obj, reduction = reduction_name)
  dims <- seq_len(min(umap_dims, ncol(emb)))
  obj <- Seurat::RunUMAP(obj, reduction = reduction_name, dims = dims,
                         reduction.name = paste0(reduction_name, ".umap"),
                         reduction.key = "trexUMAP_", seed.use = seed,
                         verbose = FALSE)
  list(object = obj, n_removed = length(bad))
}

#' Weighted nearest neighbours of RNA and receptor embeddings
#'
#' @param object Seurat object with both reductions.
#' @param rna_reduction,tcr_reduction Reduction names.
#' @param rna_dims,tcr_dims Dimensions of each modality.
#' @param resolution,algorithm Clustering of the joint graph ("wsnn").
#' @param k_nn,knn_range Neighbour search sizes (Seurat defaults; lower
#'   them only for very small objects).
#' @param seed Random seed.
#' @return The object with `weighted.nn`, `wsnn`, a `wnn.umap` and an
#'   RNA-only `rna.umap` reduction, and joint clusters in `wsnn_clusters`.
#' @export
multimodal_neighbours <- function(object, rna_reduction = "pca",
                                  tcr_reduction = "Trex", rna_dims = 20L,
                                  tcr_dims = 20L, resolution = 0.6,
                                  algorithm = 3L, k_nn = 20L,
                                  knn_range = 200L, seed = 1234L) {
  need(c("Seurat", "SeuratObject"), "multimodal_neighbours()")
  nd <- function(red, n) {
    seq_len(min(n, ncol(SeuratObject::Embeddings(object, reduction = red))))
  }
  set.seed(seed)
  object <- Seurat::FindMultiModalNeighbors(
    object, reduction.list = list(rna_reduction, tcr_reduction),
    dims.list = list(nd(rna_reduction, rna_dims), nd(tcr_reduction, tcr_dims)),
    modality.weight.name = c("RNA.weight", "TCR.weight"), k.nn = k_nn,
    knn.range = knn_range,
    verbose = FALSE)
  object <- Seurat::RunUMAP(object, nn.name = "weighted.nn",
                            reduction.name = "wnn.umap",
                            reduction.key = "wnnUMAP_", seed.use = seed,
                            verbose = FALSE)
  object <- Seurat::FindClusters(object, graph.name = "wsnn",
                                 resolution = resolution,
                                 algorithm = algorithm, random.seed = seed,
                                 verbose = FALSE)
  object$wsnn_clusters <- object$seurat_clusters
  object <- Seurat::RunUMAP(object, reduction = rna_reduction,
                            dims = nd(rna_reduction, rna_dims),
                            reduction.name = "rna.umap",
                            reduction.key = "rnaUMAP_", seed.use = seed,
                            verbose = FALSE)
  object
}

#' One representative cell per clonotype (CoNGA-style) with its embedding
#'
#' @param object Seurat object with `CTaa` and `CTgene` metadata.
#' @param method "dist" (cell nearest the clonotype centre) or "mean".
#' @param chain,theta Passed to `Trex::runTrex()` (geometric method).
#' @param neighbour_dims,resolution,algorithm Clustering of the
#'   representatives on the receptor embedding.
#' @param umap_dims Dimensions of the receptor UMAP.
#' @param seed Random seed.
#' @return Seurat object of representatives with `Trex` and `Trex.umap`.
#' @export
conga_reduction <- function(object, method = "dist", chain = "TRB",
                            theta = pi, neighbour_dims = 10L,
                            resolution = 0.8, algorithm = 3L,
                            umap_dims = 20L, seed = 1234L) {
  need(c("Seurat", "SeuratObject", "Trex"), "conga_reduction()")
  set.seed(seed)
  rep_obj <- Trex::CoNGAfy(object, method = method)
  rep_obj <- Trex::runTrex(rep_obj, chains = chain, method = "geometric",
                           theta = theta, reduction.name = "Trex")
  emb <- SeuratObject::Embeddings(rep_obj, reduction = "Trex")
  bad <- rownames(emb)[apply(emb, 1, function(x) any(is.na(x)))]
  bad <- union(bad, setdiff(colnames(rep_obj), rownames(emb)))
  # Incomplete embeddings are removed before the neighbour search rather
  # than after it, so that the search never sees a missing coordinate.
  if (length(bad)) {
    rep_obj <- subset(rep_obj, cells = setdiff(colnames(rep_obj), bad))
  }
  emb <- SeuratObject::Embeddings(rep_obj, reduction = "Trex")
  rep_obj <- Seurat::FindNeighbors(
    rep_obj, reduction = "Trex",
    dims = seq_len(min(neighbour_dims, ncol(emb))), verbose = FALSE)
  rep_obj <- Seurat::FindClusters(rep_obj, resolution = resolution,
                                  algorithm = algorithm, random.seed = seed,
                                  verbose = FALSE)
  rep_obj <- Seurat::RunUMAP(rep_obj, reduction = "Trex",
                             dims = seq_len(min(umap_dims, ncol(emb))),
                             reduction.name = "Trex.umap",
                             reduction.key = "trexUMAP_", seed.use = seed,
                             verbose = FALSE)
  rep_obj
}

#' Receptor embedding and joint RNA + TCR neighbourhoods
#'
#' @param object Seurat object of the single-cell transcriptome, with a PCA,
#'   to which the receptors are attached here.
#' @param cells Output of `cell_receptor_summary()` for the same cells.
#' @param run_dir,cfg,report,prefix As for the other sections.
#' @return The embedded object, invisibly (NULL when skipped).
#' @export
section_tcr_embedding <- function(object, cells, run_dir, cfg, report,
                                  prefix = "16") {
  report$h1("Receptor embedding and joint RNA + TCR neighbourhoods")
  if (is.null(object)) {
    report$p("No Seurat object; section skipped.")
    return(invisible(NULL))
  }
  if (!requireNamespace("Trex", quietly = TRUE)) {
    report$p(paste0(
      "The Trex package is not installed (remotes::install_github(",
      "\"BorchLab/Trex\")); section skipped."))
    return(invisible(NULL))
  }
  need(c("Seurat", "SeuratObject", "ggplot2"), "section_tcr_embedding()")
  em <- embedding_config(cfg)
  if (!em$rna_reduction %in% names(object@reductions)) {
    report$p(sprintf("The object has no '%s' reduction; section skipped.",
                     em$rna_reduction))
    return(invisible(NULL))
  }
  if (!all(c("CTaa", "CTgene") %in% colnames(object[[]]))) {
    object <- add_receptor_metadata(object, cells)
  }
  report$p(sprintf(paste0(
    "Trex embedding of %s CDR3s, method '%s' (theta = %.4f). Cells with a ",
    "receptor: %d of %d."), em$chain, em$method, em$theta,
    sum(!is.na(object$CTaa)), ncol(object)))

  tx <- run_trex_embedding(object, chain = em$chain, method = em$method,
                           theta = em$theta, umap_dims = em$umap_dims,
                           seed = em$seed)
  obj <- tx$object
  report$p(sprintf("Cells removed for an incomplete embedding: %d; kept: %d.",
                   tx$n_removed, ncol(obj)))

  obj <- multimodal_neighbours(
    obj, rna_reduction = em$rna_reduction, tcr_reduction = "Trex",
    rna_dims = em$wnn_rna_dims, tcr_dims = em$wnn_tcr_dims,
    resolution = em$wnn_resolution, algorithm = em$cluster_algorithm,
    k_nn = min(em$k_nn, ncol(obj) - 1L), knn_range = em$knn_range,
    seed = em$seed)
  w <- obj$RNA.weight
  report$p(sprintf(paste0(
    "WNN graph on %s (%d dims) and Trex (%d dims); median RNA weight %.2f ",
    "(the share of each cell's neighbourhood defined by expression). ",
    "Joint clusters at resolution %.2f: %d."),
    em$rna_reduction, em$wnn_rna_dims, em$wnn_tcr_dims, stats::median(w),
    em$wnn_resolution, length(unique(obj$wsnn_clusters))))

  ct <- function(o, red, title) {
    Seurat::DimPlot(o, reduction = red, group.by = "CTaa") +
      ggplot2::scale_colour_viridis_d(option = "B") +
      ggplot2::labs(title = title) + theme_pub(grid = "none") +
      Seurat::NoLegend()
  }
  cl <- function(o, red, by, title) {
    Seurat::DimPlot(o, reduction = red, group.by = by) +
      ggplot2::labs(title = title) + theme_pub(grid = "none")
  }
  figs <- list(
    trex_clusters = cl(obj, "Trex.umap", "wsnn_clusters",
                       "Receptor embedding (joint clusters)"),
    trex_clonotypes = ct(obj, "Trex.umap", "Receptor embedding (clonotypes)"),
    wnn_clusters = cl(obj, "wnn.umap", "wsnn_clusters", "WNN RNA + TCR"),
    wnn_clonotypes = ct(obj, "wnn.umap", "WNN RNA + TCR (clonotypes)"),
    rna_clusters = cl(obj, "rna.umap", "wsnn_clusters", "RNA only")
  )
  legends <- c(
    trex_clusters = "UMAP of the Trex receptor embedding, by joint cluster.",
    trex_clonotypes = "UMAP of the Trex embedding coloured by clonotype.",
    wnn_clusters = "UMAP of the WNN graph (RNA + TCR), by joint cluster.",
    wnn_clonotypes = "UMAP of the WNN graph coloured by clonotype.",
    rna_clusters = paste("UMAP of the RNA PCA alone, coloured by the joint",
                         "clusters, for comparison with the WNN UMAP."))
  for (nm in names(figs)) {
    save_figure(figs[[nm]], paste0(prefix, "_", nm), run_dir, cfg,
                width = 6, height = 5, legend = legends[[nm]])
  }

  if (em$conga) {
    rep_obj <- tryCatch(
      conga_reduction(obj, method = em$conga_method, chain = em$chain,
                      theta = em$theta,
                      neighbour_dims = em$conga_neighbour_dims,
                      resolution = em$conga_resolution,
                      algorithm = em$cluster_algorithm,
                      umap_dims = em$umap_dims, seed = em$seed),
      error = function(e) {
        report$p(paste("CoNGA reduction failed:", conditionMessage(e)))
        NULL
      })
    if (!is.null(rep_obj)) {
      report$p(sprintf(paste0(
        "CoNGA ('%s'): %d representative cells, one per clonotype; %d ",
        "clusters on the receptor embedding."), em$conga_method,
        ncol(rep_obj), length(unique(rep_obj$seurat_clusters))))
      sp <- specificity_config(cfg)
      ref <- cell_reference_lookup(sp, em$chain, report)
      if (!is.null(ref) && "CTaa" %in% colnames(rep_obj[[]])) {
        part <- if (em$chain == "TRA") 1L else 2L
        cdr3 <- vapply(strsplit(as.character(rep_obj$CTaa), "_", fixed = TRUE),
                       function(z) if (length(z) >= part) z[part] else NA_character_,
                       character(1))
        ann <- annotate_cells_by_cdr3(cdr3, ref$lookup, max_distance = 0)
        rep_obj$species <- ann$Epitope.species
        p <- Seurat::DimPlot(rep_obj, reduction = "Trex.umap",
                             group.by = "species") +
          ggplot2::labs(title = sprintf("CoNGA representatives by %s species",
                                        em$chain)) +
          theme_pub(grid = "none")
        save_figure(p, paste0(prefix, "_conga_species"), run_dir, cfg,
                    width = 7, height = 5,
                    legend = paste("One representative cell per clonotype on",
                                   "the receptor UMAP, coloured by the",
                                   "antigen species of an exact CDR3 match."))
      }
      saveRDS(rep_obj, file.path(run_dir, "seurat_conga.rds"))
    }
  }
  saveRDS(obj, file.path(run_dir, "seurat_trex_wnn.rds"))
  invisible(obj)
}
