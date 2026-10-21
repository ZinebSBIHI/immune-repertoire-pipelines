# ---------------------------------------------------------------------------
# Single-cell transcriptome: marker panels, cluster markers and contrasts
#
# Three things a transcriptome analysis has to do after clustering, written
# once here so that the single-cell-gex pipeline and any reanalysis produce
# the same tables:
#
#   1. show a panel of known markers across the clusters, which is what
#      cluster identity is argued from;
#   2. find the markers of each cluster without a prior;
#   3. compare named groups of cells against each other.
#
# WHAT IS DELIBERATELY NOT HERE. Naming a cluster is a judgement about one
# dataset: it depends on that experiment's design, on which markers came
# out, and on the biology under study. So is the order clusters are shown
# in, and so is which contrast matters. All three are read from the
# configuration or passed in; none is written into the code. The functions
# produce the evidence, the project supplies the interpretation.
#
# GENE PANELS ARE DATA, NOT CODE. A panel hard-coded in a function hides
# which publication it came from and cannot be changed without editing the
# analysis. Panels therefore come from `gex$marker_panels` in the
# configuration, as named lists of genes, and every panel is checked against
# the genes actually present in the object: a panel scored or plotted on a
# fraction of its genes is not that panel, and the coverage is reported so
# that a thin panel is not read as a negative result.
#
# WHY COVERAGE IS CHECKED BEFORE PLOTTING. Seurat's DotPlot aborts on a
# feature the object does not carry, so a single gene absent from one
# dataset would take down the whole section. Filtering first, and reporting
# what was dropped, keeps the step informative instead of fragile.
# ---------------------------------------------------------------------------

#' Genes of a panel that the object actually carries
#'
#' @param object A Seurat object.
#' @param panels Named list of gene vectors.
#' @return List with `kept` (the filtered panels, panels left with no gene
#'   removed) and `coverage` (one row per panel: genes requested, genes
#'   present, percentage, and the missing genes).
#' @export
panel_coverage <- function(object, panels) {
  need("SeuratObject", "panel_coverage()")
  present <- rownames(object)
  kept <- lapply(panels, function(g) intersect(as.character(g), present))

  coverage <- data.frame(
    panel = names(panels),
    n_genes = vapply(panels, length, integer(1)),
    n_present = vapply(kept, length, integer(1)),
    stringsAsFactors = FALSE)
  coverage$pct_present <- 100 * coverage$n_present / coverage$n_genes
  coverage$missing <- vapply(seq_along(panels), function(i) {
    miss <- setdiff(as.character(panels[[i]]), present)
    if (!length(miss)) "" else paste(miss, collapse = ",")
  }, character(1))
  rownames(coverage) <- NULL

  list(kept = kept[vapply(kept, length, integer(1)) > 0L],
       coverage = coverage)
}

#' Dot plot of a marker panel across groups of cells
#'
#' Scaled expression on a diverging colour scale with a neutral midpoint, so
#' that "above the mean for this gene" and "below it" read as opposite
#' directions rather than as two shades of one colour. Dot size is the
#' fraction of cells in the group expressing the gene, which is the part of
#' a dot plot that distinguishes a marker expressed by a few cells very
#' strongly from one expressed by most of the group.
#'
#' @param object A Seurat object.
#' @param features Gene vector, or a named list of gene vectors, which
#'   facets the panel by group of genes.
#' @param group_by Metadata column defining the groups on the y axis.
#' @param title Plot title.
#' @param y_levels Optional explicit order for the y axis, which the
#'   configuration supplies; the order of clusters is a presentation choice
#'   and is not decided here.
#' @param midpoint Value of scaled expression shown as neutral.
#' @param dot_scale Dot size scaling, passed to Seurat.
#' @return A ggplot object.
#' @export
marker_dotplot <- function(object, features, group_by = "seurat_clusters",
                           title = NULL, y_levels = NULL, midpoint = 0.1,
                           dot_scale = 8) {
  need("Seurat", "marker_dotplot()")
  div <- palette_diverging(3)

  p <- Seurat::DotPlot(object, features = features, group.by = group_by,
                       col.min = -1, col.max = 1, dot.scale = dot_scale,
                       scale = TRUE) +
    ggplot2::scale_colour_gradient2(low = div[1], mid = div[2],
                                    high = div[3], midpoint = midpoint,
                                    space = "Lab",
                                    name = "Scaled\nexpression") +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    theme_pub(grid = "none") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = ggplot2::element_text(size = 9),
      strip.text = ggplot2::element_text(size = 8))

  if (!is.null(y_levels)) {
    p <- p + ggplot2::scale_y_discrete(limits = rev(y_levels))
  }
  p
}

#' Markers of every cluster
#'
#' Wilcoxon test of each gene in one cluster against all other cells, which
#' is Seurat's default and what cluster identity is usually argued from.
#' Positive markers only by default: the genes a cluster expresses more than
#' the rest are what name it, and the negative side of the same test is the
#' mirror image of the other clusters' positive markers.
#'
#' @param object A Seurat object.
#' @param group_by Metadata column defining the clusters.
#' @param assay Assay to test. The RNA assay is the right one for
#'   differential expression even when the embedding came from an
#'   integrated assay: integration alters expression values to make samples
#'   resemble each other, which is not what should be tested.
#' @param only_pos Keep only genes enriched in the cluster.
#' @param min_pct Minimum fraction of cells expressing the gene in either
#'   group.
#' @param logfc_threshold Minimum log fold change to test.
#' @return Data frame of markers with a `gene` column, or an empty data
#'   frame when no gene passed.
#' @export
cluster_markers <- function(object, group_by = "seurat_clusters",
                            assay = "RNA", only_pos = TRUE, min_pct = 0.25,
                            logfc_threshold = 0.25) {
  need(c("Seurat", "SeuratObject"), "cluster_markers()")
  SeuratObject::DefaultAssay(object) <- assay
  if (!group_by %in% names(object@meta.data)) {
    stop("cluster_markers() needs metadata column '", group_by, "'.",
         call. = FALSE)
  }
  SeuratObject::Idents(object) <- group_by

  # Seurat v5 keeps counts and data split across layers after a merge, and
  # the differential-expression tests need them joined. Joining is a no-op
  # on an object that is already joined, so it is attempted and ignored on
  # failure rather than made a precondition the caller has to remember.
  object <- tryCatch(SeuratObject::JoinLayers(object, assay = assay),
                     error = function(e) object)

  out <- Seurat::FindAllMarkers(object, only.pos = only_pos,
                                min.pct = min_pct,
                                logfc.threshold = logfc_threshold,
                                verbose = FALSE)
  if (!nrow(out)) {
    return(data.frame(p_val = numeric(0), avg_log2FC = numeric(0),
                      pct.1 = numeric(0), pct.2 = numeric(0),
                      p_val_adj = numeric(0), cluster = character(0),
                      gene = character(0), stringsAsFactors = FALSE))
  }
  if (!"gene" %in% names(out)) out$gene <- rownames(out)
  rownames(out) <- NULL
  out
}

#' Strongest markers per cluster
#'
#' @param markers Output of [cluster_markers()].
#' @param n Markers to keep per cluster.
#' @param min_log2fc Optional minimum log fold change applied before the
#'   ranking, which is how a heatmap gets genes that separate the clusters
#'   rather than the n best of a weak set.
#' @param by Ranking column: effect size (`avg_log2FC`) or adjusted p-value.
#' @return Data frame, at most `n` rows per cluster.
#' @export
top_markers <- function(markers, n = 10L, min_log2fc = NULL,
                        by = c("avg_log2FC", "p_val_adj")) {
  need("dplyr", "top_markers()")
  by <- match.arg(by)
  if (!nrow(markers)) return(markers)

  if (!is.null(min_log2fc)) {
    markers <- markers[!is.na(markers$avg_log2FC) &
                         markers$avg_log2FC > min_log2fc, , drop = FALSE]
    if (!nrow(markers)) return(markers)
  }
  grouped <- dplyr::group_by(markers, .data$cluster)
  out <- if (by == "avg_log2FC") {
    dplyr::slice_max(grouped, .data$avg_log2FC, n = n, with_ties = FALSE)
  } else {
    dplyr::slice_min(grouped, .data$p_val_adj, n = n, with_ties = FALSE)
  }
  dplyr::ungroup(out)
}

#' Differential expression between two named groups of cells
#'
#' @param object A Seurat object.
#' @param group1,group2 Levels of `group_by` to compare; the fold change is
#'   positive for genes higher in `group1`.
#' @param group_by Metadata column holding the grouping.
#' @param assay Assay to test.
#' @param min_pct,logfc_threshold Test thresholds. Looser than the
#'   cluster-marker defaults on purpose: a contrast between two conditions
#'   of the same cell type is expected to be subtler than the difference
#'   between two clusters, and a threshold tuned for the latter would return
#'   nothing.
#' @return Data frame with `gene`, the test columns, and `comparison`,
#'   `group1`, `group2`, ordered by adjusted p-value.
#' @export
compare_conditions <- function(object, group1, group2,
                               group_by = "condition", assay = "RNA",
                               min_pct = 0.10, logfc_threshold = 0.10) {
  need(c("Seurat", "SeuratObject"), "compare_conditions()")
  if (!group_by %in% names(object@meta.data)) {
    stop("compare_conditions() needs metadata column '", group_by, "'.",
         call. = FALSE)
  }
  levels_present <- unique(as.character(object@meta.data[[group_by]]))
  missing_groups <- setdiff(c(group1, group2), levels_present)
  if (length(missing_groups)) {
    stop("Group(s) not present in '", group_by, "': ",
         paste(missing_groups, collapse = ", "),
         "\nPresent: ", paste(levels_present, collapse = ", "),
         call. = FALSE)
  }

  SeuratObject::DefaultAssay(object) <- assay
  SeuratObject::Idents(object) <- group_by
  object <- tryCatch(SeuratObject::JoinLayers(object, assay = assay),
                     error = function(e) object)

  res <- Seurat::FindMarkers(object, ident.1 = group1, ident.2 = group2,
                             only.pos = FALSE, min.pct = min_pct,
                             logfc.threshold = logfc_threshold,
                             verbose = FALSE)
  if (!nrow(res)) {
    return(data.frame(gene = character(0), p_val = numeric(0),
                      avg_log2FC = numeric(0), p_val_adj = numeric(0),
                      comparison = character(0), group1 = character(0),
                      group2 = character(0), stringsAsFactors = FALSE))
  }
  res$gene <- rownames(res)
  rownames(res) <- NULL
  res$comparison <- paste0(group1, "_vs_", group2)
  res$group1 <- group1
  res$group2 <- group2
  res[order(res$p_val_adj, res$p_val), , drop = FALSE]
}

#' Run a list of configured contrasts
#'
#' @param object A Seurat object.
#' @param contrasts List of contrasts, each a list or character vector
#'   giving the two groups, optionally named.
#' @param group_by Metadata column holding the grouping.
#' @param ... Passed to [compare_conditions()].
#' @return Data frame of every contrast's results stacked, with a
#'   `comparison` column. A contrast whose groups are absent is reported as
#'   a warning and skipped rather than aborting the set.
#' @export
compare_condition_pairs <- function(object, contrasts,
                                    group_by = "condition", ...) {
  need("dplyr", "compare_condition_pairs()")
  rows <- list()

  for (i in seq_along(contrasts)) {
    spec <- contrasts[[i]]
    pair <- if (!is.null(spec$group1) && !is.null(spec$group2)) {
      c(as.character(spec$group1), as.character(spec$group2))
    } else if (length(unlist(spec)) >= 2L) {
      as.character(unlist(spec))[1:2]
    } else {
      warning("Contrast ", i, " does not name two groups; skipped.",
              call. = FALSE)
      next
    }

    res <- try(compare_conditions(object, group1 = pair[1], group2 = pair[2],
                                  group_by = group_by, ...), silent = TRUE)
    if (inherits(res, "try-error")) {
      warning("Contrast ", paste(pair, collapse = " vs "), " skipped: ",
              conditionMessage(attr(res, "condition")), call. = FALSE)
      next
    }
    rows[[length(rows) + 1L]] <- res
  }

  if (!length(rows)) {
    return(data.frame(gene = character(0), comparison = character(0),
                      stringsAsFactors = FALSE))
  }
  dplyr::bind_rows(rows)
}

#' Score gene modules on every cell
#'
#' One score per panel, each the mean expression of the panel's genes minus
#' the mean of a matched set of control genes, which is what makes scores of
#' different panels comparable (Tirosh et al. Science. 2016;352:189-196).
#'
#' Each panel is scored in its own call so that the resulting column can be
#' named after the panel. Scored together, `AddModuleScore()` names its
#' output columns positionally and ignores the names of the list it was
#' given, so a panel dropped for insufficient coverage silently shifts every
#' later panel onto the wrong column.
#'
#' @param object A Seurat object.
#' @param panels Named list of gene vectors.
#' @param assay,layer Assay and layer to score.
#' @param min_genes Minimum genes present for a panel to be scored. A score
#'   over one or two genes is that gene's expression, not a module.
#' @param nbin Expression bins the control genes are drawn from.
#' @param ctrl Control genes per panel gene. By default derived from the
#'   number of genes in the object: Seurat's default of 100 cannot be
#'   sampled from an expression bin of a small or heavily filtered object,
#'   and the failure there is an error about sampling rather than anything
#'   that names the cause.
#' @param seed Random seed, for the control gene sampling.
#' @return List with `object` (carrying one `<panel>_score` column per
#'   scored panel), `coverage` and `scored` (the panel names scored).
#' @export
module_scores <- function(object, panels, assay = "RNA", layer = "data",
                          min_genes = 3L, nbin = 24L, ctrl = NULL,
                          seed = 1234L) {
  need(c("Seurat", "SeuratObject"), "module_scores()")
  SeuratObject::DefaultAssay(object) <- assay

  cov <- panel_coverage(object, panels)
  coverage <- cov$coverage
  coverage$scored <- coverage$n_present >= min_genes
  coverage$score_column <- NA_character_

  # Control genes are sampled without replacement from the expression bin of
  # each panel gene, so the number asked for has to fit inside a bin.
  if (is.null(ctrl)) {
    ctrl <- min(100L, max(1L, floor(nrow(object) / nbin) - 1L))
  }

  scored <- character(0)
  for (i in which(coverage$scored)) {
    nm <- coverage$panel[i]
    genes <- intersect(as.character(panels[[nm]]), rownames(object))
    col <- paste0(make.names(nm), "_score")
    out <- try(Seurat::AddModuleScore(object, features = list(genes),
                                      name = col, assay = assay,
                                      layer = layer, nbin = nbin,
                                      ctrl = ctrl, seed = seed),
               silent = TRUE)
    if (inherits(out, "try-error")) {
      coverage$scored[i] <- FALSE
      warning("Panel '", nm, "' could not be scored: ",
              conditionMessage(attr(out, "condition")), call. = FALSE)
      next
    }
    object <- out
    # AddModuleScore appends the cluster index to `name`, so the column it
    # created is "<col>1"; it is renamed to the panel's own name so that
    # downstream code never has to know that.
    produced <- paste0(col, "1")
    if (produced %in% names(object@meta.data)) {
      names(object@meta.data)[names(object@meta.data) == produced] <- col
    }
    coverage$score_column[i] <- col
    scored <- c(scored, nm)
  }

  list(object = object, coverage = coverage, scored = scored)
}

#' Heatmap of selected genes across groups of cells
#'
#' @param object A Seurat object.
#' @param features Genes to show.
#' @param group_by Metadata column defining the column blocks.
#' @param assay Assay to display.
#' @param downsample Cells kept per group. A heatmap of every cell is a
#'   figure whose columns are narrower than a pixel and a file too large to
#'   open; downsampling changes nothing a reader can see.
#' @param title Plot title.
#' @return A ggplot object, or NULL when no feature is present.
#' @export
marker_heatmap <- function(object, features, group_by = "seurat_clusters",
                           assay = "RNA", downsample = 200L, title = NULL) {
  need(c("Seurat", "SeuratObject"), "marker_heatmap()")
  features <- unique(as.character(features))
  features <- intersect(features, rownames(object))
  if (!length(features)) return(NULL)

  SeuratObject::DefaultAssay(object) <- assay
  SeuratObject::Idents(object) <- group_by
  object <- tryCatch(SeuratObject::JoinLayers(object, assay = assay),
                     error = function(e) object)
  # Scaling is restricted to the genes shown: scaling the whole matrix to
  # display a few dozen genes costs time and memory for nothing.
  object <- Seurat::ScaleData(object, features = features, assay = assay,
                              verbose = FALSE)
  sub <- if (is.null(downsample)) object else
    subset(object, downsample = downsample)

  Seurat::DoHeatmap(sub, features = features, group.by = group_by,
                    assay = assay) +
    ggplot2::scale_fill_gradientn(colours = palette_diverging(64),
                                  na.value = "white",
                                  name = "Scaled\nexpression") +
    ggplot2::labs(title = title) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6))
}

# ---------------------------------------------------------------------------
# Sections
# ---------------------------------------------------------------------------

#' Marker panels across clusters and across conditions
#'
#' @param object A clustered Seurat object.
#' @param run_dir,cfg,report Run directory, configuration and report writer.
#' @param prefix File-name prefix.
#' @return The coverage table, invisibly.
#' @export
section_marker_panels <- function(object, run_dir, cfg, report,
                                  prefix = "06b") {
  need("Seurat", "section_marker_panels()")
  report$h1("Marker panels")

  panels <- cfg$gex$marker_panels %||% list()
  if (!length(panels)) {
    report$p(paste0(
      "No marker panel in the configuration (gex$marker_panels), so no ",
      "panel figure is drawn. Panels are configuration rather than code so ",
      "that the publication each came from stays visible and a panel can ",
      "be changed without editing the analysis."))
    return(invisible(NULL))
  }

  group_var <- cfg$gex$comparisons$group_variable %||% "condition"
  cluster_levels <- cfg$gex$cluster_order %||% NULL

  all_cov <- list()
  for (nm in names(panels)) {
    genes <- panels[[nm]]
    # A panel may itself be a named list of gene groups, which facets the
    # dot plot by group; coverage is reported on the flattened panel.
    groups <- if (is.list(genes)) genes else list(genes)
    names(groups) <- if (is.list(genes)) names(genes) else nm

    cov <- panel_coverage(object, groups)
    cov$coverage$panel_set <- nm
    all_cov[[nm]] <- cov$coverage
    if (!length(cov$kept)) {
      report$p(sprintf(paste0("Panel '%s': none of its genes is present in ",
                              "the object; not drawn."), nm))
      next
    }

    features <- if (is.list(genes)) cov$kept else unlist(cov$kept,
                                                         use.names = FALSE)
    for (grp in unique(c("seurat_clusters", group_var))) {
      if (!grp %in% names(object@meta.data)) next
      if (all(is.na(object@meta.data[[grp]]))) next
      p <- try(marker_dotplot(
        object, features = features, group_by = grp,
        title = sprintf("%s by %s", nm, grp),
        y_levels = if (identical(grp, "seurat_clusters")) cluster_levels
                   else NULL), silent = TRUE)
      if (inherits(p, "try-error")) {
        report$p(sprintf("Panel '%s' by %s could not be drawn: %s", nm, grp,
                         conditionMessage(attr(p, "condition"))))
        next
      }
      save_figure(p, sprintf("%s_%s_by_%s", prefix, make.names(nm), grp),
                  run_dir, cfg,
                  width = max(7, 0.22 * length(unlist(features))),
                  height = 5,
                  legend = sprintf(paste0("Expression of the '%s' panel ",
                                          "across %s. Colour is scaled mean ",
                                          "expression, dot size the fraction ",
                                          "of cells in the group expressing ",
                                          "the gene."), nm, grp))
    }
  }

  coverage <- do.call(rbind, all_cov)
  if (!is.null(coverage)) {
    report$h2("Panel coverage")
    report$tbl(as.data.frame(coverage[, c("panel_set", "panel", "n_genes",
                                          "n_present", "pct_present")]),
               digits = 1)
    save_table(coverage, paste0(prefix, "_panel_coverage"), run_dir,
               "genes of each marker panel present in the object")
    report$p("")
    report$p(paste0(
      "Genes absent from the object are listed in the coverage table and ",
      "were dropped before plotting. A panel shown on a fraction of its ",
      "genes is not that panel, so the fraction is reported rather than ",
      "left for the reader to assume."))
  }
  invisible(coverage)
}

#' Cluster markers, the strongest per cluster, and their heatmap
#'
#' @param object A clustered Seurat object.
#' @param run_dir,cfg,report Run directory, configuration and report writer.
#' @param prefix File-name prefix.
#' @return The full marker table, invisibly.
#' @export
section_cluster_markers <- function(object, run_dir, cfg, report,
                                    prefix = "06") {
  need(c("Seurat", "dplyr"), "section_cluster_markers()")
  report$h1("Cluster markers")

  mk <- cfg$gex$markers %||% list()
  min_pct <- as.numeric(mk$min_pct %||% 0.25)
  logfc   <- as.numeric(mk$logfc_threshold %||% 0.25)
  top_n   <- as.integer(mk$top_n %||% 10L)
  heat_fc <- mk$heatmap_min_log2fc %||% 1
  downsample <- as.integer(mk$heatmap_downsample %||% 200L)

  markers <- cluster_markers(object, min_pct = min_pct,
                             logfc_threshold = logfc)
  if (!nrow(markers)) {
    report$p(sprintf(paste0("No gene passed min.pct = %.2f and a log fold ",
                            "change of %.2f. Either the clusters are not ",
                            "separated at this resolution, or the ",
                            "thresholds are too strict for this dataset."),
                     min_pct, logfc))
    return(invisible(markers))
  }

  save_table(markers, paste0(prefix, "_cluster_markers"), run_dir,
             sprintf("positive markers per cluster (Wilcoxon, min.pct %.2f, logFC %.2f)",
                     min_pct, logfc))
  top <- top_markers(markers, n = top_n)
  save_table(top, paste0(prefix, "_cluster_markers_top"), run_dir,
             sprintf("%d strongest markers per cluster", top_n))

  report$p(sprintf("%s markers across %d clusters; the %d strongest per cluster are saved.",
                   format(nrow(markers), big.mark = " "),
                   length(unique(as.character(markers$cluster))), top_n))
  report$p("")
  report$p(paste0(
    "Tested on the RNA assay, one cluster against all other cells. Where ",
    "the embedding came from an integrated assay the test is still run on ",
    "RNA: integration alters expression to make samples resemble each ",
    "other, which is not what should be tested."))

  heat <- top_markers(markers, n = top_n, min_log2fc = heat_fc)
  if (nrow(heat)) {
    p <- marker_heatmap(object, heat$gene, group_by = "seurat_clusters",
                        downsample = downsample,
                        title = sprintf("Top %d markers per cluster", top_n))
    if (!is.null(p)) {
      save_figure(p, paste0(prefix, "a_cluster_marker_heatmap"), run_dir, cfg,
                  width = 10, height = max(5, 0.12 * nrow(heat)),
                  legend = sprintf(paste0("Scaled expression of the %d ",
                                          "strongest markers per cluster ",
                                          "above a log fold change of %g, ",
                                          "on up to %d cells per cluster."),
                                   top_n, heat_fc, downsample))
    }
  } else {
    report$p("")
    report$p(sprintf(paste0("No marker exceeds a log fold change of %g, so ",
                            "the heatmap was not drawn; the full table ",
                            "still holds the weaker markers."), heat_fc))
  }
  invisible(markers)
}

#' Configured contrasts between groups of cells
#'
#' @param object A clustered Seurat object.
#' @param run_dir,cfg,report Run directory, configuration and report writer.
#' @param prefix File-name prefix.
#' @return The stacked contrast results, invisibly.
#' @export
section_condition_contrasts <- function(object, run_dir, cfg, report,
                                        prefix = "06c") {
  need(c("Seurat", "dplyr"), "section_condition_contrasts()")
  report$h1("Differential expression between groups")

  cmp <- cfg$gex$comparisons %||% list()
  group_var <- cmp$group_variable %||% "condition"
  contrasts <- cmp$contrasts %||% list()
  min_pct <- as.numeric(cmp$min_pct %||% 0.10)
  logfc   <- as.numeric(cmp$logfc_threshold %||% 0.10)
  top_n   <- as.integer(cmp$top_n %||% 10L)

  if (!length(contrasts)) {
    report$p(paste0(
      "No contrast in the configuration (gex$comparisons$contrasts), so no ",
      "comparison was run. Which contrast answers the question is a design ",
      "decision and is not inferred from the data: with k groups there are ",
      "k(k-1)/2 possible pairs, and running all of them turns the analysis ",
      "into a search."))
    return(invisible(NULL))
  }
  if (!group_var %in% names(object@meta.data)) {
    report$p(sprintf(paste0("The grouping variable '%s' is not in the ",
                            "object's metadata, so the contrasts were not ",
                            "run."), group_var))
    return(invisible(NULL))
  }

  res <- compare_condition_pairs(object, contrasts, group_by = group_var,
                                 min_pct = min_pct, logfc_threshold = logfc)
  if (!nrow(res)) {
    report$p("No contrast returned a result; see the warnings for which were skipped.")
    return(invisible(res))
  }

  save_table(res, paste0(prefix, "_condition_contrasts"), run_dir,
             sprintf("differential expression for each configured contrast on '%s'",
                     group_var))

  report$p(sprintf(paste0("Contrasts run: %s. Thresholds: min.pct %.2f, log ",
                          "fold change %.2f, looser than the cluster-marker ",
                          "thresholds because a contrast within one cell ",
                          "type is expected to be subtler than the ",
                          "difference between two clusters."),
                   paste(unique(res$comparison), collapse = ", "),
                   min_pct, logfc))
  report$p("")

  counts <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(res, .data$comparison),
    n_tested = dplyr::n(),
    n_significant = sum(.data$p_val_adj < (cfg$statistics$alpha %||% 0.05),
                        na.rm = TRUE),
    n_up_in_group1 = sum(.data$p_val_adj < (cfg$statistics$alpha %||% 0.05) &
                           .data$avg_log2FC > 0, na.rm = TRUE),
    .groups = "drop"))
  report$h2("Genes per contrast")
  report$tbl(as.data.frame(counts))

  top <- dplyr::ungroup(dplyr::slice_min(
    dplyr::group_by(res, .data$comparison), .data$p_val_adj, n = top_n,
    with_ties = FALSE))
  save_table(top, paste0(prefix, "_condition_contrasts_top"), run_dir,
             sprintf("%d most significant genes per contrast", top_n))

  p <- marker_heatmap(object, unique(top$gene), group_by = group_var,
                      downsample = as.integer(cmp$heatmap_downsample %||%
                                                200L),
                      title = "Top differential genes by group")
  if (!is.null(p)) {
    save_figure(p, paste0(prefix, "a_contrast_heatmap"), run_dir, cfg,
                width = 10, height = max(5, 0.12 * length(unique(top$gene))),
                legend = sprintf(paste0("Scaled expression of the %d most ",
                                        "significant genes of each contrast, ",
                                        "with cells grouped by %s."),
                                 top_n, group_var))
  }
  invisible(res)
}

#' Module scores for the configured gene signatures
#'
#' The default panels are an activation and an exhaustion signature, both
#' small and widely used; any project signature replaces them through
#' `gex$module_scores`, and the coverage of each is reported.
#'
#' @param object A clustered Seurat object.
#' @param run_dir,cfg,report Run directory, configuration and report writer.
#' @param prefix File-name prefix.
#' @return List with the scored `object` and the `coverage` table,
#'   invisibly.
#' @export
section_module_scores <- function(object, run_dir, cfg, report,
                                  prefix = "07b") {
  need(c("Seurat", "dplyr"), "section_module_scores()")
  report$h1("Module scores")

  # Configured panels (activation and exhaustion by default) plus any
  # signature file; see gex_signature_panels().
  panels <- gex_signature_panels(cfg)
  min_genes <- as.integer(cfg$gex$module_score_min_genes %||% 3L)

  res <- module_scores(object, panels, min_genes = min_genes,
                       nbin = as.integer(cfg$gex$module_score_nbin %||% 24L),
                       ctrl = cfg$gex$module_score_ctrl,
                       seed = cfg$project$seed %||% 1234L)
  report$tbl(as.data.frame(res$coverage[, c("panel", "n_genes", "n_present",
                                            "pct_present", "scored")]),
             digits = 1)
  save_table(res$coverage, paste0(prefix, "_module_score_coverage"), run_dir,
             "genes of each signature present, and the column holding its score")
  report$p("")
  report$p(sprintf(paste0(
    "Each score is the mean expression of the signature's genes minus the ",
    "mean of a matched control set, so scores of different signatures are ",
    "comparable. A signature with fewer than %d genes present is not ",
    "scored: a score over one or two genes is that gene's expression. ",
    "Coverage is reported so that a low-coverage score is not read as a ",
    "negative result."), min_genes))

  if (!length(res$scored)) {
    report$p("")
    report$p("No signature had enough genes present to be scored.")
    return(invisible(list(object = object, coverage = res$coverage)))
  }

  group_var <- cfg$gex$comparisons$group_variable %||% "condition"
  for (nm in res$scored) {
    col <- res$coverage$score_column[res$coverage$panel == nm]
    for (grp in unique(c("seurat_clusters", group_var))) {
      if (!grp %in% names(res$object@meta.data)) next
      if (all(is.na(res$object@meta.data[[grp]]))) next
      n_lv <- length(unique(stats::na.omit(res$object@meta.data[[grp]])))
      p <- Seurat::VlnPlot(res$object, features = col, group.by = grp,
                           pt.size = 0) +
        ggplot2::scale_fill_manual(values = palette_categorical(n_lv),
                                   guide = "none") +
        ggplot2::labs(x = NULL, y = sprintf("%s score", nm),
                      title = sprintf("%s score by %s", nm, grp)) +
        theme_pub() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                           hjust = 1))
      save_figure(p, sprintf("%s_%s_by_%s", prefix, make.names(nm), grp),
                  run_dir, cfg, width = 7, height = 4,
                  legend = sprintf(paste0("Distribution of the %s module ",
                                          "score across %s. Individual cells ",
                                          "are not drawn."), nm, grp))
    }
  }

  score_cols <- res$coverage$score_column[!is.na(res$coverage$score_column)]
  if ("seurat_clusters" %in% names(res$object@meta.data)) {
    summ <- res$object@meta.data[, c("seurat_clusters", score_cols),
                                 drop = FALSE]
    by_cluster <- dplyr::ungroup(dplyr::summarise(
      dplyr::group_by(summ, .data$seurat_clusters),
      dplyr::across(dplyr::all_of(score_cols),
                    function(x) mean(x, na.rm = TRUE)),
      n_cells = dplyr::n(), .groups = "drop"))
    report$h2("Mean score per cluster")
    report$tbl(as.data.frame(by_cluster), digits = 3)
    save_table(by_cluster, paste0(prefix, "_module_scores_by_cluster"),
               run_dir, "mean module score per cluster")
  }

  invisible(list(object = res$object, coverage = res$coverage))
}
