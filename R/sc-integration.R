# ---------------------------------------------------------------------------
# Joining a single-cell repertoire to a transcriptome object
# ---------------------------------------------------------------------------

#' Attach clonotypes and clone sizes to a Seurat object
#'
#' Barcode conventions differ between Cell Ranger contigs and a merged
#' Seurat object (prefixes added by `merge(add.cell.ids = )`, `-1` suffixes),
#' and a silent mismatch attaches the repertoire to nothing while appearing
#' to succeed. The overlap is therefore checked and reported before anything
#' is written, and a zero overlap is reported as a failure, not as an empty
#' result.
#'
#' @param seurat_path Path to an .rds Seurat object (produced by the
#'   single-cell-gex pipeline).
#' @param cell_rows One row per cell with `cell_id`, `sample_id`,
#'   `clonotype` (one receptor chosen per cell).
#' @param cells_per_clone Clone sizes: `sample_id`, `clonotype`, `abundance`.
#' @inheritParams section_arguments
#' @return The annotated object, invisibly, or NULL when nothing matched.
#' @export
section_seurat_attachment <- function(seurat_path, cell_rows, cells_per_clone,
                                      run_dir, cfg, report, prefix = "09") {
  report$h1("Transcriptome integration")

  if (is.null(seurat_path) || !file.exists(seurat_path)) {
    report$p(sprintf(paste0(
      "No Seurat object at %s; integration skipped. The single-cell-gex ",
      "pipeline writes a suitable object (seurat_annotated.rds)."),
      seurat_path %||% "(unset)"))
    return(invisible(NULL))
  }
  need(c("Seurat", "SeuratObject"), "section_seurat_attachment()")
  obj <- readRDS(seurat_path)

  one <- cell_rows[!duplicated(cell_rows$cell_id), c("cell_id", "sample_id",
                                                    "clonotype")]
  one <- merge(one, cells_per_clone, by = c("sample_id", "clonotype"),
               all.x = TRUE)

  matched <- sum(colnames(obj) %in% one$cell_id)
  report$p(sprintf("Cells in the object: %d; matched to a receptor: %d (%.1f%%)",
                   ncol(obj), matched, 100 * matched / max(1, ncol(obj))))

  if (matched == 0L) {
    report$p("")
    report$p(paste0(
      "No barcode matched, so the repertoire was NOT attached. The cell ",
      "names in the object and the contig barcodes follow different ",
      "conventions (library prefix, '-1' suffix); reconcile them before ",
      "integrating."))
    return(invisible(NULL))
  }

  idx <- match(colnames(obj), one$cell_id)
  obj$clonotype   <- one$clonotype[idx]
  obj$clone_cells <- one$abundance[idx]
  obj$clone_size  <- cut(obj$clone_cells, breaks = c(0, 1, 2, 5, 10, Inf),
                         labels = c("1", "2", "3-5", "6-10", ">10"))
  saveRDS(obj, file.path(run_dir, "seurat_with_repertoire.rds"))
  report$p("Repertoire attached; object saved as seurat_with_repertoire.rds.")

  if ("umap" %in% names(obj@reductions)) {
    p <- Seurat::DimPlot(obj, group.by = "clone_size", order = TRUE) +
      ggplot2::scale_colour_manual(values = palette_sequential(5, direction = -1),
                                   na.value = "grey90") +
      ggplot2::labs(title = "Clone size on the transcriptome embedding") +
      theme_pub(grid = "none")
    save_figure(p, paste0(prefix, "_umap_clone_size"), run_dir, cfg,
                width = 6, height = 5,
                legend = paste("Cells on the transcriptome UMAP coloured by",
                               "the number of cells sharing their",
                               "clonotype. Grey: no receptor recovered."))
  }
  invisible(obj)
}
