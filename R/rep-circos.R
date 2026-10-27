# ---------------------------------------------------------------------------
# V-J pairing as chord diagrams
#
# Two figures share this file.
#
#   1. PER-SAMPLE V-J PAIRING. For each sample, the ten most used V genes and
#      every J gene they pair with, as a directional chord diagram (V to J)
#      under a title band; then a table of each sample's top V genes, a bar
#      chart comparing the samples and a row-scaled heatmap.
#
#   2. V-J PAIRING OF SIGNIFICANTLY EXPANDED CLONOTYPES (bulk, longitudinal).
#      One chord diagram per post-baseline timepoint pooled over subjects,
#      one per subject, and the whole baseline repertoire for reference.
#
# WHY THE COLOURS ARE FIXED BEFORE ANY PANEL IS DRAWN. A chord diagram is
# read by colour. If each panel assigns colours to its own genes, the same V
# gene is blue in one sample and orange in the next, and the reader ends up
# comparing colours that mean different things. The mapping is therefore
# computed once over the union of the genes shown in any panel and written
# to a table beside the figures. In the second figure the sector ORDER is
# fixed as well, so a gene also keeps its position on the circle.
#
# RIBBON WIDTH is a count, and what is counted differs between the figures:
# rows of the input in the first (chains, or cells when the input has one
# chain per cell); distinct expanded clonotypes in the second, except for
# its baseline panel, which is weighted by abundance because no clonotype
# is significant at the baseline the test is relative to. Every legend
# states which.
#
# The diagrams are drawn with circlize (Gu et al. 2014, Bioinformatics
# 30:2811-2812, PMID 24930139). circlize draws with base graphics, so
# save_base_figure() writes them in the configured formats and records the
# legend, as save_figure() does for ggplot figures.
# ---------------------------------------------------------------------------

#' V gene palette of the per-sample chord diagrams (twenty colours)
#' @noRd
VJ_V_PALETTE <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
  "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5"
)

#' J gene palette of the per-sample chord diagrams (fifteen colours)
#' @noRd
VJ_J_PALETTE <- c(
  "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
  "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
  "#ccebc5", "#ffed6f", "#ff7f00", "#1b9e77", "#d95f02"
)

#' V palette of the expanded-clonotype chord diagrams (twenty colours)
#' @noRd
VJ_FIXED_V_PALETTE <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
  "#CC79A7", "#000000", "#8C564B", "#17BECF", "#BCBD22", "#9467BD",
  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#E377C2", "#7F7F7F",
  "#AEC7E8", "#FFBB78"
)

#' J palette of the expanded-clonotype chord diagrams (thirteen colours)
#' @noRd
VJ_FIXED_J_PALETTE <- c(
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
  "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928",
  "#999999"
)

# --- internal helpers -------------------------------------------------------

#' Short V or J gene name for sector labels
#'
#' TRBV20-1 -> V20-1, IGHJ4 -> J4. Allele and ambiguity are dropped first.
#' @noRd
vj_short_name <- function(x, segment = c("V", "J")) {
  segment <- match.arg(segment)
  sub(paste0("^(TR[ABGD]|IG[HKL])", segment), segment, strip_allele(x))
}

#' Circos parameters from the configuration, with their defaults
#'
#' Read from cfg$repertoire$circos. The defaults are the values of the
#' original analyses.
#' @noRd
circos_config <- function(cfg) {
  cc <- cfg$repertoire$circos %||% list()
  ov <- cc$v_gene_colours
  out <- list(
    top_v               = as.integer(cc$top_v %||% 10L),
    productive_only     = isTRUE(cc$productive_only %||% FALSE),
    v_gene_colours      = if (length(ov)) unlist(ov) else NULL,
    expanded_top_v      = as.integer(cc$expanded_top_v %||% 15L),
    expanded_top_j      = as.integer(cc$expanded_top_j %||% 13L),
    expanded_coloured_v = as.integer(cc$expanded_coloured_v %||% 10L),
    per_subject         = isTRUE(cc$per_subject %||% TRUE)
  )
  ints <- c(out$top_v, out$expanded_top_v, out$expanded_top_j,
            out$expanded_coloured_v)
  if (anyNA(ints) || any(ints < 1L)) {
    stop("repertoire$circos: top_v, expanded_top_v, expanded_top_j and ",
         "expanded_coloured_v must be positive integers.", call. = FALSE)
  }
  if (length(out$v_gene_colours) && is.null(names(out$v_gene_colours))) {
    stop("repertoire$circos$v_gene_colours must map gene names to colours.",
         call. = FALSE)
  }
  out
}

#' Genes with the largest total, keeping ties at the boundary
#' @noRd
vj_top_with_ties <- function(pairs, col, k) {
  tot <- tapply(pairs$n, pairs[[col]], sum)
  if (!length(tot)) return(character(0))
  cutoff <- sort(tot, decreasing = TRUE)[min(k, length(tot))]
  names(tot)[tot >= cutoff]
}

#' V and J genes of each clonotype
#'
#' For a clonotype definition whose key does not carry the genes ("aa",
#' "nt"), the pair carried by the largest abundance across the table.
#' @noRd
clonotype_vj_lookup <- function(db, definition, abundance_col = NULL) {
  need("dplyr", "clonotype_vj_lookup()")
  d <- data.frame(clonotype = clonotype_key(db, definition),
                  v_gene = strip_allele(db$v_call),
                  j_gene = strip_allele(db$j_call),
                  w = if (is.null(abundance_col)) 1 else
                    as.numeric(db[[abundance_col]]),
                  stringsAsFactors = FALSE)
  d <- d[!is.na(d$v_gene) & !is.na(d$j_gene) & nzchar(d$v_gene) &
           nzchar(d$j_gene), , drop = FALSE]
  agg <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(d, .data$clonotype, .data$v_gene, .data$j_gene),
    w = sum(.data$w, na.rm = TRUE), .groups = "drop"))
  agg <- agg[order(agg$clonotype, -agg$w, agg$v_gene, agg$j_gene,
                   method = "radix"), , drop = FALSE]
  as.data.frame(agg[!duplicated(agg$clonotype),
                    c("clonotype", "v_gene", "j_gene")])
}

#' Count V-J pairs
#' @noRd
vj_count_pairs <- function(v, j, w = NULL) {
  if (!length(v)) {
    return(data.frame(v_gene = character(0), j_gene = character(0),
                      n = numeric(0), stringsAsFactors = FALSE))
  }
  if (is.null(w)) w <- rep(1, length(v))
  agg <- stats::aggregate(list(n = w), by = list(v_gene = v, j_gene = j),
                          FUN = sum, na.rm = TRUE)
  agg[order(-agg$n, agg$v_gene, agg$j_gene, method = "radix"), , drop = FALSE]
}

# --- base-graphics output -----------------------------------------------------

#' Save a base-graphics figure in every configured format
#'
#' The base-graphics counterpart of `save_figure()`, for figures drawn with
#' circlize or other base-graphics code: opens each configured device in
#' the run's figures directory, calls `draw`, closes the device, and appends
#' the legend to FIGURE_LEGENDS.txt. A format that fails is reported as a
#' warning and the others are still written.
#'
#' @param draw A function drawing the figure on the current device. If it
#'   takes an argument it receives the format ("pdf", "png", "svg", "tiff"),
#'   for layouts that must differ between vector and raster output.
#' @param name File stem, without extension.
#' @param run_dir Run directory from `init_run()`.
#' @param cfg Configuration list (`output$formats`, `output$dpi`).
#' @param width,height Size in inches.
#' @param legend Optional figure legend, appended to FIGURE_LEGENDS.txt.
#' @param px_width,px_height Raster size in pixels; by default width and
#'   height times `res`.
#' @param res Raster resolution in pixels per inch; by default
#'   `cfg$output$dpi`.
#' @return The paths written, invisibly.
#' @export
save_base_figure <- function(draw, name, run_dir, cfg, width = 7, height = 7,
                             legend = NULL, px_width = NULL, px_height = NULL,
                             res = NULL) {
  stopifnot(is.function(draw))
  res <- res %||% cfg$output$dpi %||% 300
  px_width  <- px_width %||% round(width * res)
  px_height <- px_height %||% round(height * res)
  fig_dir <- file.path(run_dir, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  written <- character(0)
  for (fmt in cfg$output$formats %||% c("pdf", "png")) {
    path <- file.path(fig_dir, paste0(name, ".", fmt))
    ok <- tryCatch({
      switch(fmt,
             pdf  = grDevices::pdf(path, width = width, height = height,
                                   useDingbats = FALSE),
             svg  = grDevices::svg(path, width = width, height = height),
             png  = grDevices::png(path, width = px_width, height = px_height,
                                   res = res),
             tiff = grDevices::tiff(path, width = px_width, height = px_height,
                                    res = res, compression = "lzw"),
             stop("unsupported format '", fmt, "'", call. = FALSE))
      dev_id <- grDevices::dev.cur()
      tryCatch(if (length(formals(draw))) draw(fmt) else draw(),
               finally = grDevices::dev.off(dev_id))
      TRUE
    }, error = function(e) {
      warning("Could not write ", path, ": ", conditionMessage(e),
              call. = FALSE)
      if (file.exists(path)) unlink(path)
      FALSE
    })
    if (isTRUE(ok)) written <- c(written, path)
  }

  if (!is.null(legend)) {
    append_utf8(c(name, legend, ""), file.path(run_dir, "FIGURE_LEGENDS.txt"))
  }
  invisible(written)
}

# --- 1. per-sample V-J pairing ----------------------------------------------

#' V-J pairing matrix of one sample: its top V genes and every J they use
#'
#' Rows of the locus with a V and a J call are counted; the `top_v` most
#' used V genes are kept with every J gene they pair with. Sector names
#' carry a "V_" or "J_" prefix so a V and a J sector can never share a name;
#' the labels drawn drop it.
#'
#' @param db AIRR table, one row per chain (contig) or per cell.
#' @param sample Sample to use; NULL uses every row.
#' @param top_v Number of V genes kept, by number of rows.
#' @param locus Locus kept (e.g. "TRB"); NULL applies no locus filter.
#' @param productive_only Keep productive rows only. FALSE by default, as in
#'   the original analysis, which counted every chain of the locus that had
#'   a V and a J call.
#' @param v_col,j_col,sample_col Column names.
#' @return NULL when no row is left, else a list: `matrix` (V x J counts, V
#'   rows by decreasing use, J columns by decreasing use), `top_v`, `all_j`,
#'   `sample_id`, `total` (rows of the locus with V and J calls) and
#'   `v_gene_freq` (every V gene of the sample with its count).
#' @export
vj_circos_data <- function(db, sample = NULL, top_v = 10, locus = NULL,
                           productive_only = FALSE, v_col = "v_call",
                           j_col = "j_call", sample_col = "sample_id") {
  tracking_check_columns(db, c(v_col, j_col, if (!is.null(sample)) sample_col),
                         "vj_circos_data()")
  x <- db
  if (!is.null(sample)) {
    x <- x[as.character(x[[sample_col]]) %in% as.character(sample), , drop = FALSE]
  }
  if (!is.null(locus) && "locus" %in% names(x)) {
    x <- x[which(x$locus %in% locus), , drop = FALSE]
  }
  if (productive_only && "productive" %in% names(x)) {
    x <- x[which(parse_productive(x$productive)), , drop = FALSE]
  }
  v <- strip_allele(x[[v_col]])
  j <- strip_allele(x[[j_col]])
  ok <- !is.na(v) & !is.na(j) & nzchar(v) & nzchar(j) &
    v != "None" & j != "None"
  if (!any(ok)) return(NULL)

  v_s <- paste0("V_", vj_short_name(v[ok], "V"))
  j_s <- paste0("J_", vj_short_name(j[ok], "J"))

  tab <- table(v_s)
  v_gene_freq <- data.frame(v_gene = names(tab), n = as.integer(tab),
                            stringsAsFactors = FALSE)
  v_gene_freq <- v_gene_freq[order(-v_gene_freq$n, v_gene_freq$v_gene,
                                   method = "radix"), , drop = FALSE]
  rownames(v_gene_freq) <- NULL
  top <- utils::head(v_gene_freq$v_gene, top_v)

  keep <- v_s %in% top
  m <- table(v_s[keep], j_s[keep])
  m <- matrix(as.numeric(m), nrow = nrow(m),
              dimnames = list(rownames(m), colnames(m)))
  m <- m[order(-rowSums(m), rownames(m), method = "radix"), , drop = FALSE]
  if (ncol(m) > 1L) {
    m <- m[, order(-colSums(m), colnames(m), method = "radix"), drop = FALSE]
  }

  list(matrix = m, top_v = rownames(m), all_j = colnames(m),
       sample_id = if (is.null(sample)) "all" else as.character(sample)[1],
       total = sum(ok), v_gene_freq = v_gene_freq)
}

#' Fixed V and J colours shared by every chord diagram of an analysis
#'
#' V genes take the twenty-colour V palette and J genes the fifteen-colour J
#' palette, in the order given; beyond those, colours continue along viridis
#' (V) and plasma (J). Pass the genes sorted, so the mapping does not depend
#' on which sample happened to be processed first.
#'
#' @param v_genes,j_genes Sector names, as in `vj_circos_data()$matrix`.
#' @param override_v Optional named vector fixing the colour of some V genes;
#'   names may be full gene names (TRBV20-1), short names (V20-1) or sector
#'   names (V_V20-1).
#' @return A list: `v_colors`, `j_colors` (named character vectors).
#' @export
vj_colour_map <- function(v_genes, j_genes, override_v = NULL) {
  v_genes <- unique(as.character(v_genes))
  j_genes <- unique(as.character(j_genes))
  v_pal <- VJ_V_PALETTE
  if (length(v_genes) > length(v_pal)) {
    v_pal <- c(v_pal, scales::viridis_pal()(length(v_genes) - length(v_pal)))
  }
  j_pal <- VJ_J_PALETTE
  if (length(j_genes) > length(j_pal)) {
    j_pal <- c(j_pal, scales::viridis_pal(option = "C")(length(j_genes) -
                                                          length(j_pal)))
  }
  v_map <- stats::setNames(v_pal[seq_along(v_genes)], v_genes)
  j_map <- stats::setNames(j_pal[seq_along(j_genes)], j_genes)

  if (length(override_v)) {
    nm <- names(override_v)
    nm <- ifelse(grepl("^V_", nm), nm, paste0("V_", vj_short_name(nm, "V")))
    common <- intersect(nm, names(v_map))
    v_map[common] <- as.character(override_v)[match(common, nm)]
  }
  list(v_colors = v_map, j_colors = j_map)
}

#' Draw one sample's V-J chord diagram
#'
#' A title band above a directional chord diagram from V to J (height
#' difference and big arrows, links sorted with the largest on top, black
#' link borders), V labels facing clockwise and J labels facing downward,
#' and the totals in the bottom-right corner. Draws on the current device;
#' `save_base_figure()` opens and closes the devices.
#'
#' @param circos_data Output of `vj_circos_data()`.
#' @param colour_map Output of `vj_colour_map()`.
#' @param title Title band, first line.
#' @param subtitle Optional second line, printed in brackets.
#' @param unit What one counted row is, for the corner legend ("chains",
#'   "cells").
#' @param heights Relative heights of the title band and the diagram.
#' @return The matrix drawn, invisibly.
#' @export
draw_vj_circos <- function(circos_data, colour_map, title, subtitle = NULL,
                           unit = "chains", heights = c(1.5, 13.5)) {
  need("circlize", "draw_vj_circos()")
  mat <- circos_data$matrix
  grid_col <- c(colour_map$v_colors[rownames(mat)],
                colour_map$j_colors[colnames(mat)])
  names(grid_col) <- c(rownames(mat), colnames(mat))
  grid_col[is.na(grid_col)] <- "grey85"
  labels <- c(stats::setNames(sub("^V_", "", rownames(mat)), rownames(mat)),
              stats::setNames(sub("^J_", "", colnames(mat)), colnames(mat)))
  v_names <- rownames(mat)

  on.exit(circlize::circos.clear(), add = TRUE)
  graphics::layout(matrix(1:2, nrow = 2), heights = heights)
  graphics::par(mar = c(0, 0, 0.5, 0))
  graphics::plot.new()
  graphics::text(0.5, 0.70, title, cex = 1.7, font = 2)
  if (!is.null(subtitle) && !is.na(subtitle) && nzchar(subtitle)) {
    graphics::text(0.5, 0.20, paste0("(", subtitle, ")"), cex = 1.0, font = 3)
  }

  graphics::par(mar = c(2, 2, 2, 2))
  circlize::circos.clear()
  circlize::circos.par(
    start.degree = 90,
    gap.after = c(rep(3, nrow(mat)), rep(10, ncol(mat))),
    track.margin = c(0.01, 0.01),
    points.overflow.warning = FALSE,
    canvas.xlim = c(-1.3, 1.3),
    canvas.ylim = c(-1.3, 1.3)
  )
  circlize::chordDiagram(
    mat,
    grid.col = grid_col,
    transparency = 0.3,
    directional = 1,
    direction.type = c("diffHeight", "arrows"),
    diffHeight = -circlize::mm_h(3),
    link.arr.type = "big.arrow",
    link.sort = TRUE,
    link.largest.ontop = TRUE,
    link.border = "black",
    link.lwd = 0.7,
    annotationTrack = "grid",
    annotationTrackHeight = circlize::mm_h(8),
    preAllocateTracks = list(track.height = 0.15)
  )
  circlize::circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      sector <- circlize::get.cell.meta.data("sector.index")
      if (sector %in% v_names) {
        circlize::circos.text(circlize::CELL_META$xcenter,
                              circlize::CELL_META$ycenter, labels[[sector]],
                              facing = "clockwise", niceFacing = TRUE,
                              adj = c(0.5, 0.5), cex = 1.1, font = 2)
      } else {
        circlize::circos.text(circlize::CELL_META$xcenter,
                              circlize::CELL_META$ycenter, labels[[sector]],
                              facing = "downward", niceFacing = TRUE,
                              adj = c(0.5, 0.5), cex = 0.9, font = 2)
      }
    },
    bg.border = NA,
    track.height = 0.15
  )
  graphics::legend(
    "bottomright",
    legend = c(sprintf("Total %s: %s", unit,
                       format(circos_data$total, big.mark = ",")),
               paste("V genes shown:", nrow(mat)),
               paste("J genes shown:", ncol(mat))),
    bty = "n", cex = 0.8, text.font = 2)
  invisible(mat)
}

#' Use of each sample's top V genes, compared across samples
#'
#' For every V gene in the top list of at least one sample: its count in
#' every sample, its percentage of the sample's total, its rank in the
#' sample and whether it is in that sample's own top list. A gene outside a
#' sample's top list carries its true count there, not zero; zero means the
#' sample has no row with that gene.
#'
#' @param circos_list Named list of `vj_circos_data()` results, in display
#'   order.
#' @return Data frame: sample_id (factor, display order), v_gene (short
#'   name), count, percentage, rank and in_top.
#' @export
top_v_usage <- function(circos_list) {
  circos_list <- circos_list[!vapply(circos_list, is.null, logical(1))]
  if (!length(circos_list)) return(NULL)
  genes <- unique(unlist(lapply(circos_list, `[[`, "top_v")))
  rows <- lapply(names(circos_list), function(s) {
    cd <- circos_list[[s]]
    vf <- cd$v_gene_freq
    idx <- match(genes, vf$v_gene)
    cnt <- vf$n[idx]
    cnt[is.na(cnt)] <- 0
    data.frame(sample_id = s, v_gene = sub("^V_", "", genes),
               count = cnt, percentage = round(100 * cnt / cd$total, 2),
               rank = idx, in_top = genes %in% cd$top_v,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$sample_id <- factor(out$sample_id, levels = names(circos_list))
  out[order(out$sample_id, -out$count, out$v_gene, method = "radix"), ,
      drop = FALSE]
}

#' Bar chart of the top V genes, one bar per sample
#'
#' Dodged bars with thin black outlines; V genes ordered by their total
#' count over the samples shown.
#'
#' @param usage Output of `top_v_usage()`.
#' @param title,subtitle Plot title and subtitle.
#' @param unit What a counted row is ("chains", "cells"), for the axis.
#' @return A ggplot object.
#' @export
plot_top_v_usage <- function(usage, title = "Top V gene usage",
                             subtitle = "Comparison between samples",
                             unit = "chains") {
  samples <- levels(droplevels(as.factor(usage$sample_id)))
  tot <- tapply(usage$count, usage$v_gene, sum)
  v_order <- names(tot)[order(-tot, names(tot), method = "radix")]
  usage$v_gene <- factor(usage$v_gene, levels = v_order)
  usage$sample_id <- factor(as.character(usage$sample_id), levels = samples)

  ggplot2::ggplot(usage, ggplot2::aes(x = .data$v_gene, y = .data$count,
                                      fill = .data$sample_id)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8),
                      width = 0.7, colour = "black", linewidth = 0.2) +
    ggplot2::scale_fill_manual(
      values = stats::setNames(palette_categorical(length(samples)), samples)) +
    ggplot2::scale_y_continuous(labels = scales::label_comma(),
                                expand = ggplot2::expansion(mult = c(0, 0.1))) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "V gene",
                  y = sprintf("Number of %s", unit), fill = "Sample") +
    theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "bottom")
}

#' Heatmap of the top V genes across samples, scaled by row
#'
#' Counts are centred and scaled within each V gene (a row z-score, as
#' `pheatmap(scale = "row")` does), so the colour shows in which samples a
#' gene is used more or less, not how much it is used overall. Rows and
#' columns are ordered by complete-linkage hierarchical clustering on the
#' scaled values, the pheatmap defaults; the dendrograms are not drawn.
#'
#' @param usage Output of `top_v_usage()`.
#' @param title Plot title.
#' @return A ggplot object.
#' @export
plot_top_v_heatmap <- function(usage, title = "V gene usage heatmap") {
  samples <- levels(droplevels(as.factor(usage$sample_id)))
  genes <- unique(as.character(usage$v_gene))
  m <- matrix(0, nrow = length(genes), ncol = length(samples),
              dimnames = list(genes, samples))
  m[cbind(match(as.character(usage$v_gene), genes),
          match(as.character(usage$sample_id), samples))] <- usage$count

  mu <- rowMeans(m)
  s <- apply(m, 1, stats::sd)
  z <- (m - mu) / s
  z[!is.finite(z)] <- 0

  row_ord <- if (nrow(z) > 1L) stats::hclust(stats::dist(z), "complete")$order else 1L
  col_ord <- if (ncol(z) > 1L) stats::hclust(stats::dist(t(z)), "complete")$order else 1L
  long <- data.frame(
    v_gene = factor(rep(genes, times = length(samples)),
                    levels = rev(genes[row_ord])),
    sample_id = factor(rep(samples, each = length(genes)),
                       levels = samples[col_ord]),
    z = as.vector(z), stringsAsFactors = FALSE)

  ggplot2::ggplot(long, ggplot2::aes(.data$sample_id, .data$v_gene,
                                     fill = .data$z)) +
    ggplot2::geom_tile(colour = "grey60", linewidth = 0.2) +
    ggplot2::scale_fill_viridis_c(name = "Row z-score") +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    theme_pub(grid = "none") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

# --- 2. V-J pairing of significantly expanded clonotypes ---------------------

#' Draw a V-J chord diagram with fixed colours and a fixed sector order
#'
#' The `top_v` V genes and `top_j` J genes with the largest total are kept
#' (ties kept). Sectors follow the order of the names of `v_cols` then
#' `j_cols`, so a gene keeps its position from panel to panel, not only its
#' colour. The title is written left-aligned in a margin reserved above the
#' circle. Draws on the current device.
#'
#' @param pairs Data frame with v_gene, j_gene and n.
#' @param v_cols,j_cols Named colour vectors over every gene of the analysis.
#' @param title Panel title.
#' @param top_v,top_j Number of V and J genes shown.
#' @param label_cex Size of the gene labels.
#' @return The matrix drawn, invisibly (NULL when nothing is left).
#' @export
draw_vj_chord <- function(pairs, v_cols, j_cols, title = NULL, top_v = 15,
                          top_j = 13, label_cex = 0.55) {
  need("circlize", "draw_vj_chord()")
  op <- graphics::par(mar = c(0.5, 0.5, 2.5, 0.5))
  on.exit(graphics::par(op), add = TRUE)

  keep_v <- vj_top_with_ties(pairs, "v_gene", top_v)
  keep_j <- vj_top_with_ties(pairs, "j_gene", top_j)
  d <- pairs[pairs$v_gene %in% keep_v & pairs$j_gene %in% keep_j &
               pairs$n > 0, , drop = FALSE]
  if (!nrow(d)) {
    graphics::plot.new()
    graphics::mtext(paste0(title %||% "", " (no data)"), side = 3, line = 0.7,
                    adj = 0, cex = 0.9, font = 2)
    return(invisible(NULL))
  }

  vs <- sort(unique(d$v_gene), method = "radix")
  js <- sort(unique(d$j_gene), method = "radix")
  m <- matrix(0, nrow = length(vs), ncol = length(js), dimnames = list(vs, js))
  m[cbind(match(d$v_gene, vs), match(d$j_gene, js))] <- d$n

  grid_col <- c(v_cols[rownames(m)], j_cols[colnames(m)])
  names(grid_col) <- c(rownames(m), colnames(m))
  grid_col[is.na(grid_col)] <- "grey85"
  ord <- c(intersect(names(v_cols), rownames(m)),
           intersect(names(j_cols), colnames(m)))
  ord <- c(ord, setdiff(c(rownames(m), colnames(m)), ord))

  on.exit(circlize::circos.clear(), add = TRUE)
  circlize::circos.clear()
  circlize::circos.par(gap.degree = 2, start.degree = 90,
                       points.overflow.warning = FALSE)
  circlize::chordDiagram(m, grid.col = grid_col, order = ord,
                         transparency = 0.35, annotationTrack = "grid",
                         preAllocateTracks = list(track.height = 0.12))
  circlize::circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    circlize::circos.text(circlize::CELL_META$xcenter,
                          circlize::CELL_META$ylim[1] + 0.3,
                          circlize::CELL_META$sector.index,
                          facing = "clockwise", niceFacing = TRUE,
                          adj = c(0, 0.5), cex = label_cex)
  }, bg.border = NA)
  if (!is.null(title)) {
    graphics::mtext(title, side = 3, line = 0.7, adj = 0, cex = 0.9, font = 2)
  }
  invisible(m)
}

# --- sections -----------------------------------------------------------------

#' Per-sample V-J pairing: chord diagrams, top V genes, comparison
#'
#' For each sample, the `repertoire$circos$top_v` (10) most used V genes and
#' every J gene they pair with, as a chord diagram, with one colour mapping
#' over all samples; then the table of top V genes, a bar chart comparing
#' samples and, with two samples or more, a row-scaled heatmap. Written for
#' the single-cell pipeline; on bulk data it counts rearrangements, that is,
#' clonotype-level usage.
#'
#' @param db AIRR table, one row per chain (contig) or per cell.
#' @param locus Locus to use, as in the configuration.
#' @param run_dir Run directory from `init_run()`.
#' @param cfg Configuration list.
#' @param report Report writer from `open_report()`.
#' @param sample_meta Optional sample metadata (sample_id, condition). Samples
#'   follow the configured condition order when a condition is available,
#'   and the condition is printed under each diagram's title.
#' @param unit What one row of `db` is: "chains" for contigs, "cells" when
#'   the input has one chain per cell.
#' @param prefix File name prefix.
#' @return Invisibly, a list: `usage`, `colour_map` and `circos` (per-sample
#'   data).
#' @export
section_vj_circos <- function(db, locus, run_dir, cfg, report,
                              sample_meta = NULL, unit = "chains",
                              prefix = "14") {
  need("circlize", "section_vj_circos()")
  report$h1("V-J pairing per sample")
  cc <- circos_config(cfg)
  locus_lab <- paste(locus, collapse = "/")

  samples <- tracking_levels(db$sample_id)
  cond <- NULL
  if (!is.null(sample_meta) && all(c("sample_id", "condition") %in% names(sample_meta))) {
    cond <- stats::setNames(as.character(sample_meta$condition),
                            as.character(sample_meta$sample_id))
    cond_lv <- tracking_levels(sample_meta$condition, cfg$cohort$conditions)
    samples <- samples[order(match(cond[samples], cond_lv), na.last = TRUE)]
  }

  circos_list <- stats::setNames(lapply(samples, function(s) {
    vj_circos_data(db, sample = s, top_v = cc$top_v, locus = locus,
                   productive_only = cc$productive_only)
  }), samples)
  empty <- names(circos_list)[vapply(circos_list, is.null, logical(1))]
  circos_list <- circos_list[!vapply(circos_list, is.null, logical(1))]
  if (length(empty)) {
    report$p(sprintf("No %s chain with V and J calls in: %s.", locus_lab,
                     paste(empty, collapse = ", ")))
  }
  if (!length(circos_list)) {
    report$p("Nothing to draw.")
    return(invisible(NULL))
  }

  all_v <- sort(unique(unlist(lapply(circos_list, `[[`, "top_v"))),
                method = "radix")
  all_j <- sort(unique(unlist(lapply(circos_list, `[[`, "all_j"))),
                method = "radix")
  colour_map <- vj_colour_map(all_v, all_j, override_v = cc$v_gene_colours)
  colours <- data.frame(
    gene = c(sub("^V_", "", names(colour_map$v_colors)),
             sub("^J_", "", names(colour_map$j_colors))),
    sector = c(names(colour_map$v_colors), names(colour_map$j_colors)),
    type = c(rep("V", length(colour_map$v_colors)),
             rep("J", length(colour_map$j_colors))),
    colour = unname(c(colour_map$v_colors, colour_map$j_colors)),
    stringsAsFactors = FALSE)
  save_table(colours, paste0(prefix, "_vj_colour_mapping"), run_dir,
             "colour of every V and J gene, shared by all chord diagrams")

  report$p(sprintf(paste0(
    "For each sample, the %d most used %s V genes and every J gene they pair ",
    "with. Ribbon width: number of %s carrying the pair. One colour mapping ",
    "is used for every sample (table %s_vj_colour_mapping), so a colour ",
    "means the same gene in every diagram.%s"), cc$top_v, locus_lab, unit,
    prefix, if (cc$productive_only) " Productive chains only." else
      " All chains with a V and a J call are counted, productive or not."))

  for (s in names(circos_list)) {
    cd <- circos_list[[s]]
    sub_title <- if (!is.null(cond)) cond[[s]] else NULL
    title <- sprintf("%s repertoire - %s", locus_lab, s)
    save_base_figure(
      function(fmt) draw_vj_circos(cd, colour_map, title = title,
                                   subtitle = sub_title, unit = unit,
                                   heights = if (fmt %in% c("png", "tiff"))
                                     c(2.5, 10.7) else c(1.5, 13.5)),
      name = sprintf("%sa_vj_circos_%s", prefix, tracking_stem(s)),
      run_dir = run_dir, cfg = cfg, width = 12, height = 12.5,
      px_width = 1400, px_height = 1550, res = 150,
      legend = sprintf(paste0(
        "%s V-J pairing in %s: the %d most used V genes and every J gene ",
        "they pair with. Ribbon width: number of %s carrying the pair (total ",
        "%s); arrows run from V to J. Colours are fixed across samples."),
        locus_lab, s, length(cd$top_v), unit,
        format(cd$total, big.mark = ",")))
  }

  usage <- top_v_usage(circos_list)
  if (!is.null(cond)) usage$condition <- unname(cond[as.character(usage$sample_id)])
  save_table(usage, paste0(prefix, "_top_v_genes"), run_dir, sprintf(paste0(
    "top V genes of each sample: count of %s, percentage, rank, whether in ",
    "the sample's own top %d"), unit, cc$top_v))

  report$h2("Per sample")
  per <- do.call(rbind, lapply(names(circos_list), function(s) {
    cd <- circos_list[[s]]
    data.frame(sample_id = s, total = cd$total, v_shown = nrow(cd$matrix),
               j_shown = ncol(cd$matrix),
               top_v_gene = sub("^V_", "", cd$top_v[1]),
               top_v_pct = round(100 * sum(cd$matrix[1, ]) / cd$total, 2),
               stringsAsFactors = FALSE)
  }))
  report$tbl(per)

  n_s <- length(circos_list)
  p_bar <- plot_top_v_usage(
    usage, unit = unit,
    title = sprintf("Top %d %s V gene usage", cc$top_v, locus_lab))
  save_figure(p_bar, paste0(prefix, "b_top_v_usage"), run_dir, cfg,
              width = if (n_s > 4) 14 else 12, height = 8,
              legend = sprintf(paste0(
                "Number of %s using each V gene that is among the %d most used ",
                "of at least one sample, one bar per sample; genes ordered by ",
                "their total over samples. A missing bar means the sample has ",
                "no %s with that gene."), unit, cc$top_v, unit))
  if (n_s >= 2L) {
    p_hm <- plot_top_v_heatmap(usage, title = sprintf("%s V gene usage heatmap",
                                                      locus_lab))
    save_figure(p_hm, paste0(prefix, "c_top_v_heatmap"), run_dir, cfg,
                width = 8, height = 10,
                legend = sprintf(paste0(
                  "Use of the same V genes across samples, scaled within each ",
                  "gene (row z-score of the count of %s). Rows and columns ",
                  "ordered by complete-linkage clustering on the scaled values."),
                  unit))
  }
  invisible(list(usage = usage, colour_map = colour_map, circos = circos_list))
}

#' V-J pairing of significantly expanded clonotypes, against the baseline
#'
#' Bulk pipeline, after `test_clonal_abundance()`. One chord diagram per
#' post-baseline timepoint (expanded clonotypes pooled over subjects), one
#' per subject (all its post-baseline timepoints pooled), and one of the
#' whole baseline repertoire, all with the same colours and the same sector
#' order, plus a combined figure. Colours go to the
#' `repertoire$circos$expanded_coloured_v` (10) V genes most carried by
#' expanded clonotypes together with the `expanded_top_v` (15) most abundant
#' baseline V genes; other V genes are grey.
#'
#' @param db AIRR table with sample_id, subject_id, timepoint, v_call,
#'   j_call and the abundance column.
#' @param calls Output of `test_clonal_abundance()`.
#' @param abundance_col Abundance column of `db` (weights the baseline panel).
#' @param locus Locus, as in the configuration.
#' @param run_dir Run directory from `init_run()`.
#' @param cfg Configuration list.
#' @param report Report writer from `open_report()`.
#' @param unit Unit of the abundance ("molecules").
#' @param prefix File name prefix.
#' @return Invisibly, the V-J pairs of every panel.
#' @export
section_expanded_vj_circos <- function(db, calls, abundance_col, locus,
                                       run_dir, cfg, report,
                                       unit = "molecules", prefix = "13") {
  need(c("circlize", "dplyr"), "section_expanded_vj_circos()")
  report$h1("V-J pairing of significantly expanded clonotypes")
  if (is.null(calls) || !any(calls$sig_expansion, na.rm = TRUE)) {
    report$p("No clonotype called significantly expanded; chord diagrams not drawn.")
    return(invisible(NULL))
  }
  cc <- circos_config(cfg)
  def <- cfg$repertoire$clonotype_definition %||% "strict"
  ref <- cfg$differential_abundance$reference_timepoint
  locus_lab <- paste(locus, collapse = "/")

  dbl <- filter_locus_productive(db, locus)
  genes <- clonotype_vj_lookup(dbl, def, abundance_col)

  sig <- calls[which(calls$sig_expansion), c("subject_id", "sample_id",
                                             "timepoint", "clonotype")]
  sig$subject_id <- as.character(sig$subject_id)
  sig$timepoint_chr <- as.character(sig$timepoint)
  idx <- match(as.character(sig$clonotype), genes$clonotype)
  sig$v_gene <- genes$v_gene[idx]
  sig$j_gene <- genes$j_gene[idx]
  n_no_gene <- sum(is.na(idx))
  sig <- sig[!is.na(idx), , drop = FALSE]
  if (!nrow(sig)) {
    report$p("No V and J call found for the expanded clonotypes.")
    return(invisible(NULL))
  }

  has_base <- !is.null(ref) && "timepoint" %in% names(dbl) &&
    ref %in% as.character(dbl$timepoint)
  vj_base <- if (has_base) {
    b <- dbl[as.character(dbl$timepoint) == ref, , drop = FALSE]
    vj_count_pairs(strip_allele(b$v_call), strip_allele(b$j_call),
                   as.numeric(b[[abundance_col]]))
  } else NULL
  if (!is.null(vj_base)) {
    vj_base <- vj_base[!is.na(vj_base$v_gene) & !is.na(vj_base$j_gene) &
                         nzchar(vj_base$v_gene) & nzchar(vj_base$j_gene), ,
                       drop = FALSE]
  }

  # --- one palette and one sector order for every panel -----------------
  all_v <- sort(unique(c(sig$v_gene, vj_base$v_gene)), method = "radix")
  all_j <- sort(unique(c(sig$j_gene, vj_base$j_gene)), method = "radix")
  occ <- unique(sig[, c("subject_id", "timepoint_chr", "clonotype", "v_gene")])
  v_tab <- table(occ$v_gene)
  v_rank_sig <- names(v_tab)[order(-as.integer(v_tab), names(v_tab),
                                   method = "radix")]
  v_rank_base <- if (!is.null(vj_base) && nrow(vj_base)) {
    tot <- tapply(vj_base$n, vj_base$v_gene, sum)
    utils::head(names(tot)[order(-tot, names(tot), method = "radix")],
                cc$expanded_top_v)
  } else character(0)
  v_top <- unique(c(utils::head(v_rank_sig, cc$expanded_coloured_v), v_rank_base))
  if (length(v_top) > length(VJ_FIXED_V_PALETTE)) {
    report$p(sprintf("%d V genes to colour for %d colours: the palette is recycled.",
                     length(v_top), length(VJ_FIXED_V_PALETTE)))
  }
  v_cols <- stats::setNames(rep("grey85", length(all_v)), all_v)
  v_cols[v_top] <- VJ_FIXED_V_PALETTE[((seq_along(v_top) - 1L) %%
                                         length(VJ_FIXED_V_PALETTE)) + 1L]
  j_cols <- stats::setNames(VJ_FIXED_J_PALETTE[((seq_along(all_j) - 1L) %%
                                                  length(VJ_FIXED_J_PALETTE)) + 1L],
                            all_j)
  save_table(data.frame(gene = c(names(v_cols), names(j_cols)),
                        type = c(rep("V", length(v_cols)), rep("J", length(j_cols))),
                        colour = unname(c(v_cols, j_cols)),
                        stringsAsFactors = FALSE),
             paste0(prefix, "_vj_chord_colours"), run_dir,
             paste0("colour of every V and J gene, shared by all ",
                    "expanded-clonotype chord diagrams"))

  # --- panels ---------------------------------------------------------------
  panels <- list()
  tp_lv <- tracking_levels(calls$timepoint, cfg$cohort$timepoints)
  for (tp in setdiff(tp_lv, ref)) {
    d <- unique(sig[sig$timepoint_chr == tp,
                    c("subject_id", "clonotype", "v_gene", "j_gene")])
    if (!nrow(d)) next
    panels[[length(panels) + 1L]] <- list(
      kind = "timepoint", label = tp,
      pairs = vj_count_pairs(d$v_gene, d$j_gene),
      title = sprintf("V-J pairing - significantly expanded - %s (all subjects)", tp),
      name = sprintf("%sa_vj_chord_%s", prefix, tracking_stem(tp)),
      size = 8,
      legend = sprintf(paste0(
        "V-J pairing of the clonotypes called significantly expanded at %s ",
        "against each subject's baseline, pooled over subjects. Ribbon width: ",
        "number of distinct expanded clonotypes (per subject) carrying the ",
        "pair. The %d most used V and %d most used J genes are shown; colours ",
        "and sector positions are the same in every panel, and V genes outside ",
        "the coloured set are grey."), tp, cc$expanded_top_v, cc$expanded_top_j))
  }
  if (!is.null(vj_base) && nrow(vj_base)) {
    panels[[length(panels) + 1L]] <- list(
      kind = "baseline", label = ref, pairs = vj_base,
      title = sprintf("V-J pairing at baseline (%s, top %d V genes by abundance)",
                      ref, cc$expanded_top_v),
      name = sprintf("%sb_vj_chord_baseline", prefix), size = 8,
      legend = sprintf(paste0(
        "V-J pairing of the whole baseline (%s) repertoire, pooled over ",
        "subjects. Ribbon width: abundance (%s), since no clonotype is ",
        "significant at the baseline the test is relative to. Same genes ",
        "shown, colours and sector positions as the expanded-clonotype panels."),
        ref, unit))
  }
  if (cc$per_subject) {
    for (subj in tracking_levels(calls$subject_id)) {
      d <- unique(sig[sig$subject_id == subj, c("clonotype", "v_gene", "j_gene")])
      if (!nrow(d)) next
      panels[[length(panels) + 1L]] <- list(
        kind = "subject", label = subj,
        pairs = vj_count_pairs(d$v_gene, d$j_gene),
        title = sprintf("V-J pairing - significantly expanded - %s", subj),
        name = sprintf("%sc_vj_chord_%s", prefix, tracking_stem(subj)), size = 6,
        legend = sprintf(paste0(
          "V-J pairing of the clonotypes of %s called significantly expanded ",
          "at any post-baseline timepoint. Ribbon width: number of distinct ",
          "expanded clonotypes carrying the pair. Colours and sector positions ",
          "as in the other panels."), subj))
    }
  }

  for (pn in panels) {
    save_base_figure(
      function() draw_vj_chord(pn$pairs, v_cols, j_cols, title = pn$title,
                               top_v = cc$expanded_top_v,
                               top_j = cc$expanded_top_j),
      name = pn$name, run_dir = run_dir, cfg = cfg, width = pn$size,
      height = pn$size, legend = pn$legend)
  }

  if (length(panels) >= 2L) {
    n_col <- 2L
    n_row <- ceiling(length(panels) / n_col)
    save_base_figure(function() {
      graphics::layout(matrix(seq_len(n_row * n_col), nrow = n_row,
                              byrow = TRUE))
      for (i in seq_along(panels)) {
        draw_vj_chord(panels[[i]]$pairs, v_cols, j_cols,
                      title = sprintf("%s. %s", LETTERS[(i - 1L) %% 26L + 1L],
                                      panels[[i]]$title),
                      top_v = cc$expanded_top_v, top_j = cc$expanded_top_j)
      }
    }, name = paste0(prefix, "d_vj_chord_combined"), run_dir = run_dir,
    cfg = cfg, width = 12, height = max(5, 5 * n_row),
    legend = paste0("All V-J chord panels side by side: expanded clonotypes ",
                    "per timepoint, the baseline repertoire, then expanded ",
                    "clonotypes per subject. Same colours and sector order ",
                    "throughout."))
  }

  all_pairs <- do.call(rbind, lapply(panels, function(pn) {
    if (!nrow(pn$pairs)) return(NULL)
    data.frame(panel = pn$kind, label = pn$label, pn$pairs,
               stringsAsFactors = FALSE)
  }))
  if (!is.null(all_pairs)) {
    save_table(all_pairs, paste0(prefix, "_vj_chord_pairs"), run_dir,
               paste0("V-J pairs drawn in each chord panel (n: clonotypes, ",
                      "or abundance at baseline)"))
  }

  occ_all <- unique(sig[, c("subject_id", "timepoint_chr", "clonotype",
                            "v_gene", "j_gene")])
  vj_sig_all <- vj_count_pairs(occ_all$v_gene, occ_all$j_gene)
  save_table(vj_sig_all, paste0(prefix, "_vj_pairs_expanded"), run_dir,
             paste0("V-J pairs of expanded clonotypes, all timepoints ",
                    "(n: subject x timepoint x clonotype)"))

  report$p(sprintf(paste0(
    "%d significantly expanded clonotype calls with V and J genes%s. ",
    "Panels: %d post-baseline timepoint(s)%s%s. Colours and sector positions ",
    "are fixed across panels (table %s_vj_chord_colours)."),
    nrow(sig),
    if (n_no_gene) sprintf(" (%d without a gene call left out)", n_no_gene) else "",
    sum(vapply(panels, function(p) p$kind == "timepoint", logical(1))),
    if (has_base) ", the baseline repertoire" else ", no baseline repertoire",
    if (cc$per_subject) sprintf(", %d subject(s)",
                                sum(vapply(panels, function(p) p$kind == "subject",
                                           logical(1)))) else "",
    prefix))
  report$h2("Most frequent V-J pairs among expanded clonotypes")
  report$tbl(utils::head(vj_sig_all, 10))
  invisible(all_pairs)
}
