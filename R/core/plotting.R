# ---------------------------------------------------------------------------
# Figure theme and palettes
#
# One theme and one palette source for every figure in every pipeline, so a
# manuscript's panels match without per-script adjustment.
#
# The categorical palette is Okabe–Ito, which stays distinguishable under the
# common forms of colour vision deficiency and in greyscale. Series are
# therefore never distinguished by colour alone: shape or linetype carries
# the same information wherever a plot has more than one series.
#   Okabe M, Ito K. Color universal design (CUD): how to make figures and
#   presentations that are friendly to colorblind people. 2008.
#   https://jfly.uni-koeln.de/color/
# ---------------------------------------------------------------------------

#' Okabe–Ito categorical palette
#'
#' @param n Number of colours required.
#' @return Character vector of hex colours.
palette_categorical <- function(n = 8) {
  base <- c(
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#009E73",  # bluish green
    "#CC79A7",  # reddish purple
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#F0E442",  # yellow
    "#000000"   # black
  )
  if (n <= length(base)) return(base[seq_len(n)])

  # Beyond eight categories, colour stops separating groups reliably. Extend
  # by interpolation but warn: the figure almost certainly needs faceting or
  # grouping instead.
  warning("Requested ", n, " categorical colours. Beyond 8 the categories ",
          "are no longer reliably distinguishable; consider faceting or ",
          "collapsing groups.", call. = FALSE)
  grDevices::colorRampPalette(base)(n)
}

#' Sequential palette for continuous quantities
palette_sequential <- function(n = 256, direction = 1) {
  cols <- grDevices::colorRampPalette(
    c("#FFF7EC", "#FDD49E", "#FC8D59", "#D7301F", "#7F0000"))(n)
  if (direction < 0) rev(cols) else cols
}

#' Diverging palette for signed quantities (log fold change, z-scores)
#'
#' Symmetric about a neutral midpoint, so zero reads as zero.
palette_diverging <- function(n = 256) {
  grDevices::colorRampPalette(
    c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"))(n)
}

#' Resolve condition colours from the configuration
#'
#' Uses explicit colours where the configuration gives them, and the
#' categorical palette for the rest, so a cohort can fix the colour of one
#' arm without having to name all of them.
condition_colours <- function(cfg) {
  conds <- cfg$cohort$conditions
  if (!length(conds)) return(NULL)

  explicit <- cfg$cohort$condition_colours %||% list()
  cols <- setNames(palette_categorical(length(conds)), conds)
  for (nm in intersect(names(explicit), conds)) {
    cols[[nm]] <- explicit[[nm]]
  }
  cols
}

#' Publication theme
#'
#' @param base_size Base font size in points.
#' @param grid Which grid lines to keep: "y" (default), "x", "both", "none".
theme_pub <- function(base_size = 11, grid = "y") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }
  g <- ggplot2::theme_minimal(base_size = base_size)

  th <- g + ggplot2::theme(
    panel.border     = ggplot2::element_rect(fill = NA, colour = "grey25",
                                             linewidth = 0.4),
    axis.ticks       = ggplot2::element_line(colour = "grey25",
                                             linewidth = 0.3),
    axis.text        = ggplot2::element_text(colour = "grey15"),
    axis.title       = ggplot2::element_text(colour = "grey15"),
    # Panel tags (A, B, C) are set in bold at a fixed size so a composite
    # figure assembled with patchwork does not rescale them per panel.
    plot.title       = ggplot2::element_text(face = "bold",
                                             size = base_size * 1.1,
                                             hjust = 0),
    plot.subtitle    = ggplot2::element_text(colour = "grey35",
                                             size = base_size * 0.95,
                                             hjust = 0),
    plot.caption     = ggplot2::element_text(colour = "grey45",
                                             size = base_size * 0.8,
                                             hjust = 0),
    strip.text       = ggplot2::element_text(face = "bold",
                                             size = base_size * 0.95),
    strip.background = ggplot2::element_blank(),
    legend.key       = ggplot2::element_blank(),
    legend.title     = ggplot2::element_text(size = base_size * 0.95),
    plot.margin      = ggplot2::margin(6, 8, 6, 6)
  )

  th + switch(
    grid,
    "y"    = ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                            panel.grid.minor   = ggplot2::element_blank()),
    "x"    = ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                            panel.grid.minor   = ggplot2::element_blank()),
    "both" = ggplot2::theme(panel.grid.minor = ggplot2::element_blank()),
    "none" = ggplot2::theme(panel.grid = ggplot2::element_blank()),
    stop("grid must be one of 'y', 'x', 'both', 'none'.", call. = FALSE)
  )
}

#' Apply the publication theme to every subsequent plot
use_pub_theme <- function(base_size = 11) {
  ggplot2::theme_set(theme_pub(base_size))
  invisible(TRUE)
}

#' Annotate a panel with the number of observations behind it
#'
#' A repertoire figure is only interpretable next to the depth it was
#' computed at. This writes n into the caption rather than leaving the
#' reader to look it up in a supplementary table.
caption_n <- function(n, unit = "cells") {
  ggplot2::labs(caption = sprintf("n = %s %s",
                                  format(n, big.mark = " "), unit))
}

#' Thin a large point cloud or curve for plotting
#'
#' Repertoire scatterplots routinely hold hundreds of thousands of points,
#' most of them overplotted in the dense corner. Thinning the background
#' while keeping every point above a frequency cut-off produces a file that
#' opens in a vector editor, without changing what the figure shows.
#'
#' @param df Data frame.
#' @param keep Logical vector, or name of a logical column: rows always kept.
#' @param max_background Maximum number of background rows retained.
thin_points <- function(df, keep = NULL, max_background = 20000) {
  stopifnot(is.data.frame(df))
  if (is.character(keep) && length(keep) == 1L) keep <- df[[keep]]
  if (is.null(keep)) keep <- rep(FALSE, nrow(df))
  keep <- as.logical(keep) & !is.na(keep)

  fg <- df[keep, , drop = FALSE]
  bg <- df[!keep, , drop = FALSE]
  if (nrow(bg) > max_background) {
    bg <- bg[sample.int(nrow(bg), max_background), , drop = FALSE]
  }
  rbind(bg, fg)
}
