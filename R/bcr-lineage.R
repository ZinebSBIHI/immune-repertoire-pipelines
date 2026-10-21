# ---------------------------------------------------------------------------
# B-cell receptor lineage trees
#
# One maximum-likelihood tree per clone, built from its heavy chains and
# rooted on the clone's germline (dowser, with phangorn as the tree
# engine). Trees are built for the clones that link the sets of cells an
# analysis compares -- for example clones found both in a sorted
# compartment and among cells that divided in culture -- because those are
# the clones whose genealogy says how the sets are related.
#
# Which sets are compared is a property of the study and is passed in as
# data (column = accepted values); nothing here knows a compartment name.
#
# References
#   Hoehn KB, Pybus OG, Kleinstein SH. Phylogenetic analysis of migration,
#   differentiation, and class switching in B cells. PLoS Comput Biol.
#   2022;18(4):e1009885. doi:10.1371/journal.pcbi.1009885
#   Schliep KP. phangorn: phylogenetic analysis in R. Bioinformatics.
#   2011;27(4):592-593. doi:10.1093/bioinformatics/btq706
# ---------------------------------------------------------------------------

#' Rows belonging to one set: every listed column holds an accepted value
#' @noRd
lineage_in_set <- function(db, set) {
  if (!is.list(set) || !length(set) || is.null(names(set)) ||
      any(!nzchar(names(set)))) {
    stop("Each set must be a named list mapping columns to accepted values, ",
         "e.g. list(compartment = \"memory\").", call. = FALSE)
  }
  require_columns(db, names(set), "A lineage set")
  keep <- rep(TRUE, nrow(db))
  for (col in names(set)) {
    keep <- keep & as.character(db[[col]]) %in% as.character(unlist(set[[col]]))
  }
  keep
}

#' Clones present in every one of several sets of cells
#'
#' A clone is selected when it has at least one heavy chain in each set and
#' at least `min_cells` heavy chains in total. Size is counted over the whole
#' input, not only within the sets, as in the original analysis.
#'
#' @param db Heavy chains with clone assignments (one per cell in
#'   single-cell data).
#' @param sets List of sets. Each set is a named list mapping columns to the
#'   values accepted; a row belongs to the set when every listed column holds
#'   one of its values, for example
#'   `list(list(compartment = "memory"),
#'   list(proliferation = "divided", condition = c("A", "B")))`.
#'   NULL or an empty list selects on size alone.
#' @param min_cells Minimum number of heavy chains in the clone.
#' @param clone_col,subject_col Clone and subject columns; clones are
#'   qualified by subject.
#' @return One row per selected clone, largest first: `clone_key`, the
#'   subject and clone columns, `n_sequences`, and one logical column
#'   `in_set_<k>` per set.
#' @export
select_lineage_clones <- function(db, sets = NULL, min_cells = 3L,
                                  clone_col = "clone_id",
                                  subject_col = "subject_id") {
  heavy <- igh_rows(db)
  require_columns(heavy, clone_col, "select_lineage_clones()")
  heavy <- heavy[!is.na(heavy[[clone_col]]), , drop = FALSE]
  heavy$clone_key <- subject_clone_key(heavy, clone_col, subject_col)

  set_cols <- character(0)
  for (k in seq_along(sets)) {
    col <- paste0("in_set_", k)
    heavy[[col]] <- lineage_in_set(heavy, sets[[k]])
    set_cols <- c(set_cols, col)
  }
  id_cols <- c(intersect(subject_col, names(heavy)), clone_col)

  out <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(heavy, .data$clone_key),
    dplyr::across(dplyr::all_of(id_cols), dplyr::first),
    n_sequences = dplyr::n(),
    dplyr::across(dplyr::all_of(set_cols), any),
    .groups = "drop"))

  keep <- out$n_sequences >= min_cells
  for (col in set_cols) keep <- keep & out[[col]]
  out <- out[which(keep), , drop = FALSE]
  out[order(-out$n_sequences, out$clone_key), , drop = FALSE]
}

#' Dominant V gene, J gene and CDR3 of each clone
#'
#' The most frequent value of each among the clone's heavy chains, ties
#' going to the first in alphabetical order, joined as
#' "IGHV | IGHJ | CDR3 (aa)": the subtitle of each tree in the original
#' figures.
#'
#' @param db Heavy chains with clone assignments.
#' @param clone_col,subject_col Clone and subject columns.
#' @return Data frame with `clone_key` and `label`.
#' @export
clone_dominant_labels <- function(db, clone_col = "clone_id",
                                  subject_col = "subject_id") {
  heavy <- with_junction_aa(with_gene_cols(igh_rows(db)))
  require_columns(heavy, clone_col, "clone_dominant_labels()")
  heavy <- heavy[!is.na(heavy[[clone_col]]), , drop = FALSE]
  key <- subject_clone_key(heavy, clone_col, subject_col)
  mode_of <- function(x) {
    t <- table(x)
    if (!length(t)) NA_character_ else names(which.max(t))
  }
  idx <- split(seq_len(nrow(heavy)), key)
  data.frame(
    clone_key = names(idx),
    label = vapply(idx, function(i) paste(
      mode_of(heavy$v_gene[i]), mode_of(heavy$j_gene[i]),
      mode_of(heavy$junction_aa[i]), sep = " | "), character(1)),
    row.names = NULL, stringsAsFactors = FALSE)
}

#' Maximum-likelihood lineage trees for selected clones
#'
#' As in the original analyses: `dowser::formatClones()` on heavy chains
#' (chain "H"), keeping the `traits` distinct when identical sequences are
#' collapsed, then `dowser::getTrees()` with `build = "pml"`, a
#' maximum-likelihood tree estimated with phangorn and rooted on the
#' germline.
#'
#' @param db Heavy chains with clone assignments and alignments.
#' @param clone_keys Subject-qualified clone keys, as `clone_key` of
#'   `select_lineage_clones()`.
#' @param mode "single_cell" (cells identified by `cell_id`) or "bulk"
#'   (each sequence its own unit, `cell = id = "sequence_id"`).
#' @param traits Columns kept distinct during collapsing; they are also the
#'   columns available to colour the tips.
#' @param germline_col See `observed_mutation_profile()`.
#' @param minseq Minimum number of sequences per clone; by default 3 in
#'   single-cell and 2 in bulk mode, as in the original analyses.
#' @param build Tree-building method passed to `dowser::getTrees()`.
#' @param clone_col,subject_col Clone and subject columns.
#' @param nproc Cores passed to dowser.
#' @return The dowser tree table, one row per clone; its `clone_id` holds
#'   the clone key.
#' @export
build_lineage_trees <- function(db, clone_keys, mode = c("single_cell", "bulk"),
                                traits = NULL, germline_col = NULL,
                                minseq = NULL, build = "pml",
                                clone_col = "clone_id",
                                subject_col = "subject_id", nproc = 1L) {
  mode <- match.arg(mode)
  need("dowser", "build_lineage_trees()")
  germline_col <- selection_germline_col(db, germline_col)
  minseq <- as.integer(minseq %||% if (mode == "single_cell") 3L else 2L)

  heavy <- selection_heavy(db, germline_col)
  require_columns(heavy, clone_col, "build_lineage_trees()")
  heavy$.clone_key <- subject_clone_key(heavy, clone_col, subject_col)
  tree_db <- heavy[heavy$.clone_key %in% clone_keys, , drop = FALSE]
  if (!nrow(tree_db)) {
    stop("None of the requested clones has heavy chains with both a sequence ",
         "and a germline alignment.", call. = FALSE)
  }
  # dowser refuses factor traits.
  traits <- intersect(traits, names(tree_db))
  for (tr in traits) tree_db[[tr]] <- as.character(tree_db[[tr]])

  args <- list(tree_db, clone = ".clone_key", chain = "H", germ = germline_col,
               traits = if (length(traits)) traits else NULL, minseq = minseq,
               nproc = nproc)
  if (mode == "bulk") {
    args$id <- "sequence_id"
    args$cell <- "sequence_id"
  }
  clones <- do.call(dowser::formatClones, args)
  dowser::getTrees(clones, build = build, nproc = nproc)
}

#' What is missing to draw trees, if anything
#'
#' dowser::plotTrees() needs ggtree; some dowser builds ship without the
#' plotting functions.
#' @noRd
lineage_plot_missing <- function() {
  miss <- character(0)
  has_plot <- !is.null(tryCatch(getExportedValue("dowser", "plotTrees"),
                                error = function(e) NULL))
  if (!has_plot) miss <- c(miss, "dowser::plotTrees()")
  for (p in c("ggtree", "patchwork")) {
    if (!requireNamespace(p, quietly = TRUE)) miss <- c(miss, p)
  }
  miss
}

#' Draw the trees, `per_page` to a page, subtitled with each clone's label
#' @noRd
lineage_plot_pages <- function(trees, labels, tip_col = NULL, tip_size = 2.5,
                               per_page = 6L, ncol = 3L,
                               title = "Lineage trees") {
  plot_trees <- getExportedValue("dowser", "plotTrees")
  plist <- if (is.null(tip_col)) plot_trees(trees, tipsize = tip_size) else
    plot_trees(trees, tips = tip_col, tipsize = tip_size)
  ids <- trees$clone_id
  plist <- lapply(seq_along(plist), function(i) {
    lab <- labels$label[match(ids[i], labels$clone_key)]
    plist[[i]] + ggplot2::labs(subtitle = ifelse(is.na(lab), "", lab))
  })
  chunks <- split(seq_along(plist), ceiling(seq_along(plist) / per_page))
  lapply(seq_along(chunks), function(k) {
    patchwork::wrap_plots(plist[chunks[[k]]], ncol = ncol) +
      patchwork::plot_annotation(
        title = sprintf("%s (page %d)", title, k),
        theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold")))
  })
}

#' Lineage trees of clones shared between sets of cells
#'
#' Selects the clones present in every set and of the minimum size, builds
#' one maximum-likelihood heavy-chain tree per clone, saves the trees (.rds)
#' and a table of the clones, and draws the trees six to a page, each
#' subtitled with the clone's dominant IGHV | IGHJ | CDR3. Settings under
#' `bcr$lineage`: `min_cells` (3), `minseq_single_cell` (3), `minseq_bulk`
#' (2), `build` ("pml"), `trees_per_page` (6), `page_columns` (3),
#' `tip_size` (2.5), `max_clones` (none).
#'
#' Single-cell use, as in the original analysis: clones of at least 3 cells
#' found both in a sorted compartment and among cells that divided in the
#' stimulated conditions; tips coloured by compartment; compartment,
#' condition and subject kept as traits. Bulk use: the original bulk
#' analyses prepared every heavy-chain clone of at least 2 sequences, with
#' subject and isotype as traits, without selecting on sets.
#'
#' If the plotting functions are unavailable, the trees are still built and
#' saved, and the report says what is missing.
#'
#' @inheritParams section_mutation_load
#' @param sets See `select_lineage_clones()`.
#' @param tip_col Column colouring the tips.
#' @param traits Columns kept distinct when identical sequences are
#'   collapsed. By default the tip column, the columns used in `sets` and
#'   the subject, plus `c_call` in bulk mode.
#' @param mode "single_cell" or "bulk".
#' @param clone_col,subject_col Clone and subject columns.
#' @param title Title of the tree pages.
#' @return The dowser tree table, invisibly; NULL when no clone qualifies.
#' @export
section_lineage_trees <- function(db, sets = NULL, tip_col = NULL,
                                  traits = NULL,
                                  mode = c("single_cell", "bulk"),
                                  clone_col = "clone_id",
                                  subject_col = "subject_id",
                                  title = "Lineage trees of selected clones",
                                  name = "lineage_trees",
                                  run_dir, cfg, report, prefix = "14") {
  mode <- match.arg(mode)
  need("dowser", "section_lineage_trees()")
  s <- bcr_advanced_settings(cfg)
  report$h1("Lineage trees")

  germline_col <- selection_germline_col(db, NULL)
  min_cells <- if (mode == "single_cell") s$lineage_min_cells else
    s$lineage_minseq_bulk
  # Selection counts every heavy chain of a clone, as in the original;
  # only the trees need the alignments.
  sel <- select_lineage_clones(db, sets = sets, min_cells = min_cells,
                               clone_col = clone_col, subject_col = subject_col)
  heavy <- selection_heavy(db, germline_col)
  set_text <- if (length(sets)) {
    paste(vapply(sets, function(st) paste(
      sprintf("%s in {%s}", names(st),
              vapply(st, function(v) paste(unlist(v), collapse = ", "),
                     character(1))),
      collapse = " and "), character(1)), collapse = "; ")
  } else {
    "no set restriction"
  }
  report$p(sprintf(paste0(
    "Clones with at least %d heavy chains and present in every set (%s): %d."),
    min_cells, set_text, nrow(sel)))
  if (!nrow(sel)) {
    report$p("No clone qualifies; no tree built.")
    return(invisible(NULL))
  }
  if (!is.null(s$lineage_max_clones) && nrow(sel) > s$lineage_max_clones) {
    report$p(sprintf("Trees limited to the %d largest of them (bcr$lineage$max_clones).",
                     s$lineage_max_clones))
    sel <- utils::head(sel, s$lineage_max_clones)
  }

  if (is.null(traits)) {
    traits <- c(tip_col, unlist(lapply(sets, names)), subject_col,
                if (mode == "bulk") "c_call")
  }
  traits <- intersect(unique(traits), names(heavy))
  tree_db <- heavy[subject_clone_key(heavy, clone_col, subject_col) %in%
                     sel$clone_key, , drop = FALSE]
  labels <- clone_dominant_labels(tree_db, clone_col = clone_col,
                                  subject_col = subject_col)

  minseq <- if (mode == "single_cell") s$lineage_minseq_sc else
    s$lineage_minseq_bulk
  trees <- tryCatch(
    build_lineage_trees(tree_db, sel$clone_key, mode = mode, traits = traits,
                        germline_col = germline_col, minseq = minseq,
                        build = s$lineage_build, clone_col = clone_col,
                        subject_col = subject_col, nproc = s$nproc),
    error = function(e) e)
  if (inherits(trees, "error")) {
    report$p(paste("Tree building failed:", conditionMessage(trees)))
    return(invisible(NULL))
  }

  file_stem <- paste0(prefix, "_", name)
  rds <- file.path(run_dir, paste0(file_stem, ".rds"))
  saveRDS(trees, rds)

  tab <- merge(sel, labels, by = "clone_key", all.x = TRUE, sort = FALSE)
  tab$n_tree_sequences <- trees$seqs[match(tab$clone_key, trees$clone_id)]
  tab <- tab[order(-tab$n_sequences, tab$clone_key), , drop = FALSE]
  report$p(sprintf(paste0(
    "%d trees built (%s, minimum %d sequences per clone; traits kept ",
    "distinct: %s). Trees saved as %s."),
    nrow(trees), s$lineage_build, minseq,
    if (length(traits)) paste(traits, collapse = ", ") else "none",
    basename(rds)))
  report$tbl(as.data.frame(tab[, c("clone_key", "n_sequences",
                                   "n_tree_sequences", "label")]))
  save_table(tab, file_stem, run_dir,
             "clones with lineage trees: size, sequences in the tree, dominant genes and CDR3")

  missing_parts <- lineage_plot_missing()
  if (length(missing_parts)) {
    report$p(sprintf(paste0(
      "Tree figures not drawn: %s not available. The trees in %s can be ",
      "drawn later with dowser::plotTrees()."),
      paste(missing_parts, collapse = ", "), basename(rds)))
    return(invisible(trees))
  }
  pages <- lineage_plot_pages(trees, labels, tip_col = tip_col,
                              tip_size = s$lineage_tip_size,
                              per_page = s$lineage_per_page,
                              ncol = s$lineage_ncol, title = title)
  for (k in seq_along(pages)) {
    save_figure(pages[[k]], sprintf("%s_p%02d", file_stem, k), run_dir, cfg,
                width = 12, height = 7,
                legend = paste0(
                  "Maximum-likelihood heavy-chain trees (dowser). Tips ",
                  "coloured by ", tip_col %||% "(none)", "; subtitle = ",
                  "dominant IGHV | IGHJ | CDR3 (aa); root = germline."))
  }
  invisible(trees)
}
