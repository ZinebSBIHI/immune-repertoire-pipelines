# ---------------------------------------------------------------------------
# Antigen specificity of single cells
#
# Each cell's TRA and TRB CDR3s are annotated against a reference of
# receptors with a known specificity: the reference bundled with Trex by
# default, as in the original analysis, or a VDJdb / McPAS-TCR file.
#
# Annotation of a cell, at edit distance k (the rule of Trex::annotateDB):
#   1. exact: its CDR3 is in the reference;
#   2. otherwise, within k edits (0 < d <= k) of a CDR3 of the SAME dataset
#      that matched exactly, whose annotation it takes. The anchors are the
#      exactly matched CDR3s in the order of the cells; the first within
#      reach wins.
# The approximate step therefore spreads an annotation to the near
# neighbours of an exact match, and never links a cell to a reference CDR3
# that no cell carries exactly. Reproducing the rule here, rather than
# calling annotateDB(), lets the same code run on any reference and keeps
# every cell barcode intact: annotateDB() appends rows under duplicated
# barcodes, which R renames, so a cell within reach of two anchors keeps the
# first annotation and the second is silently discarded. The outcome is the
# same; the test suite checks the equivalence against annotateDB() itself.
#
# Clone size categories follow the scRepertoire convention used in the
# original analysis, counted in cells per clonotype within each sample:
# Single (1), Small (2), Medium (3-5), Large (6-10), Hyperexpanded (> 10),
# with Medium and above counted as expanded.
#   Borcherding N, Bormann NL, Kraus G. F1000Res. 2020;9:47.
#   doi:10.12688/f1000research.22139.2. PMID 32789006.
# ---------------------------------------------------------------------------

#' One row per cell with its receptor chains, in the scRepertoire layout
#'
#' Several chains of the same locus in one cell are joined with ";" in
#' decreasing order of abundance, and the two loci with "_", which is the
#' layout of scRepertoire's `CTaa` and `CTgene` columns expected by Trex.
#' A missing chain is written "NA", a missing D gene "None".
#'
#' @param db AIRR data frame with `cell_id`, `locus`, `junction_aa`,
#'   `v_call`, `j_call` and optionally `d_call`, `c_call`.
#' @param abundance_col Column used to order chains within a cell.
#' @param productive_only Keep productive chains only.
#' @return Data frame, one row per cell: cell_id, sample metadata, TRA_cdr3,
#'   TRB_cdr3, TRA_v, TRB_v, TRA_j, TRB_j, CTaa, CTgene.
#' @export
cell_receptor_summary <- function(db, abundance_col = "umi_count",
                                  productive_only = TRUE) {
  for (cl in c("cell_id", "locus", "junction_aa", "v_call", "j_call")) {
    if (!cl %in% names(db)) {
      stop("cell_receptor_summary() needs column '", cl, "'.", call. = FALSE)
    }
  }
  d <- db[db$locus %in% c("TRA", "TRB"), , drop = FALSE]
  if (productive_only && "productive" %in% names(d)) {
    keep <- parse_productive(d$productive)
    d <- d[which(keep), , drop = FALSE]
  }
  d <- d[!is.na(d$junction_aa) & nzchar(d$junction_aa), , drop = FALSE]
  if (!nrow(d)) stop("No productive TRA or TRB chain.", call. = FALSE)
  # Cells keep the order in which they first appear in the input: the
  # propagation step of annotate_cells_by_cdr3() depends on cell order.
  input_order <- unique(as.character(d$cell_id))

  ab <- if (abundance_col %in% names(d)) as.numeric(d[[abundance_col]]) else 0
  ab[is.na(ab)] <- 0
  d <- d[order(d$cell_id, d$locus, -ab, d$junction_aa), , drop = FALSE]
  gene <- function(x) {
    x <- strip_allele(x)
    x[is.na(x) | !nzchar(x)] <- NA_character_
    x
  }
  d$v <- gene(d$v_call)
  d$j <- gene(d$j_call)
  d$dg <- if ("d_call" %in% names(d)) gene(d$d_call) else NA_character_
  d$cg <- if ("c_call" %in% names(d)) gene(d$c_call) else NA_character_

  join <- function(v) paste(ifelse(is.na(v), "NA", v), collapse = ";")
  out <- data.frame(cell_id = input_order, stringsAsFactors = FALSE)
  meta_cols <- intersect(c("sample_id", "subject_id", "condition",
                           "timepoint"), names(d))
  first <- d[!duplicated(d$cell_id), c("cell_id", meta_cols), drop = FALSE]
  for (m in meta_cols) {
    out[[m]] <- first[[m]][match(out$cell_id, first$cell_id)]
  }

  for (ch in c("TRA", "TRB")) {
    sub <- d[d$locus == ch, , drop = FALSE]
    idx <- split(seq_len(nrow(sub)), sub$cell_id)
    pick <- function(col) {
      v <- vapply(idx, function(i) join(sub[[col]][i]), character(1))
      v[match(out$cell_id, names(v))]
    }
    out[[paste0(ch, "_cdr3")]] <- pick("junction_aa")
    out[[paste0(ch, "_v")]] <- pick("v")
    out[[paste0(ch, "_j")]] <- pick("j")
    genes <- vapply(idx, function(i) {
      parts <- if (ch == "TRA") {
        paste(ifelse(is.na(sub$v[i]), "NA", sub$v[i]),
              ifelse(is.na(sub$j[i]), "NA", sub$j[i]),
              ifelse(is.na(sub$cg[i]), "NA", sub$cg[i]), sep = ".")
      } else {
        paste(ifelse(is.na(sub$v[i]), "NA", sub$v[i]),
              ifelse(is.na(sub$dg[i]), "None", sub$dg[i]),
              ifelse(is.na(sub$j[i]), "NA", sub$j[i]),
              ifelse(is.na(sub$cg[i]), "NA", sub$cg[i]), sep = ".")
      }
      paste(parts, collapse = ";")
    }, character(1))
    out[[paste0(ch, "_genes")]] <- genes[match(out$cell_id, names(genes))]
  }
  na_str <- function(x) ifelse(is.na(x), "NA", x)
  out$CTaa <- paste(na_str(out$TRA_cdr3), na_str(out$TRB_cdr3), sep = "_")
  out$CTgene <- paste(na_str(out$TRA_genes), na_str(out$TRB_genes),
                      sep = "_")
  out$TRA_genes <- NULL
  out$TRB_genes <- NULL
  out
}

#' Annotate cells by CDR3, exact match then propagation within the data
#'
#' @param cdr3 Character vector, one CDR3 per cell in cell order (several
#'   chains joined with ";" are matched as one string, as in Trex).
#' @param lookup Output of `collapse_reference()`.
#' @param max_distance Edit distance for the propagation step; 0 for exact
#'   matching only.
#' @return Data frame, one row per cell: match ("exact", "distance" or NA),
#'   matched_cdr3, and the annotation columns named as in Trex:
#'   Epitope.target, Epitope.sequence, Epitope.species, Tissue, Cell.type,
#'   Database.
#' @export
annotate_cells_by_cdr3 <- function(cdr3, lookup, max_distance = 0) {
  need("stringdist", "annotate_cells_by_cdr3()")
  cdr3 <- as.character(cdr3)
  cdr3[cdr3 %in% c("", "NA")] <- NA_character_
  n <- length(cdr3)
  ref_idx <- rep(NA_integer_, n)
  matched <- rep(NA_character_, n)
  how <- rep(NA_character_, n)

  exact <- match(cdr3, lookup$cdr3)
  exact[is.na(cdr3)] <- NA_integer_
  hit <- !is.na(exact)
  ref_idx[hit] <- exact[hit]
  matched[hit] <- cdr3[hit]
  how[hit] <- "exact"

  if (max_distance > 0 && any(hit)) {
    anchors <- unique(cdr3[hit])
    todo <- which(!hit & !is.na(cdr3))
    if (length(todo)) {
      u <- unique(cdr3[todo])
      assigned <- rep(NA_integer_, length(u))
      for (a in seq_along(anchors)) {
        free <- which(is.na(assigned))
        if (!length(free)) break
        dd <- stringdist::stringdist(anchors[a], u[free], method = "lv")
        near <- free[!is.na(dd) & dd > 0 & dd <= max_distance]
        assigned[near] <- a
      }
      a_of_cell <- assigned[match(cdr3[todo], u)]
      ok <- !is.na(a_of_cell)
      cells <- todo[ok]
      matched[cells] <- anchors[a_of_cell[ok]]
      ref_idx[cells] <- match(matched[cells], lookup$cdr3)
      how[cells] <- "distance"
    }
  }

  data.frame(
    match = how,
    matched_cdr3 = matched,
    Epitope.target = lookup$epitope_genes[ref_idx],
    Epitope.sequence = lookup$epitopes[ref_idx],
    Epitope.species = lookup$specificities[ref_idx],
    Tissue = lookup$tissue[ref_idx],
    Cell.type = lookup$cell_type[ref_idx],
    Database = lookup$sources[ref_idx],
    stringsAsFactors = FALSE)
}

#' Clone size category from the number of cells per clonotype
#'
#' @param n_cells Cells sharing the clonotype within the sample.
#' @param breaks Named upper bounds of the categories below the largest;
#'   the defaults are those of the original analysis.
#' @param top Name of the open-ended largest category.
#' @return Factor of category names, smallest first. The attribute
#'   `labels` holds display labels with the intervals, e.g.
#'   "Medium (2 < X <= 5)".
#' @export
clone_size_category <- function(n_cells,
                                breaks = c(Single = 1, Small = 2, Medium = 5,
                                           Large = 10),
                                top = "Hyperexpanded") {
  breaks <- sort(breaks)
  cuts <- c(0, unname(breaks), Inf)
  lv <- c(names(breaks), top)
  out <- cut(as.numeric(n_cells), breaks = cuts, labels = lv, right = TRUE)
  lo <- cuts[-length(cuts)]
  top_max <- suppressWarnings(max(as.numeric(n_cells), na.rm = TRUE))
  labels <- sprintf("%s (%g < X <= %g)", lv[-length(lv)], lo[-length(lo)],
                    unname(breaks))
  # The open category is labelled with the largest clone observed, as
  # scRepertoire does, or as open-ended when no clone reaches it.
  labels <- c(labels, if (is.finite(top_max) && top_max > max(breaks))
    sprintf("%s (%g < X <= %g)", top, max(breaks), top_max) else
      sprintf("%s (X > %g)", top, max(breaks)))
  attr(out, "labels") <- stats::setNames(labels, lv)
  out
}

#' Count the V genes of a set of cells
#'
#' @param genes Per-cell V genes of one chain, several joined with ";".
#' @param top_n Number of genes kept, most frequent first.
#' @return Data frame: v_gene, n, freq (over all genes counted).
#' @export
cell_v_usage <- function(genes, top_n = 30) {
  g <- unlist(strsplit(as.character(genes[!is.na(genes)]), ";", fixed = TRUE))
  g <- g[!is.na(g) & nzchar(g) & g != "NA"]
  if (!length(g)) {
    return(data.frame(v_gene = character(0), n = integer(0),
                      freq = numeric(0)))
  }
  t <- sort(table(g), decreasing = TRUE)
  out <- data.frame(v_gene = names(t), n = as.integer(t),
                    stringsAsFactors = FALSE)
  out$freq <- out$n / sum(out$n)
  utils::head(out, top_n)
}

#' Most frequent CDR3 sequences of a set of cells
#'
#' @param cdr3 Per-cell CDR3 of one chain.
#' @param top_n Number kept.
#' @return Data frame: cdr3, n, freq_within_top (the share among the
#'   sequences kept, as computed in the original analysis).
#' @export
top_cell_cdr3 <- function(cdr3, top_n = 20) {
  x <- as.character(cdr3[!is.na(cdr3) & nzchar(cdr3) & cdr3 != "NA"])
  if (!length(x)) {
    return(data.frame(cdr3 = character(0), n = integer(0),
                      freq_within_top = numeric(0)))
  }
  t <- sort(table(x), decreasing = TRUE)
  out <- data.frame(cdr3 = names(t), n = as.integer(t),
                    stringsAsFactors = FALSE)
  out <- utils::head(out, top_n)
  out$freq_within_top <- out$n / sum(out$n)
  out
}

#' Expansion and specificity groups of cells
#'
#' @param size Clone size category (output of `clone_size_category()`).
#' @param target Logical, cell flagged for the target species.
#' @param target_label Display label of the target species.
#' @param expanded Categories counted as expanded.
#' @param top Name of the largest category.
#' @return Data frame: clone_expansion_group, expanded_clone,
#'   expanded_target_group (factors with the original level order).
#' @export
cell_specificity_groups <- function(size, target, target_label = "Target",
                                    expanded = c("Medium", "Large",
                                                 "Hyperexpanded"),
                                    top = "Hyperexpanded") {
  size <- as.character(size)
  target <- !is.na(target) & target
  g <- rep("Other", length(size))
  g[size %in% "Medium"] <- "Medium clones"
  g[size %in% "Large"] <- "Large clones"
  g[size %in% top] <- paste(top, "clones")
  hyper_target <- paste(top, target_label)
  g[size %in% top & target] <- hyper_target
  clone_group <- factor(g, levels = c("Other", "Medium clones",
                                      "Large clones", paste(top, "clones"),
                                      hyper_target))

  is_exp <- size %in% expanded
  e <- ifelse(is_exp, "Expanded clone", "Other")
  et <- rep("Other", length(size))
  et[is_exp] <- "Expanded clone"
  et[target] <- target_label
  both <- paste("Expanded +", target_label)
  et[is_exp & target] <- both
  data.frame(
    clone_expansion_group = clone_group,
    expanded_clone = factor(e, levels = c("Other", "Expanded clone")),
    expanded_target_group = factor(et, levels = c("Other", "Expanded clone",
                                                  target_label, both)),
    stringsAsFactors = FALSE)
}

#' Reference lookup for single-cell annotation, by chain
#' @noRd
cell_reference_lookup <- function(sp, chain, report) {
  src <- sp$single_cell_reference
  if (identical(src, "trex")) {
    if (requireNamespace("Trex", quietly = TRUE)) {
      return(list(lookup = collapse_reference(trex_reference(chain)),
                  label = "Trex reference"))
    }
    report$p(paste0(
      "The Trex package is not installed; falling back to the configured ",
      "VDJdb / McPAS-TCR files."))
    src <- "all"
  }
  which <- if (identical(src, "all")) c("vdjdb", "mcpas") else src
  refs <- load_tcr_references(sp, chain, which)
  if (!length(refs)) return(NULL)
  list(lookup = collapse_reference(do.call(rbind, refs)),
       label = paste(names(refs), collapse = " + "))
}

#' Antigen specificity of single T cells
#'
#' Annotates the TRA and TRB CDR3 of every cell at each configured edit
#' distance, flags the cells of the target species, combines the two chains,
#' describes the V genes and CDR3s of the flagged cells, crosses
#' specificity with clonal expansion, and counts the expanded flagged cells
#' per condition. Figures on the transcriptome embedding are drawn when a
#' Seurat object is given.
#'
#' @param db AIRR data frame of the single-cell TCR (TRA and TRB rows).
#' @param abundance_col Abundance column (UMIs), used to order chains.
#' @param run_dir,cfg,report,prefix As for the other sections.
#' @param object Optional Seurat object whose cells carry the same
#'   `cell_id`s (from `section_seurat_attachment()`).
#' @return Per-cell annotation table, invisibly.
#' @export
section_cell_specificity <- function(db, abundance_col, run_dir, cfg, report,
                                     object = NULL, prefix = "15") {
  need(c("dplyr", "ggplot2", "stringdist"), "section_cell_specificity()")
  report$h1("Antigen specificity of single cells")
  sp <- specificity_config(cfg)
  report_annotation_caveat(report)

  cells <- cell_receptor_summary(db, abundance_col = abundance_col)
  chains <- c("TRA", "TRB")[c(any(!is.na(cells$TRA_cdr3)),
                              any(!is.na(cells$TRB_cdr3)))]
  distances <- sort(unique(sp$edit_distances))
  report$p(sprintf(paste0(
    "%d cells with a productive TRA or TRB chain. Annotation rule: exact ",
    "CDR3 identity, then propagation within the dataset to CDR3s within an ",
    "edit distance of %s of an exactly matched one."),
    nrow(cells), paste(distances, collapse = " or ")))

  # --- clone size per sample, on the configured clonotype -----------------
  locus <- cfg$repertoire$locus %||% "TRB"
  cl <- filter_locus_productive(db, locus)
  cl <- cl[cl$cell_id %in% cells$cell_id, , drop = FALSE]
  ab <- if (abundance_col %in% names(cl)) as.numeric(cl[[abundance_col]]) else 0
  cl <- cl[order(cl$cell_id, -ab), , drop = FALSE]
  cl$clonotype <- clonotype_key(cl, cfg$repertoire$clonotype_definition)
  cl <- cl[!duplicated(cl$cell_id), c("cell_id", "sample_id", "clonotype"),
           drop = FALSE]
  sizes <- dplyr::count(cl, .data$sample_id, .data$clonotype,
                        name = "clone_cells")
  cl$clone_cells <- sizes$clone_cells[match(
    paste(cl$sample_id, cl$clonotype), paste(sizes$sample_id, sizes$clonotype))]
  size_cat <- clone_size_category(cl$clone_cells, sp$clone_size_breaks)
  m <- match(cells$cell_id, cl$cell_id)
  cells$clonotype <- cl$clonotype[m]
  cells$clone_cells <- cl$clone_cells[m]
  cells$clone_size <- size_cat[m]
  report$p(paste0("Clone size categories (cells per clonotype within a ",
                  "sample): ", paste(attr(size_cat, "labels"),
                                     collapse = ", "), "."))

  # Annotation propagates from anchors in cell order; with a Seurat object
  # that order is the object's, as in Trex::annotateDB().
  if (!is.null(object)) {
    ord <- c(intersect(colnames(object), cells$cell_id),
             setdiff(cells$cell_id, colnames(object)))
    cells <- cells[match(ord, cells$cell_id), , drop = FALSE]
    rownames(cells) <- NULL
  }

  # --- annotation per chain and distance ----------------------------------
  species_tab <- list()
  for (ch in chains) {
    ref <- cell_reference_lookup(sp, ch, report)
    if (is.null(ref)) {
      report$p(sprintf("No reference available for %s; skipped.", ch))
      next
    }
    for (k in distances) {
      ann <- annotate_cells_by_cdr3(cells[[paste0(ch, "_cdr3")]],
                                    ref$lookup, max_distance = k)
      sfx <- paste0("_d", k)
      names(ann) <- paste0(ch, "_", names(ann), sfx)
      cells <- cbind(cells, ann)
      spc <- cells[[paste0(ch, "_Epitope.species", sfx)]]
      if (any(!is.na(spc))) {
        t <- as.data.frame(table(species = spc), stringsAsFactors = FALSE)
        t$chain <- ch
        t$edit_distance <- k
        species_tab[[length(species_tab) + 1L]] <- t
      }
      if (!is.null(sp$target_species)) {
        cells[[paste0(ch, "_target", sfx)]] <-
          grepl(sp$target_species, spc, ignore.case = TRUE)
      }
    }
    report$p(sprintf("%s annotated against: %s.", ch, ref$label))
  }
  if (!length(species_tab)) {
    report$p("No reference could be loaded; section stopped.")
    return(invisible(NULL))
  }
  species_tab <- dplyr::bind_rows(species_tab)
  names(species_tab)[names(species_tab) == "Freq"] <- "n_cells"
  save_table(species_tab, paste0(prefix, "_cell_species_counts"), run_dir,
             "annotated cells per antigen species, chain and edit distance")
  report$h2("Annotated cells per chain and edit distance")
  report$tbl(as.data.frame(dplyr::summarise(
    dplyr::group_by(species_tab, .data$chain, .data$edit_distance),
    n_cells_annotated = sum(.data$n_cells), .groups = "drop")))

  if (is.null(sp$target_species)) {
    report$p(paste0(
      "No target species configured (specificity$target_species), so no ",
      "cell is flagged and the expansion-by-specificity analysis is not run."))
    save_table(cells, paste0(prefix, "_cell_annotation"), run_dir,
               "per-cell receptors, clone size and annotation")
    return(invisible(cells))
  }

  # --- target flag: per chain, then combined at the largest distance -------
  lab <- sp$target_label
  kmax <- max(distances)
  flag <- function(ch, k) {
    col <- paste0(ch, "_target_d", k)
    if (col %in% names(cells)) cells[[col]] else rep(FALSE, nrow(cells))
  }
  group_col <- if ("condition" %in% names(cells) &&
                   !all(is.na(cells$condition))) "condition" else "sample_id"
  per_group <- list()
  for (ch in chains) for (k in distances) {
    f <- flag(ch, k)
    t <- dplyr::count(data.frame(group = cells[[group_col]], flag = f),
                      .data$group, .data$flag, name = "n_cells")
    t$chain <- ch
    t$edit_distance <- k
    per_group[[length(per_group) + 1L]] <- t
  }
  per_group <- dplyr::bind_rows(per_group)
  names(per_group)[names(per_group) == "group"] <- group_col
  save_table(per_group, paste0(prefix, "_target_cells_by_", group_col),
             run_dir, sprintf("cells flagged %s per chain, distance and %s",
                              lab, group_col))

  tra <- flag("TRA", kmax)
  trb <- flag("TRB", kmax)
  cells$target_any <- tra | trb
  cells$target_chain <- factor(
    ifelse(tra & trb, "Both_TRA_TRB", ifelse(tra, "TRA_only",
                                             ifelse(trb, "TRB_only",
                                                    "Other_NA"))),
    levels = c("TRA_only", "TRB_only", "Both_TRA_TRB", "Other_NA"))
  report$h2(sprintf("Cells flagged %s (TRA or TRB, edit distance %d)",
                    lab, kmax))
  report$tbl(as.data.frame(table(category = cells$target_chain)))

  # --- V genes and CDR3s of the flagged cells ------------------------------
  for (ch in chains) for (k in distances) {
    f <- flag(ch, k)
    if (!any(f)) next
    vu <- cell_v_usage(cells[[paste0(ch, "_v")]][f], sp$top_v_genes)
    save_table(vu, sprintf("%s_target_%sV_usage_d%d", prefix, ch, k),
               run_dir, sprintf("V genes of the %s cells flagged on %s, d=%d",
                                lab, ch, k))
    if (k == kmax && nrow(vu)) {
      p <- ggplot2::ggplot(vu, ggplot2::aes(stats::reorder(.data$v_gene,
                                                           .data$n), .data$n)) +
        ggplot2::geom_col(fill = "grey35") + ggplot2::coord_flip() +
        ggplot2::labs(x = "V gene", y = sprintf("%s cells", lab),
                      title = sprintf("Top %d %sV, %s cells (d = %d)",
                                      sp$top_v_genes, ch, lab, k)) +
        theme_pub(grid = "x")
      save_figure(p, sprintf("%s_target_%sV_usage", prefix, ch), run_dir, cfg,
                  width = 5, height = 6,
                  legend = sprintf(paste("The %d most frequent V genes of %s",
                                         "among cells flagged %s."),
                                   sp$top_v_genes, ch, lab))
    }
  }
  if ("TRB" %in% chains) {
    for (k in distances) {
      f <- flag("TRB", k)
      if (!any(f)) next
      tc <- top_cell_cdr3(cells$TRB_cdr3[f], sp$top_cdr3)
      save_table(tc, sprintf("%s_target_TRB_top_cdr3_d%d", prefix, k),
                 run_dir, sprintf("most frequent CDR3b of the %s cells, d=%d",
                                  lab, k))
      if (k == kmax && nrow(tc)) {
        p <- ggplot2::ggplot(tc, ggplot2::aes(stats::reorder(.data$cdr3,
                                                             .data$n),
                                              .data$n)) +
          ggplot2::geom_col(fill = "grey35") + ggplot2::coord_flip() +
          ggplot2::labs(x = "CDR3 beta", y = sprintf("%s cells", lab),
                        title = sprintf("Top %d TRB CDR3, %s cells (d = %d)",
                                        sp$top_cdr3, lab, k)) +
          theme_pub(grid = "x")
        save_figure(p, paste0(prefix, "_target_TRB_top_cdr3"), run_dir, cfg,
                    width = 6, height = 6,
                    legend = sprintf(paste("The %d most frequent TRB CDR3s",
                                           "among cells flagged %s."),
                                     sp$top_cdr3, lab))
      }
    }
  }

  # --- expansion x specificity ---------------------------------------------
  grp <- cell_specificity_groups(cells$clone_size, cells$target_any, lab,
                                 expanded = sp$expanded_categories)
  cells <- cbind(cells, grp)
  report$h2("Clonal expansion and specificity")
  report$tbl(as.data.frame(table(group = cells$clone_expansion_group)))
  report$tbl(as.data.frame(table(group = cells$expanded_target_group)))

  focus <- cells$expanded_clone == "Expanded clone" & cells$target_any
  lv <- if (identical(group_col, "condition") && length(cfg$cohort$conditions))
    cfg$cohort$conditions else sort(unique(as.character(cells[[group_col]])))
  totals <- table(factor(as.character(cells[[group_col]]), levels = lv))
  counts <- table(factor(as.character(cells[[group_col]][focus]), levels = lv))
  bars <- data.frame(group = lv, count = as.integer(counts),
                     total = as.integer(totals), stringsAsFactors = FALSE)
  bars$percentage <- ifelse(bars$total > 0, 100 * bars$count / bars$total,
                            NA_real_)
  names(bars)[1] <- group_col
  save_table(bars, paste0(prefix, "_expanded_target_by_", group_col), run_dir,
             sprintf("expanded %s cells per %s (count and percentage)",
                     lab, group_col))
  report$tbl(bars)

  cols <- if (identical(group_col, "condition")) condition_colours(cfg) else
    NULL
  for (what in c("count", "percentage")) {
    d <- bars
    d$.g <- factor(d[[group_col]], levels = rev(lv))
    p <- ggplot2::ggplot(d, ggplot2::aes(.data$.g, .data[[what]],
                                         fill = .data$.g)) +
      ggplot2::geom_col(width = 0.7) + ggplot2::coord_flip() +
      ggplot2::labs(x = NULL,
                    y = if (what == "count") "Number of cells" else
                      "Percentage of cells (%)",
                    title = sprintf("Expanded + %s cells (n = %d)", lab,
                                    sum(bars$count))) +
      theme_pub(grid = "x") + ggplot2::theme(legend.position = "none")
    if (length(cols)) {
      p <- p + ggplot2::scale_fill_manual(values = cols)
    } else {
      p <- p + ggplot2::scale_fill_manual(
        values = rep(palette_categorical(8), length.out = length(lv)))
    }
    save_figure(p, sprintf("%s_expanded_target_%s", prefix, what), run_dir,
                cfg, width = 5.5, height = 3.5,
                legend = sprintf(paste("Cells belonging to expanded",
                                       "clonotypes (%s) and flagged %s, per",
                                       "%s (%s)."),
                                 paste(sp$expanded_categories, collapse = ", "),
                                 lab, group_col, what))
  }

  save_table(cells, paste0(prefix, "_cell_annotation"), run_dir,
             "per-cell receptors, clone size, annotation and groups")

  if (!is.null(object)) {
    cell_specificity_umaps(object, cells, chains, distances, lab, run_dir,
                           cfg, report, prefix)
  }
  invisible(cells)
}

#' Embedding figures of the single-cell specificity section
#' @noRd
cell_specificity_umaps <- function(object, cells, chains, distances, lab,
                                   run_dir, cfg, report, prefix) {
  need(c("Seurat", "SeuratObject"), "cell_specificity_umaps()")
  if (!"umap" %in% names(object@reductions)) {
    report$p("The Seurat object has no 'umap' reduction; figures skipped.")
    return(invisible(NULL))
  }
  idx <- match(colnames(object), cells$cell_id)
  if (!any(!is.na(idx))) {
    report$p("No cell of the Seurat object matched; figures skipped.")
    return(invisible(NULL))
  }
  keep <- c(grep("_target_d[0-9]+$|_Epitope\\.", names(cells), value = TRUE),
            "target_any", "target_chain", "clone_expansion_group",
            "expanded_clone", "expanded_target_group", "clone_size")
  md <- cells[idx, intersect(keep, names(cells)), drop = FALSE]
  rownames(md) <- colnames(object)
  object <- Seurat::AddMetaData(object, md)
  saveRDS(object, file.path(run_dir, "seurat_with_specificity.rds"))

  flag_plot <- function(col, title) {
    v <- object[[col]][, 1]
    object$.flag <- factor(ifelse(!is.na(v) & v, lab, "Other / NA"),
                           levels = c("Other / NA", lab))
    Seurat::DimPlot(object, reduction = "umap", group.by = ".flag",
                    cols = c("grey85", "red"), pt.size = 1.5,
                    order = TRUE) +
      ggplot2::labs(title = title) + theme_pub(grid = "none")
  }
  panels <- list()
  for (ch in chains) for (k in distances) {
    col <- paste0(ch, "_target_d", k)
    if (col %in% colnames(object[[]])) {
      panels[[paste(ch, k)]] <- flag_plot(col, sprintf("%s, d = %d", ch, k))
    }
  }
  for (nm in names(panels)) {
    save_figure(panels[[nm]], sprintf("%s_umap_target_%s", prefix,
                                      gsub(" ", "_d", nm)),
                run_dir, cfg, width = 6, height = 5,
                legend = sprintf("Cells flagged %s on %s.", lab, nm))
  }

  p <- Seurat::DimPlot(object, reduction = "umap",
                       group.by = "clone_expansion_group",
                       cols = c("grey85", "#9ecae1", "#fdae6b", "#fb6a4a",
                                "darkred"),
                       pt.size = 1.4, order = TRUE) +
    ggplot2::labs(title = "Clonal expansion") + theme_pub(grid = "none")
  save_figure(p, paste0(prefix, "_umap_expansion_groups"), run_dir, cfg,
              width = 7, height = 5,
              legend = paste("Cells by clone size category, with the",
                             "largest category split by the target flag."))

  focus <- colnames(object)[!is.na(object$expanded_clone) &
                              object$expanded_clone == "Expanded clone" &
                              !is.na(object$target_any) & object$target_any]
  p <- Seurat::DimPlot(object, reduction = "umap", cells.highlight = focus,
                       cols = "grey85", cols.highlight = "red",
                       pt.size = 1.2, sizes.highlight = 2) +
    ggplot2::labs(title = sprintf("Expanded + %s (n = %d)", lab,
                                  length(focus))) +
    theme_pub(grid = "none") + ggplot2::theme(legend.position = "none")
  save_figure(p, paste0(prefix, "_umap_expanded_target"), run_dir, cfg,
              width = 6, height = 5,
              legend = sprintf("Cells of expanded clonotypes flagged %s.",
                               lab))

  p <- Seurat::DimPlot(object, reduction = "umap",
                       group.by = "expanded_target_group",
                       cols = c("grey85", "orange", "royalblue", "red"),
                       pt.size = 1.5, order = TRUE) +
    ggplot2::labs(title = sprintf("Clonal expansion and %s specificity",
                                  lab)) +
    theme_pub(grid = "none")
  save_figure(p, paste0(prefix, "_umap_expanded_by_target"), run_dir, cfg,
              width = 7, height = 5,
              legend = "Cells by expansion and by the target flag.")
  invisible(object)
}
