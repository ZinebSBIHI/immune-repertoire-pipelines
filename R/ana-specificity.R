# ---------------------------------------------------------------------------
# Annotation of repertoires against receptors of known specificity
#
# Matching rules, as in the original analyses:
#
#   * TCR, bulk: exact identity of the CDR3 amino-acid sequence (the junction,
#     conserved C and F/W included, as recorded by VDJdb and McPAS-TCR).
#   * TCR, single cell: exact identity, then propagation to receptors of the
#     same dataset within an edit distance of 1 or 2 of an exactly matched
#     CDR3 (see annotate_cells_by_cdr3()).
#   * BCR: same IGHV gene and CDRH3 identity of at least 0.80, normalised by
#     the longer sequence, against a reference of antibodies with a known
#     target (see match_antibody_reference()).
#
# An annotation is a statement about sequence identity with a published
# receptor, not a demonstration of specificity: the same CDR3 can pair with
# other chains, and identical CDR3s arise by convergent recombination in
# receptors of unrelated specificity. It is reported as a hypothesis to test.
#
# Bystander categories. Among clonotypes that changed significantly, a match
# to a pathogen unrelated to the antigen under study (CMV, EBV, influenza,
# SARS-CoV-2 by default) points to a change that is not antigen-driven by
# the exposure being studied, which is the reason to separate them. The
# categories are checked in a fixed order and the first match wins, so a
# CDR3 recorded against several species takes the first listed.
#
#   Annotation logic adapted from immunarch::dbAnnotate (ImmunoMind,
#   https://immunarch.com) and from Trex::annotateDB (Borcherding N et al.
#   Nat Immunol. 2024;25(9):1731-1741. doi:10.1038/s41590-024-01888-9.
#   PMID 39164479).
# ---------------------------------------------------------------------------

#' Specificity parameters from the configuration, with their defaults
#'
#' Relative file paths are resolved against `paths$reference_dir`, where
#' the other downloaded references live.
#' @noRd
specificity_config <- function(cfg) {
  sp <- cfg$specificity %||% list()
  ref_dir <- cfg$paths$reference_dir %||% "reference"
  resolve <- function(f) {
    if (is.null(f) || !length(f) || is.na(f) || !nzchar(f)) return(NULL)
    if (grepl("^(/|~|[A-Za-z]:)", f)) path.expand(f) else file.path(ref_dir, f)
  }
  byst <- sp$bystanders %||% as.list(default_bystanders())
  out <- list(
    vdjdb_file = resolve(sp$vdjdb_file),
    mcpas_file = resolve(sp$mcpas_file),
    antibody_file = resolve(sp$antibody_reference_file),
    vdjdb_host = sp$vdjdb_host_species %||% "HomoSapiens",
    mcpas_host = sp$mcpas_host_species %||% "Human",
    columns = sp$columns %||% list(),
    single_cell_reference = as.character(sp$single_cell_reference %||% "trex"),
    edit_distances = as.integer(unlist(sp$edit_distances %||% c(1L, 2L))),
    antibody_min_identity = as.numeric(sp$antibody_min_identity %||% 0.80),
    target_species = if (length(sp$target_species))
      as.character(sp$target_species)[1] else NULL,
    target_label = if (length(sp$target_label))
      as.character(sp$target_label)[1] else NULL,
    bystanders = vapply(byst, function(x) as.character(x)[1], character(1)),
    top_v_genes = as.integer(sp$top_v_genes %||% 30L),
    top_cdr3 = as.integer(sp$top_cdr3 %||% 20L),
    clone_size_breaks = unlist(sp$clone_size_breaks %||%
                                 list(Single = 1, Small = 2, Medium = 5,
                                      Large = 10)),
    expanded_categories = as.character(unlist(
      sp$expanded_categories %||% c("Medium", "Large", "Hyperexpanded")))
  )
  if (is.null(names(out$bystanders)) || any(!nzchar(names(out$bystanders)))) {
    stop("specificity$bystanders must be a named list of patterns, e.g. ",
         "{CMV: 'CMV|cytomegalo'}.", call. = FALSE)
  }
  if (anyNA(out$edit_distances) || any(out$edit_distances < 0L)) {
    stop("specificity$edit_distances must be non-negative integers.",
         call. = FALSE)
  }
  if (out$antibody_min_identity <= 0 || out$antibody_min_identity > 1) {
    stop("specificity$antibody_min_identity must lie in (0, 1].",
         call. = FALSE)
  }
  if (!out$single_cell_reference %in% c("trex", "vdjdb", "mcpas", "all")) {
    stop("specificity$single_cell_reference must be one of trex, vdjdb, ",
         "mcpas, all.", call. = FALSE)
  }
  if (!is.null(out$target_species)) {
    out$target_label <- out$target_label %||% out$target_species
    # A pathogen cannot be both the target of the study and a bystander:
    # its category is removed from the bystander list.
    clash <- grepl(out$target_species, names(out$bystanders),
                   ignore.case = TRUE)
    out$bystanders <- out$bystanders[!clash]
  }
  out
}

#' Default bystander categories
#'
#' Named regular expressions matched, case-insensitively and in this order,
#' against the antigen species of an annotated CDR3: the patterns of the
#' original analysis.
#'
#' @return Named character vector.
#' @export
default_bystanders <- function() {
  c(CMV = "CMV|cytomegalo",
    EBV = "EBV|Epstein",
    Influenza = "Influenza|flu",
    `SARS-CoV-2` = "SARS|CoV")
}

#' Annotate sequences by exact CDR3 identity
#'
#' @param x Data frame with a CDR3 column.
#' @param lookup Output of `collapse_reference()`.
#' @param cdr3_col Name of the CDR3 amino-acid column in `x`.
#' @return `x` with `annotated`, `specificities`, `epitopes`,
#'   `epitope_genes` and `sources` added; one row per input row.
#' @export
annotate_cdr3_exact <- function(x, lookup, cdr3_col = "junction_aa") {
  if (!cdr3_col %in% names(x)) {
    stop("Column '", cdr3_col, "' not found.", call. = FALSE)
  }
  idx <- match(as.character(x[[cdr3_col]]), lookup$cdr3)
  idx[is.na(x[[cdr3_col]])] <- NA_integer_
  x$annotated <- !is.na(idx)
  x$specificities <- lookup$specificities[idx]
  x$epitopes <- lookup$epitopes[idx]
  x$epitope_genes <- lookup$epitope_genes[idx]
  x$sources <- lookup$sources[idx]
  x
}

#' Classify annotated sequences into target, bystander, other, unannotated
#'
#' @param specificities Collapsed antigen species (NA when unannotated).
#' @param bystanders Named patterns, checked in order.
#' @param target_species Optional pattern for the antigen under study; its
#'   matches are classified first.
#' @param target_label Display label for the target category.
#' @return Factor with levels in display order.
#' @export
specificity_status <- function(specificities,
                               bystanders = default_bystanders(),
                               target_species = NULL,
                               target_label = target_species) {
  s <- as.character(specificities)
  annotated <- !is.na(s) & nzchar(s)
  status <- rep(NA_character_, length(s))
  levels_out <- character(0)

  if (!is.null(target_species)) {
    lab <- sprintf("Target (%s)", target_label)
    hit <- annotated & is.na(status) &
      grepl(target_species, s, ignore.case = TRUE)
    status[hit] <- lab
    levels_out <- c(levels_out, lab)
  }
  for (nm in names(bystanders)) {
    lab <- sprintf("Bystander (%s)", nm)
    hit <- annotated & is.na(status) &
      grepl(bystanders[[nm]], s, ignore.case = TRUE)
    status[hit] <- lab
    levels_out <- c(levels_out, lab)
  }
  other <- "Known specificity (other)"
  status[annotated & is.na(status)] <- other
  status[!annotated] <- "Unannotated"
  factor(status, levels = c(levels_out, other, "Unannotated"))
}

#' Colours for the specificity status categories
#' @noRd
status_colours <- function(levels) {
  fixed <- c(`Known specificity (other)` = "#E69F00",
             Unannotated = "#D55E00")
  pool <- c("#0072B2", "#009E73", "#CC79A7", "#56B4E9", "#000000",
            "#F0E442")
  other_lv <- setdiff(levels, names(fixed))
  tgt <- grepl("^Target", other_lv)
  cols <- character(0)
  cols[other_lv[tgt]] <- "#B2182B"
  byst <- other_lv[!tgt]
  cols[byst] <- rep(pool, length.out = length(byst))
  c(cols, fixed[intersect(names(fixed), levels)])[levels]
}

#' Map clonotype keys to their CDR3 amino-acid sequence
#'
#' A clonotype key built on the junction (the "nt", "strict" and "aa"
#' definitions) translates to exactly one CDR3; a key built on genes only
#' does not, and those keys map to NA.
#'
#' @param db AIRR data frame.
#' @param definition Clonotype definition used to build the keys.
#' @return Data frame: clonotype, junction_aa.
#' @export
clonotype_cdr3_map <- function(db, definition = "strict") {
  d <- db[!is.na(db$junction_aa), , drop = FALSE]
  d$clonotype <- clonotype_key(d, definition)
  m <- unique(data.frame(clonotype = d$clonotype,
                         junction_aa = as.character(d$junction_aa),
                         stringsAsFactors = FALSE))
  amb <- m$clonotype[duplicated(m$clonotype)]
  m$junction_aa[m$clonotype %in% amb] <- NA_character_
  unique(m)
}

#' Remove the conserved junction anchors from a CDR3
#'
#' Junction to CDR3 in the IMGT sense: the leading cysteine and the trailing
#' tryptophan or phenylalanine are removed, which is how CoV-AbDab records
#' CDRH3 sequences.
#' @noRd
strip_junction_anchors <- function(x) {
  sub("[WF]$", "", sub("^C", "", as.character(x)))
}

#' Match heavy chains to a reference of antibodies with a known target
#'
#' For every query sequence, the reference antibodies with the same IGHV
#' gene are scored by CDRH3 identity, one minus the Levenshtein distance
#' divided by the longer length, and the best one is kept when it reaches
#' `min_identity`. Requiring the same V gene first follows the observation
#' that public antibody classes are defined by V gene and CDRH3 together.
#'
#' @param query Data frame with V gene and junction columns, typically the
#'   members of public groups from `section_convergence()`.
#' @param reference Output of `read_antibody_reference()`.
#' @param min_identity Minimum identity of the best hit.
#' @param cdr3_col,v_col Columns of `query`; the junction is converted to a
#'   CDR3 by removing its conserved anchors.
#' @return Data frame of the query rows with a hit: the query columns,
#'   `cdr3`, `ab_hit`, `ab_cdr3`, `identity`, sorted by identity.
#' @export
match_antibody_reference <- function(query, reference, min_identity = 0.80,
                                     cdr3_col = "junction_aa",
                                     v_col = "v_gene") {
  need("stringdist", "match_antibody_reference()")
  for (cl in c(cdr3_col, v_col)) {
    if (!cl %in% names(query)) {
      stop("Column '", cl, "' not found in the query.", call. = FALSE)
    }
  }
  q <- query[!is.na(query[[cdr3_col]]), , drop = FALSE]
  q$cdr3 <- strip_junction_anchors(q[[cdr3_col]])
  q_v <- ighv_gene(q[[v_col]])

  hit <- rep(NA_character_, nrow(q))
  hit_cdr3 <- rep(NA_character_, nrow(q))
  ident <- rep(NA_real_, nrow(q))
  by_v <- split(seq_len(nrow(reference)), reference$ab_vgene)
  for (i in seq_len(nrow(q))) {
    cand <- by_v[[q_v[i]]]
    if (is.na(q_v[i]) || !length(cand)) next
    s <- 1 - stringdist::stringdist(q$cdr3[i], reference$ab_cdr3[cand],
                                    method = "lv") /
      pmax(nchar(q$cdr3[i]), nchar(reference$ab_cdr3[cand]))
    b <- which.max(s)
    hit[i] <- reference$ab_name[cand[b]]
    hit_cdr3[i] <- reference$ab_cdr3[cand[b]]
    ident[i] <- s[b]
  }
  q$ab_hit <- hit
  q$ab_cdr3 <- hit_cdr3
  q$identity <- ident
  out <- q[!is.na(q$identity) & q$identity >= min_identity, , drop = FALSE]
  out <- out[order(-out$identity), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Load the configured TCR references for one chain
#' @noRd
load_tcr_references <- function(sp, chain, which = c("vdjdb", "mcpas")) {
  refs <- list()
  if ("vdjdb" %in% which && !is.null(sp$vdjdb_file) &&
      file.exists(sp$vdjdb_file)) {
    refs$VDJdb <- read_vdjdb(sp$vdjdb_file, chain = chain,
                             host_species = sp$vdjdb_host,
                             columns = sp$columns$vdjdb)
  }
  if ("mcpas" %in% which && !is.null(sp$mcpas_file) &&
      file.exists(sp$mcpas_file)) {
    refs$McPAS <- read_mcpas(sp$mcpas_file, chain = chain,
                             host_species = sp$mcpas_host,
                             columns = sp$columns$mcpas)
  }
  refs
}

#' Report text on what an annotation can and cannot say
#' @noRd
report_annotation_caveat <- function(report) {
  report$p(paste0(
    "Public receptor databases are dominated by a few intensively studied ",
    "pathogens (CMV, EBV, influenza, SARS-CoV-2). An unannotated clonotype ",
    "is therefore not evidence of specificity for the antigen under study, ",
    "and an annotated one is a sequence match to a published receptor, a ",
    "hypothesis to test rather than a demonstration."))
  report$p("")
}

#' Antigen specificity of a bulk TCR repertoire
#'
#' Annotates the observed CDR3s against VDJdb and McPAS-TCR by exact
#' identity, describes the antigen species found per subject and timepoint,
#' and, when differential-abundance calls are given, classifies the
#' significantly expanded clonotypes as target, bystander, other annotated
#' or unannotated.
#'
#' @param db AIRR data frame (bulk), with `junction_aa`.
#' @param abundance_col Abundance column.
#' @param locus Analysed locus; only TRA and TRB are annotated.
#' @param calls Optional output of `test_clonal_abundance()`.
#' @param run_dir,cfg,report,prefix As for the other sections.
#' @return List with `per_sample`, `composition` and `status`, invisibly.
#' @export
section_specificity <- function(db, abundance_col, locus, run_dir, cfg,
                                report, calls = NULL, prefix = "15") {
  need(c("dplyr", "ggplot2"), "section_specificity()")
  report$h1("Antigen specificity against reference databases")
  sp <- specificity_config(cfg)
  chain <- as.character(locus)[1]
  if (!chain %in% c("TRA", "TRB")) {
    report$p(sprintf(paste0("The TCR reference databases do not apply to ",
                            "locus %s; section skipped."), chain))
    return(invisible(NULL))
  }
  refs <- load_tcr_references(sp, chain)
  if (!length(refs)) {
    report$p(paste0(
      "No reference database found (specificity$vdjdb_file, ",
      "specificity$mcpas_file); section skipped. See data/README.md for ",
      "where to download them."))
    return(invisible(NULL))
  }
  report_annotation_caveat(report)

  # --- what the references contain --------------------------------------
  overview <- reference_species_overview(dplyr::bind_rows(refs))
  save_table(overview, paste0(prefix, "_reference_species"), run_dir,
             "records and distinct CDR3s per antigen species, per database")
  report$h2("Antigen species in the references (largest ten per database)")
  report$tbl(as.data.frame(dplyr::ungroup(dplyr::slice_head(
    dplyr::group_by(overview, .data$source), n = 10))))

  obs <- filter_locus_productive(db, locus)
  obs <- obs[!is.na(obs$junction_aa) & nzchar(obs$junction_aa), ,
             drop = FALSE]
  obs$.abund <- as.numeric(obs[[abundance_col]])
  has_tp <- "timepoint" %in% names(obs) && !all(is.na(obs$timepoint))
  grp <- if (has_tp) c("subject_id", "timepoint") else
    c("subject_id", "sample_id")

  # --- exact annotation per database, per sample -------------------------
  per_sample <- list()
  keys <- unique(c("sample_id", grp))
  for (nm in names(refs)) {
    lk <- collapse_reference(refs[[nm]])
    hit <- obs$junction_aa %in% lk$cdr3
    tmp <- data.frame(obs[, keys, drop = FALSE], .hit = hit,
                      .abund = obs$.abund, .cdr3 = obs$junction_aa,
                      stringsAsFactors = FALSE)
    s <- dplyr::summarise(
      dplyr::group_by(tmp, dplyr::across(dplyr::all_of(keys))),
      n_cdr3 = dplyr::n_distinct(.data$.cdr3),
      n_cdr3_annotated = dplyr::n_distinct(.data$.cdr3[.data$.hit]),
      abundance_fraction_annotated = sum(.data$.abund[.data$.hit]) /
        sum(.data$.abund),
      .groups = "drop")
    s$database <- nm
    per_sample[[nm]] <- s
    report$p(sprintf("%-6s distinct CDR3s annotated across samples: %d",
                     nm, length(unique(obs$junction_aa[hit]))))
  }
  per_sample <- dplyr::bind_rows(per_sample)
  save_table(per_sample, paste0(prefix, "_annotation_per_sample"), run_dir,
             "distinct CDR3s annotated by exact identity, per database")
  report$p("")

  # --- detailed join and species composition (first database) ------------
  main <- names(refs)[1]
  ref <- refs[[main]]
  detail <- merge(
    obs[, intersect(c("sample_id", "subject_id", "timepoint", "condition",
                      "junction_aa", "v_call", "j_call", abundance_col),
                    names(obs)), drop = FALSE],
    unique(ref[, c("cdr3", "epitope", "epitope_gene", "epitope_species",
                   "mhc")]),
    by.x = "junction_aa", by.y = "cdr3")
  save_table(detail, paste0(prefix, "_annotated_rearrangements_", main),
             run_dir, paste("observed rearrangements with their", main,
                            "epitope, antigen gene, species and MHC"))

  # Each (CDR3, species) pair contributes the clonotype's abundance, as in
  # the original figure: a CDR3 recorded against k species is counted under
  # each of them, so fractions are of annotations rather than of molecules.
  pairs <- unique(ref[, c("cdr3", "epitope_species")])
  ann <- merge(obs[, c(grp, "junction_aa", ".abund"), drop = FALSE], pairs,
               by.x = "junction_aa", by.y = "cdr3")
  composition <- NULL
  if (nrow(ann)) {
    composition <- dplyr::mutate(
      dplyr::group_by(
        dplyr::summarise(
          dplyr::group_by(ann, dplyr::across(dplyr::all_of(grp)),
                          .data$epitope_species),
          abundance = sum(.data$.abund), .groups = "drop"),
        dplyr::across(dplyr::all_of(grp))),
      fraction = .data$abundance / sum(.data$abundance))
    composition <- dplyr::ungroup(composition)
    save_table(composition, paste0(prefix, "_species_composition"), run_dir,
               paste("abundance-weighted antigen species among annotated",
                     "clonotypes (", main, ")"))
    n_sp <- length(unique(composition$epitope_species))
    p <- ggplot2::ggplot(composition,
                         ggplot2::aes(.data[[grp[2]]], .data$fraction,
                                      fill = .data$epitope_species)) +
      ggplot2::geom_col(width = 0.8) +
      ggplot2::facet_wrap(~ subject_id, nrow = 1) +
      ggplot2::scale_fill_manual(
        values = grDevices::hcl.colors(max(n_sp, 2), "viridis")[seq_len(n_sp)]) +
      ggplot2::labs(x = NULL, y = "Fraction of annotated abundance",
                    fill = "Antigen species",
                    title = sprintf("Annotated specificities (%s)", main)) +
      theme_pub() +
      ggplot2::theme(legend.position = "bottom",
                     axis.text.x = ggplot2::element_text(angle = 90,
                                                         hjust = 1))
    save_figure(p, paste0(prefix, "_species_composition"), run_dir, cfg,
                width = max(6, 2 + 2 * length(unique(composition$subject_id))),
                height = 5,
                legend = paste("Antigen species of the clonotypes annotated by",
                               "exact CDR3 identity, weighted by abundance,",
                               "per subject and", grp[2], "."))
  } else {
    report$p(sprintf("No observed CDR3 matched %s exactly.", main))
  }

  # --- significantly expanded clonotypes ----------------------------------
  status <- NULL
  if (!is.null(calls) && nrow(calls)) {
    status <- specificity_of_expanded(
      calls, db, refs, sp, run_dir, cfg, report, prefix)
  }
  invisible(list(per_sample = per_sample, composition = composition,
                 status = status))
}

#' Status of significantly expanded clonotypes (part of section_specificity)
#' @noRd
specificity_of_expanded <- function(calls, db, refs, sp, run_dir, cfg,
                                    report, prefix) {
  report$h2("Status of significantly expanded clonotypes")
  da <- cfg$differential_abundance
  map <- clonotype_cdr3_map(db, cfg$repertoire$clonotype_definition)
  if (all(is.na(map$junction_aa))) {
    report$p(paste0(
      "The clonotype definition does not identify a CDR3 sequence, so the ",
      "expanded clonotypes cannot be annotated."))
    return(invisible(NULL))
  }
  # Both databases are merged into one dictionary, so a CDR3 recorded in
  # either counts as annotated.
  lk <- collapse_reference(dplyr::bind_rows(refs))
  calls$junction_aa <- map$junction_aa[match(calls$clonotype, map$clonotype)]
  calls <- annotate_cdr3_exact(calls, lk)
  calls$specificity_status <- specificity_status(
    calls$specificities, sp$bystanders, sp$target_species, sp$target_label)

  sig <- calls[calls$sig_expansion, , drop = FALSE]
  report$p(sprintf(paste0(
    "Significantly expanded: fold change >= %g and FDR < %g against each ",
    "subject's baseline. %d clonotype-timepoint pairs; %d annotated, %d ",
    "unannotated (%.1f%%)."),
    da$fold_change_threshold, da$fdr_threshold, nrow(sig),
    sum(sig$annotated), sum(!sig$annotated),
    if (nrow(sig)) 100 * mean(!sig$annotated) else 0))
  save_table(sig, paste0(prefix, "_expanded_annotated"), run_dir,
             "significantly expanded clonotypes with their annotation status")

  if (nrow(sig)) {
    counts <- dplyr::count(sig, .data$subject_id, .data$timepoint,
                           .data$specificity_status, name = "n_clonotypes")
    save_table(counts, paste0(prefix, "_expanded_status_counts"), run_dir,
               "significantly expanded clonotypes per status")
    report$tbl(as.data.frame(dplyr::count(sig, .data$specificity_status,
                                          name = "n_clonotypes")))
    lv <- levels(sig$specificity_status)
    p <- ggplot2::ggplot(counts, ggplot2::aes(.data$timepoint,
                                              .data$n_clonotypes,
                                              fill = .data$specificity_status)) +
      ggplot2::geom_col(width = 0.8) +
      ggplot2::facet_wrap(~ subject_id, nrow = 1) +
      ggplot2::scale_fill_manual(values = status_colours(lv), drop = FALSE) +
      ggplot2::labs(x = NULL, y = "Significantly expanded clonotypes",
                    fill = NULL,
                    title = "Annotation status of expanded clonotypes") +
      theme_pub() +
      ggplot2::theme(legend.position = "bottom",
                     axis.text.x = ggplot2::element_text(angle = 90,
                                                         hjust = 1))
    save_figure(p, paste0(prefix, "_expanded_status"), run_dir, cfg,
                width = max(6, 2 + 2 * length(unique(counts$subject_id))),
                height = 5,
                legend = paste("Significantly expanded clonotypes by",
                               "annotation status: bystander species,",
                               "other known specificity, or no match in the",
                               "references."))
  }

  # Volcano with the thresholds drawn and the annotation of the expanded
  # clonotypes shown; capped at 50 as in the differential-abundance volcano.
  v <- calls[is.finite(calls$log2_fold_change) & !is.na(calls$fdr), ,
             drop = FALSE]
  if (nrow(v)) {
    v$neglog_fdr <- pmin(-log10(pmax(v$fdr, 1e-300)), 50)
    v$group <- ifelse(v$sig_expansion & v$annotated,
                      "Significant, annotated",
                      ifelse(v$sig_expansion, "Significant, unannotated",
                             "Not significant"))
    v$group <- factor(v$group, levels = c("Not significant",
                                          "Significant, annotated",
                                          "Significant, unannotated"))
    v <- thin_points(v, keep = "sig_expansion", max_background = 15000)
    v <- v[order(v$group), , drop = FALSE]
    p <- ggplot2::ggplot(v, ggplot2::aes(.data$log2_fold_change,
                                         .data$neglog_fdr,
                                         colour = .data$group)) +
      ggplot2::geom_point(size = 0.8, alpha = 0.6) +
      ggplot2::geom_vline(xintercept = log2(da$fold_change_threshold),
                          linetype = "dashed", linewidth = 0.3) +
      ggplot2::geom_hline(yintercept = -log10(da$fdr_threshold),
                          linetype = "dashed", linewidth = 0.3) +
      ggplot2::facet_grid(subject_id ~ timepoint) +
      ggplot2::scale_colour_manual(
        values = c(`Not significant` = "grey75",
                   `Significant, annotated` = "#0072B2",
                   `Significant, unannotated` = "#D55E00"), drop = FALSE) +
      ggplot2::labs(x = expression(log[2] ~ "fold change vs baseline"),
                    y = expression(-log[10] ~ "FDR (capped at 50)"),
                    colour = NULL,
                    title = "Expanded clonotypes and database annotation") +
      theme_pub(grid = "both") +
      ggplot2::theme(legend.position = "bottom")
    save_figure(p, paste0(prefix, "_volcano_annotation"), run_dir, cfg,
                width = 9, height = 7,
                legend = paste("Fold change against each subject's baseline",
                               "versus FDR, coloured by significance and",
                               "database annotation; dashed lines mark the",
                               "configured thresholds."))
  }
  invisible(sig)
}

#' Public heavy chains matched to a reference of known antibodies
#'
#' @param groups Members of public groups, as returned in
#'   `section_convergence()$groups$groups` (columns subject_id, v_gene,
#'   junction_aa, public_group).
#' @param label Label of the sequence set, used in file names.
#' @param run_dir,cfg,report,prefix As for the other sections.
#' @return The hits, invisibly.
#' @export
section_antibody_lookup <- function(groups, run_dir, cfg, report,
                                    label = "public", prefix = "07c") {
  report$h1("Public heavy chains against a reference of known antibodies")
  sp <- specificity_config(cfg)
  if (is.null(sp$antibody_file) || !file.exists(sp$antibody_file)) {
    report$p(paste0(
      "No antibody reference configured (specificity$",
      "antibody_reference_file); lookup skipped."))
    return(invisible(NULL))
  }
  if (is.null(groups) || !nrow(groups)) {
    report$p("No public group to look up.")
    return(invisible(NULL))
  }
  ref <- read_antibody_reference(sp$antibody_file,
                                 columns = sp$columns$antibody)
  q <- unique(groups[, intersect(c("public_group", "v_gene", "subject_id",
                                   "junction_aa"), names(groups)),
                     drop = FALSE])
  hits <- match_antibody_reference(q, ref,
                                   min_identity = sp$antibody_min_identity)
  report$p(sprintf(paste0(
    "%d reference antibodies read. Rule: same IGHV gene, CDRH3 identity >= ",
    "%.2f normalised by the longer sequence (the junction anchors C and W/F ",
    "are removed from the query). Query sequences: %d; with a hit: %d."),
    nrow(ref), sp$antibody_min_identity, nrow(q), nrow(hits)))
  if (nrow(hits)) {
    save_table(hits, paste0(prefix, "_antibody_hits_", label), run_dir,
               "public heavy chains with a close reference antibody")
    report$tbl(as.data.frame(utils::head(hits, 15)))
  }
  invisible(hits)
}

#' Entries of a reference for one antigen species, and their composition
#'
#' The species is selected by exact equality, as in the original export,
#' so a CDR3 recorded against several species (joined with ";") is not
#' included; the per-cell flag of `section_cell_specificity()` uses a
#' pattern instead and does include them.
#'
#' @param reference Standard reference table.
#' @param species Antigen species label.
#' @return The matching rows.
#' @export
reference_species_subset <- function(reference, species) {
  reference[!is.na(reference$epitope_species) &
              reference$epitope_species == species, , drop = FALSE]
}

#' Composition of the reference entries of the target species
#'
#' @param reference Standard reference table.
#' @param species Antigen species label (exact).
#' @param run_dir,cfg,report,prefix As for the other sections.
#' @return The subset, invisibly.
#' @export
section_reference_composition <- function(reference, species, run_dir, cfg,
                                          report, prefix = "15b") {
  need(c("dplyr", "ggplot2"), "section_reference_composition()")
  sub <- reference_species_subset(reference, species)
  report$h2(sprintf("Reference entries for %s", species))
  report$p(sprintf("%d entries, %d distinct epitopes.", nrow(sub),
                   length(unique(stats::na.omit(sub$epitope)))))
  if (!nrow(sub)) return(invisible(sub))
  save_table(sub, paste0(prefix, "_reference_", tracking_stem(species)),
             run_dir, paste("reference entries recorded against", species))

  top_count <- function(v, n) {
    v <- v[!is.na(v)]
    if (!length(v)) return(NULL)
    t <- sort(table(v), decreasing = TRUE)
    t <- t[seq_len(min(n, length(t)))]
    data.frame(label = names(t), n = as.integer(t), stringsAsFactors = FALSE)
  }
  bar <- function(d, title, ylab) {
    ggplot2::ggplot(d, ggplot2::aes(stats::reorder(.data$label, .data$n),
                                    .data$n, fill = .data$n)) +
      ggplot2::geom_col() + ggplot2::coord_flip() +
      ggplot2::scale_fill_gradient(low = "#9FE1CB", high = "#0F6E56",
                                   guide = "none") +
      ggplot2::labs(x = NULL, y = ylab, title = title,
                    subtitle = sprintf("%s entries of the reference", species)) +
      theme_pub(grid = "x")
  }
  v <- top_count(sub$v_gene, 20)
  if (!is.null(v)) {
    save_figure(bar(v, "V gene usage", "CDR3 entries"),
                paste0(prefix, "_reference_v_genes"), run_dir, cfg,
                width = 5, height = 5,
                legend = paste("The 20 most frequent V genes among the",
                               "reference entries of the species."))
  }
  tg <- top_count(sub$epitope_gene, 20)
  if (!is.null(tg)) {
    save_figure(bar(tg, "Antigen targets", "CDR3 entries"),
                paste0(prefix, "_reference_targets"), run_dir, cfg,
                width = 5, height = 5,
                legend = paste("The 20 most frequent antigen targets among",
                               "the reference entries of the species."))
    tg15 <- utils::head(tg, 15)
    p <- ggplot2::ggplot(tg15, ggplot2::aes(stats::reorder(.data$label,
                                                           .data$n),
                                            .data$n)) +
      ggplot2::geom_segment(ggplot2::aes(xend = .data$label, y = 0,
                                         yend = .data$n),
                            colour = "#BA7517") +
      ggplot2::geom_point(size = 3, colour = "#BA7517") +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "Unique CDR3", title = "Targets recognised",
                    subtitle = sprintf("Top 15 - %s", species)) +
      theme_pub(grid = "x")
    save_figure(p, paste0(prefix, "_reference_targets_top15"), run_dir, cfg,
                width = 5, height = 4.5,
                legend = "The 15 most frequent antigen targets (lollipop).")
  }
  invisible(sub)
}
