# ---------------------------------------------------------------------------
# Shared analysis sections
#
# Each function here is one section of a repertoire analysis: it computes,
# writes its tables and figures into the run directory, appends its numbers
# to the report, and returns what the next section needs.
#
# WHY THESE EXIST. The bulk and the single-cell pipeline ask the same
# questions of their data — how clonal is the repertoire, how diverse at
# common depth, which V genes, how much overlap, is there convergence
# between subjects. Written twice, the two pipelines drifted apart: one grew
# a convergence section the other never got, the other grew a diversity
# section the first was missing, and neither reader could tell which
# analyses a given pipeline actually offered. Calling the same function from
# both makes "the same analyses everywhere" a property of the code rather
# than something to remember.
#
# THE UNIT OF OBSERVATION differs between them and that is the one thing
# these functions must stay honest about. In bulk data the unit is the
# molecule (UMI); in single-cell data it is the cell. Every section takes a
# `unit` argument used in axis labels, legends and report text, so a figure
# never silently implies cells when it counted molecules.
#
# CONTRACT. `clono` is the output of collapse_clonotypes(): one row per
# clonotype per sample, with sample_id, subject_id, clonotype, abundance and
# frequency. `report` is the writer list from open_report().
# ---------------------------------------------------------------------------

#' Composition and depth QC
#'
#' @param db AIRR data frame before collapsing.
#' @param locus Locus the analysis will use.
#' @return The QC table, invisibly.
section_composition <- function(db, locus, run_dir, cfg, report,
                                prefix = "01") {
  need("dplyr", "section_composition()")
  report$h1("Input and composition")

  has_locus <- "locus" %in% names(db)
  qc <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(db, .data$sample_id, .data$subject_id),
    n_rearrangements = dplyr::n(),
    pct_productive = 100 * mean(parse_productive(.data$productive), na.rm = TRUE),
    pct_target_locus = if (has_locus)
      100 * mean(.data$locus %in% locus, na.rm = TRUE) else NA_real_,
    .groups = "drop"))

  report$h2("Per sample")
  report$tbl(as.data.frame(qc))
  save_table(qc, paste0(prefix, "_composition_qc"), run_dir,
             "productive and locus fractions per sample")
  report$p("")
  report$p(paste0(
    "A sample whose productive fraction or locus composition is out of line ",
    "with the rest is a library problem. It is better caught here than ",
    "inferred later from an odd diversity value."))
  invisible(qc)
}

#' Clonotype definition sensitivity
#'
#' The same sample yields a different clonotype count under each of the four
#' definitions in common use, and the spread is often several-fold. Printing
#' all four makes the configured choice visible to a reader instead of
#' implicit, and is worth carrying as a supplementary table.
section_clonotype_definitions <- function(db, abundance_col, locus,
                                          run_dir, cfg, report,
                                          prefix = "02") {
  need(c("dplyr", "tidyr"), "section_clonotype_definitions()")
  report$h1("Clonotype definition")

  defs <- c("aa", "nt", "gene", "strict")
  counts <- dplyr::bind_rows(lapply(defs, function(d) {
    cl <- collapse_clonotypes(db, abundance_col = abundance_col,
                              definition = d, min_abundance = 1,
                              locus = locus)
    out <- dplyr::ungroup(dplyr::summarise(
      dplyr::group_by(cl, .data$sample_id),
      n_clonotypes = dplyr::n(), .groups = "drop"))
    out$definition <- d
    out
  }))
  tab <- tidyr::pivot_wider(counts, names_from = "definition",
                            values_from = "n_clonotypes")

  report$p(paste0(
    "Clonotype counts under each definition, at no abundance floor. The ",
    "configured definition ('", cfg$repertoire$clonotype_definition,
    "') is used everywhere below."))
  report$p("")
  report$tbl(as.data.frame(tab))
  save_table(tab, paste0(prefix, "_clonotype_definition_sensitivity"), run_dir,
             "clonotype count under each of the four definitions")
  invisible(tab)
}

#' Repertoire structure: depth, richness, clonality, size classes
section_structure <- function(clono, sample_meta, run_dir, cfg, report,
                              unit = "molecules", prefix = "03") {
  need("dplyr", "section_structure()")
  report$h1("Repertoire structure")

  rep_summary <- repertoire_summary(clono)
  join_keys <- intersect(c("sample_id", "subject_id"), names(sample_meta))
  rep_summary <- dplyr::left_join(rep_summary, sample_meta, by = join_keys)

  report$h2(sprintf("Per sample (depth in %s)", unit))
  report$tbl(as.data.frame(rep_summary))
  save_table(rep_summary, paste0(prefix, "_repertoire_summary"), run_dir,
             sprintf("depth, richness and clonality per sample (%s)", unit))

  dr <- range(rep_summary$depth, na.rm = TRUE)
  report$p("")
  report$p(sprintf(
    paste0("Depth spans %s to %s %s, a %.1f-fold range. Richness, clone ",
           "counts and clonal concentration are all monotone in depth, so ",
           "those comparisons are made at common depth in the diversity ",
           "section."),
    format(dr[1], big.mark = " "), format(dr[2], big.mark = " "), unit,
    dr[2] / max(1, dr[1])))

  homeo <- clonal_homeostasis(clono)
  save_table(homeo, paste0(prefix, "_clonal_homeostasis"), run_dir,
             "repertoire fraction by clonotype size class")

  p <- ggplot2::ggplot(homeo, ggplot2::aes(.data$sample_id,
                                           .data$repertoire_fraction,
                                           fill = .data$size_class)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::scale_fill_manual(
      values = palette_sequential(length(levels(homeo$size_class)),
                                  direction = -1)) +
    ggplot2::scale_y_continuous(labels = scales::label_percent()) +
    ggplot2::labs(x = NULL, y = "Repertoire", fill = "Clonotype frequency",
                  title = "Clonal space homeostasis") +
    theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_figure(p, paste0(prefix, "_clonal_homeostasis"), run_dir, cfg,
              width = 8, height = 4,
              legend = sprintf(paste0("Fraction of the repertoire held by ",
                                      "clonotypes in each frequency class, ",
                                      "per sample. Depth is counted in %s."),
                               unit))
  invisible(rep_summary)
}

#' Diversity at native and at common depth
section_diversity <- function(clono, sample_meta, run_dir, cfg, report,
                              within = "subject_id", unit = "molecules",
                              prefix = "04") {
  need("dplyr", "section_diversity()")
  report$h1("Diversity")

  if (!isTRUE(cfg$repertoire$rarefaction$enabled)) {
    report$p(paste0(
      "Rarefaction disabled in the configuration. Any comparison of ",
      "richness, clone counts or clonal concentration between samples of ",
      "unequal depth is confounded by depth."))
    return(invisible(NULL))
  }

  rar <- rarefy_to_common_depth(
    clono, within = within,
    n_bootstrap = cfg$repertoire$rarefaction$n_bootstrap)

  join_keys <- intersect(names(rar), names(sample_meta))
  join_keys <- intersect(join_keys, c("sample_id", "subject_id"))
  rar <- dplyr::left_join(rar, sample_meta, by = join_keys)

  save_table(rar, paste0(prefix, "_rarefied_metrics"), run_dir,
             sprintf("metrics at common depth (%s), with 95%% intervals", unit))

  report$h2(sprintf("Common depth (%s)", unit))
  depth_cols <- intersect(c(within, "target_depth"), names(rar))
  report$tbl(as.data.frame(unique(rar[, depth_cols, drop = FALSE])))
  report$p("")
  report$p(sprintf(
    paste0("The common depth is taken within %s rather than across the ",
           "cohort: degrading every subject to the depth of the shallowest ",
           "sample elsewhere in the cohort would discard information for ",
           "nothing."),
    if (is.null(within)) "the whole cohort" else within))

  report$h2("Rarefied metrics")
  keep <- intersect(c("subject_id", "sample_id", "statistic", "mean",
                      "q025", "q975"), names(rar))
  report$tbl(as.data.frame(rar[, keep, drop = FALSE]))

  p <- ggplot2::ggplot(rar, ggplot2::aes(.data$sample_id, .data$mean)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$q025, ymax = .data$q975),
                           width = 0.15, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_point(ggplot2::aes(colour = .data$subject_id), size = 2.2) +
    ggplot2::facet_wrap(~statistic, scales = "free_y") +
    ggplot2::labs(x = NULL, y = sprintf("Value at common depth (%s)", unit),
                  colour = "Subject",
                  title = "Repertoire metrics at common depth",
                  subtitle = sprintf("%d resamples; bars are 95%% intervals",
                                     cfg$repertoire$rarefaction$n_bootstrap)) +
    theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_figure(p, paste0(prefix, "_rarefied_metrics"), run_dir, cfg,
              width = 9, height = 5,
              legend = sprintf(paste0("Repertoire metrics after resampling ",
                                      "every sample to a common depth in %s. ",
                                      "Points are the resampling mean, bars ",
                                      "the 2.5th and 97.5th percentiles."),
                               unit))
  invisible(rar)
}

#' CDR3 length and V gene usage, at clonotype and at observation level
section_gene_usage <- function(db, abundance_col, locus, sample_meta,
                               run_dir, cfg, report, unit = "molecules",
                               prefix = "05") {
  need("dplyr", "section_gene_usage()")
  report$h1("CDR3 length and gene usage")

  db_locus <- filter_locus_productive(db, locus)

  len <- unique(db_locus[, intersect(c("sample_id", "junction_aa"),
                                     names(db_locus)), drop = FALSE])
  len$cdr3_length <- nchar(as.character(len$junction_aa))
  p_len <- ggplot2::ggplot(len, ggplot2::aes(.data$cdr3_length)) +
    ggplot2::geom_histogram(binwidth = 1, fill = "grey70", colour = "white",
                            linewidth = 0.2) +
    ggplot2::facet_wrap(~sample_id, scales = "free_y") +
    ggplot2::labs(x = "CDR3 length (amino acids)", y = "Clonotypes",
                  title = sprintf("%s CDR3 length", paste(locus, collapse = "/"))) +
    theme_pub()
  save_figure(p_len, paste0(prefix, "a_cdr3_length"), run_dir, cfg,
              width = 9, height = 5,
              legend = "CDR3 length distribution per sample, one count per distinct clonotype.")

  usage <- dplyr::bind_rows(
    gene_usage(db_locus, "v_call", level = "clonotype",
               definition = cfg$repertoire$clonotype_definition),
    gene_usage(db_locus, "v_call", level = "molecule",
               abundance_col = abundance_col,
               definition = cfg$repertoire$clonotype_definition))
  usage$level[usage$level == "molecule"] <- unit

  save_table(usage, paste0(prefix, "_v_gene_usage"), run_dir,
             sprintf("V gene usage at clonotype and %s level", unit))

  report$p(paste0(
    "Usage is reported at clonotype level and at ", unit, " level. At ", unit,
    " level a single large clone dominates its sample's profile, so an ",
    "apparent gene bias can be one expanded lineage. Where the two levels ",
    "disagree, the clone is the explanation and the finding is about that ",
    "clone, not about the repertoire."))
  report$p("")

  top_v <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(usage[usage$level == "clonotype", ], .data$gene),
    m = mean(.data$frequency), .groups = "drop"))
  top_v <- utils::head(top_v$gene[order(-top_v$m)], 15)

  p_use <- ggplot2::ggplot(usage[usage$gene %in% top_v, ],
                           ggplot2::aes(stats::reorder(.data$gene,
                                                       .data$frequency),
                                        .data$frequency,
                                        colour = .data$level)) +
    ggplot2::stat_summary(fun = mean, geom = "point", size = 1.8,
                          position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar",
                          width = 0.2, linewidth = 0.3,
                          position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::coord_flip() +
    ggplot2::scale_colour_manual(values = palette_categorical(2)) +
    ggplot2::scale_y_continuous(labels = scales::label_percent()) +
    ggplot2::labs(x = NULL, y = "Usage", colour = "Level",
                  title = "V gene usage, 15 most used genes",
                  subtitle = "Clonotype level separates repertoire-wide bias from clonal expansion") +
    theme_pub(grid = "x")
  save_figure(p_use, paste0(prefix, "b_v_gene_usage"), run_dir, cfg,
              width = 6.5, height = 5,
              legend = sprintf(paste0("Mean V gene usage across samples at ",
                                      "clonotype and %s level, with standard ",
                                      "error."), unit))
  invisible(usage)
}

#' Pairwise repertoire overlap
section_overlap <- function(clono, run_dir, cfg, report,
                            within = NULL, prefix = "06") {
  need("dplyr", "section_overlap()")
  report$h1("Repertoire overlap")

  if (length(unique(as.character(clono$sample_id))) < 2L) {
    report$p("Fewer than two samples; overlap not computed.")
    return(invisible(NULL))
  }

  compute <- function(sub, label) {
    out <- dplyr::bind_rows(
      repertoire_overlap(sub, method = "morisita"),
      repertoire_overlap(sub, method = "jaccard"))
    if (!is.null(label)) out$group <- label
    out
  }

  ov <- if (is.null(within) || !within %in% names(clono)) {
    compute(clono, NULL)
  } else {
    groups <- unique(as.character(clono[[within]]))
    dplyr::bind_rows(lapply(groups, function(g) {
      sub <- clono[as.character(clono[[within]]) == g, , drop = FALSE]
      if (length(unique(as.character(sub$sample_id))) < 2L) return(NULL)
      compute(sub, g)
    }))
  }

  if (!nrow(ov)) {
    report$p("No group had two or more samples to compare.")
    return(invisible(NULL))
  }

  report$tbl(as.data.frame(ov))
  save_table(ov, paste0(prefix, "_overlap"), run_dir,
             "pairwise overlap, Morisita-Horn and Jaccard")
  report$p("")
  report$p(paste0(
    "Morisita-Horn is abundance-weighted and comparatively robust to depth. ",
    "Jaccard counts clonotypes and is the more depth-dependent of the two, ",
    "so where depths differ it should be read at rarefied depth or not at ",
    "all."))

  p <- ggplot2::ggplot(ov, ggplot2::aes(.data$sample_a, .data$sample_b,
                                        fill = .data$value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$value)),
                       size = 2.6, colour = "grey15") +
    ggplot2::facet_wrap(~method) +
    ggplot2::scale_fill_gradientn(colours = palette_sequential(64)) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Overlap",
                  title = "Repertoire overlap between samples") +
    theme_pub(grid = "none") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_figure(p, paste0(prefix, "_overlap"), run_dir, cfg,
              width = 9, height = 4.5,
              legend = paste("Pairwise repertoire overlap. Morisita-Horn is",
                             "abundance-weighted; Jaccard counts clonotypes."))
  invisible(ov)
}

#' Convergence between subjects, against a permutation null
#'
#' @param shm_col,isotype_col Optional columns enabling the
#'   germline-proximity filter; supply them for BCR, leave NULL for TCR.
section_convergence <- function(db, locus, run_dir, cfg, report,
                                shm_col = NULL, isotype_col = NULL,
                                prefix = "07") {
  need("dplyr", "section_convergence()")
  report$h1("Convergence between subjects")

  conv <- filter_locus_productive(db, locus)
  conv$v_gene <- strip_allele(conv$v_call)
  conv <- conv[!is.na(conv$junction_aa) & nzchar(conv$junction_aa), ,
               drop = FALSE]

  if (isTRUE(cfg$convergence$require_mutated_or_switched) &&
      (!is.null(shm_col) || !is.null(isotype_col))) {
    conv <- filter_informative_sequences(
      as.data.frame(conv),
      shm_col = shm_col %||% ".missing",
      isotype_col = isotype_col %||% ".missing",
      shm_threshold = cfg$convergence$shm_threshold)
    report$p(sprintf(
      paste0("Restricted to sequences mutated above %.0f%% or ",
             "class-switched: %s sequences. Unmutated unswitched sequences ",
             "sit close to germline and therefore close to each other, so ",
             "including them measures shared germline rather than shared ",
             "selection."),
      100 * cfg$convergence$shm_threshold,
      format(nrow(conv), big.mark = " ")))
    report$p("")
  } else if (is.null(shm_col) && is.null(isotype_col)) {
    report$p(paste0(
      "No mutation or isotype information, so the germline-proximity filter ",
      "does not apply. For receptors that do not hypermutate this is ",
      "correct; convergence is assessed on the sequences as observed."))
    report$p("")
  }

  conv <- unique(as.data.frame(
    conv[, c("subject_id", "v_gene", "junction_aa")]))

  if (nrow(conv) < 2L ||
      length(unique(as.character(conv$subject_id))) < 2L) {
    report$p("Too few subjects or sequences for a convergence analysis.")
    return(invisible(NULL))
  }

  args <- list(
    cdr3_col = "junction_aa", v_col = "v_gene", subject_col = "subject_id",
    identity_threshold = cfg$convergence$cdr3_identity,
    min_subjects = cfg$convergence$min_subjects,
    normalise = cfg$convergence$distance_normalisation)

  pub  <- do.call(detect_public_groups, c(list(conv), args))
  perm <- do.call(permute_public_groups,
                  c(list(conv), args,
                    list(n_permutations = cfg$convergence$n_permutations)))

  report$h2("Public groups against a label-permutation null")
  report$p(sprintf("observed            %d", perm$observed))
  report$p(sprintf("expected by chance  %.1f", perm$null_mean))
  report$p(sprintf("fold enrichment     %.2f", perm$fold_enrichment))
  report$p(sprintf("empirical p         %.4f  (%d permutations)",
                   perm$p_empirical, perm$n_permutations))
  report$p("")
  report$p(paste0(
    "Subject labels are permuted with the sequences held fixed, so sequence ",
    "composition and V gene structure are preserved and only the assignment ",
    "to subjects is randomised. Without this null the number of public ",
    "groups cannot be read: it rises with depth, with the number of ",
    "subjects and with a looser threshold."))

  if (nrow(pub$summary)) {
    save_table(pub$summary, paste0(prefix, "_public_groups"), run_dir,
               "CDR3 groups shared between subjects")
    save_table(pub$groups, paste0(prefix, "_public_group_members"), run_dir,
               "sequences belonging to each public group")
    report$h2("Largest public groups")
    report$tbl(as.data.frame(utils::head(pub$summary, 15)))
  }

  p <- ggplot2::ggplot(data.frame(n = perm$null), ggplot2::aes(.data$n)) +
    ggplot2::geom_histogram(bins = 30, fill = "grey75", colour = "white",
                            linewidth = 0.2) +
    ggplot2::geom_vline(xintercept = perm$observed, colour = "#D55E00",
                        linewidth = 0.7) +
    ggplot2::annotate("text", x = perm$observed, y = Inf, label = " observed",
                      hjust = 0, vjust = 1.8, colour = "#D55E00", size = 3.2) +
    ggplot2::labs(x = "Public groups per permutation", y = "Permutations",
                  title = "Public groups against a label-permutation null",
                  subtitle = sprintf("empirical p = %.4f over %d permutations",
                                     perm$p_empirical, perm$n_permutations)) +
    theme_pub()
  save_figure(p, paste0(prefix, "_convergence_permutation"), run_dir, cfg,
              width = 5.5, height = 3.6,
              legend = paste("Distribution of the number of public CDR3",
                             "groups after permuting subject labels, with",
                             "the observed value marked."))
  invisible(list(groups = pub, permutation = perm))
}

#' Paired comparison of repertoire metrics between conditions
section_paired_comparisons <- function(rep_summary, run_dir, cfg, report,
                                       group_col = "condition",
                                       metrics = c("n_clonotypes", "top1_pct",
                                                   "shannon", "inv_simpson",
                                                   "gini"),
                                       prefix = "08") {
  need("dplyr", "section_paired_comparisons()")
  report$h1("Paired comparisons between conditions")

  if (!group_col %in% names(rep_summary) ||
      all(is.na(rep_summary[[group_col]]))) {
    report$p(sprintf("No '%s' column; paired comparisons not computed.",
                     group_col))
    return(invisible(NULL))
  }

  input <- rep_summary[!is.na(rep_summary[[group_col]]), , drop = FALSE]
  input$subject_id <- as.character(input$subject_id)
  input[[group_col]] <- as.character(input[[group_col]])
  metrics <- intersect(metrics, names(input))

  report$p(paste0(
    "Each metric is compared between conditions, paired within subject. The ",
    "attainable bound is reported beside every p-value: with few subjects a ",
    "paired signed-rank test has a floor it cannot go below whatever the ",
    "effect size, and in that regime the consistency of the direction ",
    "across subjects is the evidence that does not depend on a test the ",
    "design cannot support."))
  report$p("")

  rows <- list()
  for (m in metrics) {
    res <- try(paired_comparison_all(
      as.data.frame(input), value_col = m, group_col = group_col,
      subject_col = "subject_id", alpha = cfg$statistics$alpha,
      fdr_method = cfg$differential_abundance$fdr_method), silent = TRUE)
    if (!inherits(res, "try-error")) rows[[m]] <- res
  }

  if (!length(rows)) {
    report$p("No metric had two conditions with paired observations.")
    return(invisible(NULL))
  }

  paired <- dplyr::bind_rows(rows)
  save_table(paired, paste0(prefix, "_paired_comparisons"), run_dir,
             "paired comparisons with attainable p and direction concordance")

  keep <- intersect(c("metric", "group_a", "group_b", "n_pairs",
                      "n_effective", "mean_difference", "p_value",
                      "exact_test", "attainable_p", "powered", "direction",
                      "direction_concordant"), names(paired))
  report$tbl(as.data.frame(paired[, keep, drop = FALSE]))

  report$h2("Written out")
  for (ln in format_paired_result(paired)) report$p("  ", ln)

  if (!any(paired$powered, na.rm = TRUE)) {
    report$p("")
    report$p(sprintf(
      paste0("No comparison in this table can reach alpha = %.2f: at this ",
             "sample size the exact paired test bottoms out above it. ",
             "Report direction concordance, not significance."),
      cfg$statistics$alpha))
  }
  invisible(paired)
}

# --- small shared helpers ---------------------------------------------------

#' Parse an AIRR productive flag
#'
#' Files in the wild carry TRUE/FALSE, T/F, True/False and t/f. Parsing them
#' two different ways in two places is how one section ends up analysing
#' rows another discarded, so every caller goes through this.
parse_productive <- function(x) {
  if (is.logical(x)) return(x)
  v <- toupper(trimws(as.character(x)))
  out <- rep(NA, length(v))
  out[v %in% c("T", "TRUE", "1")]  <- TRUE
  out[v %in% c("F", "FALSE", "0")] <- FALSE
  out
}

#' Strip allele and ambiguity from a gene call
#'
#' IGHV1-69*01 -> IGHV1-69, and "TRBV20-1,TRBV20-2" -> TRBV20-1. Allele
#' resolution is not reliable at typical repertoire sequencing depth, and an
#' allele-level key splits genuine clonotypes.
strip_allele <- function(x) {
  sub("\\*.*$", "", sub(",.*$", "", as.character(x)))
}

#' Restrict to productive rearrangements of the requested locus
#'
#' which() rather than a bare logical subscript: a single NA in the
#' productive flag would otherwise inject a row of all-NA instead of
#' dropping the row.
filter_locus_productive <- function(db, locus) {
  keep <- parse_productive(db$productive)
  keep[is.na(keep)] <- FALSE
  if (!is.null(locus) && "locus" %in% names(db)) {
    keep <- keep & db$locus %in% locus
  }
  out <- db[which(keep), , drop = FALSE]
  if (!nrow(out)) {
    stop("No productive ", paste(locus, collapse = "/"), " rearrangements ",
         "after filtering. Check repertoire$locus and the productive column.",
         call. = FALSE)
  }
  out
}
