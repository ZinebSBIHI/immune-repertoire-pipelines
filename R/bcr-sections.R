# ---------------------------------------------------------------------------
# B-cell receptor analysis sections
#
# Shared by the bulk pipeline (locus IGH) and the single-cell pipeline
# (locus IGH). The difference between the two is confined to `mode`:
#
#   "single_cell"  one heavy chain per cell first, doublets removed; light
#                  chains resolve heavy-chain clones into subclones.
#   "bulk"         each row is a sequence; no cell resolution, no light
#                  chains.
#
# Clonal inference follows the order of the original analysis: clones are
# assigned on heavy chains within subject (scoper::hierarchicalClones), and,
# in single-cell data, light chains then split them into subclones
# (dowser::resolveLightChains). scoper can also do both at once from paired
# chains; the two-step order is kept because it is the one the published
# analyses used.
# ---------------------------------------------------------------------------

#' Clonal inference, germlines and somatic hypermutation
#'
#' @param db AIRR data frame (productive filtering is applied here).
#' @param mode "single_cell" or "bulk".
#' @param umi_col Abundance column (UMI count).
#' @inheritParams section_arguments
#' @return The heavy chains (plus light chains in single-cell mode) with
#'   `clone_id`, the germline column used, and `mu_freq_<region>` columns; the
#'   name of the germline column and of the primary mutation column are
#'   attached as attributes `germline_col` and `shm_col`.
#' @export
section_bcr_clones <- function(db, mode = c("single_cell", "bulk"), umi_col,
                               run_dir, cfg, report, prefix = "00") {
  mode <- match.arg(mode)
  need(c("shazam", "scoper"), "section_bcr_clones()")
  report$h1(sprintf("B-cell clonal inference (%s)", sub("_", "-", mode)))

  if (mode == "single_cell") {
    res <- resolve_one_heavy_per_cell(
      db, umi_col = umi_col,
      prefer_switched = cfg$bcr$cell_resolution$prefer_switched_isotype,
      aggregate_umi = cfg$bcr$cell_resolution$aggregate_umi,
      doublet_min_umi = cfg$bcr$cell_resolution$doublet_min_umi)
    cells <- res$db
    save_table(res$filter_log, paste0(prefix, "a_cell_filter_log"), run_dir,
               "rows and cells retained at each filtering step")
    report$h2("One heavy chain per cell")
    report$tbl(res$filter_log)
    report$p("")
    report$p(sprintf(paste0(
      "A barcode carrying two or more heavy chains each at or above %g UMI ",
      "is a doublet and is removed (%d cells) rather than resolved to ",
      "whichever chain had more support."),
      cfg$bcr$cell_resolution$doublet_min_umi, length(res$doublets)))
    heavy <- cells[cells$locus == "IGH", , drop = FALSE]
    light <- cells[cells$locus %in% c("IGK", "IGL"), , drop = FALSE]
  } else {
    heavy <- filter_locus_productive(db, "IGH")
    light <- heavy[0, , drop = FALSE]
  }

  if (!"junction_length" %in% names(heavy)) {
    heavy$junction_length <- nchar(as.character(heavy$junction))
  }
  if (nrow(light) && !"junction_length" %in% names(light)) {
    light$junction_length <- nchar(as.character(light$junction))
  }

  # --- clonal threshold --------------------------------------------------
  threshold <- cfg$bcr$clone_threshold_fixed
  if (identical(cfg$bcr$clone_threshold_mode, "estimate")) {
    dn <- shazam::distToNearest(heavy, sequenceColumn = "junction",
                                vCallColumn = "v_call", jCallColumn = "j_call",
                                model = "ham", normalize = "len",
                                fields = cfg$bcr$clone_fields)
    thr <- tryCatch(
      shazam::findThreshold(dn$dist_nearest, method = "gmm",
                            model = "gamma-norm", cutoff = "optimal"),
      error = function(e) NULL)
    if (is.null(thr) || is.na(thr@threshold)) {
      # The mixture fit needs a bimodal distance distribution, which small
      # or clonally dominated datasets do not always give. Falling back to
      # the fixed value is stated rather than silent.
      report$p(sprintf(paste0(
        "The distance-to-nearest distribution could not be fitted, so the ",
        "fixed threshold %.4f is used. This happens when the data are too ",
        "few or too clonal for a bimodal fit; state it in Methods."),
        threshold))
    } else {
      threshold <- thr@threshold
      report$p(sprintf("Clonal threshold estimated from the data: %.4f",
                       threshold))
    }

    dn_plot <- dn[!is.na(dn$dist_nearest), , drop = FALSE]
    if (nrow(dn_plot)) {
      p <- ggplot2::ggplot(dn_plot, ggplot2::aes(.data$dist_nearest)) +
        ggplot2::geom_histogram(bins = 60, fill = "grey75", colour = "white",
                                linewidth = 0.2) +
        ggplot2::geom_vline(xintercept = threshold, colour = "#D55E00",
                            linewidth = 0.7) +
        ggplot2::facet_wrap(~subject_id, scales = "free_y") +
        ggplot2::labs(x = "Normalised Hamming distance to nearest neighbour",
                      y = "Sequences", title = "Clonal assignment threshold",
                      subtitle = sprintf("threshold %.4f", threshold)) +
        theme_pub()
      save_figure(p, paste0(prefix, "b_clone_threshold"), run_dir, cfg,
                  width = 8, height = 3.4,
                  legend = paste("Distribution of the normalised distance",
                                 "from each heavy chain to its nearest",
                                 "neighbour, per subject, with the clonal",
                                 "assignment threshold."))
    }
  } else {
    report$p(sprintf(paste0(
      "Clonal threshold fixed by configuration: %.4f. A fixed threshold ",
      "makes runs comparable across batches at the cost of not adapting to ",
      "each dataset; either choice must be stated in Methods."), threshold))
  }

  # scoper refuses input that already carries its working columns. Data that
  # went through an earlier Change-O DefineClones step, or a previous run,
  # arrives with a clone_id; it is discarded and recomputed here, and the
  # report says so, because two clone definitions in one table would be
  # worse than either.
  reserved <- c("clone_id", "vj_group", "vjl_group", "junction_l", "cdr3_col",
                "clone_temp", "cell_id_temp")
  pre <- intersect(reserved, names(heavy))
  if (length(pre)) {
    heavy <- heavy[, setdiff(names(heavy), pre), drop = FALSE]
    report$p(sprintf(paste0(
      "The input already carried %s; existing clone assignments are ",
      "discarded and recomputed with the threshold above."),
      paste(pre, collapse = ", ")))
  }
  if (nrow(light)) light <- light[, setdiff(names(light), pre), drop = FALSE]

  clones <- scoper::hierarchicalClones(
    heavy, threshold = threshold, method = cfg$bcr$clone_method,
    linkage = cfg$bcr$clone_linkage, normalize = "len",
    fields = cfg$bcr$clone_fields, summarize_clones = FALSE)

  n_clones <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(clones, .data$subject_id),
    n_sequences = dplyr::n(),
    n_clones = dplyr::n_distinct(.data$clone_id), .groups = "drop"))
  report$h2("Clones per subject")
  report$tbl(as.data.frame(n_clones))
  save_table(n_clones, paste0(prefix, "c_clones_per_subject"), run_dir,
             "sequences and clones per subject")

  # --- germlines -----------------------------------------------------------
  # Three sources, in order of preference. Reconstruction masks the D
  # segment, because D assignment is unreliable at the junction and a wrong
  # call would be scored as mutation exactly where it is hardest to detect.
  imgt_dir <- file.path(cfg$paths$reference_dir, "imgt", "human", "vdj")
  has_col <- function(col) col %in% names(clones) &&
    any(!is.na(clones[[col]]) & nzchar(as.character(clones[[col]])))

  if (has_col("germline_alignment_d_mask")) {
    germline_col <- "germline_alignment_d_mask"
    report$p("Germlines: D-masked germlines already present in the input are used.")
  } else if (dir.exists(imgt_dir)) {
    need("dowser", "Germline reconstruction")
    references <- dowser::readIMGT(dir = imgt_dir, quiet = TRUE)
    clones <- dowser::createGermlines(clones, references, nproc = 1)
    germline_col <- "germline_alignment_d_mask"
    report$p(paste0("Germlines: reconstructed from the IMGT references with ",
                    "the D segment masked."))
  } else if (has_col("germline_alignment")) {
    germline_col <- "germline_alignment"
    report$p(paste0(
      "Germlines: the input's germline_alignment is used. It is NOT ",
      "D-masked, so mutation counts include the D region, where germline ",
      "assignment is least reliable. Provide IMGT references to reconstruct ",
      "D-masked germlines."))
  } else {
    stop("No germline available: the input has no germline_alignment(_d_mask) ",
         "column and no IMGT references were found at ", imgt_dir,
         ".\nSee the pipeline's data/README.md.", call. = FALSE)
  }

  # --- light chains (single-cell) -----------------------------------------
  if (mode == "single_cell" && nrow(light)) {
    need("dowser", "Light-chain subclones")
    both <- dplyr::bind_rows(clones, light)
    subclones <- dowser::resolveLightChains(both, cell = "cell_id",
                                            locus = "locus")
    pct <- 100 * length(unique(light$cell_id)) /
      max(1L, length(unique(clones$cell_id)))
    report$p(sprintf(paste0(
      "Light chains resolve heavy-chain clones into subclones; a light chain ",
      "is available for %.1f%% of cells."), pct))
  } else {
    subclones <- clones
  }

  # --- somatic hypermutation ----------------------------------------------
  # Counted on heavy chains only: they are the chains carrying a germline.
  heavy_rows <- subclones$locus == "IGH"
  shm_heavy <- subclones[heavy_rows, , drop = FALSE]
  shm_other <- subclones[!heavy_rows, , drop = FALSE]

  for (region in cfg$bcr$shm$regions) {
    # IMGT_V and IMGT_VDJ are lazy-loaded *datasets* of shazam, not objects
    # of its namespace: get(envir = asNamespace()) does not see them, while
    # getExportedValue() -- what `::` uses -- looks in the lazy-data table.
    region_def <- getExportedValue("shazam", region)
    shm_heavy <- shazam::observedMutations(
      shm_heavy, sequenceColumn = "sequence_alignment",
      germlineColumn = germline_col, regionDefinition = region_def,
      frequency = TRUE, combine = TRUE, nproc = 1)
    names(shm_heavy)[names(shm_heavy) == "mu_freq"] <-
      paste0("mu_freq_", tolower(region))
  }
  out <- dplyr::bind_rows(shm_heavy, shm_other)
  shm_col <- paste0("mu_freq_", tolower(cfg$bcr$shm$regions[1]))
  out$mu_freq <- out[[shm_col]]
  out$c_gene  <- strip_allele(out$c_call)

  shm_by_sample <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(out[out$locus == "IGH", ], .data$subject_id,
                    .data$sample_id),
    n = dplyr::n(),
    median_mu_freq = stats::median(.data$mu_freq, na.rm = TRUE),
    pct_unmutated = 100 * mean(.data$mu_freq == 0, na.rm = TRUE),
    .groups = "drop"))
  report$h2(sprintf("Somatic hypermutation (%s)", cfg$bcr$shm$regions[1]))
  report$tbl(as.data.frame(shm_by_sample))
  save_table(shm_by_sample, paste0(prefix, "d_shm_by_sample"), run_dir,
             "median mutation frequency and unmutated fraction per sample")

  attr(out, "germline_col") <- germline_col
  attr(out, "shm_col") <- "mu_freq"
  attr(out, "threshold") <- threshold
  out
}

#' Isotype composition, with the unswitched compartment combined
#' @inheritParams section_arguments
#' @export
section_isotypes <- function(db, run_dir, cfg, report, unit = "cells",
                             prefix = "00e") {
  need("dplyr", "section_isotypes()")
  report$h1("Isotypes")
  report$p(paste0(
    "One isotype is retained per cell (or sequence), so an IgM versus IgD ",
    "split would describe which transcript dominated rather than two ",
    "populations: naive and marginal-zone B cells co-transcribe IGHM and ",
    "IGHD. The combined unswitched fraction is reported instead."))
  report$p("")

  iso <- summarise_isotypes(db, group_cols = c("subject_id", "sample_id"),
                            combine_unswitched = cfg$bcr$report_unswitched_combined)
  sw <- switched_fraction(db, group_cols = c("subject_id", "sample_id"))
  save_table(iso, paste0(prefix, "_isotype_composition"), run_dir,
             "isotype composition per sample")
  save_table(sw, paste0(prefix, "_switched_fraction"), run_dir,
             "class-switched fraction per sample")
  report$tbl(as.data.frame(sw))

  p <- ggplot2::ggplot(iso, ggplot2::aes(.data$sample_id, .data$fraction,
                                         fill = .data$isotype_class)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::facet_wrap(~subject_id, scales = "free_x") +
    ggplot2::scale_fill_manual(
      values = palette_categorical(length(unique(iso$isotype_class)))) +
    ggplot2::scale_y_continuous(labels = scales::label_percent()) +
    ggplot2::labs(x = NULL, y = tools_title(unit), fill = "Isotype",
                  title = "Isotype composition") +
    theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_figure(p, paste0(prefix, "_isotype_composition"), run_dir, cfg,
              width = 9, height = 4,
              legend = paste("Isotype composition per sample. IGHM and IGHD",
                             "are reported as a single unswitched class."))
  invisible(list(isotypes = iso, switched = sw))
}

#' Clonal expansion, with a threshold matched to the material
#'
#' `stimulated` must be a column of the data, declared by the cohort. It is
#' never inferred from a sample label, because guessing it would silently
#' apply the wrong threshold.
#' @inheritParams section_arguments
#' @export
section_clonal_expansion <- function(db, run_dir, cfg, report,
                                     unit = "cells", prefix = "00f") {
  need("dplyr", "section_clonal_expansion()")
  report$h1("Clonal expansion")
  report$p(sprintf(paste0(
    "Expansion thresholds differ by material: %d %s per clone for resting ",
    "and %d for stimulated. A multi-day stimulation culture produces ",
    "identical daughter cells, so the cut-off that marks genuine expansion ",
    "in resting material is met trivially after stimulation."),
    cfg$bcr$shm$expansion_min_cells_resting, unit,
    cfg$bcr$shm$expansion_min_cells_stimulated))
  report$p("")

  heavy <- db[db$locus == "IGH", , drop = FALSE]
  if (!"stimulated" %in% names(heavy)) {
    heavy$stimulated <- FALSE
    report$p(paste0(
      "No 'stimulated' column; the resting threshold is applied throughout. ",
      "Add the column upstream if the cohort contains stimulated material."))
    report$p("")
  }

  sizes <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(heavy, .data$subject_id, .data$sample_id,
                    .data$stimulated, .data$clone_id),
    n = dplyr::n(), .groups = "drop"))
  sizes$threshold <- expansion_threshold(
    sizes$stimulated, cfg$bcr$shm$expansion_min_cells_resting,
    cfg$bcr$shm$expansion_min_cells_stimulated)
  sizes$expanded <- sizes$n >= sizes$threshold

  expansion <- dplyr::ungroup(dplyr::summarise(
    dplyr::group_by(sizes, .data$subject_id, .data$sample_id),
    n_clones = dplyr::n(),
    n_expanded = sum(.data$expanded),
    pct_in_expanded = 100 * sum(.data$n[.data$expanded]) / sum(.data$n),
    largest_clone_pct = 100 * max(.data$n) / sum(.data$n),
    .groups = "drop"))
  report$tbl(as.data.frame(expansion))
  save_table(expansion, paste0(prefix, "_clonal_expansion"), run_dir,
             "expanded clones and the share they hold")
  invisible(expansion)
}

#' Dominant clones: profiled, never silently excluded
#' @inheritParams section_arguments
#' @param min_fraction Smallest share of a subject's cells for a clone to be profiled.
#' @export
section_dominant_clones <- function(db, umi_col, run_dir, cfg, report,
                                    min_fraction = 0.01, prefix = "00g") {
  report$h1("Dominant clones")
  report$p(paste0(
    "Clones above ", 100 * min_fraction, "% of a subject are profiled, ",
    "because any repertoire-level statistic on that compartment follows ",
    "them: gene usage, CDR3 length and mutation load all do. Two criteria ",
    "separate real expansion from artefact: UMI support at the level of the ",
    "subject's other cells, and a paired light chain. Intraclonal diversity ",
    "is deliberately not a criterion, since a genuine expansion produces ",
    "identical daughter cells. Nothing is excluded here; the verdict ",
    "records which criteria each clone meets, and any statistic an ",
    "exclusion would change should be reported both ways."))
  report$p("")

  if (!"cell_id" %in% names(db) || !umi_col %in% names(db)) {
    report$p("No per-cell barcode or UMI column; dominant-clone profiling skipped.")
    return(invisible(NULL))
  }
  dom <- profile_dominant_clones(db, clone_col = "clone_id",
                                 subject_col = "subject_id",
                                 umi_col = umi_col, min_fraction = min_fraction)
  if (!nrow(dom)) {
    report$p(sprintf("No clone exceeds %.0f%% of any subject.",
                     100 * min_fraction))
    return(invisible(dom))
  }
  keep <- intersect(c("subject_id", "clone_id", "n_cells",
                      "fraction_of_subject", "v_gene", "isotype", "umi_ratio",
                      "pct_with_light_chain", "verdict"), names(dom))
  report$tbl(as.data.frame(dom[, keep, drop = FALSE]))
  save_table(dom, paste0(prefix, "_dominant_clones"), run_dir,
             "clones above the threshold, with supporting criteria")
  invisible(dom)
}

#' Title-case a unit for an axis label
#' @noRd
tools_title <- function(unit) {
  paste0(toupper(substr(unit, 1, 1)), substr(unit, 2, nchar(unit)))
}
