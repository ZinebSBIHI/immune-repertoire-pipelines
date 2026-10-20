#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Pipeline 03 — Bulk TCR/BCR repertoire
#
# Input   AIRR Rearrangement TSV files, one per sample, as produced by MiXCR,
#         IgBLAST + MakeDb, or any AIRR-compliant assembler.
# Output  A timestamped run directory holding figures, tables, a numeric
#         report, the resolved configuration and sessionInfo().
#
# Run:
#   Rscript examples/simulate_airr.R          # once, to create example data
#   cp config/config.example.yml config/config.yml
#   Rscript pipelines/03_bulk_repertoire.R
#
# Steps
#   1  Read inputs, parse metadata, QC the composition of each sample
#   2  Sweep the abundance floor and report what each value costs
#   3  Repertoire structure: depth, richness, clonality, size classes
#   4  Diversity at native and at common (rarefied) depth
#   5  Repertoire overlap against the subject's own baseline
#   6  Gene usage, at clonotype and at molecule level
#   7  Differential clonal abundance: category and significance, kept apart
#   8  Convergence between subjects, against a permutation null
#   9  Numeric report
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

# Locate the repository root, whether the pipeline is launched from the
# repository directory, from pipelines/, or by absolute path.
ROOT <- local({
  candidates <- getwd()
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) candidates <- c(dirname(dirname(normalizePath(f[1]))), candidates)
  for (start in candidates) {
    d <- normalizePath(start, mustWork = FALSE)
    for (i in 1:6) {
      if (dir.exists(file.path(d, "R")) &&
          file.exists(file.path(d, "config", "config.example.yml"))) return(d)
      if (identical(dirname(d), d)) break
      d <- dirname(d)
    }
  }
  stop("Could not locate the repository root. Run from the repository ",
       "directory:\n  Rscript pipelines/03_bulk_repertoire.R", call. = FALSE)
})
for (.f in list.files(file.path(ROOT, "R"), pattern = "[.]R$",
                      full.names = TRUE)) source(.f)

cfg     <- load_config(file.path(ROOT, "config/config.yml"))
run_dir <- init_run(cfg, "03_bulk_repertoire")
use_pub_theme()

report <- open_report(file.path(run_dir, "REPORT.txt"),
                      sprintf("Bulk repertoire analysis - %s",
                              cfg$project$name))

# ===========================================================================
# 1  Input, metadata, composition QC
# ===========================================================================
report$h1("1. Input and composition")

files <- list.files(cfg$paths$data_dir, pattern = "[.](airr[.])?tsv$",
                    full.names = TRUE, recursive = TRUE)
if (!length(files)) {
  stop("No AIRR .tsv files found under ", cfg$paths$data_dir,
       "\nRun  Rscript examples/simulate_airr.R  to create an example dataset.",
       call. = FALSE)
}
message(sprintf("Found %d AIRR file(s).", length(files)))

meta <- parse_filename_metadata(files, cfg$cohort$filename_pattern,
                                required = c("subject", "timepoint"))

db <- bind_rows(lapply(seq_len(nrow(meta)), function(i) {
  x <- utils::read.delim(meta$file[i], stringsAsFactors = FALSE,
                         check.names = FALSE)
  x$sample_id  <- sub("[.](airr[.])?tsv$", "", meta$basename[i])
  x$subject_id <- meta$subject[i]
  x$timepoint  <- meta$timepoint[i]
  x$condition  <- meta$condition[i] %||% NA_character_
  x
}))

db$subject_id <- as_configured_factor(db$subject_id, cfg$cohort$subjects,
                                      "subject")
db$timepoint  <- as_configured_factor(db$timepoint, cfg$cohort$timepoints,
                                      "timepoint")
if (!all(is.na(db$condition))) {
  db$condition <- as_configured_factor(db$condition, cfg$cohort$conditions,
                                       "condition")
}

# The abundance column comes from the configuration and is never guessed:
# UMI counts and read counts are not interchangeable, and `duplicate_count`
# means molecules in some assemblers' output and reads in others'.
abundance_col <- cfg$repertoire$abundance_column
if (!abundance_col %in% names(db)) {
  stop("repertoire$abundance_column is '", abundance_col, "' but no such ",
       "column is present. Available: ", paste(names(db), collapse = ", "),
       call. = FALSE)
}
report$p(sprintf("Abundance column: %s", abundance_col))
report$p(sprintf("Locus analysed:   %s", cfg$repertoire$locus))
report$p(sprintf("Clonotype:        %s", cfg$repertoire$clonotype_definition))
report$p("")

# Composition QC. A sample whose productive fraction or locus composition is
# out of line with the rest is a library problem, and is better caught here
# than inferred later from an odd diversity value.
qc <- db %>%
  group_by(sample_id, subject_id, timepoint) %>%
  summarise(
    n_rearrangements = n(),
    pct_productive = 100 * mean(as.logical(productive), na.rm = TRUE),
    pct_target_locus = 100 * mean(locus == cfg$repertoire$locus, na.rm = TRUE),
    .groups = "drop"
  )
report$h2("Composition per sample")
report$tbl(as.data.frame(qc))
save_table(qc, "01_composition_qc", run_dir, "productive and locus fractions")

# ===========================================================================
# 2  Abundance floor sweep
# ===========================================================================
report$h1("2. Abundance floor")
report$p(
  "The abundance floor applies to both timepoints, and it acts on de novo",
  "\ncalls in two opposing directions. Raising it removes single-molecule",
  "\nclonotypes, which reduces the calls; it also reclassifies low-abundance",
  "\nbaseline clonotypes as not detected, which creates calls, because",
  "\nabsence at baseline is the definition of de novo. The net direction is",
  "\na property of the dataset, not of the method, which is why the sweep is",
  "\nrun rather than a default being assumed. Report the chosen value in",
  "\nMethods and the sweep as a supplementary figure.")
report$p("")

sweep <- sweep_abundance_floor(
  db, abundance_col = abundance_col,
  floors = cfg$repertoire$min_abundance_sweep,
  definition = cfg$repertoire$clonotype_definition,
  reference_timepoint = cfg$differential_abundance$reference_timepoint,
  fold_change_threshold = cfg$differential_abundance$fold_change_threshold,
  fdr_threshold = cfg$differential_abundance$fdr_threshold,
  fdr_method = cfg$differential_abundance$fdr_method
)
save_table(sweep, "02_abundance_floor_sweep", run_dir,
           "cost of each abundance floor, per sample")

sweep_summary <- sweep %>%
  group_by(min_abundance) %>%
  summarise(
    median_pct_clonotypes = stats::median(pct_clonotypes_retained, na.rm = TRUE),
    median_pct_abundance  = stats::median(pct_abundance_retained,  na.rm = TRUE),
    median_pct_de_novo    = stats::median(pct_de_novo_retained,    na.rm = TRUE),
    .groups = "drop"
  )
report$h2("Retained at each floor, relative to the most permissive value")
report$tbl(as.data.frame(sweep_summary), digits = 1)

p_sweep <- sweep %>%
  select(sample_id, subject_id, min_abundance,
         Clonotypes = pct_clonotypes_retained,
         Molecules  = pct_abundance_retained,
         `De novo calls` = pct_de_novo_retained) %>%
  pivot_longer(c(Clonotypes, Molecules, `De novo calls`),
               names_to = "quantity", values_to = "pct") %>%
  # Baseline samples have no de novo calls by construction, so they
  # contribute no point to that panel. Dropping them here is explicit;
  # leaving it to the geom would produce a missing-value warning that looks
  # like a defect.
  filter(!is.na(pct)) %>%
  ggplot(aes(factor(min_abundance), pct, group = sample_id)) +
  geom_line(alpha = 0.35, linewidth = 0.3) +
  geom_point(aes(shape = subject_id), size = 1.6, alpha = 0.8) +
  facet_wrap(~quantity) +
  labs(x = "Minimum abundance (UMI)", y = "Retained (%)", shape = "Subject",
       title = "Cost of the abundance floor") +
  theme_pub()
save_figure(p_sweep, "02_abundance_floor_sweep", run_dir, cfg,
            width = 9, height = 3.4,
            legend = paste("Percentage of clonotypes, molecules and de novo",
                           "calls retained at each abundance floor, relative",
                           "to the most permissive floor, one line per sample."))

# ===========================================================================
# 3  Repertoire structure
# ===========================================================================
report$h1("3. Repertoire structure")

clono <- collapse_clonotypes(
  db, abundance_col = abundance_col,
  definition = cfg$repertoire$clonotype_definition,
  min_abundance = cfg$repertoire$min_abundance,
  locus = cfg$repertoire$locus
)
sample_meta <- unique(db[, c("sample_id", "subject_id", "timepoint",
                             "condition")])
clono <- left_join(clono %>% select(-subject_id), sample_meta, by = "sample_id")

rep_summary <- repertoire_summary(clono) %>%
  left_join(sample_meta, by = c("sample_id", "subject_id"))
report$h2("Per-sample summary")
report$tbl(as.data.frame(rep_summary))
save_table(rep_summary, "03_repertoire_summary", run_dir,
           "depth, richness and clonality per sample")

depth_range <- range(rep_summary$depth)
report$p("")
report$p(sprintf(
  paste0("Sequencing depth spans %s to %s molecules, a %.1f-fold range. ",
         "Richness, clone counts and clonal concentration are all monotone ",
         "in depth, so comparisons of those quantities are made at common ",
         "depth in section 4."),
  format(depth_range[1], big.mark = " "),
  format(depth_range[2], big.mark = " "),
  depth_range[2] / depth_range[1]))

p_depth <- ggplot(rep_summary, aes(timepoint, depth, group = subject_id)) +
  geom_line(alpha = 0.5) +
  geom_point(aes(colour = subject_id), size = 2.4) +
  scale_y_continuous(labels = scales::label_number(big.mark = " ")) +
  labs(x = NULL, y = "Productive molecules", colour = "Subject",
       title = "Sequencing depth") +
  theme_pub()

p_rich <- ggplot(rep_summary, aes(timepoint, n_clonotypes, group = subject_id)) +
  geom_line(alpha = 0.5) +
  geom_point(aes(colour = subject_id), size = 2.4) +
  scale_y_continuous(labels = scales::label_number(big.mark = " ")) +
  labs(x = NULL, y = "Distinct clonotypes", colour = "Subject",
       title = "Observed richness",
       subtitle = "Depth-dependent; see rarefied values") +
  theme_pub()

save_figure(p_depth, "03a_depth", run_dir, cfg, width = 5, height = 3.4)
save_figure(p_rich,  "03b_richness_observed", run_dir, cfg, width = 5,
            height = 3.4)

homeo <- clonal_homeostasis(clono)
save_table(homeo, "03_clonal_homeostasis", run_dir,
           "repertoire fraction by clonotype size class")

p_homeo <- homeo %>%
  left_join(sample_meta, by = c("sample_id", "subject_id")) %>%
  ggplot(aes(timepoint, repertoire_fraction, fill = size_class)) +
  geom_col(width = 0.75) +
  facet_wrap(~subject_id) +
  scale_fill_manual(values = palette_sequential(
    length(levels(homeo$size_class)), direction = -1)) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(x = NULL, y = "Repertoire", fill = "Clonotype frequency",
       title = "Clonal space homeostasis") +
  theme_pub()
save_figure(p_homeo, "03c_clonal_homeostasis", run_dir, cfg,
            width = 8, height = 3.4,
            legend = paste("Fraction of the repertoire held by clonotypes in",
                           "each frequency class, per sample."))

# ===========================================================================
# 4  Diversity at native and common depth
# ===========================================================================
report$h1("4. Diversity")

if (isTRUE(cfg$repertoire$rarefaction$enabled)) {
  # The common depth is taken within subject, not across the cohort: a
  # within-subject comparison of timepoints should not be degraded to
  # accommodate the shallowest sample in another subject.
  rar <- rarefy_to_common_depth(
    clono, within = "subject_id",
    n_bootstrap = cfg$repertoire$rarefaction$n_bootstrap)
  rar <- left_join(rar, sample_meta, by = c("sample_id", "subject_id"))
  save_table(rar, "04_rarefied_metrics", run_dir,
             "metrics at within-subject common depth, with 95% intervals")

  report$h2("Common depth per subject")
  report$tbl(as.data.frame(unique(rar[, c("subject_id", "target_depth")])))

  report$h2("Rarefied metrics")
  report$tbl(as.data.frame(rar[, c("subject_id", "sample_id", "statistic",
                                   "mean", "q025", "q975")]))

  p_rar <- ggplot(rar, aes(timepoint, mean, group = subject_id)) +
    geom_errorbar(aes(ymin = q025, ymax = q975), width = 0.12,
                  linewidth = 0.3, colour = "grey50") +
    geom_line(alpha = 0.5) +
    geom_point(aes(colour = subject_id), size = 2.2) +
    facet_wrap(~statistic, scales = "free_y") +
    labs(x = NULL, y = "Value at common depth", colour = "Subject",
         title = "Repertoire metrics, rarefied to within-subject common depth",
         subtitle = sprintf("%d resamples; bars are 95%% intervals",
                            cfg$repertoire$rarefaction$n_bootstrap)) +
    theme_pub()
  save_figure(p_rar, "04_rarefied_metrics", run_dir, cfg,
              width = 9, height = 5,
              legend = paste("Repertoire metrics after resampling every",
                             "sample to the common depth of its subject.",
                             "Points are the resampling mean, bars the 2.5th",
                             "and 97.5th percentiles."))
} else {
  report$p("Rarefaction disabled in the configuration. Any comparison of",
           "\nrichness, clone counts or clonal concentration between samples",
           "\nof unequal depth is confounded by depth.")
}

# ===========================================================================
# 5  Overlap with the subject's own baseline
# ===========================================================================
report$h1("5. Overlap with baseline")

ref_tp <- cfg$differential_abundance$reference_timepoint
ov_rows <- list()
for (subj in levels(droplevels(clono$subject_id))) {
  sub_clono <- clono[clono$subject_id == subj, , drop = FALSE]
  if (length(unique(sub_clono$sample_id)) < 2L) next
  ov <- repertoire_overlap(sub_clono, method = "morisita")
  ov$subject_id <- subj
  ov_rows[[length(ov_rows) + 1L]] <- ov
}
if (length(ov_rows)) {
  ov_all <- bind_rows(ov_rows)
  save_table(ov_all, "05_overlap_morisita", run_dir,
             "pairwise Morisita-Horn overlap within subject")
  report$h2("Morisita-Horn overlap within subject")
  report$tbl(as.data.frame(ov_all))
  report$p("")
  report$p("Morisita-Horn is abundance-weighted and comparatively robust to",
           "\ndepth. Jaccard counts clonotypes and would need rarefied input.")
}

# ===========================================================================
# 6  Gene usage
# ===========================================================================
report$h1("6. V and J gene usage")
report$p(
  "Usage is reported at clonotype level and at molecule level. At molecule",
  "\nlevel a single large clone dominates its sample's profile, so an",
  "\napparent gene bias can be one expanded lineage. Where the two levels",
  "\ndisagree, the clone is the explanation.")
report$p("")

# which() rather than a bare logical subscript: a single NA in `productive`
# would otherwise inject a row of all-NA rather than drop the row. The
# productive flag is parsed the same way as in collapse_clonotypes(), so
# the two cannot disagree on a file that writes "t"/"f" or "True"/"False".
prod_flag <- db$productive
if (is.character(prod_flag)) prod_flag <- toupper(prod_flag) %in% c("T", "TRUE")
db_target <- db[which(db$locus == cfg$repertoire$locus &
                        as.logical(prod_flag)), , drop = FALSE]
if (!nrow(db_target)) {
  stop("No productive ", cfg$repertoire$locus, " rearrangements after ",
       "filtering. Check repertoire$locus and the productive column.",
       call. = FALSE)
}

usage_clone <- gene_usage(db_target, "v_call", level = "clonotype",
                          definition = cfg$repertoire$clonotype_definition)
usage_mol   <- gene_usage(db_target, "v_call", level = "molecule",
                          abundance_col = abundance_col,
                          definition = cfg$repertoire$clonotype_definition)
usage <- bind_rows(usage_clone, usage_mol) %>%
  left_join(sample_meta, by = "sample_id")
save_table(usage, "06_v_gene_usage", run_dir,
           "V gene usage at clonotype and molecule level")

top_v <- usage_clone %>%
  group_by(gene) %>% summarise(m = mean(frequency), .groups = "drop") %>%
  slice_max(m, n = 15) %>% pull(gene)

p_usage <- usage %>%
  filter(gene %in% top_v) %>%
  ggplot(aes(reorder(gene, frequency), frequency, colour = level)) +
  stat_summary(fun = mean, geom = "point", size = 1.8,
               position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2,
               linewidth = 0.3, position = position_dodge(width = 0.5)) +
  coord_flip() +
  scale_colour_manual(values = palette_categorical(2)) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(x = NULL, y = "Usage", colour = "Level",
       title = "V gene usage, 15 most used genes",
       subtitle = "Clonotype level versus molecule level") +
  theme_pub(grid = "x")
save_figure(p_usage, "06_v_gene_usage", run_dir, cfg, width = 6, height = 5,
            legend = paste("Mean V gene usage across samples at clonotype",
                           "and molecule level, with standard error."))

# ===========================================================================
# 7  Differential clonal abundance
# ===========================================================================
report$h1("7. Differential clonal abundance")
report$p(
  "Two axes, reported separately:",
  "\n  category     descriptive, relative to the subject's own baseline",
  "\n  significance Fisher exact test on counts, with FDR control",
  "\nA clonotype can be de novo and non-significant, or pre-existing and",
  "\nsignificantly expanded. No single label asserting induction is produced.")
report$p("")

calls <- test_clonal_abundance(
  clono,
  reference_timepoint = ref_tp,
  fold_change_threshold = cfg$differential_abundance$fold_change_threshold,
  detection_floor = cfg$differential_abundance$detection_floor,
  fdr_threshold = cfg$differential_abundance$fdr_threshold,
  fdr_method = cfg$differential_abundance$fdr_method
)
calls <- left_join(calls, sample_meta %>% select(sample_id, condition),
                   by = "sample_id")
save_table(calls %>% filter(sig_expansion | category == "de_novo"),
           "07_differential_abundance_calls", run_dir,
           "clonotypes significantly expanded or called de novo")

call_summary <- summarise_abundance_calls(calls)
report$h2("Calls per sample")
report$tbl(as.data.frame(call_summary))
save_table(call_summary, "07_differential_abundance_summary", run_dir,
           "counts by category and significance, per sample")

persist <- clonotype_persistence(calls)
if (nrow(persist)) {
  save_table(persist, "07_clonotype_persistence", run_dir,
             "timepoints at which each significantly expanded clonotype is called")
  report$h2("Persistence of significantly expanded clonotypes")
  report$tbl(as.data.frame(
    persist %>% count(subject_id, n_timepoints, name = "n_clonotypes")))
  report$p("")
  report$p(sprintf(
    "%.0f%% of significantly expanded clonotypes are called at a single timepoint.",
    100 * mean(persist$n_timepoints == 1)))
}

p_volcano <- calls %>%
  filter(is.finite(log2_fold_change), !is.na(fdr)) %>%
  mutate(neglog_fdr = pmin(-log10(fdr), 50)) %>%
  thin_points(keep = "sig_expansion", max_background = 15000) %>%
  ggplot(aes(log2_fold_change, neglog_fdr)) +
  geom_point(aes(colour = sig_expansion, shape = category), size = 0.8,
             alpha = 0.6) +
  geom_hline(yintercept = -log10(cfg$differential_abundance$fdr_threshold),
             linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = log2(cfg$differential_abundance$fold_change_threshold),
             linetype = "dashed", linewidth = 0.3) +
  facet_grid(subject_id ~ timepoint) +
  scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#D55E00"),
                      labels = c("not significant", "significant")) +
  labs(x = expression(log[2]~"fold change vs baseline"),
       y = expression(-log[10]~"FDR (capped at 50)"),
       colour = NULL, shape = "Category",
       title = "Differential clonal abundance") +
  theme_pub(grid = "both")
save_figure(p_volcano, "07a_volcano", run_dir, cfg, width = 9, height = 6,
            legend = paste("Fold change against the subject's own baseline",
                           "versus FDR. Dashed lines mark the configured",
                           "thresholds. Colour gives statistical",
                           "significance, shape the descriptive category;",
                           "the two are independent."))

p_counts <- call_summary %>%
  select(subject_id, timepoint,
         `De novo` = n_sig_expanded_de_novo,
         `Pre-existing` = n_sig_expanded_pre_existing) %>%
  pivot_longer(c(`De novo`, `Pre-existing`), names_to = "origin",
               values_to = "n") %>%
  ggplot(aes(timepoint, n, fill = origin)) +
  geom_col(width = 0.7) +
  facet_wrap(~subject_id) +
  scale_fill_manual(values = palette_categorical(2)) +
  labs(x = NULL, y = "Significantly expanded clonotypes", fill = "Origin",
       title = "Significantly expanded clonotypes by origin") +
  theme_pub()
save_figure(p_counts, "07b_sig_expanded_by_origin", run_dir, cfg,
            width = 8, height = 3.4,
            legend = paste("Clonotypes meeting both the fold-change and FDR",
                           "criteria, split by whether they were detected at",
                           "baseline."))

# ===========================================================================
# 8  Paired comparison of repertoire metrics across conditions
# ===========================================================================
report$h1("8. Paired comparisons across conditions")
report$p(
  "Each metric is compared between conditions, paired within subject. The",
  "\nattainable bound is reported beside every p-value: with few subjects a",
  "\npaired signed-rank test has a floor it cannot go below whatever the",
  "\neffect size, and in that regime the consistency of the direction across",
  "\nsubjects is the evidence that does not depend on a test the design",
  "\ncannot support.")
report$p("")

metric_input <- rep_summary %>%
  filter(!is.na(condition)) %>%
  select(subject_id, condition, n_clonotypes, top1_pct, shannon,
         inv_simpson, gini) %>%
  mutate(subject_id = as.character(subject_id),
         condition = as.character(condition))

metrics <- c("n_clonotypes", "top1_pct", "shannon", "inv_simpson", "gini")
paired_rows <- list()
for (m in metrics) {
  res <- try(paired_comparison_all(
    as.data.frame(metric_input), value_col = m, group_col = "condition",
    subject_col = "subject_id", alpha = cfg$statistics$alpha,
    fdr_method = cfg$differential_abundance$fdr_method), silent = TRUE)
  if (!inherits(res, "try-error")) paired_rows[[m]] <- res
}

if (length(paired_rows)) {
  paired <- bind_rows(paired_rows)
  save_table(paired, "08_paired_comparisons", run_dir,
             "paired comparisons with attainable p and direction concordance")

  report$tbl(as.data.frame(
    paired[, c("metric", "group_a", "group_b", "n_pairs", "n_effective",
               "mean_difference", "p_value", "exact_test", "attainable_p",
               "powered", "direction", "direction_concordant")]))

  report$h2("Written out")
  for (ln in format_paired_result(paired)) report$p("  ", ln)

  if (!any(paired$powered, na.rm = TRUE)) {
    report$p("")
    report$p(sprintf(
      paste0("No comparison in this table can reach alpha = %.2f: with %d ",
             "subjects the exact paired test bottoms out above it. Report ",
             "direction concordance, not significance."),
      cfg$statistics$alpha, max(paired$n_pairs, na.rm = TRUE)))
  }
} else {
  report$p("No metric had two conditions with paired observations.")
}

# ===========================================================================
# 9  Convergence between subjects
# ===========================================================================
report$h1("9. Convergence between subjects")

conv_input <- db_target %>%
  mutate(v_gene = sub("\\*.*$", "", sub(",.*$", "", v_call))) %>%
  filter(!is.na(junction_aa), nchar(junction_aa) > 0) %>%
  distinct(subject_id, v_gene, junction_aa)

if (nrow(conv_input) > 1L && n_distinct(conv_input$subject_id) >= 2L) {
  pub <- detect_public_groups(
    as.data.frame(conv_input),
    cdr3_col = "junction_aa", v_col = "v_gene", subject_col = "subject_id",
    identity_threshold = cfg$convergence$cdr3_identity,
    min_subjects = cfg$convergence$min_subjects,
    normalise = cfg$convergence$distance_normalisation)

  report$p(sprintf("Public groups observed: %d", nrow(pub$summary)))

  if (nrow(pub$summary)) {
    save_table(pub$summary, "09_public_groups", run_dir,
               "CDR3 groups shared between subjects")
    report$h2("Largest public groups")
    report$tbl(as.data.frame(utils::head(pub$summary, 15)))
  }

  perm <- permute_public_groups(
    as.data.frame(conv_input),
    cdr3_col = "junction_aa", v_col = "v_gene", subject_col = "subject_id",
    identity_threshold = cfg$convergence$cdr3_identity,
    min_subjects = cfg$convergence$min_subjects,
    normalise = cfg$convergence$distance_normalisation,
    n_permutations = cfg$convergence$n_permutations)

  report$h2("Permutation null")
  report$p(sprintf("observed            %d", perm$observed))
  report$p(sprintf("expected by chance  %.1f", perm$null_mean))
  report$p(sprintf("fold enrichment     %.2f", perm$fold_enrichment))
  report$p(sprintf("empirical p         %.4f  (%d permutations)",
                   perm$p_empirical, perm$n_permutations))
  report$p("")
  report$p("Subject labels are permuted with the sequences held fixed, so",
           "\nsequence composition and V gene structure are preserved and only",
           "\nthe assignment to subjects is randomised. Without this null the",
           "\nnumber of public groups cannot be interpreted: it rises with",
           "\ndepth, with the number of subjects and with a looser threshold.")

  p_perm <- ggplot(data.frame(n = perm$null), aes(n)) +
    geom_histogram(bins = 30, fill = "grey75", colour = "white",
                   linewidth = 0.2) +
    geom_vline(xintercept = perm$observed, colour = "#D55E00",
               linewidth = 0.7) +
    annotate("text", x = perm$observed, y = Inf, label = " observed",
             hjust = 0, vjust = 1.8, colour = "#D55E00", size = 3.2) +
    labs(x = "Public groups per permutation", y = "Permutations",
         title = "Public groups against a label-permutation null",
         subtitle = sprintf("empirical p = %.4f over %d permutations",
                            perm$p_empirical, perm$n_permutations)) +
    theme_pub()
  save_figure(p_perm, "09_convergence_permutation", run_dir, cfg,
              width = 5.5, height = 3.6,
              legend = paste("Distribution of the number of public CDR3",
                             "groups after permuting subject labels, with the",
                             "observed value marked."))
} else {
  report$p("Too few subjects or sequences for a convergence analysis.")
}

# ===========================================================================
# 10  Close
# ===========================================================================
report$h1("10. Notes")
report$p("Parameters are recorded in resolved_config.yml; package versions in",
         "\nsessionInfo.txt. Figure legends are in FIGURE_LEGENDS.txt and the",
         "\ntable inventory in TABLE_MANIFEST.txt.")
report$close()

finalise_run(cfg, run_dir)
