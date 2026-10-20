# ---------------------------------------------------------------------------
# Run directories, figure and table output
#
# Two rules are enforced here rather than left to each pipeline:
#
#   1. No run overwrites another. Every pipeline call creates a fresh
#      timestamped directory. Re-running an analysis to answer a reviewer
#      three months later must not destroy the figures that went into the
#      submission.
#
#   2. Every output directory records how it was produced: the resolved
#      configuration and sessionInfo(). A figure whose parameters cannot be
#      recovered is not reproducible, however clean the code that made it.
# ---------------------------------------------------------------------------

#' Create a timestamped run directory
#'
#' @param cfg Configuration list from `load_config()`.
#' @param pipeline Short pipeline identifier, e.g. "03_bulk_repertoire".
#' @return Path to the created directory, with `figures/` and `tables/`
#'   subdirectories.
init_run <- function(cfg, pipeline) {
  # A run that aborts on a missing input or an absent package leaves behind
  # a directory holding only scaffolding. Those accumulate and make it hard
  # to see which runs actually produced something, so they are swept here.
  prune_empty_runs(cfg$paths$output_dir)

  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  run_dir <- file.path(cfg$paths$output_dir,
                       sprintf("%s_%s_%s", pipeline, cfg$project$name, stamp))

  dir.create(file.path(run_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(run_dir, "tables"),  recursive = TRUE, showWarnings = FALSE)

  set.seed(cfg$project$seed)

  if (isTRUE(cfg$output$write_resolved_config)) {
    if (requireNamespace("yaml", quietly = TRUE)) {
      to_write <- cfg
      to_write$.meta <- NULL
      yaml::write_yaml(to_write, file.path(run_dir, "resolved_config.yml"))
    }
  }

  writeLines(
    c(RUN_MARKER,
      sprintf("pipeline:  %s", pipeline),
      sprintf("project:   %s", cfg$project$name),
      sprintf("started:   %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
      sprintf("seed:      %d", cfg$project$seed),
      sprintf("config:    %s", cfg$.meta$config_path %||% "(unknown)"),
      sprintf("R:         %s", paste(R.version$major, R.version$minor, sep = ".")),
      sprintf("platform:  %s", R.version$platform)),
    file.path(run_dir, "RUN_INFO.txt")
  )

  message("Run directory: ", run_dir)
  run_dir
}

#' Remove run directories that hold no results
#'
#' A directory is pruned only when BOTH conditions hold:
#'   - it is one of ours, proven by a readable RUN_INFO.txt carrying the
#'     marker line that `init_run()` writes;
#'   - it contains no figure, no table and no saved object, only the
#'     scaffolding.
#'
#' The ownership check is the important one. The output directory is a path
#' the user chose, and it may well hold their own folders: a directory of
#' manuscript panels, a notes folder. Deleting on emptiness alone would
#' silently remove those, irreversibly. Anything without our marker is left
#' untouched however empty it looks.
prune_empty_runs <- function(output_dir) {
  if (!dir.exists(output_dir)) return(invisible(0L))

  scaffolding <- c("RUN_INFO.txt", "resolved_config.yml", "sessionInfo.txt",
                   "REPORT.txt", "FIGURE_LEGENDS.txt", "TABLE_MANIFEST.txt")
  removed <- 0L

  for (d in list.dirs(output_dir, recursive = FALSE)) {
    info <- file.path(d, "RUN_INFO.txt")
    if (!file.exists(info)) next
    first <- tryCatch(readLines(info, n = 1L, warn = FALSE),
                      error = function(e) character(0))
    if (!length(first) || !identical(first[1], RUN_MARKER)) next

    files <- list.files(d, recursive = TRUE)
    results <- setdiff(files, scaffolding)
    if (!length(results)) {
      unlink(d, recursive = TRUE)
      removed <- removed + 1L
    }
  }
  if (removed) {
    message(sprintf("Removed %d earlier run director%s that produced no output.",
                    removed, if (removed == 1L) "y" else "ies"))
  }
  invisible(removed)
}

#' First line of RUN_INFO.txt, marking a directory as created by this code
RUN_MARKER <- "# immune-repertoire-pipelines run directory"

#' Write sessionInfo() to the run directory
#'
#' Called at the end of every pipeline. Package versions move; a figure that
#' cannot be tied to the versions that produced it cannot be reproduced.
finalise_run <- function(cfg, run_dir) {
  if (isTRUE(cfg$output$write_session_info)) {
    si <- utils::capture.output(utils::sessionInfo())
    writeLines(si, file.path(run_dir, "sessionInfo.txt"))
  }
  cat(sprintf("\nFinished: %s\nOutputs in: %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), run_dir))
  invisible(run_dir)
}

#' Save a ggplot in every configured format
#'
#' @param plot A ggplot object.
#' @param name File stem, without extension.
#' @param run_dir Run directory from `init_run()`.
#' @param cfg Configuration list.
#' @param width,height Size in inches.
#' @param legend Optional figure legend text, appended to FIGURE_LEGENDS.txt.
#'   Writing the legend beside the figure, at the moment the figure is made,
#'   is the only reliable way to keep the two consistent.
save_figure <- function(plot, name, run_dir, cfg,
                        width = 7, height = 5, legend = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }
  fig_dir <- file.path(run_dir, "figures")

  for (fmt in cfg$output$formats) {
    path <- file.path(fig_dir, paste0(name, ".", fmt))
    ok <- try(
      ggplot2::ggsave(path, plot = plot, width = width, height = height,
                      dpi = cfg$output$dpi, units = "in", device = fmt,
                      limitsize = FALSE),
      silent = TRUE
    )
    # PDF devices are the usual failure point on headless machines; a figure
    # that cannot be written as PDF should still be written as PNG rather
    # than aborting the run.
    if (inherits(ok, "try-error")) {
      warning("Could not write ", path, ": ",
              conditionMessage(attr(ok, "condition")), call. = FALSE)
    }
  }

  if (!is.null(legend)) {
    cat(sprintf("%s\n%s\n\n", name, legend),
        file = file.path(run_dir, "FIGURE_LEGENDS.txt"), append = TRUE)
  }
  invisible(plot)
}

#' Write a table to the run directory
#'
#' @param x A data frame.
#' @param name File stem, without extension.
#' @param run_dir Run directory.
#' @param description Optional one-line description, recorded in
#'   TABLE_MANIFEST.txt so a directory of CSVs stays readable.
save_table <- function(x, name, run_dir, description = NULL) {
  stopifnot(is.data.frame(x))
  path <- file.path(run_dir, "tables", paste0(name, ".csv"))
  utils::write.csv(x, path, row.names = FALSE, na = "")

  cat(sprintf("%-50s %6d rows  %s\n", paste0(name, ".csv"), nrow(x),
              description %||% ""),
      file = file.path(run_dir, "TABLE_MANIFEST.txt"), append = TRUE)
  invisible(path)
}

#' Open a plain-text report and return an appender
#'
#' Numeric results quoted in a manuscript should come from a file that the
#' analysis wrote, not from a console session that has since been lost. This
#' returns a small set of writers used by the pipelines to build such a file.
#'
#' @return A list of functions: `h1`, `h2`, `p`, `tbl`, `close`.
open_report <- function(path, title) {
  con <- file(path, open = "wt", encoding = "UTF-8")

  # writeLines(useBytes = TRUE) rather than cat(): under a C/POSIX locale,
  # cat() re-encodes and warns on any non-ASCII character, and reports carry
  # gene names and free text. Note that cat() has no useBytes argument, so
  # it cannot be made locale-safe this way.
  w <- function(...) writeLines(paste0(...), con = con, useBytes = TRUE)

  w(strrep("=", 78))
  w(title)
  w(sprintf("generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
  w(strrep("=", 78))

  list(
    h1 = function(x) { w(""); w(strrep("=", 78)); w(toupper(x)); w(strrep("=", 78)) },
    h2 = function(x) { w(""); w(x); w(strrep("-", nchar(x))) },
    p  = function(...) w(...),
    tbl = function(df, digits = 3) {
      if (!is.data.frame(df) || !nrow(df)) { w("(no rows)"); return(invisible()) }
      num <- vapply(df, is.numeric, logical(1))
      df[num] <- lapply(df[num], round, digits = digits)
      out <- utils::capture.output(print(df, row.names = FALSE))
      for (ln in out) w(ln)
    },
    close = function() { w(""); w(strrep("=", 78)); w("end of report"); close(con) }
  )
}

#' Require a package, with an informative message naming the pipeline step
need <- function(pkg, what) {
  missing_pkgs <- pkg[!vapply(pkg, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) {
    stop(what, " requires package(s) not installed: ",
         paste(missing_pkgs, collapse = ", "),
         "\nSee README.md for installation, including the Bioconductor and ",
         "GitHub sources.", call. = FALSE)
  }
  invisible(TRUE)
}
