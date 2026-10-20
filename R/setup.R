# ---------------------------------------------------------------------------
# Shared pipeline preamble: locate the repository root and source R/
#
# Every pipeline begins with `source("R/setup.R")` after setting ROOT, or by
# copying the root-finding block below. Kept in one place so the four
# pipelines cannot drift in how they bootstrap.
# ---------------------------------------------------------------------------

#' Locate the repository root
#'
#' Works whether the pipeline is launched from the repository directory, from
#' pipelines/, or by absolute path (Rscript sets --file=).
find_repo_root <- function() {
  candidates <- getwd()
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) {
    candidates <- c(dirname(dirname(normalizePath(f[1]))), candidates)
  }
  for (start in candidates) {
    d <- normalizePath(start, mustWork = FALSE)
    for (i in 1:6) {
      if (dir.exists(file.path(d, "R")) &&
          file.exists(file.path(d, "config", "config.example.yml"))) {
        return(d)
      }
      if (identical(dirname(d), d)) break
      d <- dirname(d)
    }
  }
  stop("Could not locate the repository root. Run from the repository ",
       "directory, e.g.  Rscript pipelines/01_bulk_repertoire.R", call. = FALSE)
}

#' Source every module under R/, in dependency order
source_modules <- function(root) {
  order <- c("R/core", "R/repertoire", "R/analysis", "R/bcr")
  for (dir in order) {
    files <- list.files(file.path(root, dir), pattern = "[.]R$",
                        full.names = TRUE)
    for (f in files) source(f)
  }
}
