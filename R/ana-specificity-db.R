# ---------------------------------------------------------------------------
# Reference tables of receptors with a known specificity
#
# Public receptor databases are read into one standard layout, so that the
# matching code downstream never depends on which database, or which export
# of it, the user downloaded:
#
#   cdr3, chain, v_gene, j_gene, epitope, epitope_gene, epitope_species,
#   mhc, tissue, cell_type, source
#
# The files are downloaded by the user once (see the pipelines' data
# README); nothing here reaches the network at run time, so an annotation
# is reproducible from the archived file and its release date.
#
#   VDJdb      Shugay M et al. Nucleic Acids Res. 2018;46(D1):D419-D427.
#              doi:10.1093/nar/gkx760. PMID 28977646.
#              Goncharov M et al. Nat Methods. 2022;19(9):1017-1019.
#              doi:10.1038/s41592-022-01578-0. PMID 35970936.
#   McPAS-TCR  Tickotsky N et al. Bioinformatics. 2017;33(18):2924-2929.
#              doi:10.1093/bioinformatics/btx286. PMID 28481982.
#   Trex       the reference bundled with the Trex package (compiled from
#              VDJdb, McPAS-TCR, IEDB and PIRD), one row per CDR3.
#
# COVERAGE. These databases are dominated by a few intensively studied
# pathogens (CMV, EBV, influenza, SARS-CoV-2). An unannotated receptor is
# therefore not evidence of anything: most receptors specific for a less
# studied antigen are simply absent from the references.
# ---------------------------------------------------------------------------

#' Standard columns of a reference table
#' @noRd
reference_columns <- function() {
  c("cdr3", "chain", "v_gene", "j_gene", "epitope", "epitope_gene",
    "epitope_species", "mhc", "tissue", "cell_type", "source")
}

#' Map a raw table onto the standard reference columns
#'
#' @param x Raw data frame.
#' @param layouts List of named character vectors, standard name -> raw
#'   column, tried in order; the first whose `cdr3` and `epitope_species`
#'   columns are both present is used.
#' @param columns Optional named overrides, applied on top of the layout.
#' @param source Database label written in the `source` column.
#' @noRd
standardise_reference <- function(x, layouts, columns = NULL, source,
                                  what = "reference") {
  chosen <- NULL
  for (lay in layouts) {
    if (all(c(lay[["cdr3"]], lay[["epitope_species"]]) %in% names(x))) {
      chosen <- lay
      break
    }
  }
  if (is.null(chosen)) chosen <- layouts[[1]]
  if (length(columns)) chosen[names(columns)] <- unlist(columns)

  miss <- setdiff(c(chosen[["cdr3"]], chosen[["epitope_species"]]), names(x))
  if (length(miss)) {
    stop(what, ": column(s) not found: ", paste(miss, collapse = ", "),
         "\nColumns present: ", paste(utils::head(names(x), 30),
                                      collapse = ", "),
         "\nName them with the `columns` argument (or the configuration) ",
         "if the export uses another layout.", call. = FALSE)
  }

  out <- data.frame(row.names = seq_len(nrow(x)))
  for (std in reference_columns()) {
    raw <- chosen[std]
    out[[std]] <- if (!is.na(raw) && raw %in% names(x)) {
      as.character(x[[raw]])
    } else NA_character_
  }
  # A layout may map `source` onto a column (Trex records the contributing
  # database per row); otherwise the database label is written.
  if (!is.na(source)) out$source <- source
  for (std in reference_columns()) {
    v <- trimws(out[[std]])
    v[v %in% c("", "NA")] <- NA_character_
    out[[std]] <- v
  }
  out <- out[!is.na(out$cdr3), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Read a VDJdb file into the standard reference layout
#'
#' Two layouts are in circulation and both are recognised: the release files
#' (`vdjdb.txt`, `vdjdb.slim.txt`, with lower-case dotted names such as
#' `cdr3` and `antigen.species`) and the export of the web search page
#' (`CDR3`, `Epitope species`, ...). In the release layout the column
#' `species` is the HOST species of the receptor (HomoSapiens), not the
#' species of the antigen; it is used to filter hosts and is never read as
#' the antigen species.
#'
#' The file is read without quote processing: the web export carries JSON
#' in its `Method` and `Meta` columns, whose double quotes would otherwise
#' be taken as field delimiters and shift every following column.
#'
#' @param path Path to the TSV file.
#' @param chain Receptor chain to keep ("TRB", "TRA"); NULL keeps all.
#' @param host_species Host species to keep; NULL keeps all.
#' @param columns Optional named character vector overriding column names
#'   (names among cdr3, chain, v_gene, j_gene, epitope, epitope_gene,
#'   epitope_species, mhc, host).
#' @return Data frame in the standard reference layout, `source = "VDJdb"`.
#' @export
read_vdjdb <- function(path, chain = "TRB", host_species = "HomoSapiens",
                       columns = NULL) {
  if (!file.exists(path)) stop("VDJdb file not found: ", path, call. = FALSE)
  x <- utils::read.delim(path, quote = "", comment.char = "",
                         check.names = FALSE, stringsAsFactors = FALSE,
                         na.strings = c("", "NA"), colClasses = "character")
  layouts <- list(
    release = c(cdr3 = "cdr3", chain = "gene", host = "species",
                v_gene = "v.segm", j_gene = "j.segm",
                epitope = "antigen.epitope", epitope_gene = "antigen.gene",
                epitope_species = "antigen.species", mhc = "mhc.a"),
    web = c(cdr3 = "CDR3", chain = "Gene", host = "Species", v_gene = "V",
            j_gene = "J", epitope = "Epitope", epitope_gene = "Epitope gene",
            epitope_species = "Epitope species", mhc = "MHC A")
  )
  filter_reference_rows(x, layouts, columns, chain, host_species,
                        source = "VDJdb", what = "VDJdb")
}

#' Read a McPAS-TCR file into the standard reference layout
#'
#' The pathology column is taken as the specificity label, as in the
#' original analysis: McPAS-TCR records the pathology a receptor is
#' associated with, which for infections names the pathogen.
#'
#' The published CSV is not valid UTF-8 throughout (some gene names carry
#' Windows-1252 dashes); invalid lines are converted from that encoding so
#' that the file reads under any locale.
#'
#' @inheritParams read_vdjdb
#' @param host_species Host species to keep ("Human"); NULL keeps all.
#' @return Data frame in the standard reference layout, `source = "McPAS"`.
#' @export
read_mcpas <- function(path, chain = "TRB", host_species = "Human",
                       columns = NULL) {
  if (!file.exists(path)) stop("McPAS-TCR file not found: ", path,
                               call. = FALSE)
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  bad <- !validUTF8(lines)
  if (any(bad)) {
    lines[bad] <- iconv(lines[bad], from = "CP1252", to = "UTF-8", sub = "?")
  }
  x <- utils::read.csv(text = lines, check.names = FALSE,
                       stringsAsFactors = FALSE, na.strings = c("", "NA"),
                       colClasses = "character")
  # McPAS-TCR stores both chains of a receptor on one row, so the chain is
  # chosen by the column read rather than by a filter.
  ch <- chain %||% "TRB"
  lay <- if (identical(ch, "TRA")) {
    c(cdr3 = "CDR3.alpha.aa", v_gene = "TRAV", j_gene = "TRAJ")
  } else {
    c(cdr3 = "CDR3.beta.aa", v_gene = "TRBV", j_gene = "TRBJ")
  }
  layouts <- list(c(lay, host = "Species", epitope = "Epitope.peptide",
                    epitope_gene = "Antigen.protein",
                    epitope_species = "Pathology", mhc = "MHC",
                    tissue = "Tissue", cell_type = "T.Cell.Type"))
  out <- filter_reference_rows(x, layouts, columns, chain = NULL,
                               host_species = host_species,
                               source = "McPAS", what = "McPAS-TCR")
  out$chain <- ch
  out
}

#' Shared host and chain filtering for the readers
#' @noRd
filter_reference_rows <- function(x, layouts, columns, chain, host_species,
                                  source, what) {
  lay <- NULL
  for (l in layouts) {
    if (all(c(l[["cdr3"]], l[["epitope_species"]]) %in% names(x))) {
      lay <- l
      break
    }
  }
  if (is.null(lay)) lay <- layouts[[1]]
  if (length(columns)) lay[names(columns)] <- unlist(columns)

  keep <- rep(TRUE, nrow(x))
  if (!is.null(host_species) && !is.na(lay["host"]) &&
      lay[["host"]] %in% names(x)) {
    keep <- keep & x[[lay[["host"]]]] %in% host_species
  }
  if (!is.null(chain) && !is.na(lay["chain"]) &&
      lay[["chain"]] %in% names(x)) {
    keep <- keep & x[[lay[["chain"]]]] %in% chain
  }
  x <- x[which(keep), , drop = FALSE]
  standardise_reference(x, list(lay[setdiff(names(lay), "host")]),
                        columns = NULL, source = source, what = what)
}

#' The reference bundled with the Trex package, in the standard layout
#'
#' Trex ships one table per chain with one row per CDR3; multi-valued
#' fields are joined with ";" by the package and are kept as they are.
#'
#' @param chain "TRB" or "TRA".
#' @return Data frame in the standard reference layout; `source` holds the
#'   contributing database recorded by Trex.
#' @export
trex_reference <- function(chain = "TRB") {
  need("Trex", "trex_reference()")
  env <- new.env()
  utils::data("Trex.database", package = "Trex", envir = env)
  db <- as.data.frame(env$Trex.database[[chain]], stringsAsFactors = FALSE)
  if (!nrow(db)) stop("No ", chain, " table in Trex.database.", call. = FALSE)
  out <- standardise_reference(
    db,
    list(c(cdr3 = "CDR3", v_gene = "V.gene", j_gene = "J.gene",
           epitope = "Epitope.sequence", epitope_gene = "Epitope.target",
           epitope_species = "Epitope.species", mhc = "HLA1",
           tissue = "Tissue", cell_type = "Cell.type", source = "Database")),
    source = NA_character_, what = "Trex.database")
  out$chain <- chain
  out
}

#' Collapse a reference table to one row per CDR3
#'
#' A CDR3 can be recorded against several epitopes, species and studies. The
#' lookup keeps one row per CDR3 with every distinct value joined, so that
#' joining it onto a repertoire never duplicates a clonotype.
#'
#' @param reference Standard reference table (or several, row-bound).
#' @param sep Separator between distinct values; " | " as in the original
#'   analysis.
#' @return Data frame: cdr3, specificities, epitopes, epitope_genes, sources,
#'   tissue, cell_type, n_records.
#' @export
collapse_reference <- function(reference, sep = " | ") {
  need("dplyr", "collapse_reference()")
  if (!nrow(reference)) {
    return(data.frame(cdr3 = character(0), specificities = character(0),
                      epitopes = character(0), epitope_genes = character(0),
                      sources = character(0), tissue = character(0),
                      cell_type = character(0), n_records = integer(0),
                      stringsAsFactors = FALSE))
  }
  join <- function(v, s = sep) {
    v <- unique(v[!is.na(v) & nzchar(v)])
    if (length(v)) paste(v, collapse = s) else NA_character_
  }
  out <- dplyr::summarise(
    dplyr::group_by(reference, .data$cdr3),
    specificities = join(.data$epitope_species),
    epitopes = join(.data$epitope),
    epitope_genes = join(.data$epitope_gene),
    sources = join(.data$source, ","),
    tissue = join(.data$tissue),
    cell_type = join(.data$cell_type),
    n_records = dplyr::n(),
    .groups = "drop")
  as.data.frame(out, stringsAsFactors = FALSE)
}

#' Count reference records per antigen species
#'
#' The first thing to look at before annotating: which species the
#' reference can recognise at all, and how unevenly they are represented.
#'
#' @param reference Standard reference table.
#' @return Data frame: source, epitope_species, n_records, n_cdr3.
#' @export
reference_species_overview <- function(reference) {
  need("dplyr", "reference_species_overview()")
  out <- dplyr::summarise(
    dplyr::group_by(reference, .data$source, .data$epitope_species),
    n_records = dplyr::n(),
    n_cdr3 = dplyr::n_distinct(.data$cdr3),
    .groups = "drop")
  out <- out[order(out$source, -out$n_records), , drop = FALSE]
  as.data.frame(out, stringsAsFactors = FALSE)
}

#' Read a reference table of antibodies with a known target
#'
#' Columns are found as in the original analysis: the first column is the
#' antibody name, the heavy-chain V gene is the first column matching
#' "Heavy.*V.*Gene" and the CDRH3 the first matching "CDRH3", each
#' overridable. This fits the CSV export of CoV-AbDab (Raybould MIJ et al.
#' Bioinformatics. 2021;37(5):734-735. doi:10.1093/bioinformatics/btaa739.
#' PMID 32805021).
#'
#' The V gene is reduced to "IGHV<family>-<number>" by the original
#' pattern. Entries whose CDRH3 is missing or recorded as not determined
#' ("ND") are dropped.
#'
#' @param path Path to the CSV file.
#' @param columns Optional named overrides: name, v_gene, cdr3.
#' @return Data frame: ab_name, ab_vgene, ab_cdr3.
#' @export
read_antibody_reference <- function(path, columns = NULL) {
  if (!file.exists(path)) stop("Antibody reference not found: ", path,
                               call. = FALSE)
  ab <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                        colClasses = "character", na.strings = c("", "NA"))
  pick <- function(key, pattern) {
    if (!is.null(columns[[key]])) return(columns[[key]])
    hit <- grep(pattern, names(ab), value = TRUE, ignore.case = TRUE)
    if (length(hit)) hit[1] else NA_character_
  }
  name_col <- columns[["name"]] %||% names(ab)[1]
  v_col <- pick("v_gene", "Heavy.*V.*Gene")
  cdr3_col <- pick("cdr3", "CDRH3")
  if (is.na(v_col) || is.na(cdr3_col) ||
      !all(c(name_col, v_col, cdr3_col) %in% names(ab))) {
    stop("Antibody reference: heavy-chain V gene or CDRH3 column not ",
         "recognised. Columns present: ",
         paste(utils::head(names(ab), 30), collapse = ", "),
         "\nName them with the `columns` argument.", call. = FALSE)
  }
  out <- data.frame(ab_name = ab[[name_col]],
                    ab_vgene = ighv_gene(ab[[v_col]]),
                    ab_cdr3 = toupper(trimws(ab[[cdr3_col]])),
                    stringsAsFactors = FALSE)
  out[!is.na(out$ab_cdr3) & nzchar(out$ab_cdr3) & out$ab_cdr3 != "ND", ,
      drop = FALSE]
}

#' Reduce a heavy-chain V call to IGHV<family>-<number>
#'
#' The pattern of the original analysis. It is applied to the query and to
#' the reference alike, so that a call such as IGHV3-30-3 and a reference
#' entry recorded as IGHV3-30 are compared at the same resolution.
#' @noRd
ighv_gene <- function(x) {
  x <- as.character(x)
  pat <- "IGHV[0-9]+-[0-9A-Za-z]+"
  out <- rep(NA_character_, length(x))
  hit <- !is.na(x) & regexpr(pat, x) > 0
  out[hit] <- regmatches(x[hit], regexpr(pat, x[hit]))
  out
}
