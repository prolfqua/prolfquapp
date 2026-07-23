#' Drop proteins supported by fewer than `nr_peptides` distinct peptides
#'
#' Reader-local minimum-peptides-per-protein filter. Counts the number of
#' distinct `peptide_col` values per `protein_col` in `data` and keeps only the
#' rows whose parent protein reaches `nr_peptides`. The reader supplies the
#' columns because only the reader knows which column holds the stripped
#' (unmodified) peptide sequence.
#'
#' Filtering the quant/peptide table is sufficient: `ProteinAnnotation` is
#' right-joined onto the `LFQData` protein set (it seeds `row_annot` from the
#' `LFQData` proteins and left-joins annotation onto it), so the annotation,
#' the logged protein summary, and IBAQ all follow the filtered quant
#' automatically -- see the "ProteinAnnotation is annotation, not a filter"
#' invariant in `AGENTS.md`. Apply this to the long table *before*
#' `setup_analysis()` / `LFQData$new()`.
#'
#' A no-op when `nr_peptides <= 1` (or `NULL`), so default runs are unaffected.
#'
#' @param data long-format quant table (before `setup_analysis`)
#' @param protein_col name of the parent protein column (single column)
#' @param peptide_col name of the (stripped) peptide column to count (single column)
#' @param nr_peptides minimum number of distinct peptides per protein (>= 1)
#' @return `data` with the rows of under-supported proteins removed
#' @export
#' @family preprocessing
#' @examples
#' d <- data.frame(
#'   prot = c("A", "A", "A", "B", "C", "C"),
#'   pep = c("p1", "p2", "p2", "p1", "p1", "p2")
#' )
#' # A has 2 distinct peptides, B has 1, C has 2
#' nrow(filter_by_peptide_count(d, "prot", "pep", 2)) # drops B
#' nrow(filter_by_peptide_count(d, "prot", "pep", 1)) # keeps all
filter_by_peptide_count <- function(data, protein_col, peptide_col, nr_peptides = 1) {
  if (is.null(nr_peptides) || nr_peptides <= 1) {
    return(data)
  }
  stopifnot(
    length(protein_col) == 1L,
    length(peptide_col) == 1L,
    protein_col %in% colnames(data),
    peptide_col %in% colnames(data)
  )
  counts <- data |>
    dplyr::distinct(dplyr::across(dplyr::all_of(c(protein_col, peptide_col)))) |>
    dplyr::count(dplyr::across(dplyr::all_of(protein_col)), name = ".nr_peptides")
  keep <- counts[[protein_col]][counts[[".nr_peptides"]] >= nr_peptides]

  n_total <- nrow(counts)
  n_kept <- length(keep)
  logger::log_info(
    "nr_peptides filter (>= ",
    nr_peptides,
    " distinct ",
    peptide_col,
    " per ",
    protein_col,
    "): kept ",
    n_kept,
    " / ",
    n_total,
    " proteins (dropped ",
    n_total - n_kept,
    ")"
  )

  data[data[[protein_col]] %in% keep, , drop = FALSE]
}


#' Validate and normalize an `nr_peptides` threshold
#'
#' Central guard for the minimum-peptides-per-protein threshold, applied where
#' the value enters a `ProlfquAppConfig` (native YAML, programmatic config, and
#' CLI override). Accepts clean YAML numerics (`2`, `2.0`) and numeric-looking
#' strings (`"2"`); rejects values that are not a single whole number `>= 1`.
#'
#' @param x candidate value (numeric, integer, or character scalar); `NULL`
#'   defaults to `1L`
#' @return the value coerced to a positive whole-number integer
#' @keywords internal
#' @noRd
.validate_nr_peptides <- function(x) {
  if (is.null(x)) {
    return(1L)
  }
  if (length(x) != 1L) {
    stop("nr_peptides must be a single whole number >= 1", call. = FALSE)
  }
  num <- suppressWarnings(as.numeric(x))
  if (
    !is.finite(num) ||
      num < 1 ||
      num != round(num) ||
      num > .Machine$integer.max
  ) {
    stop(
      "nr_peptides must be a whole number >= 1, got: ",
      format(x),
      call. = FALSE
    )
  }
  as.integer(round(num))
}


#' Add `nr_peptides` to reader args only for readers that declare it
#'
#' Each reader owns the minimum-peptides-per-protein count against its own
#' stripped-peptide column, so the threshold is forwarded only to readers whose
#' formals declare an `nr_peptides` argument. A reader without it is left
#' unfiltered; a warning is raised when the caller asked for filtering
#' (`nr_peptides > 1`) so the skip is explicit rather than a silent no-op.
#'
#' @param base_args named list of arguments passed to the reader
#' @param preprocess_fn the resolved reader function
#' @param nr_peptides minimum distinct peptides per protein
#' @return `base_args`, with `nr_peptides` added iff the reader supports it
#' @keywords internal
#' @noRd
.forward_nr_peptides <- function(base_args, preprocess_fn, nr_peptides) {
  if ("nr_peptides" %in% names(formals(preprocess_fn))) {
    base_args$nr_peptides <- nr_peptides
  } else if (!is.null(nr_peptides) && nr_peptides > 1) {
    warning(
      "reader does not support nr_peptides filtering; ignoring nr_peptides = ",
      nr_peptides,
      call. = FALSE
    )
  }
  base_args
}
