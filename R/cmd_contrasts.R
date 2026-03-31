#' Generate contrasts for a single-factor design
#'
#' Reads an annotation file, validates the control level exists,
#' and adds a CONTROL column (C for control, T for treatment).
#'
#' @param annotation_file path to annotation CSV/TSV/XLSX
#' @param control reference level name (e.g. "WT")
#' @param group group column name, or NULL to auto-detect
#' @return data.frame with CONTROL column added
#' @export
#' @examples
#' csv <- system.file("application/contrasts/scenario1_single_factor.csv",
#'   package = "prolfquapp")
#' result <- run_contrasts_single(csv, control = "WT")
#' table(result$group, result$CONTROL)
#'
run_contrasts_single <- function(annotation_file, control, group = NULL) {
  stopifnot(file.exists(annotation_file))

  res <- prolfquapp::read_annotation(annotation_file, QC = TRUE)
  annot <- res$annot
  group_col <- if (!is.null(group)) group else res$atable$factors[["G_"]]

  if (!control %in% annot[[group_col]]) {
    stop(
      "control '", control, "' not found in column '", group_col, "'.",
      call. = FALSE
    )
  }

  annot$CONTROL <- ifelse(annot[[group_col]] == control, "C", "T")
  annot
}

#' Generate contrasts for a two-factor design
#'
#' Reads an annotation file and adds ContrastName/Contrast columns
#' using \code{\link[prolfqua]{annotation_add_contrasts}}.
#'
#' @param annotation_file path to annotation CSV/TSV/XLSX
#' @param f1 primary factor column name
#' @param f2 secondary factor column name
#' @param interactions logical; include interaction contrasts? Default TRUE.
#' @return data.frame with ContrastName and Contrast columns added
#' @export
#' @examples
#' csv <- system.file("application/contrasts/scenario2_two_factor.csv",
#'   package = "prolfquapp")
#' result <- run_contrasts_twofactor(csv, f1 = "treatment", f2 = "time")
#' unique(result[!is.na(result$ContrastName), c("ContrastName", "Contrast")])
#'
run_contrasts_twofactor <- function(
  annotation_file,
  f1,
  f2,
  interactions = TRUE
) {
  stopifnot(file.exists(annotation_file))

  df <- prolfquapp::read_table_data(annotation_file)
  missing_cols <- setdiff(c(f1, f2), colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      "Column(s) not found: ", paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  res <- prolfqua::annotation_add_contrasts(
    df,
    primary_col = f1,
    secondary_col = f2,
    interactions = interactions
  )
  res$annot
}
