#' Convert AnnData back to LFQData + ProteinAnnotation
#'
#' Reconstructs the \code{list(lfqdata, protein_annotation)} pair from
#' an AnnData object that was created by \code{\link{preprocess_DIANN_anndata}}
#' (or any producer that writes the \code{prolfquapp} uns namespace).
#'
#' @param adata an \code{anndataR::AnnData} object with \code{uns$prolfquapp}
#' @return list with \code{lfqdata} (LFQData) and \code{protein_annotation}
#'   (ProteinAnnotation)
#' @export
#'
#' @examples
#' \dontrun{
#' res <- sim_data_protAnnot()
#' adata <- preprocess_DIANN_anndata_from_lfq(res$lfqdata, res$pannot)
#' back <- LFQData_from_anndata(adata)
#' }
LFQData_from_anndata <- function(adata) {
  validate_prolfquapp_anndata(adata)

  pmeta <- adata$uns[["prolfquapp"]]
  ac <- pmeta$analysis_configuration
  pa <- pmeta$protein_annotation

  # --- rebuild AnalysisConfiguration ---
  config <- prolfqua::AnalysisConfiguration$new()
  config$sep <- ac$sep
  config$fileName <- ac$fileName
  config$sampleName <- ac$sampleName
  config$isotopeLabel <- ac$isotopeLabel
  config$ident_qValue <- ac$ident_qValue
  config$ident_Score <- ac$ident_Score %||% character()
  config$nr_children <- ac$nr_children
  config$is_response_transformed <- ac$is_response_transformed
  config$factors <- ac$factors
  config$factorDepth <- ac$factorDepth
  config$hierarchy <- ac$hierarchy
  config$hierarchyDepth <- ac$hierarchyDepth
  config$min_peptides_protein <- ac$min_peptides_protein

  for (wi in ac$workIntensity) {
    config$set_response(wi)
  }

  # --- convert wide → long ---
  long_data <- anndata_to_long(adata, config)

  lfqdata <- prolfqua::LFQData$new(long_data, config)

  # --- rebuild ProteinAnnotation ---
  var_df <- as.data.frame(adata$var)

  protAnnot <- prolfquapp::ProteinAnnotation$new(
    lfqdata,
    var_df,
    description = pa$description,
    cleaned_ids = pa$cleaned_ids,
    full_id = pa$full_id,
    exp_nr_children = pa$exp_nr_children,
    pattern_contaminants = pa$pattern_contaminants,
    pattern_decoys = pa$pattern_decoys
  )

  return(list(lfqdata = lfqdata, protein_annotation = protAnnot))
}


#' Validate that an AnnData has the prolfquapp uns namespace
#'
#' @param adata an AnnData object
#' @export
validate_prolfquapp_anndata <- function(adata) {
  if (is.null(adata$uns)) {
    stop("AnnData has no 'uns' slot.")
  }
  if (is.null(adata$uns[["prolfquapp"]])) {
    stop(
      "AnnData uns is missing the 'prolfquapp' namespace. ",
      "This AnnData was not created by prolfquapp."
    )
  }
  pmeta <- adata$uns[["prolfquapp"]]
  required <- c(
    "schema_version",
    "source_software",
    "analysis_configuration",
    "protein_annotation"
  )
  missing <- setdiff(required, names(pmeta))
  if (length(missing) > 0) {
    stop(
      "prolfquapp uns is missing required keys: ",
      paste(missing, collapse = ", ")
    )
  }
  invisible(TRUE)
}


#' Convert AnnData wide matrix to long-format tibble
#'
#' Internal helper. Melts the X matrix and any additional layers,
#' joins with obs (sample factors) and var (hierarchy columns).
#'
#' @param adata AnnData object
#' @param config AnalysisConfiguration
#' @return tibble in long format
#' @keywords internal
anndata_to_long <- function(adata, config) {
  obs_df <- as.data.frame(adata$obs)
  var_df <- as.data.frame(adata$var)

  # Determine which columns from var are hierarchy keys
  hierarchy_cols <- config$hierarchy_keys()
  # isotopeLabel column
  iso_col <- config$isotopeLabel

  # Get sample names (row names of X = obs rownames)
  sample_names <- rownames(obs_df)
  # Get feature IDs (column names of X = var rownames)
  feature_ids <- rownames(var_df)

  # Primary intensity layer
  X <- adata$X
  rownames(X) <- sample_names
  colnames(X) <- feature_ids

  # Melt X to long format: sampleName, featureID, value
  X_long <- as.data.frame(X, check.names = FALSE)
  X_long[[config$sampleName]] <- sample_names
  X_long <- tidyr::pivot_longer(
    X_long,
    cols = -dplyr::all_of(config$sampleName),
    names_to = ".feature_id",
    values_to = config$get_response()
  )

  # Add additional layers as columns
  layer_names <- adata$uns[["prolfquapp"]]$layer_names
  if (!is.null(layer_names)) {
    for (lname in layer_names) {
      if (lname == config$get_response()) {
        next
      }
      layer_mat <- adata$layers[[lname]]
      if (is.null(layer_mat)) {
        next
      }
      rownames(layer_mat) <- sample_names
      colnames(layer_mat) <- feature_ids
      layer_long <- as.data.frame(layer_mat, check.names = FALSE)
      layer_long[[config$sampleName]] <- sample_names
      layer_long <- tidyr::pivot_longer(
        layer_long,
        cols = -dplyr::all_of(config$sampleName),
        names_to = ".feature_id",
        values_to = lname
      )
      X_long <- dplyr::left_join(
        X_long,
        layer_long,
        by = c(config$sampleName, ".feature_id")
      )
    }
  }

  # Build var lookup: feature_id → hierarchy columns + metadata columns
  var_lookup <- var_df
  var_lookup$.feature_id <- feature_ids

  # Select hierarchy keys + isotopeLabel + identification columns from var
  var_cols_to_join <- intersect(
    c(hierarchy_cols, iso_col, config$ident_qValue, config$nr_children),
    colnames(var_lookup)
  )
  var_join <- var_lookup[, c(".feature_id", var_cols_to_join), drop = FALSE]

  # Join var metadata onto long data
  long_data <- dplyr::left_join(X_long, var_join, by = ".feature_id")
  long_data$.feature_id <- NULL

  # Join obs metadata (factors, fileName)
  obs_cols_to_join <- intersect(
    c(config$fileName, config$factor_keys(), iso_col),
    colnames(obs_df)
  )
  obs_join <- obs_df[, c(config$sampleName, obs_cols_to_join), drop = FALSE]
  # Deduplicate in case sampleName is already in obs_cols_to_join
  obs_join <- obs_join[, !duplicated(colnames(obs_join)), drop = FALSE]

  long_data <- dplyr::left_join(long_data, obs_join, by = config$sampleName)

  # Ensure isotopeLabel exists

  if (!iso_col %in% colnames(long_data)) {
    long_data[[iso_col]] <- "light"
  }

  tibble::as_tibble(long_data)
}
