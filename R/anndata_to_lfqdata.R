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
  artifact_type <- pmeta$artifact_type %||% "lfqdata"
  if (!identical(artifact_type, "lfqdata")) {
    stop("LFQData_from_anndata() requires a prolfquapp LFQData artifact; got '", artifact_type, "'.")
  }
  ac <- pmeta$analysis_configuration
  pa <- pmeta$protein_annotation

  config <- prolfqua::AnalysisConfiguration$new()
  config$sep <- ac$sep
  config$file_name <- ac$file_name %||% ac$fileName
  config$sample_name <- ac$sample_name %||% ac$sampleName
  config$isotope_label <- ac$isotope_label %||% ac$isotopeLabel
  config$ident_q_value <- ac$ident_q_value %||% ac$ident_qValue
  config$ident_score <- (ac$ident_score %||% ac$ident_Score) %||% character()
  config$nr_children <- ac$nr_children
  config$is_response_transformed <- ac$is_response_transformed
  config$factors <- ac$factors
  config$factor_depth <- ac$factor_depth %||% ac$factorDepth
  config$hierarchy <- ac$hierarchy
  config$hierarchy_depth <- ac$hierarchy_depth %||% ac$hierarchyDepth
  config$min_peptides_protein <- ac$min_peptides_protein
  for (wi in (ac$work_intensity %||% ac$workIntensity)) {
    config$set_response(wi)
  }

  lfqdata <- prolfqua::LFQData$new(anndata_to_long(adata, config), config)
  protAnnot <- prolfquapp::ProteinAnnotation$new(
    lfqdata,
    as.data.frame(adata$var),
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
  pmeta <- adata$uns[["prolfquapp"]]
  if (is.null(pmeta)) {
    stop("AnnData uns is missing the 'prolfquapp' namespace. This AnnData was not created by prolfquapp.")
  }
  artifact_type <- pmeta$artifact_type %||% "lfqdata"
  common <- c("schema_version", "source_software", "analysis_configuration")
  required <- list(
    lfqdata = c(common, "protein_annotation"),
    dea_results = c("artifact_type", common, "layer_names", "contrasts", "formula", "provenance")
  )[[artifact_type]]
  if (is.null(required)) {
    stop("Unsupported prolfquapp AnnData artifact type: ", artifact_type)
  }
  missing <- setdiff(required, names(pmeta))
  if (length(missing) > 0) {
    stop("prolfquapp uns is missing required keys: ", paste(missing, collapse = ", "))
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
#' @noRd
anndata_to_long <- function(adata, config) {
  obs_df <- as.data.frame(adata$obs)
  var_df <- as.data.frame(adata$var)
  iso_col <- config$isotope_label
  melt <- function(mat, value_name) {
    dimnames(mat) <- list(rownames(obs_df), rownames(var_df))
    long <- as.data.frame(mat, check.names = FALSE)
    long[[config$sample_name]] <- rownames(obs_df)
    tidyr::pivot_longer(
      long,
      cols = -dplyr::all_of(config$sample_name),
      names_to = ".feature_id",
      values_to = value_name
    )
  }

  long_data <- melt(adata$X, config$get_response())
  for (lname in adata$uns[["prolfquapp"]]$layer_names) {
    if (lname != config$get_response() && !is.null(adata$layers[[lname]])) {
      layer_long <- melt(adata$layers[[lname]], lname)
      long_data <- dplyr::left_join(long_data, layer_long, by = c(config$sample_name, ".feature_id"))
    }
  }

  # Join hierarchy keys, isotope label and identification columns from var
  var_cols <- intersect(c(config$hierarchy_keys(), iso_col, config$ident_q_value, config$nr_children), colnames(var_df))
  var_df$.feature_id <- rownames(var_df)
  long_data <- dplyr::left_join(long_data, var_df[, c(".feature_id", var_cols), drop = FALSE], by = ".feature_id")
  long_data$.feature_id <- NULL

  # Join obs metadata (factors, fileName); drop a duplicated sampleName column
  obs_cols <- intersect(c(config$file_name, config$factor_keys(), iso_col), colnames(obs_df))
  obs_join <- obs_df[, c(config$sample_name, obs_cols), drop = FALSE]
  obs_join <- obs_join[, !duplicated(colnames(obs_join)), drop = FALSE]
  long_data <- dplyr::left_join(long_data, obs_join, by = config$sample_name)

  if (!iso_col %in% colnames(long_data)) {
    long_data[[iso_col]] <- "light"
  }
  tibble::as_tibble(long_data)
}
