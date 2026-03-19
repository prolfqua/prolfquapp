#' Preprocess DIANN output and return AnnData
#'
#' Same interface as \code{\link{preprocess_DIANN}} but returns an
#' \code{anndataR::AnnData} object instead of
#' \code{list(lfqdata, protein_annotation)}.
#'
#' The AnnData \code{uns} slot contains three namespaces:
#' \describe{
#'   \item{X_layer_name}{Name of the primary intensity column stored in X}
#'   \item{exploreDE}{Column role metadata compatible with anndata_omics_bridge}
#'   \item{prolfquapp}{Round-trip reconstruction metadata (config, protein annotation)}
#' }
#'
#' @inheritParams preprocess_DIANN
#' @return \code{anndataR::AnnData} object
#' @export
#'
#' @examples
#' \dontrun{
#' x <- get_DIANN_files("inst/application/DIANN/2706527/")
#' annotation <- file.path("inst/application/DIANN/2706527/dataset.csv") |>
#'   readr::read_csv() |>
#'   prolfquapp::read_annotation(QC = TRUE)
#' adata <- preprocess_DIANN_anndata(x$data, x$fasta, annotation)
#' }
preprocess_DIANN_anndata <- function(
  quant_data,
  fasta_file,
  annotation,
  pattern_contaminants = "^zz|^CON|Cont_",
  pattern_decoys = "^REV_|^rev",
  q_value = 0.01,
  hierarchy_depth = 1,
  nr_peptides = 1
) {
  if (!requireNamespace("anndataR", quietly = TRUE)) {
    stop(
      "Package 'anndataR' is required. Install it with: ",
      "install.packages('anndataR')"
    )
  }

  # Delegate to existing preprocessor
  result <- preprocess_DIANN(
    quant_data = quant_data,
    fasta_file = fasta_file,
    annotation = annotation,
    pattern_contaminants = pattern_contaminants,
    pattern_decoys = pattern_decoys,
    q_value = q_value,
    hierarchy_depth = hierarchy_depth,
    nr_peptides = nr_peptides
  )

  # Convert to AnnData with round-trip metadata
  preprocess_anndata_from_lfq(
    result$lfqdata,
    result$protein_annotation,
    source_software = "DIANN"
  )
}


#' Convert LFQData + ProteinAnnotation to AnnData with round-trip metadata
#'
#' Builds an AnnData object from LFQData and ProteinAnnotation, with the
#' \code{prolfquapp} uns namespace for lossless round-trip via
#' \code{\link{LFQData_from_anndata}}.
#'
#' Handles both protein-level (hierarchyDepth=1) and peptide-level
#' (hierarchyDepth=2) data. For peptide-level data, var contains one row per
#' feature (peptide) with protein annotation joined in.
#'
#' @param lfqdata LFQData object
#' @param protAnnot ProteinAnnotation object
#' @param source_software character, name of the source software (e.g. "DIANN")
#' @return \code{anndataR::AnnData} object
#' @export
#'
#' @examples
#' if (requireNamespace("anndataR", quietly = TRUE)) {
#'   # Protein-level
#'   res <- sim_data_protAnnot(Nprot = 10, PROTEIN = TRUE)
#'   adata <- preprocess_anndata_from_lfq(res$lfqdata, res$pannot, "simulated")
#'   adata  # 12 obs x 10 var
#'
#'   # Peptide-level (multiple peptides per protein)
#'   res2 <- sim_data_protAnnot(Nprot = 10, PROTEIN = FALSE)
#'   adata2 <- preprocess_anndata_from_lfq(res2$lfqdata, res2$pannot, "simulated")
#'   adata2  # 12 obs x ~28 var (peptides)
#'
#'   # Round-trip back to LFQData
#'   back <- LFQData_from_anndata(adata2)
#'   back$lfqdata$hierarchy_counts()
#' }
preprocess_anndata_from_lfq <- function(
  lfqdata,
  protAnnot,
  source_software = "unknown"
) {
  if (!requireNamespace("anndataR", quietly = TRUE)) {
    stop(
      "Package 'anndataR' is required. Install it with: ",
      "install.packages('anndataR')"
    )
  }

  config <- lfqdata$config
  hierarchy_keys <- config$hierarchy_keys()

  # Build layers using to_wide; use rowdata for var construction
  value_cols <- config$value_vars()
  layers <- list()
  rowdata <- NULL
  for (val in value_cols) {
    wide <- lfqdata$to_wide(as.matrix = TRUE, value = val)
    mat <- wide$data
    if (is.null(rowdata)) {
      rowdata <- wide$rowdata
    }
    layers[[val]] <- t(mat) # samples × features (rownames = original matrix rownames)
  }

  X <- layers[[1]]
  sample_names <- rownames(X)
  # Use rowdata to create composite feature IDs with config$sep
  rowdata_df <- as.data.frame(rowdata)
  feature_ids <- do.call(
    paste,
    c(rowdata_df[, hierarchy_keys, drop = FALSE], list(sep = config$sep))
  )

  # Rename matrix columns to composite feature IDs
  for (lname in names(layers)) {
    colnames(layers[[lname]]) <- feature_ids
  }
  X <- layers[[1]]

  # Build obs (sample annotation)
  obs_df <- as.data.frame(lfqdata$factors())
  rownames(obs_df) <- obs_df[[config$sampleName]]
  obs_df <- obs_df[sample_names, , drop = FALSE]

  # Build var (feature annotation) from rowdata + protein annotation
  var_df <- rowdata_df[, hierarchy_keys, drop = FALSE]
  rownames(var_df) <- feature_ids

  # Add identification metadata (aggregate per feature)
  ident_cols <- intersect(
    c(config$ident_qValue, config$nr_children),
    colnames(lfqdata$data)
  )
  if (length(ident_cols) > 0) {
    ident_data <- lfqdata$data[, c(hierarchy_keys, ident_cols), drop = FALSE] |>
      as.data.frame() |>
      dplyr::group_by(dplyr::across(dplyr::all_of(hierarchy_keys))) |>
      dplyr::summarize(
        dplyr::across(
          dplyr::all_of(ident_cols),
          ~ stats::median(.x, na.rm = TRUE)
        ),
        .groups = "drop"
      ) |>
      as.data.frame()

    var_df <- dplyr::left_join(var_df, ident_data, by = hierarchy_keys)
    rownames(var_df) <- feature_ids
  }

  # Join protein annotation
  prot_annot_df <- as.data.frame(protAnnot$row_annot)
  pID <- protAnnot$pID
  var_df <- dplyr::left_join(var_df, prot_annot_df, by = pID)
  rownames(var_df) <- feature_ids

  # Remove isotopeLabel from var if present (it's constant)
  var_df[[config$isotopeLabel]] <- NULL

  # Build uns metadata
  uns <- build_prolfquapp_uns(
    config = config,
    protAnnot = protAnnot,
    layer_names = value_cols,
    source_software = source_software
  )

  adata <- anndataR::AnnData(
    X = X,
    var = var_df,
    obs = obs_df,
    layers = layers,
    uns = uns
  )

  return(adata)
}


#' Build the prolfquapp uns namespace for round-trip AnnData conversion
#'
#' Serializes \code{AnalysisConfiguration} and \code{ProteinAnnotation} fields
#' into a nested list suitable for storage in the AnnData \code{uns} slot.
#'
#' @param config AnalysisConfiguration object
#' @param protAnnot ProteinAnnotation object
#' @param layer_names character vector of layer names
#' @param source_software character, name of the source software
#' @return list to be assigned to \code{adata$uns}
#' @keywords internal
build_prolfquapp_uns <- function(
  config,
  protAnnot,
  layer_names,
  source_software = "unknown"
) {
  # Primary intensity name
  x_layer_name <- config$get_response()

  # exploreDE namespace (anndata_omics_bridge convention)
  exploreDE <- list(
    column_roles = list(
      var = list(
        description = list(protAnnot$description),
        label = list(protAnnot$cleaned_ids, protAnnot$full_id)
      ),
      obs = list(
        factor = as.list(config$factor_keys()),
        label = list(config$sampleName, config$fileName)
      )
    )
  )

  # prolfquapp namespace for round-trip reconstruction
  prolfquapp_ns <- list(
    schema_version = "1.0.0",
    source_software = source_software,
    analysis_configuration = list(
      sep = config$sep,
      fileName = config$fileName,
      sampleName = config$sampleName,
      isotopeLabel = config$isotopeLabel,
      ident_qValue = config$ident_qValue,
      ident_Score = config$ident_Score,
      nr_children = config$nr_children,
      workIntensity = config$workIntensity,
      is_response_transformed = config$is_response_transformed,
      factors = config$factors,
      factorDepth = config$factorDepth,
      hierarchy = config$hierarchy,
      hierarchyDepth = config$hierarchyDepth,
      min_peptides_protein = config$min_peptides_protein
    ),
    protein_annotation = list(
      pID = protAnnot$pID,
      full_id = protAnnot$full_id,
      description = protAnnot$description,
      cleaned_ids = protAnnot$cleaned_ids,
      exp_nr_children = protAnnot$exp_nr_children,
      pattern_contaminants = protAnnot$pattern_contaminants,
      pattern_decoys = protAnnot$pattern_decoys
    ),
    layer_names = as.list(layer_names)
  )

  list(
    X_layer_name = x_layer_name,
    exploreDE = exploreDE,
    prolfquapp = prolfquapp_ns
  )
}
