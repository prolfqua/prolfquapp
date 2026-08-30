.dea_anndata_required_assays <- c("rawData", "transformedData")

.dea_anndata_result_markers <- c(
  "modelName",
  "estimate_type",
  "contrast",
  "Bait",
  "diff",
  "statistic",
  "p.value",
  "FDR",
  "BFDR"
)

.summarized_experiment_axis_names <- function(se, axis) {
  values <- if (identical(axis, "obs")) colnames(se) else rownames(se)
  if (is.null(values) || anyNA(values) || any(!nzchar(values))) {
    stop("SummarizedExperiment ", axis, " names must be present and non-empty.")
  }
  if (anyDuplicated(values) > 0L) {
    stop("SummarizedExperiment ", axis, " names must be unique.")
  }
  unname(as.character(values))
}

.aligned_row_data_frame <- function(value, var_names, key) {
  value <- as.data.frame(value)
  if (nrow(value) != length(var_names)) {
    stop(
      "SummarizedExperiment rowData '",
      key,
      "' has ",
      nrow(value),
      " rows; expected ",
      length(var_names),
      "."
    )
  }
  if (!identical(unname(rownames(value)), var_names)) {
    stop(
      "SummarizedExperiment rowData '",
      key,
      "' is not aligned to the feature axis."
    )
  }
  value
}

.dea_anndata_feature_table <- function(row_data, contrast_keys, var_names) {
  first_key <- contrast_keys[[1]]
  first_contrast <- .aligned_row_data_frame(
    row_data[[first_key]],
    var_names,
    first_key
  )
  boundary <- match(
    TRUE,
    names(first_contrast) %in% .dea_anndata_result_markers
  )
  if (is.na(boundary)) {
    stop(
      "SummarizedExperiment rowData '",
      first_key,
      "' does not contain a recognized differential-result column."
    )
  }
  annotation_names <- names(first_contrast)[seq_len(boundary - 1L)]
  result <- first_contrast[, annotation_names, drop = FALSE]
  rownames(result) <- var_names
  result
}

.dea_anndata_hierarchy_keys <- function(se) {
  config <- S4Vectors::metadata(se)$analysis_configuration_transformed
  hierarchy <- config$hierarchy
  if (is.null(hierarchy)) character() else names(hierarchy)
}

.encode_dea_varm_key <- function(contrast_name) {
  if (!nzchar(contrast_name)) {
    stop("Differential-result contrast names must be non-empty.")
  }
  paste0(
    "dea__",
    utils::URLencode(contrast_name, reserved = TRUE, repeated = TRUE)
  )
}

.anndata_uns_value <- function(value) {
  if (is.data.frame(value)) {
    return(lapply(value, .anndata_uns_value))
  }
  if (is.list(value)) {
    return(lapply(value, .anndata_uns_value))
  }
  if (is.factor(value)) {
    return(unname(as.character(value)))
  }
  unname(value)
}

.dea_anndata_varm <- function(se, var, var_names) {
  row_data <- SummarizedExperiment::rowData(se)
  row_data_names <- names(row_data)
  contrast_keys <- grep("^constrast_", row_data_names, value = TRUE)
  hierarchy_keys <- .dea_anndata_hierarchy_keys(se)
  drop_columns <- union(names(var), hierarchy_keys)

  values <- list()
  columns <- list()
  annotations <- list()
  for (legacy_key in contrast_keys) {
    frame <- .aligned_row_data_frame(
      row_data[[legacy_key]],
      var_names,
      legacy_key
    )
    contrast_name <- sub("^constrast_", "", legacy_key)
    key <- .encode_dea_varm_key(contrast_name)
    if (key %in% names(values)) {
      stop("Differential-result names collide after AnnData key encoding.")
    }
    frame <- frame[, setdiff(names(frame), drop_columns), drop = FALSE]
    numeric_columns <- names(frame)[vapply(
      frame,
      function(column) is.numeric(column) || is.logical(column),
      logical(1)
    )]
    values[[key]] <- as.matrix(frame[, numeric_columns, drop = FALSE])
    columns[[key]] <- numeric_columns
    annotation_columns <- setdiff(names(frame), numeric_columns)
    annotations[[key]] <- lapply(
      frame[, annotation_columns, drop = FALSE],
      unname
    )
  }

  statistics_keys <- intersect(
    c("stats_normalized_wide", "stats_raw_wide"),
    row_data_names
  )
  for (key in statistics_keys) {
    frame <- .aligned_row_data_frame(row_data[[key]], var_names, key)
    frame <- frame[, setdiff(names(frame), drop_columns), drop = FALSE]
    if (!all(vapply(frame, is.numeric, logical(1)))) {
      stop("SummarizedExperiment rowData '", key, "' is not numeric.")
    }
    values[[key]] <- as.matrix(frame)
    columns[[key]] <- names(frame)
  }
  list(values = values, columns = columns, annotations = annotations)
}

.dea_anndata_uns <- function(se, layer_names, varm) {
  metadata <- S4Vectors::metadata(se)
  config <- metadata$analysis_configuration_transformed
  provenance <- metadata$report_provenance
  source_software <- provenance$software
  if (is.null(source_software) || length(source_software) != 1L) {
    source_software <- "unknown"
  }

  .anndata_uns_value(list(
    X_layer_name = "transformed",
    prolfquapp = list(
      artifact_type = "dea_results",
      schema_version = "1.0.0",
      source_software = as.character(source_software),
      analysis_configuration = config,
      layer_names = layer_names,
      feature_keys = .dea_anndata_hierarchy_keys(se),
      sample_key = config$sample_name,
      varm_columns = varm$columns,
      varm_annotations = varm$annotations,
      contrasts = metadata$contrasts,
      formula = metadata$formula,
      default_model = metadata$default_model,
      provenance = provenance,
      bfabric_urls = metadata$bfabric_urls,
      analysis_configuration_raw = metadata$analysis_configuration_raw
    )
  ))
}

.validate_dea_result_anndata <- function(
  adata,
  expected_obs_names = NULL,
  expected_var_names = NULL
) {
  validate_prolfquapp_anndata(adata)
  metadata <- adata$uns[["prolfquapp"]]
  if (!identical(metadata$artifact_type, "dea_results")) {
    stop("AnnData is not a prolfquapp DEA-results artifact.")
  }
  if (!is.null(expected_obs_names)) {
    observed <- rownames(as.data.frame(adata$obs))
    if (!identical(observed, expected_obs_names)) {
      stop("AnnData obs names changed during conversion.")
    }
  }
  if (!is.null(expected_var_names)) {
    observed <- rownames(as.data.frame(adata$var))
    if (!identical(observed, expected_var_names)) {
      stop("AnnData var names changed during conversion.")
    }
  }
  invisible(TRUE)
}

summarized_experiment_to_anndata <- function(se) {
  if (!inherits(se, "SummarizedExperiment")) {
    stop("Expected a SummarizedExperiment object.")
  }
  missing_assays <- setdiff(
    .dea_anndata_required_assays,
    SummarizedExperiment::assayNames(se)
  )
  if (length(missing_assays) > 0L) {
    stop(
      "SummarizedExperiment is missing required assay(s): ",
      paste(missing_assays, collapse = ", ")
    )
  }

  obs_names <- .summarized_experiment_axis_names(se, "obs")
  var_names <- .summarized_experiment_axis_names(se, "var")
  row_data <- SummarizedExperiment::rowData(se)
  contrast_keys <- grep("^constrast_", names(row_data), value = TRUE)
  if (length(contrast_keys) == 0L) {
    stop("SummarizedExperiment has no differential-result rowData entries.")
  }

  obs <- as.data.frame(SummarizedExperiment::colData(se))
  rownames(obs) <- obs_names
  var <- .dea_anndata_feature_table(row_data, contrast_keys, var_names)

  layers <- list(
    raw = t(as.matrix(SummarizedExperiment::assay(se, "rawData"))),
    transformed = t(
      as.matrix(SummarizedExperiment::assay(se, "transformedData"))
    )
  )
  if ("nr_children" %in% SummarizedExperiment::assayNames(se)) {
    layers$nr_children <- t(
      as.matrix(SummarizedExperiment::assay(se, "nr_children"))
    )
  }
  varm <- .dea_anndata_varm(se, var, var_names)

  adata <- anndataR::AnnData(
    X = layers$transformed,
    obs = obs,
    var = var,
    layers = layers,
    varm = varm$values,
    uns = .dea_anndata_uns(se, names(layers), varm)
  )
  .validate_dea_result_anndata(adata, obs_names, var_names)
  adata
}

write_summarized_experiment_h5ad <- function(se, path) {
  destination_dir <- dirname(path)
  if (!dir.exists(destination_dir)) {
    stop("AnnData output directory does not exist: ", destination_dir)
  }

  obs_names <- .summarized_experiment_axis_names(se, "obs")
  var_names <- .summarized_experiment_axis_names(se, "var")
  adata <- summarized_experiment_to_anndata(se)
  temporary <- tempfile(
    pattern = ".AnnData-",
    tmpdir = destination_dir,
    fileext = ".h5ad"
  )
  on.exit(unlink(temporary), add = TRUE)

  invisible(rhdf5::H5get_libversion())
  adata$write_h5ad(temporary, compression = "gzip", mode = "w")
  restored <- anndataR::read_h5ad(temporary)
  .validate_dea_result_anndata(restored, obs_names, var_names)
  if (!file.rename(temporary, path)) {
    stop("Could not move validated AnnData file to: ", path)
  }
  normalizePath(path, mustWork = TRUE)
}
