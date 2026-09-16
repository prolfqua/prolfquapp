# as_AnnData.SummarizedExperiment ----
#
# A mechanical SummarizedExperiment -> AnnData conversion. anndataR 1.1.0 only
# ships methods for SingleCellExperiment and Seurat, so prolfquapp registers its
# own method on the generic.
#
# The conversion names no domain field. Every slot is carried across by shape:
#
#   assays become layers, and assay_name selects which one is X
#   colData becomes obs
#   the annotation rowData frame becomes var, with any atomic rowData columns
#   every other nested rowData frame becomes a varm matrix of its numeric
#     columns, its remaining columns kept in uns
#   metadata becomes uns$prolfquapp, wholesale
#
# Only keys that describe the reshape itself (layer_names, uns_table_columns,
# varm_columns, varm_annotations, varm_column_order, varm_key_order) are
# derived here. Everything a consumer needs to interpret
# the result -- artifact type, schema version, column roles, provenance -- is
# read from the SummarizedExperiment, which is the single source of truth.

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

#' Coerce a value into something the h5ad writer accepts
#'
#' Recurses through lists and data frames, flattens factors, and drops names
#' that AnnData would not preserve.
#'
#' A metadata table is stored column by column rather than as an AnnData
#' dataframe group: anndataR writes a one-row data frame's columns as HDF5
#' scalars instead of length-one arrays, which violates the dataframe encoding
#' and makes the file unreadable for `anndata` in Python -- so any
#' single-contrast result would be unreadable. The table names are recorded in
#' `uns_table_names` so a reader can rebuild the frames.
#'
#' An unnamed list has no AnnData representation at all -- anndataR writes one
#' as an empty group, losing its contents without complaining -- so it is
#' rejected here rather than silently dropped.
#' @param value any R value taken from SummarizedExperiment metadata
#' @param path where `value` sits in the metadata, for error messages
#' @param depth nesting depth; tables are only supported at the top level
#' @return the coerced value
#' @keywords internal
#' @noRd
.anndata_uns_value <- function(value, path = "metadata", depth = 0L) {
  if (is.data.frame(value)) {
    if (depth > 1L) {
      stop(
        "AnnData cannot store the nested table at ",
        path,
        "; metadata tables are supported at the top level only."
      )
    }
    return(.anndata_uns_list(as.list(value), path, depth))
  }
  if (is.list(value)) {
    if (length(value) > 0L && is.null(names(value))) {
      stop(
        "AnnData cannot store the unnamed list at ",
        path,
        "; use a vector instead, or name its elements."
      )
    }
    return(.anndata_uns_list(value, path, depth))
  }
  if (is.factor(value)) {
    return(unname(as.character(value)))
  }
  unname(value)
}

.anndata_uns_list <- function(value, path, depth) {
  Map(
    function(element, element_path) {
      .anndata_uns_value(element, element_path, depth + 1L)
    },
    value,
    paste0(path, "$", names(value)),
    USE.NAMES = TRUE
  )
}

#' Split rowData into flat feature columns and nested result frames
#'
#' @param row_data the SummarizedExperiment rowData
#' @return list with `flat` (column names) and `nested` (column names)
#' @keywords internal
#' @noRd
.row_data_partition <- function(row_data) {
  names_all <- names(row_data)
  nested <- names_all[vapply(
    names_all,
    function(name) {
      value <- row_data[[name]]
      is.data.frame(value) || methods::is(value, "DataFrame")
    },
    logical(1)
  )]
  list(flat = setdiff(names_all, nested), nested = nested)
}

#' Find per-feature columns duplicated across every nested frame
#'
#' Encode a rowData frame name for use as an AnnData key
#'
#' HDF5 treats "/" as a group separator, so a name containing one cannot be
#' written as a key. Percent-encoding is applied with `repeated = TRUE` so that
#' an already-encoded name cannot collide with a newly encoded one.
#' @param name rowData frame name
#' @return the encoded key
#' @keywords internal
#' @noRd
.encode_varm_key <- function(name) {
  if (!nzchar(name)) {
    stop("rowData frame names must be non-empty.")
  }
  utils::URLencode(name, reserved = TRUE, repeated = TRUE)
}

#' Name of the rowData frame holding per-feature annotation
#'
#' prolfquapp stores feature annotation once, in a rowData frame under this
#' name. It becomes `var`; every other nested frame becomes a `varm` entry.
#' @keywords internal
#' @noRd
.anndata_var_frame_name <- "annotation"

#' Build var from the annotation frame plus any flat rowData columns
#'
#' @param row_data the SummarizedExperiment rowData
#' @param flat_names flat column names
#' @param var_names the feature axis
#' @return data.frame indexed by `var_names`
#' @keywords internal
#' @noRd
.anndata_var_table <- function(row_data, flat_names, var_names) {
  parts <- list()
  if (.anndata_var_frame_name %in% names(row_data)) {
    annotation <- .aligned_row_data_frame(
      row_data[[.anndata_var_frame_name]],
      var_names,
      .anndata_var_frame_name
    )
    rownames(annotation) <- NULL
    parts[[length(parts) + 1L]] <- annotation
  }
  if (length(flat_names) > 0L) {
    flat <- as.data.frame(row_data[, flat_names, drop = FALSE])
    rownames(flat) <- NULL
    parts[[length(parts) + 1L]] <- flat
  }
  if (length(parts) == 0L) {
    return(data.frame(row.names = var_names))
  }
  result <- do.call(cbind, parts)
  rownames(result) <- var_names
  result
}

#' Reshape nested rowData frames into varm matrices
#'
#' Numeric and logical columns become the varm matrix; remaining columns are
#' returned separately so they can be stored in `uns`.
#' @param row_data the SummarizedExperiment rowData
#' @param nested_names nested column names
#' @param var_names the feature axis
#' @return list with `values`, `columns`, `annotations` and `order`
#' @keywords internal
#' @noRd
.anndata_varm <- function(row_data, nested_names, var_names) {
  values <- list()
  columns <- list()
  annotations <- list()
  order <- list()
  for (name in nested_names) {
    key <- .encode_varm_key(name)
    if (key %in% names(values)) {
      stop("rowData frame names collide after AnnData key encoding: ", name)
    }
    frame <- .aligned_row_data_frame(row_data[[name]], var_names, name)
    # Splitting a frame into a matrix plus annotations loses the order its
    # columns had, so record it for whoever reads the AnnData back.
    order[[key]] <- names(frame)
    numeric_columns <- names(frame)[vapply(
      frame,
      function(column) is.numeric(column) || is.logical(column),
      logical(1)
    )]
    matrix_values <- as.matrix(frame[, numeric_columns, drop = FALSE])
    rownames(matrix_values) <- NULL
    values[[key]] <- matrix_values
    columns[[key]] <- numeric_columns
    annotation_columns <- setdiff(names(frame), numeric_columns)
    if (length(annotation_columns) > 0L) {
      annotations[[key]] <- lapply(
        frame[, annotation_columns, drop = FALSE],
        unname
      )
    }
  }
  list(
    values = values,
    columns = columns,
    annotations = annotations,
    order = order
  )
}

#' Convert a SummarizedExperiment to AnnData
#'
#' Registers prolfquapp's `SummarizedExperiment` method on
#' \code{anndataR::as_AnnData}. anndataR itself only provides methods for
#' `SingleCellExperiment` and `Seurat`.
#'
#' The conversion is mechanical: assays become layers, `colData` becomes `obs`,
#' flat `rowData` columns become `var`, nested `rowData` data frames become
#' `varm` matrices, and `metadata` is carried into `uns$prolfquapp` in full.
#' Because metadata is passed through rather than enumerated, anything recorded
#' on the SummarizedExperiment reaches the AnnData without changing this code.
#'
#' @param x a \code{SummarizedExperiment}
#' @param x_mapping name of the assay to use as `X`; an alias for
#'   \code{assay_name} kept for compatibility with the generic
#' @param layers_mapping,obs_mapping,var_mapping,obsm_mapping,varm_mapping,obsp_mapping,varp_mapping,uns_mapping
#'   accepted for signature compatibility with \code{anndataR::as_AnnData};
#'   this method always maps every slot it can and only \code{TRUE} is
#'   supported
#' @param assay_name assay used as `X`; defaults to the last assay
#' @param output_class must be \code{"InMemory"}; write the result out with
#'   \code{write_h5ad()} instead of asking for another class
#' @param ... ignored
#' @return an AnnData object
#' @exportS3Method anndataR::as_AnnData
as_AnnData.SummarizedExperiment <- function(
  x,
  x_mapping = NULL,
  layers_mapping = TRUE,
  obs_mapping = TRUE,
  var_mapping = TRUE,
  obsm_mapping = TRUE,
  varm_mapping = TRUE,
  obsp_mapping = TRUE,
  varp_mapping = TRUE,
  uns_mapping = TRUE,
  assay_name = NULL,
  output_class = c("InMemory", "HDF5AnnData", "ReticulateAnnData"),
  ...
) {
  mappings <- list(
    layers_mapping = layers_mapping,
    obs_mapping = obs_mapping,
    var_mapping = var_mapping,
    obsm_mapping = obsm_mapping,
    varm_mapping = varm_mapping,
    obsp_mapping = obsp_mapping,
    varp_mapping = varp_mapping,
    uns_mapping = uns_mapping
  )
  unsupported <- names(mappings)[!vapply(mappings, isTRUE, logical(1))]
  if (length(unsupported) > 0L) {
    stop(
      "as_AnnData.SummarizedExperiment supports only TRUE for: ",
      paste(unsupported, collapse = ", ")
    )
  }
  output_class <- match.arg(output_class)
  if (!identical(output_class, "InMemory")) {
    stop(
      "as_AnnData.SummarizedExperiment builds in-memory AnnData only; ",
      "write it out with write_h5ad() instead of requesting '",
      output_class,
      "'."
    )
  }

  assay_names <- SummarizedExperiment::assayNames(x)
  if (length(assay_names) == 0L) {
    stop("SummarizedExperiment has no assays to convert.")
  }
  if (is.null(assay_name)) {
    assay_name <- x_mapping
  }
  if (is.null(assay_name)) {
    assay_name <- assay_names[[length(assay_names)]]
  }
  if (!assay_name %in% assay_names) {
    stop("SummarizedExperiment has no assay named '", assay_name, "'.")
  }

  obs_names <- .summarized_experiment_axis_names(x, "obs")
  var_names <- .summarized_experiment_axis_names(x, "var")

  layers <- lapply(assay_names, function(name) {
    values <- t(as.matrix(SummarizedExperiment::assay(x, name)))
    dimnames(values) <- NULL
    values
  })
  names(layers) <- assay_names

  obs <- as.data.frame(SummarizedExperiment::colData(x))
  rownames(obs) <- obs_names

  row_data <- SummarizedExperiment::rowData(x)
  partition <- .row_data_partition(row_data)
  var <- .anndata_var_table(row_data, partition$flat, var_names)
  varm <- .anndata_varm(
    row_data,
    setdiff(partition$nested, .anndata_var_frame_name),
    var_names
  )

  metadata_values <- as.list(S4Vectors::metadata(x))
  metadata <- .anndata_uns_value(metadata_values)
  # Column order, like every other AnnData group, is not preserved, so each
  # table's columns are recorded in the order the table had them.
  metadata$uns_table_columns <- lapply(
    metadata_values[vapply(metadata_values, is.data.frame, logical(1))],
    names
  )
  metadata$layer_names <- assay_names
  metadata$varm_columns <- varm$columns
  metadata$varm_annotations <- varm$annotations
  metadata$varm_column_order <- varm$order
  # AnnData slots are unordered groups: a reader gets varm keys and layer names
  # back in HDF5's own order, so the order they had here is recorded too.
  metadata$varm_key_order <- names(varm$values)

  anndataR::AnnData(
    X = layers[[assay_name]],
    obs = obs,
    var = var,
    layers = layers,
    varm = varm$values,
    uns = list(
      X_layer_name = assay_name,
      prolfquapp = metadata
    )
  )
}
