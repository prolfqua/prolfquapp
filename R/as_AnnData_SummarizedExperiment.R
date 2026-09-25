# as_AnnData.SummarizedExperiment ----
#
# A mechanical SummarizedExperiment -> AnnData conversion; anndataR only ships
# methods for SingleCellExperiment and Seurat. Every slot is carried by shape:
#
#   assays become layers, and assay_name selects which one is X
#   colData becomes obs
#   the annotation rowData frame, plus any atomic rowData columns, becomes var
#   every other nested rowData frame becomes a varm data frame
#   metadata becomes uns$prolfquapp, wholesale
#
# HDF5 groups are unordered, so uns$prolfquapp also records layer_names,
# uns_table_columns, varm_key_order and uns_list_order.

.summarized_experiment_axis_names <- function(se, axis) {
  values <- if (identical(axis, "obs")) colnames(se) else rownames(se)
  if (is.null(values) || anyNA(values) || !all(nzchar(values)) || anyDuplicated(values)) {
    stop("SummarizedExperiment ", axis, " names must be present, non-empty and unique.")
  }
  unname(as.character(values))
}

.aligned_row_data_frame <- function(value, var_names, key) {
  value <- as.data.frame(value)
  if (!identical(unname(rownames(value)), var_names)) {
    stop("SummarizedExperiment rowData '", key, "' is not aligned to the feature axis.")
  }
  value
}

#' Coerce a metadata value into something the h5ad writer accepts
#'
#' Flattens factors and drops names AnnData would not preserve. A table is
#' stored column by column: anndataR writes a one-row data frame's columns as
#' HDF5 scalars, which Python `anndata` cannot read. An unnamed list is
#' rejected, since anndataR would write it as an empty group.
#' @keywords internal
#' @noRd
.anndata_uns_value <- function(value, path = "metadata", depth = 0L) {
  if (is.data.frame(value) && depth > 1L) {
    stop("AnnData cannot store the nested table at ", path, "; metadata tables are supported at the top level only.")
  }
  if (is.list(value)) {
    if (!is.data.frame(value) && length(value) > 0L && is.null(names(value))) {
      stop("AnnData cannot store the unnamed list at ", path, "; use a vector instead, or name its elements.")
    }
    paths <- paste0(path, "$", names(value), recycle0 = TRUE)
    return(Map(.anndata_uns_value, value, paths, MoreArgs = list(depth = depth + 1L)))
  }
  unname(if (is.factor(value)) as.character(value) else value)
}

#' Record the names of a metadata list and of its nested lists, in order
#' @keywords internal
#' @noRd
.anndata_list_order <- function(value) {
  if (!is.list(value) || is.data.frame(value) || is.null(names(value))) {
    return(NULL)
  }
  children <- Filter(Negate(is.null), lapply(value, .anndata_list_order))
  c(list(names = names(value)), if (length(children) > 0L) list(children = children))
}

#' Convert a SummarizedExperiment to AnnData
#'
#' Registers prolfquapp's `SummarizedExperiment` method on
#' \code{anndataR::as_AnnData}. Assays become layers, `colData` becomes `obs`,
#' the `annotation` rowData frame and flat `rowData` columns become `var`,
#' other nested `rowData` data frames become `varm` data frames, and `metadata`
#' is carried into `uns$prolfquapp` in full.
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
    layers_mapping,
    obs_mapping,
    var_mapping,
    obsm_mapping,
    varm_mapping,
    obsp_mapping,
    varp_mapping,
    uns_mapping
  )
  if (!all(vapply(mappings, isTRUE, logical(1)))) {
    stop("as_AnnData.SummarizedExperiment supports only TRUE for the *_mapping arguments.")
  }
  if (!identical(match.arg(output_class), "InMemory")) {
    stop("as_AnnData.SummarizedExperiment builds in-memory AnnData only; write it out with write_h5ad().")
  }
  assay_names <- SummarizedExperiment::assayNames(x)
  assay_name <- assay_name %||% x_mapping %||% utils::tail(assay_names, 1L)
  if (length(assay_name) != 1L || !assay_name %in% assay_names) {
    stop("SummarizedExperiment has no assay named '", assay_name, "'.")
  }
  obs_names <- .summarized_experiment_axis_names(x, "obs")
  var_names <- .summarized_experiment_axis_names(x, "var")

  layers <- lapply(stats::setNames(nm = assay_names), function(name) {
    unname(t(as.matrix(SummarizedExperiment::assay(x, name))))
  })
  obs <- as.data.frame(SummarizedExperiment::colData(x))
  rownames(obs) <- obs_names

  row_data <- SummarizedExperiment::rowData(x)
  nested <- vapply(as.list(row_data), function(v) is.data.frame(v) || methods::is(v, "DataFrame"), logical(1))
  var <- as.data.frame(row_data[, !nested, drop = FALSE])
  if ("annotation" %in% names(row_data)) {
    var <- cbind(.aligned_row_data_frame(row_data$annotation, var_names, "annotation"), var)
  }
  rownames(var) <- var_names
  varm_names <- setdiff(names(row_data)[nested], "annotation")
  varm <- lapply(stats::setNames(nm = varm_names), function(name) {
    .aligned_row_data_frame(row_data[[name]], var_names, name)
  })
  # HDF5 treats "/" as a group separator; percent-encoding keeps keys writable.
  names(varm) <- utils::URLencode(varm_names, reserved = TRUE, repeated = TRUE)

  metadata_values <- as.list(S4Vectors::metadata(x))
  metadata <- .anndata_uns_value(metadata_values)
  metadata$uns_table_columns <- lapply(Filter(is.data.frame, metadata_values), names)
  metadata$layer_names <- assay_names
  metadata$varm_key_order <- names(varm)
  # List order carries meaning: the first analysis_configuration$factors entry
  # is the modelled factor.
  metadata$uns_list_order <- .anndata_list_order(metadata_values)

  anndataR::AnnData(
    X = layers[[assay_name]],
    obs = obs,
    var = var,
    layers = layers,
    varm = varm,
    uns = list(X_layer_name = assay_name, prolfquapp = metadata)
  )
}
