# AnnData -> SummarizedExperiment ----
#
# Inverse of `as_AnnData.SummarizedExperiment()`. anndataR cannot do this:
# its `as_SingleCellExperiment()` maps var -> rowData, uns -> metadata,
# obsm -> reducedDims and varp -> rowPairs, but has no destination for `varm`,
# so every contrast frame would be dropped silently.
#
#   layers become assays, transposed back
#   obs becomes colData
#   var becomes the rowData annotation frame
#   varm becomes one nested rowData frame per key, its column names and
#     non-numeric columns taken from varm_columns and varm_annotations
#   uns$prolfquapp becomes metadata
#
# `X` is ignored: it duplicates one of the layers, named by `uns$X_layer_name`.

# uns keys written by the forward conversion to describe the reshape itself.
# They are not part of the SummarizedExperiment metadata.
.anndata_reshape_keys <- c(
  "layer_names",
  "uns_table_columns",
  "varm_columns",
  "varm_annotations",
  "varm_column_order",
  "varm_key_order"
)

#' Rebuild one nested rowData frame from a varm matrix
#'
#' @param values the varm matrix
#' @param columns column names of `values`
#' @param annotations list of the frame's non-numeric columns
#' @param column_order the order the frame's columns had before the split
#' @param var_names the feature axis
#' @return data.frame indexed by `var_names`
#' @keywords internal
#' @noRd
.anndata_varm_frame <- function(
  values,
  columns,
  annotations,
  column_order,
  var_names
) {
  values <- as.matrix(values)
  if (ncol(values) != length(columns)) {
    stop(
      "AnnData varm matrix has ",
      ncol(values),
      " columns but varm_columns names ",
      length(columns),
      "."
    )
  }
  frame <- as.data.frame(values)
  names(frame) <- unlist(columns, use.names = FALSE)
  for (name in names(annotations)) {
    # HDF5 hands back a one-dimensional array; the frame wants a plain vector.
    frame[[name]] <- as.vector(unlist(annotations[[name]], use.names = FALSE))
  }
  column_order <- unlist(column_order, use.names = FALSE)
  if (!setequal(column_order, names(frame))) {
    stop(
      "AnnData varm_column_order does not match the columns of the frame: ",
      paste(setdiff(column_order, names(frame)), collapse = ", ")
    )
  }
  frame <- frame[, column_order, drop = FALSE]
  rownames(frame) <- var_names
  frame
}

#' Convert an AnnData written by prolfquapp back to a SummarizedExperiment
#'
#' @param adata an \code{anndataR::AnnData} object
#' @return a \code{SummarizedExperiment}
#' @keywords internal
#' @noRd
anndata_to_summarized_experiment <- function(adata) {
  if (!inherits(adata, "AbstractAnnData")) {
    stop("Expected an anndataR AnnData object.")
  }
  metadata <- adata$uns[["prolfquapp"]]
  if (is.null(metadata)) {
    stop(
      "AnnData uns is missing the 'prolfquapp' namespace. ",
      "This AnnData was not created by prolfquapp."
    )
  }

  obs <- as.data.frame(adata$obs)
  var <- as.data.frame(adata$var)
  obs_names <- rownames(obs)
  var_names <- rownames(var)

  # HDF5 groups are unordered, so slot order comes from the recorded names.
  layer_keys <- .anndata_ordered_keys(
    metadata$layer_names,
    adata$layers_keys(),
    "layer_names"
  )
  assays <- lapply(layer_keys, function(key) {
    values <- t(as.matrix(adata$layers[[key]]))
    dimnames(values) <- list(var_names, obs_names)
    values
  })
  names(assays) <- layer_keys

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = assays,
    colData = S4Vectors::DataFrame(obs, check.names = FALSE),
    metadata = .anndata_uns_metadata(metadata)
  )

  rownames(var) <- var_names
  SummarizedExperiment::rowData(se)[[.anndata_var_frame_name]] <- var
  varm_keys <- .anndata_ordered_keys(
    metadata$varm_key_order,
    adata$varm_keys(),
    "varm_key_order"
  )
  for (key in varm_keys) {
    SummarizedExperiment::rowData(se)[[utils::URLdecode(key)]] <-
      .anndata_varm_frame(
        adata$varm[[key]],
        metadata$varm_columns[[key]],
        metadata$varm_annotations[[key]],
        metadata$varm_column_order[[key]],
        var_names
      )
  }
  se
}

#' Slot keys in the order they were written
#'
#' @param recorded the recorded order, from `uns`
#' @param present the keys the AnnData actually carries
#' @param label which recorded key is being used, for error messages
#' @return `present`, ordered by `recorded`
#' @keywords internal
#' @noRd
.anndata_ordered_keys <- function(recorded, present, label) {
  recorded <- as.vector(unlist(recorded, use.names = FALSE))
  if (!setequal(recorded, present)) {
    stop(
      "AnnData ",
      label,
      " does not list the same keys the AnnData carries: ",
      paste(union(setdiff(recorded, present), setdiff(present, recorded)), collapse = ", ")
    )
  }
  recorded
}

#' Undo what the h5ad reader does to plain values
#'
#' HDF5 returns a vector as a one-dimensional array, which makes metadata read
#' from an `.h5ad` compare unequal to the same metadata read from an `.rds`.
#' @param value a value from `uns`
#' @return the value with one-dimensional array shapes dropped
#' @keywords internal
#' @noRd
.anndata_uns_restore <- function(value) {
  if (is.data.frame(value)) {
    value[] <- lapply(value, .anndata_uns_restore)
    return(value)
  }
  if (is.list(value)) {
    return(lapply(value, .anndata_uns_restore))
  }
  if (length(dim(value)) == 1L) {
    return(as.vector(value))
  }
  value
}

#' Rebuild SummarizedExperiment metadata from the prolfquapp uns namespace
#'
#' Drops the keys that only describe the reshape and turns the entries listed
#' in `uns_table_columns` back into data frames -- they are stored column by
#' column because anndataR cannot write a one-row dataframe group that
#' `anndata` in Python can read.
#' @param metadata the `uns$prolfquapp` list
#' @return list suitable for `S4Vectors::metadata()`
#' @keywords internal
#' @noRd
.anndata_uns_metadata <- function(metadata) {
  values <- .anndata_uns_restore(
    metadata[setdiff(names(metadata), .anndata_reshape_keys)]
  )
  tables <- metadata$uns_table_columns
  for (name in intersect(names(tables), names(values))) {
    columns <- .anndata_ordered_keys(
      tables[[name]],
      names(values[[name]]),
      paste0("uns_table_columns$", name)
    )
    values[[name]] <- as.data.frame(
      values[[name]][columns],
      stringsAsFactors = FALSE
    )
  }
  values
}
