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
.anndata_reshape_keys <- c("layer_names", "varm_columns", "varm_annotations")

#' Rebuild one nested rowData frame from a varm matrix
#'
#' @param values the varm matrix
#' @param columns column names of `values`
#' @param annotations list of the frame's non-numeric columns
#' @param var_names the feature axis
#' @return data.frame indexed by `var_names`
#' @keywords internal
#' @noRd
.anndata_varm_frame <- function(values, columns, annotations, var_names) {
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
    frame[[name]] <- unlist(annotations[[name]], use.names = FALSE)
  }
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

  layer_keys <- adata$layers_keys()
  assays <- lapply(layer_keys, function(key) {
    values <- t(as.matrix(adata$layers[[key]]))
    dimnames(values) <- list(var_names, obs_names)
    values
  })
  names(assays) <- layer_keys

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = assays,
    colData = S4Vectors::DataFrame(obs, check.names = FALSE),
    metadata = metadata[setdiff(names(metadata), .anndata_reshape_keys)]
  )

  rownames(var) <- var_names
  SummarizedExperiment::rowData(se)[[.anndata_var_frame_name]] <- var
  for (key in adata$varm_keys()) {
    SummarizedExperiment::rowData(se)[[utils::URLdecode(key)]] <-
      .anndata_varm_frame(
        adata$varm[[key]],
        metadata$varm_columns[[key]],
        metadata$varm_annotations[[key]],
        var_names
      )
  }
  se
}
