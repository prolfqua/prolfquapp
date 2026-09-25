# AnnData -> SummarizedExperiment ----
#
# Inverse of `as_AnnData.SummarizedExperiment()`. anndataR's
# `as_SingleCellExperiment()` has no destination for `varm`, so it would drop
# every contrast frame. Layers become assays, obs becomes colData, var becomes
# the rowData annotation frame, each varm data frame becomes a nested rowData
# frame, and uns$prolfquapp becomes metadata. `X` duplicates the layer named by
# `uns$X_layer_name` and is ignored.

#' Convert an AnnData written by prolfquapp back to a SummarizedExperiment
#'
#' @param adata an \code{anndataR::AnnData} object
#' @return a \code{SummarizedExperiment}
#' @keywords internal
#' @noRd
anndata_to_summarized_experiment <- function(adata) {
  validate_prolfquapp_anndata(adata)
  metadata <- adata$uns[["prolfquapp"]]
  obs <- as.data.frame(adata$obs)
  var <- as.data.frame(adata$var)

  layer_keys <- .anndata_ordered_keys(metadata$layer_names, adata$layers_keys(), "layer_names")
  assays <- lapply(stats::setNames(nm = layer_keys), function(key) {
    values <- t(as.matrix(adata$layers[[key]]))
    dimnames(values) <- list(rownames(var), rownames(obs))
    values
  })
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = assays,
    colData = S4Vectors::DataFrame(obs, check.names = FALSE),
    metadata = .anndata_uns_metadata(metadata)
  )
  SummarizedExperiment::rowData(se)$annotation <- var
  for (key in .anndata_ordered_keys(metadata$varm_key_order, adata$varm_keys(), "varm_key_order")) {
    SummarizedExperiment::rowData(se)[[utils::URLdecode(key)]] <- as.data.frame(adata$varm[[key]])
  }
  se
}

# `present`, in the order `recorded` lists them; HDF5 groups are unordered.
.anndata_ordered_keys <- function(recorded, present, label) {
  recorded <- as.vector(unlist(recorded, use.names = FALSE))
  if (!setequal(recorded, present)) {
    mismatch <- c(setdiff(recorded, present), setdiff(present, recorded))
    stop("AnnData ", label, " does not list the keys the AnnData carries: ", paste(mismatch, collapse = ", "))
  }
  recorded
}

#' Rebuild SummarizedExperiment metadata from the prolfquapp uns namespace
#'
#' Drops the reshape bookkeeping, turns the one-dimensional arrays the h5ad
#' reader returns back into vectors, rebuilds the tables listed in
#' `uns_table_columns`, and restores the recorded list order.
#' @keywords internal
#' @noRd
.anndata_uns_metadata <- function(metadata) {
  restore <- function(value) {
    if (is.list(value)) {
      lapply(value, restore)
    } else if (length(dim(value)) == 1L) {
      as.vector(value)
    } else {
      value
    }
  }
  bookkeeping <- c("layer_names", "uns_table_columns", "varm_key_order", "uns_list_order")
  values <- restore(metadata[setdiff(names(metadata), bookkeeping)])
  tables <- metadata$uns_table_columns
  for (name in intersect(names(tables), names(values))) {
    columns <- .anndata_ordered_keys(tables[[name]], names(values[[name]]), paste0("uns_table_columns$", name))
    values[[name]] <- as.data.frame(values[[name]][columns], stringsAsFactors = FALSE)
  }
  .anndata_restore_list_order(values, metadata$uns_list_order)
}

# Reorders a metadata list, and its nested lists, by the recorded
# `uns_list_order`; AnnData written without it is left unchanged.
.anndata_restore_list_order <- function(value, order) {
  if (is.null(order) || !is.list(value) || is.data.frame(value)) {
    return(value)
  }
  recorded <- as.vector(unlist(order$names, use.names = FALSE))
  value <- value[union(intersect(recorded, names(value)), names(value))]
  for (name in intersect(names(order$children), names(value))) {
    value[[name]] <- .anndata_restore_list_order(value[[name]], order$children[[name]])
  }
  value
}
