# DEA-results AnnData ----
#
# prolfquapp writes its differential-expression results as an AnnData file
# beside the SummarizedExperiment. The conversion itself lives in
# `as_AnnData.SummarizedExperiment()` and is purely mechanical; this file adds
# the DEA-specific expectations: the assays a DEA result must have, and the
# validation applied before and after the h5ad round-trip.

.dea_anndata_required_assays <- c("rawData", "transformedData")

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

#' Convert a DEA SummarizedExperiment to AnnData
#'
#' Thin wrapper over \code{\link{as_AnnData.SummarizedExperiment}} that checks
#' the DEA-specific requirements and validates the result.
#'
#' @param se a \code{SummarizedExperiment} produced by
#'   \code{DEAReportGenerator$make_SummarizedExperiment()}
#' @return an AnnData object
#' @keywords internal
#' @noRd
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
  adata <- anndataR::as_AnnData(se, assay_name = "transformedData")
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
