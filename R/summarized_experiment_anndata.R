# DEA-results AnnData ----
#
# prolfquapp writes its differential-expression results as an AnnData file
# beside the SummarizedExperiment. The conversion is
# `as_AnnData.SummarizedExperiment()`; this file adds the DEA-specific checks
# applied before and after the h5ad round-trip.

.validate_dea_result_anndata <- function(adata, obs_names, var_names) {
  validate_prolfquapp_anndata(adata)
  if (!identical(adata$uns[["prolfquapp"]]$artifact_type, "dea_results")) {
    stop("AnnData is not a prolfquapp DEA-results artifact.")
  }
  if (!identical(adata$obs_names, obs_names) || !identical(adata$var_names, var_names)) {
    stop("AnnData obs or var names changed during conversion.")
  }
  invisible(TRUE)
}

#' Convert a DEA SummarizedExperiment to AnnData, with `transformedData` as X
#'
#' @param se a \code{SummarizedExperiment} produced by
#'   \code{DEAReportGenerator$make_SummarizedExperiment()}
#' @return an AnnData object
#' @keywords internal
#' @noRd
summarized_experiment_to_anndata <- function(se) {
  missing_assays <- setdiff(c("rawData", "transformedData"), SummarizedExperiment::assayNames(se))
  if (length(missing_assays) > 0L) {
    stop("SummarizedExperiment is missing required assay(s): ", paste(missing_assays, collapse = ", "))
  }
  adata <- anndataR::as_AnnData(se, assay_name = "transformedData")
  .validate_dea_result_anndata(adata, colnames(se), rownames(se))
  adata
}

#' Write a DEA SummarizedExperiment as an h5ad file
#'
#' Converts the result of \code{DEAReportGenerator$make_SummarizedExperiment()}
#' to AnnData and writes it atomically, validating the file it wrote.
#'
#' @param se a \code{SummarizedExperiment} produced by
#'   \code{DEAReportGenerator$make_SummarizedExperiment()}
#' @param path destination \code{.h5ad} file; its directory must exist
#' @return the normalized path of the written file
#' @export
write_summarized_experiment_h5ad <- function(se, path) {
  adata <- summarized_experiment_to_anndata(se)
  write_h5ad_atomic(adata, path, validate = function(restored) {
    .validate_dea_result_anndata(restored, adata$obs_names, adata$var_names)
  })
}
