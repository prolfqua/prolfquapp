#' Write an AnnData object atomically, validating the round trip
#'
#' Writes `adata` to a temporary `.h5ad` in the destination directory, reads
#' it back with [anndataR::read_h5ad()], calls `validate` on the restored
#' object, and renames the temporary file to `path` only if validation
#' returns. The temporary file is removed on any failure.
#'
#' @param adata An `anndataR` AnnData object.
#' @param path Destination `.h5ad` path. Its directory must exist.
#' @param validate Function of one argument, the restored AnnData. Should
#'   `stop()` on a failed round trip; its return value is ignored.
#' @param compression Passed to `write_h5ad()`. Default `"gzip"`.
#' @return The normalized destination path, invisibly.
#' @export
write_h5ad_atomic <- function(
  adata,
  path,
  validate = function(restored) NULL,
  compression = "gzip"
) {
  if (!dir.exists(dirname(path))) {
    stop("AnnData output directory does not exist: ", dirname(path))
  }
  temporary <- tempfile(".h5ad-write-", tmpdir = dirname(path), fileext = ".h5ad")
  on.exit(unlink(temporary), add = TRUE)
  invisible(rhdf5::H5get_libversion())
  adata$write_h5ad(temporary, compression = compression, mode = "w")
  restored <- anndataR::read_h5ad(temporary)
  validate(restored)
  if (!file.rename(temporary, path)) {
    stop("Could not move validated AnnData file to: ", path)
  }
  invisible(normalizePath(path, mustWork = TRUE))
}
