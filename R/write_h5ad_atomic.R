# Atomic AnnData write ----
#
# One protocol for every h5ad prolfquapp or a downstream package publishes:
# write to a temporary file beside the destination, read it back, hand the
# restored object to a caller-supplied validator, and only then move it into
# place. A crash or a failed validation leaves no partial file at `path`.
# The container-level writer for multimodal (.h5mu) output builds on the same
# skeleton.

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
  destination_dir <- dirname(path)
  if (!dir.exists(destination_dir)) {
    stop("AnnData output directory does not exist: ", destination_dir)
  }
  temporary <- tempfile(
    pattern = ".h5ad-write-",
    tmpdir = destination_dir,
    fileext = ".h5ad"
  )
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
