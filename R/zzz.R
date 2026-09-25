.onLoad <- function(libname, pkgname) {
  # Load the optional (Suggests) prolfquasaint, whose .onLoad registers the "saint" facade in prolfqua's registry.
  requireNamespace("prolfquasaint", quietly = TRUE)
  invisible()
}
