.onLoad <- function(libname, pkgname) {
  # The optional "saint" modelling backend lives in prolfquasaint, which
  # registers its facade in prolfqua's facade registry from its own .onLoad.
  # prolfquasaint is a Suggests dependency, so load it on demand here (when
  # installed) to keep the "saint" key reachable via prolfqua::lookup_facade()
  # exactly like every other backend -- without making it a hard Imports
  # dependency or adding saint-specific branches elsewhere.
  if (requireNamespace("prolfquasaint", quietly = TRUE)) {
    invisible(NULL)
  }
  invisible()
}
