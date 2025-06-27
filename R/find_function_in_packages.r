find_function_packages <- function(function_name = "prolfqua_preprocess_functions", prefix = "prolfqua") {
  # Get all installed packages
  all_packages <- installed.packages()[, "Package"]
  all_packages <- grep(paste0("^", prefix), all_packages, value = TRUE)
  found_packages <- character(0)

  for (pkg in all_packages) {
    tryCatch(
      {
        # Try to get the function from the package
        func <- getFromNamespace(function_name, pkg, envir = NULL)
        if (!is.null(func)) {
          found_packages <- c(found_packages, pkg)
        }
      },
      error = function(e) {
        # Function not found in this package, continue
      }
    )
  }

  return(found_packages)
}

#' Get all processing functions from all packages
#'
#' @param function_name The name of the function to get
#' @param prefix The prefix of the package names
#' @return A list of processing functions
#'
#' @export
#' @examples
#' get_procfuncs()
#' get_procfuncs("prolfqua_preprocess_functions", "prolfqua")
#' get_procfuncs("prolfqua_preprocess_functions", "prolfquapp")
#' get_procfuncs("prolfqua_preprocess_functions", "prolfquapp")
#' get_procfuncs("prolfqua_preprocess_functions", "xdx")
#'
get_procfuncs <- function(function_name = "prolfqua_preprocess_functions", prefix = "prolfqua") {
  package_names <- find_function_packages(function_name, prefix = prefix)
  procfuncs <- list()
  for (pkg in package_names) {
    procfuncs[[pkg]] <- getFromNamespace(function_name, pkg)
  }
  combined_procfuncs <- do.call(c, procfuncs)
  return(combined_procfuncs)
}

