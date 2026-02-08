#' capture output of function to send it to log
#' @param expr expression to capture output from
#' @export
capture_output <- function(expr) {
  con <- textConnection("output", "w", local = TRUE)
  sink(con)
  on.exit({
    sink()
    close(con)
  })
  eval(expr)
  paste(output, collapse = "\n")
}

#' Synchronize opt and config
#' @param opt list of command-line options
#' @param config ProlfquAppConfig configuration object
#' @export
#' @examples
#' opt <- list()
#' config <- make_DEA_config_R6()
#' xx <- sync_opt_config(opt, config)
#' stopifnot(names(xx$opt) %in% c("software","outdir"))
#' stopifnot(xx$opt$software == "DIANN")
#' stopifnot(xx$opt$outdir == ".")
#' opt$software <- "FP_TMT"
#' opt$outdir <- "testdir"
#' xx <- sync_opt_config(opt, config)
#' stopifnot(xx$config$path == "testdir")
#' stopifnot(xx$config$software == "FP_TMT")
#'
sync_opt_config <- function(opt, config){
  if (!is.null(opt$software)) {
    config$software <-  opt$software
  } else {
    opt$software <- config$software
  }
  if (!is.null(opt$workunit)) {
    logger::log_info("Setting workunit to: " ,  opt$workunit)
    config$project_spec$workunit_Id <- opt$workunit
  }
  if (!is.null(opt$outdir)) {
    config$path <- opt$outdir
  }else {
    if (!is.null(config$path)) {
      opt$outdir <- config$path
    } else {
      opt$outdir <- "."
      config$path <- opt$outdir
    }
  }
  config$set_zipdir_name()
  return(list(opt = opt, config = config))
}


#' Function to normalize paths for both Windows and Linux
#' @export
#' @param paths character vector of file paths
#' @param os operating system type
#' @return normalized path
#' @examples
#' exp_paths <-c("E:\\projects\\p29033\\TKOiWAT\\20240123_015_S629149_iWAT_FL1.d",
#'   "E:\\projects\\p29033\\TKOiWAT\\20240123_016_S629150_iWAT_FL2.d",
#'   "E:\\projects\\p29033\\TKOiWAT\\20240123_017_S629151_iWAT_FL3.d")
#'
#' normalize_path(exp_paths)
#'
normalize_path <- function(paths, os = .Platform$OS.type) {
  # Check the operating system
  if (os == "windows") {
    # On Windows, use the native path
    normalized_paths <- normalizePath(paths, winslash = "\\", mustWork = FALSE)
  } else {
    # On Unix-like systems (Linux, macOS), replace backslashes with forward slashes
    normalized_paths <- gsub("\\\\", "/", paths)
  }
  return(normalized_paths)
}

