#' capture output of function to send it to log
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

