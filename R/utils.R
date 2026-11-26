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
sync_opt_config <- function(opt, config) {
  if (!is.null(opt$software)) {
    config$software <- opt$software
  } else {
    opt$software <- config$software
  }
  if (!is.null(opt$workunit)) {
    logger::log_info("Setting workunit to: ", opt$workunit)
    config$project_spec$workunit_Id <- opt$workunit
  }
  if (!is.null(opt$outdir)) {
    config$path <- opt$outdir
  } else {
    opt$outdir <- config$path %||% "."
    config$path <- opt$outdir
  }
  if (!is.null(opt$model)) {
    config$processing_options$model <- opt$model
  }
  if (!is.null(opt$nr_peptides)) {
    config$processing_options$nr_peptides <- .validate_nr_peptides(opt$nr_peptides)
  }
  if (!is.null(opt$flat_outdir)) {
    config$flat_outdir <- isTRUE(opt$flat_outdir)
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
  if (os == "windows") {
    normalizePath(paths, winslash = "\\", mustWork = FALSE)
  } else {
    gsub("\\\\", "/", paths)
  }
}


#' Split united hierarchy keys into their component columns
#'
#' prolfqua unites hierarchy columns that are configured as several columns (e.g. `proteinID` and `idtype2` for
#' `protein_Id`) with `config$sep`. This adds the components back as columns, after each key column, for the
#' hierarchy levels up to `config$hierarchy_depth`.
#' @param data data.frame with the united hierarchy key columns
#' @param config [prolfqua::AnalysisConfiguration]
#' @return `data` with each hierarchy key at the current depth split at `config$sep` into its component columns,
#'   which are added after the key column; keys with too few parts fill from the start
#' @export
#' @examples
#' bb <- prolfqua::sim_lfq_data_protein_config(n_proteins = 3)
#' lfq <- prolfqua::LFQData$new(bb$data, bb$config)
#' res <- separate_hierarchy(lfq$data_long(), lfq$get_config())
#' setdiff(colnames(res), colnames(lfq$data_long()))
separate_hierarchy <- function(data, config) {
  for (hkey in config$hierarchy_keys_depth()) {
    parts <- config$hierarchy[[hkey]]
    if (!(length(parts) == 1 && hkey == parts)) {
      data <- data |>
        tidyr::separate_wider_delim(
          dplyr::all_of(hkey),
          delim = config$sep,
          names = parts,
          too_few = "align_start",
          too_many = "drop",
          cols_remove = FALSE
        ) |>
        dplyr::relocate(dplyr::all_of(hkey), .before = dplyr::all_of(parts[1]))
    }
  }
  data
}

#' Floor the identification confidence of accepted observations at 0.5
#'
#' The readers turn the software's identification scores (1 - PEP, 1 - q-value, probabilities) into the
#' identification confidence that prolfqua uses as model weight. An observation the software accepted keeps a
#' confidence of at least 0.5: in prolfquabenchmark, weights near 0 for accepted match-between-runs peptides lowered
#' pAUC_10 by up to 12 points, while a floor of 0.5 removed that loss and kept or improved the peptide-level models.
#' @param x numeric vector of identification confidences in \[0, 1\]
#' @return `x` with values below 0.5 raised to 0.5; NA stays NA
#' @export
#' @examples
#' floor_identification_confidence(c(0.01, 0.4, 0.9, NA))
floor_identification_confidence <- function(x) {
  pmax(x, 0.5)
}
