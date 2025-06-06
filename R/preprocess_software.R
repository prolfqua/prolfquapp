#' Mapping of software to string representations of functions and arguments
#' Preprocess Functions Mapping
#'
#' A list containing the mappings of software to string representations of functions and arguments.
#' This is used to dynamically process different proteomics software outputs.
#'
#' @export
#'
prolfqua_preprocess_functions <- list(
  DIANN = list(
    get_files = "prolfquapp::get_DIANN_files",
    preprocess = "prolfquapp::preprocess_DIANN",
    extra_args = "list(q_value = 0.01, hierarchy_depth = 1)"
  ),
  DIANN_PEPTIDE = list(
    get_files = "prolfquapp::get_DIANN_files",
    preprocess = "prolfquapp::preprocess_DIANN",
    extra_args = "list(q_value = 0.01, hierarchy_depth = 2)"
  ),
  FP_TMT = list(
    get_files = "prolfquapp::get_FP_PSM_files",
    preprocess = "prolfquapp::preprocess_FP_PSM",
    extra_args = "list(purity_threshold = 0.5, PeptideProphetProb = 0.9, hierarchy_depth = 1)"
  ),
  FP_TMT_PEPTIDE = list(
    get_files = "prolfquapp::get_FP_PSM_files",
    preprocess = "prolfquapp::preprocess_FP_PSM",
    extra_args = "list(purity_threshold = 0.5, PeptideProphetProb = 0.9, hierarchy_depth = 2)"
  ),
  MAXQUANT = list(
    get_files = "prolfquapp::get_MQ_peptide_files",
    preprocess = "prolfquapp::preprocess_MQ_peptide",
    extra_args = "list(hierarchy_depth = 1)"
  ),
  MAXQUANT_PEPTIDE = list(
    get_files = "prolfquapp::get_MQ_peptide_files",
    preprocess = "prolfquapp::preprocess_MQ_peptide",
    extra_args = "list(hierarchy_depth = 2)"
  ),
  MSSTATS = list(
    get_files = "prolfquapp::get_MSstats_files",
    preprocess = "prolfquapp::preprocess_MSstats",
    extra_args = "list(hierarchy_depth = 1)"
  ),
  MSSTATS_PEPTIDE = list(
    get_files = "prolfquapp::get_MSstats_files",
    preprocess = "prolfquapp::preprocess_MSstats",
    extra_args = "list(hierarchy_depth = 2)"
  ),
  MSSTATS_FP_DIA = list(
    get_files = "prolfquapp::get_MSstats_files",
    preprocess = "prolfquapp::preprocess_MSstats_FPDIA",
    extra_args = "list(hierarchy_depth = 1)"
  ),
  MSSTATS_FP_DIA_PEPTIDE = list(
    get_files = "prolfquapp::get_MSstats_files",
    preprocess = "prolfquapp::preprocess_MSstats_FPDIA",
    extra_args = "list(hierarchy_depth = 2)"
  ),
  BGS = list(
    get_files = "prolfquapp::get_BGS_files",
    preprocess = "prolfquapp::preprocess_BGS",
    extra_args = "list(hierarchy_depth = 1)"
  ),
  BGS_PEPTIDE = list(
    get_files = "prolfquapp::get_BGS_files",
    preprocess = "prolfquapp::preprocess_BGS",
    extra_args = "list(hierarchy_depth = 2)"
  ),
  DUMMY = list(
    get_files = "prolfquapp::get_dummy_files",
    preprocess = "prolfquapp::preprocess_dummy",
    extra_args = "list()"
  ),
  MZMINE = list(
    get_files = "prolfquapp::get_mzMine_files",
    preprocess = "prolfquapp::preprocess_mzMine",
    extra_args = "list(annotated = FALSE)"
  ),
  MZMINEannot = list(
    get_files = "prolfquapp::get_mzMine_files",
    preprocess = "prolfquapp::preprocess_mzMine",
    extra_args = "list(annotated = TRUE)"
  )

)

#' collects preprocess methods for various software
#' @export
#' @examples
#' # example code
#' annot <- data.frame(
#' file = c("a1.raw","a2.raw","a3.raw","a4.raw"),
#' name = c("aa","ba","aa","ba"),
#' group = c("a","a","b","b"))
#' annot <- read_annotation(annot, QC = TRUE)
#' res <- preprocess_software(".",annot, software = "DUMMY" )
#'
#' xx <- prolfquapp::ExternalReader$new()
#' xx$extra_args = "list()"
#' xx$get_files = "prolfquapp::get_dummy_files"
#' xx$preprocess = "prolfquapp::preprocess_dummy"
#' res <- preprocess_software(".",annot, xx, software = "FUNNY" )
#' xx <- prolfquapp::ExternalReader$new()
#' res <- preprocess_software(".",annot, xx, software = "DUMMY" )
preprocess_software <- function(indir,
                                annotation,
                                software,
                                preprocess_functions_str = NULL,
                                pattern_contaminants = "^zz|^CON|Cont_",
                                pattern_decoys = "^rev_|^REV_",
                                extreader = NULL
                                ) {
  to_function <- function(x) {
    list(
      get_files = getFromNamespace(sub(".*::", "", x$get_files), sub("^(.*)::.*", "\\1", x$get_files)),
      preprocess = getFromNamespace(sub(".*::", "", x$preprocess), sub("^(.*)::.*", "\\1", x$preprocess)),
      extra_args = eval(parse(text = x$extra_args))
    )
  }

  if (!is.null(preprocess_functions_str) && length(preprocess_functions_str$preprocess) == 1) {
    preprocess_functions <- to_function(preprocess_functions_str)
  }else{
    preprocess_functions_str <- prolfquapp::prolfqua_preprocess_functions
      # Check if software has a corresponding preprocess function
      if (!software %in% names(preprocess_functions_str)) {
        logger::log_error("No such software: ", software)
        logger::log_error("Available readers are: ", names(preprocess_functions_str))
        stop("No such software.")
        return(NULL)
      }

    preprocess_functions <-
      to_function(preprocess_functions_str[[software]])
  }
  # Fetch files
  files_fn <- preprocess_functions$get_files

  # Preprocess the data
  preprocess_fn <- preprocess_functions$preprocess
  extra_args <- preprocess_functions$extra_args


  files <- files_fn(indir)

  # Log files information
  logger::log_info("Files data: ", paste(files$data, collapse = "; "))
  logger::log_info("Files fasta: ", paste0(files$fasta, collapse = "; "))

  xd <- do.call(preprocess_fn, c(list(
    quant_data = files$data,
    fasta_file = files$fasta,
    annotation = annotation,
    pattern_contaminants = pattern_contaminants,
    pattern_decoys = pattern_decoys
  ), extra_args))

  return(list(xd = xd, files = files))
}
