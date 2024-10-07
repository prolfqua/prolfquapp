#' Mapping of software to string representations of functions and arguments
#' Preprocess Functions Mapping
#'
#' A list containing the mappings of software to string representations of functions and arguments.
#' This is used to dynamically process different proteomics software outputs.
#'
#' @export
#'
prolfq_preprocess_functions <- list(
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
  DUMMY = list(
    get_files = "prolfquapp::get_dummy_files",
    preprocess = "prolfquapp::preprocess_dummy",
    extra_args = "list()"
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
#' res <- preprocess_software(".",annot, prolfq_preprocess_functions, software = "DUMMY" )
preprocess_software <- function(indir,
                                annotation,
                                software,
                                preprocess_functions_str,
                                pattern_contaminants = "^zz|^CON|Cont_",
                                pattern_decoys = "^rev_|^REV_"
                                ) {

  preprocess_functions <-
    lapply(preprocess_functions_str, function(x) {
      list(
        get_files = getFromNamespace(sub(".*::", "", x$get_files), sub("^(.*)::.*", "\\1", x$get_files)),
        preprocess = getFromNamespace(sub(".*::", "", x$preprocess), sub("^(.*)::.*", "\\1", x$preprocess)),
        extra_args = eval(parse(text = x$extra_args))
      )
    })

  # Check if software has a corresponding preprocess function
  if (!software %in% names(preprocess_functions)) {
    logger::log_error("No such software: ", software)
    logger::log_error("Available readers are: ", names(preprocess_functions))
    stop("No such software.")
    return(NULL)
  }

  # Fetch files
  files_fn <- preprocess_functions[[software]]$get_files
  files <- files_fn(indir)
  # Log files information
  logger::log_info("Files data: ", paste(files$data, collapse = "; "))
  logger::log_info("Files fasta: ", paste0(files$fasta, collapse = "; "))

  # Preprocess the data
  preprocess_fn <- preprocess_functions[[software]]$preprocess
  extra_args <- preprocess_functions[[software]]$extra_args

  xd <- do.call(preprocess_fn, c(list(
    quant_data = files$data,
    fasta_file = files$fasta,
    annotation = annotation,
    pattern_contaminants = pattern_contaminants,
    pattern_decoys = pattern_decoys
  ), extra_args))

  return(list(xd = xd, files = files))
}
