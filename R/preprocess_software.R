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
    extra_args = "list(q_value = 0.01, hierarchy_depth = 1)",
    dataset = "prolfquapp::dataset_template_diann"
  ),
  DIANN_PEPTIDE = list(
    get_files = "prolfquapp::get_DIANN_files",
    preprocess = "prolfquapp::preprocess_DIANN",
    extra_args = "list(q_value = 0.01, hierarchy_depth = 2)",
    dataset = "prolfquapp::dataset_template_diann"
  ),
  FP_TMT = list(
    get_files = "prolfquapp::get_FP_PSM_files",
    preprocess = "prolfquapp::preprocess_FP_PSM",
    extra_args = "list(purity_threshold = 0.5, PeptideProphetProb = 0.9, hierarchy_depth = 1)",
    dataset = "prolfquapp::dataset_template_FP_TMT"
  ),
  FP_TMT_PEPTIDE = list(
    get_files = "prolfquapp::get_FP_PSM_files",
    preprocess = "prolfquapp::preprocess_FP_PSM",
    extra_args = "list(purity_threshold = 0.5, PeptideProphetProb = 0.9, hierarchy_depth = 2)",
    dataset = "prolfquapp::dataset_template_FP_TMT"
  ),
  MAXQUANT = list(
    get_files = "prolfquapp::get_MQ_peptide_files",
    preprocess = "prolfquapp::preprocess_MQ_peptide",
    extra_args = "list(hierarchy_depth = 1)",
    dataset = "prolfquapp::dataset_template_MAXQUANT"
  ),
  MAXQUANT_PEPTIDE = list(
    get_files = "prolfquapp::get_MQ_peptide_files",
    preprocess = "prolfquapp::preprocess_MQ_peptide",
    extra_args = "list(hierarchy_depth = 2)",
    dataset = "prolfquapp::dataset_template_MAXQUANT"
  ),
  MSSTATS = list(
    get_files = "prolfquapp::get_MSstats_files",
    preprocess = "prolfquapp::preprocess_MSstats",
    extra_args = "list(hierarchy_depth = 1)",
    dataset = "prolfquapp::dataset_template_MSSTATS"

  ),
  MSSTATS_PEPTIDE = list(
    get_files = "prolfquapp::get_MSstats_files",
    preprocess = "prolfquapp::preprocess_MSstats",
    extra_args = "list(hierarchy_depth = 2)",
    dataset = "prolfquapp::dataset_template_MSSTATS"

  ),
  MSSTATS_FP_DIA = list(
    get_files = "prolfquapp::get_MSstats_files",
    preprocess = "prolfquapp::preprocess_MSstats_FPDIA",
    extra_args = "list(hierarchy_depth = 1)",
    dataset = "prolfquapp::dataset_template_MSSTATS"

  ),
  MSSTATS_FP_DIA_PEPTIDE = list(
    get_files = "prolfquapp::get_MSstats_files",
    preprocess = "prolfquapp::preprocess_MSstats_FPDIA",
    extra_args = "list(hierarchy_depth = 2)",
    dataset = "prolfquapp::dataset_template_MSSTATS"

  ),
  BGS = list(
    get_files = "prolfquapp::get_BGS_files",
    preprocess = "prolfquapp::preprocess_BGS",
    extra_args = "list(hierarchy_depth = 1)",
    dataset = "prolfquapp::dataset_template_BGS"

  ),
  BGS_PEPTIDE = list(
    get_files = "prolfquapp::get_BGS_files",
    preprocess = "prolfquapp::preprocess_BGS",
    extra_args = "list(hierarchy_depth = 2)",
    dataset = "prolfquapp::dataset_template_BGS"

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
#' @param indir input directory with quantification data
#' @param annotation annotation list from read_annotation
#' @param preprocess_functions list or R6 object with get_files, preprocess, extra_args
#' @param pattern_contaminants regex pattern for contaminants
#' @param pattern_decoys regex pattern for decoys
#' @param extreader optional external reader configuration
#' @export
#' @examples
#' # example code
#' annot <- data.frame(
#'   file = c("a1.raw", "a2.raw", "a3.raw", "a4.raw"),
#'   name = c("aa", "ba", "aa", "ba"),
#'   group = c("a", "a", "b", "b")
#' )
#'
#' annot <- read_annotation(annot, QC = TRUE)
#' preprocess_functions <- prolfquapp::prolfqua_preprocess_functions[["DUMMY"]]
#' res <- preprocess_software(".", annot, preprocess_functions)
#'
#' xx <- prolfquapp::ExternalReader$new()
#' xx$extra_args <- "list()"
#' xx$get_files <- "prolfquapp::get_dummy_files"
#' xx$preprocess <- "prolfquapp::preprocess_dummy"
#' res <- preprocess_software(".", annotation = annot, preprocess_functions = xx)
#' xx <- prolfquapp::ExternalReader$new()
#'
preprocess_software <- function(indir,
                                annotation,
                                preprocess_functions,
                                pattern_contaminants = "^zz|^CON|Cont_",
                                pattern_decoys = "^rev_|^REV_",
                                extreader = NULL) {
  to_function <- function(x) {
    list(
      get_files = getFromNamespace(sub(".*::", "", x$get_files), sub("^(.*)::.*", "\\1", x$get_files)),
      preprocess = getFromNamespace(sub(".*::", "", x$preprocess), sub("^(.*)::.*", "\\1", x$preprocess)),
      extra_args = eval(parse(text = x$extra_args))
    )
  }

  preprocess_functions <- to_function(preprocess_functions)

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


#' get functions for creating datasets
#' @param preprocess_functions list with get_files, dataset, extra_args entries
#' @export
#' @examples
#'
#' # Get dataset functions for DIANN
#' preprocess_functions <- prolfquapp::prolfqua_preprocess_functions[["DIANN"]]
#' dataset_funcs <- dataset_get_functions(preprocess_functions)
#'
#' # Use the functions
#' \dontrun{
#' files <- dataset_funcs$files_fn("path/to/data")
#' dataset <- dataset_funcs$dataset_fn(files, "output_file.csv")
#' }
dataset_get_functions <- function(preprocess_functions) {
  to_function <- function(x) {
    list(
      files_fn = getFromNamespace(sub(".*::", "", x$get_files), sub("^(.*)::.*", "\\1", x$get_files)),
      dataset_fn = getFromNamespace(sub(".*::", "", x$dataset), sub("^(.*)::.*", "\\1", x$dataset)),
      extra_args = eval(parse(text = x$extra_args))
    )
  }

  preprocess_functions <- to_function(preprocess_functions)
  return(preprocess_functions)
}
