if (!require("prolfquapp", quietly = TRUE)) {
  remotes::install_github("prolfqua/prolfquapp", dependencies = TRUE)
}
if (!require("prolfqua", quietly = TRUE)) {
  remotes::install_github("prolfqua/prolfqua", dependencies = TRUE)
}
if (!require("optparse", quietly = TRUE)) {
  install.packages("optparse", dependencies = TRUE)
}

logger::log_info(
  "LIBRARY PATHS (.libPaths()):",
  paste(.libPaths(), collapse = "\n")
)

library("optparse")
option_list <- list(
  make_option(
    c("-i", "--indir"),
    type = "character",
    default = ".",
    help = "folder containing fasta file and output of the quantification software.",
    metavar = "string"
  ),
  make_option(
    c("-d", "--dataset"),
    type = "character",
    default = "dataset.csv",
    help = "name of annotation file",
    metavar = "string"
  ),
  make_option(
    c("-s", "--software"),
    type = "character",
    default = "DIANN",
    help = paste0(
      "possible options: ",
      paste(
        names(
          prolfquapp::get_procfuncs()
        ),
        collapse = ", "
      )
    ),
    metavar = "character"
  )
)

parser <- optparse::OptionParser(
  usage = "%prog config.yaml --software DIANN --indir .",
  option_list = option_list
)
if (length(commandArgs(TRUE)) == 0) {
  optparse::print_help(parser)
}

arguments <- optparse::parse_args(parser, positional_arguments = FALSE)

opt <- arguments
logger::log_info(prolfquapp::capture_output(quote(lobstr::tree(opt))))

# Interactive debugging: set opt manually, then source from here
if (FALSE) {
  opt$indir <- "DIANN_19_all_18_50_50_MBR_v01/"
  opt$dataset <- "dataset.xlsx"
  opt$software <- "DIANN"
}

prolfqua_preprocess_functions <- prolfquapp::get_procfuncs()

if (opt$software %in% names(prolfqua_preprocess_functions)) {
  stopifnot(dir.exists(opt$indir))
  preprocess_functions <- prolfqua_preprocess_functions[[opt$software]]
  preprocess_functions <- prolfquapp::dataset_get_functions(
    preprocess_functions
  )
  files <- preprocess_functions$files_fn(opt$indir)
  datasettemplate <- preprocess_functions$dataset_fn(files)
  prolfquapp::write_annotation_file(datasettemplate, opt$dataset)
} else {
  logger::log_error("software not supported")
  stop(
    "Software '",
    opt$software,
    "' not supported. Available software: ",
    paste(names(prolfqua_preprocess_functions), collapse = ", ")
  )
}
