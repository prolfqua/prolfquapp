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
    c("-t", "--pattern_decoys"),
    type = "character",
    default = "^REV_|^rev_",
    help = "decoy pattern in fasta file",
    metavar = "string"
  ),
  make_option(
    c("-w", "--workunit"),
    type = "character",
    default = "",
    help = "workunit identifier",
    metavar = "string"
  ),
  make_option(
    c("-d", "--dataset"),
    type = "character",
    default = "dataset.csv",
    help = "annotation file",
    metavar = "string"
  ),
  make_option(
    c("-o", "--outdir"),
    type = "character",
    default = "qc_dir",
    help = "folder to write the results to.",
    metavar = "string"
  ),
  make_option(
    c("-s", "--software"),
    type = "character",
    default = "DIANN",
    help = paste0(
      "possible options: ",
      paste(
        names(prolfquapp::prolfqua_preprocess_functions)[
          !grepl("PEPTIDE", names(prolfquapp::prolfqua_preprocess_functions))
        ],
        collapse = ", "
      )
    ),
    metavar = "character"
  ),
  make_option(
    c("-p", "--project"),
    type = "character",
    default = "",
    help = "your project identifier",
    metavar = "string"
  ),
  make_option(
    c("-O", "--order"),
    type = "character",
    default = "",
    help = "order ID",
    metavar = "character"
  ),
  make_option(
    c("--libPath"),
    type = "character",
    default = NULL,
    help = " (optional) R library path",
    metavar = "string"
  ),
  optparse::make_option(
    c("-y", "--yaml"),
    type = "character",
    default = "config.yaml",
    help = "yaml configuration file",
    metavar = "character"
  )
)

parser <- OptionParser(usage = "%prog --indir . ", option_list = option_list)

if (length(commandArgs(TRUE)) == 0) {
  optparse::print_help(parser)
}

arguments <- parse_args(parser, positional_arguments = TRUE)

lobstr::tree(arguments)

opt <- arguments$options
ymlfile <- arguments$args

# Interactive debugging: set opt manually, then source from here
if (FALSE) {
  opt$indir <- "DIANN_Result_WU321864/out-2025-02-14/"
  opt$software <- "DIANN"
  opt$dataset <- "DIANN_Result_WU321864/dataset.csv"
  opt$workunit <- "helloW"
}


ymlfile <- if (length(ymlfile) == 0) {
  opt$yaml
} else {
  ymlfile
}


# set library path
if (!is.null(opt$libPath) && dir.exists(opt$libPath)) {
  prolfquapp::set_lib_path(opt$libPath)
}

library(prolfquapp)
library(logger)

logger::log_info("using : ", system.file(package = "prolfqua"))
logger::log_info("using : ", system.file(package = "prolfquapp"))


result <- tryCatch(
  prolfquapp::run_qc_preprocess(
    indir = opt$indir,
    dataset = opt$dataset,
    software = opt$software,
    yaml_file = ymlfile,
    outdir = opt$outdir,
    project = opt$project,
    order = opt$order,
    workunit = opt$workunit
  ),
  error = function(e) {
    stack_trace <- capture.output(traceback())
    logger::log_error(conditionMessage(e), "\n")
    logger::log_error("Stack trace:\n")
    logger::log_error(
      paste(stack_trace, collapse = "\n"), "\n"
    )
    if (interactive()) stop(e) else quit(save = "no", status = 1)
  }
)

xd <- result$xd
GRP2 <- result$config

dir.create(GRP2$path, showWarnings = FALSE)
if (!dir.exists(GRP2$get_zipdir())) {
  dir.create(GRP2$get_zipdir())
}

pap <- QC_generator$new(
  xd$lfqdata, xd$protein_annotation, GRP2
)
pap$copy_dataset(opt$dataset)

pap$write_xlsx()

pap$render_QC_protein_abundances()
pap$render_sample_size_QC()
pap$render_index_html()
pap$render_index_md()
