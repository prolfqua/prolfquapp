if (!require("prolfqua", quietly = TRUE)) {
  remotes::install_github("prolfqua/prolfquapp", dependencies = TRUE)
}
if (!require("prolfqua", quietly = TRUE)) {
  remotes::install_github("prolfqua/prolfquapp", dependencies = TRUE)
}
if (!require("optparse", quietly = TRUE)) {
  install.packages("optparse", dependencies = TRUE)
}

logger::log_info("LIBRARY PATHS (.libPaths()):", paste(.libPaths(), collapse = "\n"))

library("optparse")
option_list <- list(
  make_option(c("-i", "--indir"),
    type = "character", default = ".",
    help = "folder containing fasta file and output of the quantification software.",
    metavar = "string"
  ),
  make_option(c("-t", "--pattern_decoys"),
    type = "character", default = "^REV_|^rev_",
    help = "decoy pattern in fasta file",
    metavar = "string"
  ),
  make_option(c("-w", "--workunit"),
    type = "character", default = "",
    help = "workunit identifier",
    metavar = "string"
  ),
  make_option(c("-d", "--dataset"),
    type = "character", default = "dataset.csv",
    help = "annotation file",
    metavar = "string"
  ),
  make_option(c("-o", "--outdir"),
    type = "character", default = "qc_dir",
    help = "folder to write the results to.",
    metavar = "string"
  ),
  make_option(c("-s", "--software"),
    type = "character", default = "DIANN",
    help = paste0(
      "possible options: ",
      paste(names(prolfquapp::prolfqua_preprocess_functions)[!grepl("PEPTIDE", names(prolfquapp::prolfqua_preprocess_functions))], collapse = ", ")
    ),
    metavar = "character"
  ),
  make_option(c("-p", "--project"),
    type = "character", default = "",
    help = "your project identifier",
    metavar = "string"
  ),
  make_option(c("-O", "--order"),
    type = "character", default = "",
    help = "order ID",
    metavar = "character"
  ),
  make_option(c("--libPath"),
    type = "character", default = NULL,
    help = " (optional) R library path",
    metavar = "string"
  ),
  optparse::make_option(c("-y", "--yaml"),
    type = "character", default = "config.yaml",
    help = "yaml configuration file",
    metavar = "character"
  )
)

parser <- OptionParser(usage = "%prog --indir . ", option_list = option_list)

if (length(commandArgs(TRUE)) == 0) {
  optparse::print_help(parser)
  # quit(status = 1)
}

arguments <- parse_args(parser, positional_arguments = TRUE)

lobstr::tree(arguments)

opt <- arguments$options
ymlfile <- arguments$args

if (FALSE) {
  opt$indir <- "DIANN_Result_WU321864/out-2025-02-14/"
  opt$software <- "DIANN"
  opt$dataset <- "DIANN_Result_WU321864/dataset.csv"
  opt$workunit <- "helloW"
} else if (FALSE) {
  opt$indir <- "outputs-20250407T1707/mzmine/"
  opt$software <- "MZMINEannot"
  opt$dataset <- "outputs-20250407T1707/bfabric/input_dataset.tsv"
  opt$workunit <- "helloWannot"
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


if (file.exists(ymlfile)) {
  GRP2 <- prolfquapp::get_config(ymlfile)
} else {
  GRP2 <- prolfquapp::make_DEA_config_R6(
    PATH = opt$outdir,
    ORDERID = opt$order,
    PROJECTID = opt$project,
    WORKUNITID = opt$workunit,
    application = opt$software,
    prefix = "QC"
  )
}

dir.create(GRP2$path)

logger::log_info(GRP2$get_zipdir())
if (!dir.exists(GRP2$get_zipdir())) {
  dir.create(GRP2$get_zipdir())
}

output_dir <- GRP2$get_zipdir()
path <- opt$indir

if (!file.exists(opt$dataset)) {
  stop("No annotation file found : ", opt$dataset)
}

annotation <- file.path(opt$dataset) |>
  prolfquapp::read_table_data() |>
  prolfquapp::read_annotation(QC = TRUE, repeated = FALSE)


result <- tryCatch(
  {
    # Attempt to run the function
    procsoft <- preprocess_software(
      opt$indir,
      annotation,
      preprocess_functions = prolfqua_preprocess_functions[[opt$software]],
      pattern_contaminants = (GRP2$processing_options$pattern_contaminants),
      pattern_decoys = GRP2$processing_options$pattern_decoys
    )
    # Return the result if successful
    list(value = procsoft, error = NULL, stack_trace = NULL)
  },
  error = function(e) {
    # On error, capture the stack trace as text
    stack_trace <- capture.output(traceback())
    # Return the error message and stack trace
    list(
      value = NULL,
      error = conditionMessage(e),
      stack_trace = paste(stack_trace, collapse = "\n")
    )
  }
)


if (!is.null(result$error)) {
  logger::log_error(result$error, "\n")
  logger::log_error("Stack trace:\n")
  logger::log_error(result$stack_trace, "\n")
  if (interactive()) {
    stop("error occured")
  } else {
    quit(save = "no", status = 1)
  }
} else {
  xd <- result$value$xd
  files <- result$value$files
}

xd$lfqdata$config$table$hierarchyDepth <- 1

GRP2$get_zipdir()

# QC_generator$debug("get_list")
# QC_generator$debug("get_protein_per_group_small_wide")
# QC_generator$debug("get_prot_IBAQ")
# QC_generator$undebug("initialize")

# xd$lfqdata$hierarchy_counts()

pap <- QC_generator$new(xd$lfqdata, xd$protein_annotation, GRP2)
# pap$get_protein_per_group_small_wide()
# pap$lfqdata

# pap$get_list()
# pap$get_peptides_wide()
# pap$get_list
# dd <- pap$get_prot_wide()
# pap$get_prot_IBAQ_wide()

# wirte parquet
pap$write_xlsx()

# library(arrow)
# arrow::write_parquet(pap$lfqdata$data, sink = "lfqdata.parquet")
# cfg <- prolfqua::R6_extract_values(pap$lfqdata$config)
# yaml::write_yaml(prolfqua::R6_extract_values(pap$lfqdata$config), "lfqdata.yaml")
# yaml::read_yaml("lfqdata.yaml")
# arrow::write_parquet(pap$lfqdata_peptide$data, sink = "lfqdata_peptide.parquet")
# yaml::write_yaml(prolfqua::R6_extract_values(pap$lfqdata_peptide$config), "lfqdata_peptide.yaml")

pap$render_QC_protein_abundances()
pap$render_sample_size_QC()
pap$render_index_html()
pap$render_index_md()
