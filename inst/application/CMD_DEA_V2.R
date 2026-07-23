if (!require("optparse", quietly = TRUE)) {
  install.packages("optparse", dependencies = TRUE)
}

option_list <- list(
  optparse::make_option(
    c("-i", "--indir"),
    type = "character",
    default = ".",
    help = "folder containing fasta and diann-output files",
    metavar = "path"
  ),
  optparse::make_option(
    c("-d", "--dataset"),
    type = "character",
    default = "dataset.csv",
    help = "file with annotation",
    metavar = "character"
  ),
  optparse::make_option(
    c("-y", "--yaml"),
    type = "character",
    default = "config.yaml",
    help = "yaml configuration file",
    metavar = "character"
  ),
  optparse::make_option(
    c("-w", "--workunit"),
    type = "character",
    default = NULL,
    help = "yaml configuration file",
    metavar = "character"
  ),
  optparse::make_option(
    c("-s", "--software"),
    type = "character",
    default = NULL,
    help = paste0(
      "possible options: ",
      paste(names(prolfquapp::get_procfuncs()), collapse = ", ")
    ),
    metavar = "character"
  ),
  optparse::make_option(
    c("-o", "--outdir"),
    type = "character",
    default = NULL,
    help = "output directory",
    metavar = "character"
  ),
  optparse::make_option(
    c("-m", "--model"),
    type = "character",
    default = NULL,
    help = paste0(
      "contrast facade method (overrides config). Options: ",
      paste(c(names(prolfqua::FACADE_REGISTRY), "saint"), collapse = ", ")
    ),
    metavar = "character"
  ),
  optparse::make_option(
    c("--libPath"),
    type = "character",
    default = NULL,
    help = " (optional) R library path",
    metavar = "string"
  ),
  optparse::make_option(
    c("--nr_peptides"),
    type = "integer",
    default = NULL,
    help = "minimum distinct peptides per protein (overrides config); >= 1",
    metavar = "N"
  ),
  optparse::make_option(
    c("--flat_outdir"),
    action = "store_true",
    default = FALSE,
    help = paste0(
      "write outputs directly into --outdir without a dated subdir ",
      "(default: dated subdir)"
    )
  )
)

parser <- optparse::OptionParser(
  usage = "%prog config.yaml --software DIANN --indir .",
  option_list = option_list
)

if (length(commandArgs(TRUE)) == 0) {
  optparse::print_help(parser)
}

arguments <- optparse::parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
ymlfile <- arguments$args

logger::log_appender(logger::appender_console)
prolfquapp::route_messages_to_logger()
logger::log_info(
  "LIBRARY PATHS (.libPaths()):",
  paste(.libPaths(), collapse = "\n")
)

prolfquapp::set_lib_path(opt$libPath)

library(prolfquapp)
logger::log_info("using : ", system.file(package = "prolfqua"))
logger::log_info("using : ", system.file(package = "prolfquapp"))

# Interactive debugging: set opt manually, then source from here
if (FALSE) {
  ymlfile <- "config.yaml"
  opt$indir <- "ptm_example-main/data_total/FP_22"
  opt$software <- "prolfquapp.FP_TMT"
  opt$dataset <- "dataset_with_contrasts.tsv"
  opt$workunit <- "total_proteome"
}

ymlfile <- if (length(ymlfile) == 0) {
  opt$yaml
} else {
  ymlfile
}

logger::log_info("YAML file read: ", ymlfile)
stopifnot(file.exists(ymlfile))

GRP2 <- prolfquapp::get_config(ymlfile)

res <- prolfquapp::sync_opt_config(opt, GRP2)
opt <- res$opt
GRP2 <- res$config


dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(GRP2$get_zipdir(), showWarnings = FALSE, recursive = TRUE)

current_time <- Sys.time()

formatted_time <- format(current_time, "%Y%m%d%H%M")
logfile <- paste0("prolfqua_", formatted_time, ".log")
appender_combined <- logger::appender_tee(file.path(GRP2$get_zipdir(), logfile))
logger::log_appender(appender_combined)
logger::log_info(prolfquapp::capture_output(quote(lobstr::tree(opt))))
logger::log_info(
  "Writing to output directory : ",
  GRP2$get_zipdir(),
  " and file :",
  logfile
)

logger::log_info("prolfquapp paramters : ")
logger::log_info(prolfquapp::capture_output(quote(lobstr::tree(prolfqua::R6_extract_values(
  GRP2
)))))


logger::log_info("Software: ", opt$software)

result <- tryCatch(
  prolfquapp::run_dea(
    indir = opt$indir,
    dataset = opt$dataset,
    software = opt$software,
    config = GRP2
  ),
  error = function(e) {
    stack_trace <- capture.output(traceback())
    logger::log_error(conditionMessage(e), "\n")
    logger::log_error("Stack trace:\n")
    logger::log_error(
      paste(stack_trace, collapse = "\n"),
      "\n"
    )
    if (interactive()) stop(e) else quit(save = "no", status = 1)
  }
)

prolfquapp:::write_dea_run_outputs(result, GRP2, opt, ymlfile)
