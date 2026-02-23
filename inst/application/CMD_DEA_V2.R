if (!require("optparse", quietly = TRUE)) {
  install.packages("optparse", dependencies = TRUE)
}

option_list <- list(
  optparse::make_option(c("-i", "--indir"),
    type = "character", default = ".",
    help = "folder containing fasta and diann-output files",
    metavar = "path"
  ),
  optparse::make_option(c("-d", "--dataset"),
    type = "character", default = "dataset.csv",
    help = "file with annotation",
    metavar = "character"
  ),
  optparse::make_option(c("-y", "--yaml"),
    type = "character", default = "config.yaml",
    help = "yaml configuration file",
    metavar = "character"
  ),
  optparse::make_option(c("-w", "--workunit"),
    type = "character", default = NULL,
    help = "yaml configuration file",
    metavar = "character"
  ),
  optparse::make_option(c("-s", "--software"),
    type = "character", default = NULL,
    help = paste0("possible options: ", paste(names(prolfquapp::get_procfuncs()),
      collapse = ", "
    )),
    metavar = "character"
  ),
  optparse::make_option(c("-o", "--outdir"),
    type = "character", default = NULL,
    help = "output directory",
    metavar = "character"
  ),
  optparse::make_option(c("--libPath"),
    type = "character", default = NULL,
    help = " (optional) R library path",
    metavar = "string"
  )
)

parser <- optparse::OptionParser(usage = "%prog config.yaml --software DIANN --indir .", option_list = option_list)

if (length(commandArgs(TRUE)) == 0) {
  optparse::print_help(parser)
  # quit(status = 1)
}

arguments <- optparse::parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
ymlfile <- arguments$args

logger::log_appender(logger::appender_console)
logger::log_info("LIBRARY PATHS (.libPaths()):", paste(.libPaths(), collapse = "\n"))

prolfquapp::set_lib_path(opt$libPath)

library(prolfquapp)
logger::log_info("using : ", system.file(package = "prolfqua"))
logger::log_info("using : ", system.file(package = "prolfquapp"))

if (FALSE) {
  ymlfile <- "config.yaml"
  opt$indir <- "ptm_example-main/data_total/FP_22"
  opt$software <- "prolfquapp.FP_TMT"
  opt$dataset <- "dataset_with_contrasts.tsv"
  opt$workunit <- "total_proteome"
} else if (FALSE) {
  ymlfile <- "config.yaml"
  opt$indir <- "o38194_proteome_includeRepeatedSample"
  opt$software <- "prolfquapp.MSSTATS"
  opt$dataset <- "annotation.tsv"
  opt$workunit <- "total_proteome"
} else if (FALSE) {
  ymlfile <- "config.yaml"
  opt$indir <- "ptm_example-main/data_ptm/FP_22"
  opt$software <- "prolfquappPTMreaders.FP_singlesite"
  opt$dataset <- "dataset_with_contrasts.tsv"
  opt$workunit <- "singlesite_PTM"
} else if (FALSE) {
  opt$yaml <- "WU339865/config.yaml"
  opt$indir <- "WU339865"
  opt$software <- "prolfquapp.DIANN"
  opt$dataset <- "WU339865/dataset.csv"
  opt$workunit <- "xyz"
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

# if internal normalizaton set the protein identifier here
GRP2$processing_options$internal <- NULL
# GRP2$processing_options$internal <- "P01876"


dir.create(opt$outdir)
dir.create(GRP2$get_zipdir())

current_time <- Sys.time()

formatted_time <- format(current_time, "%Y%m%d%H%M")
logfile <- paste0("prolfqua_", formatted_time, ".log")
appender_combined <- logger::appender_tee(file.path(GRP2$get_zipdir(), logfile))
logger::log_appender(appender_combined)
logger::log_info(prolfquapp::capture_output(quote(lobstr::tree(opt))))
logger::log_info("Writing to output directory : ", GRP2$get_zipdir(), " and file :", logfile)

logger::log_info("prolfquapp paramters : ")
logger::log_info(prolfquapp::capture_output(quote(lobstr::tree(prolfqua::R6_extract_values(GRP2)))))


annotation <- file.path(opt$dataset) |>
  prolfquapp::read_table_data() |>
  prolfquapp::read_annotation(prefix = GRP2$group)

logger::log_info("ContrastNames: \n", paste(names(annotation$contrasts), collapse = "\n"))
logger::log_info("Contrast: \n", paste(annotation$contrasts, collapse = "\n"))

logger::log_info("Factors : ", paste(annotation$atable$factor_keys_depth(), collapse = "\n"))
prolfquapp::copy_DEA_Files()
logger::log_info("Software: ", opt$software)


prolfqua_preprocess_functions <- get_procfuncs()
if (!opt$software %in% names(prolfqua_preprocess_functions)) {
  logger::log_error(opt$software, " no in ", paste0(names(prolfqua_preprocess_functions)))
}


result <- tryCatch(
  {
    # Attempt to run the function
    procsoft <- preprocess_software(
      opt$indir,
      annotation,
      preprocess_functions = prolfqua_preprocess_functions[[opt$software]],
      pattern_contaminants = GRP2$processing_options$pattern_contaminants,
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


logger::log_info("Processing done:", opt$software)
logger::log_info(paste(c("Protein Annotation :\n", capture.output(print(xd$protein_annotation$get_summary()))), collapse = "\n"))

# ---- R6 Analysis Pipeline ----
logger::log_info("CREATING DEAnalyse R6 OBJECT")

deanalyse <- prolfquapp::DEAnalyse$new(
  xd$lfqdata,
  xd$protein_annotation,
  GRP2,
  annotation$contrasts
)

logger::log_info("CONTAMINANT/DECOY SUMMARY")
deanalyse$cont_decoy_summary()

logger::log_info("REMOVING CONTAMINANTS AND DECOYS")
deanalyse$remove_cont_decoy()

logger::log_info("AGGREGATING PEPTIDE DATA: {GRP2$processing_options$aggregate}.")
deanalyse$aggregate()

logger::log_info("TRANSFORMING DATA")
deanalyse$transform_data()

logger::log_info("BUILDING MODELS AND COMPUTING CONTRASTS")
deanalyse$get_contrasts_merged_protein()

logger::log_info("ANNOTATING CONTRASTS")
deanalyse$get_annotated_contrasts()

logger::log_info("END OF ANALYSIS")

# ---- R6 Report Pipeline ----
logger::log_info("CREATING DEAReportGenerator")

reporter <- prolfquapp::DEAReportGenerator$new(deanalyse, GRP2, name = "")

logger::log_info("Writing results to: ", GRP2$get_zipdir())

outdir <- reporter$write_DEA_all(
  boxplot = FALSE,
  markdown = "_Grp2Analysis_V2_R6.Rmd",
  markdown_qc = "_DiffExpQC_R6.Rmd"
)

# ---- Parquet + YAML export ----
arrow::write_parquet(
  deanalyse$lfq_data_transformed$data,
  sink = file.path(GRP2$get_result_dir(), "lfqdata_normalized.parquet")
)
cfg <- prolfqua::R6_extract_values(deanalyse$lfq_data_transformed$config)
yaml::write_yaml(cfg, file.path(GRP2$get_result_dir(), "lfqdata.yaml"))

# ---- IBAQ ----
lfqdataIB <- xd$lfqdata$get_subset(xd$protein_annotation$clean(
  contaminants = GRP2$processing_options$remove_cont,
  decoys = GRP2$processing_options$remove_decoys
))

# do not write when peptide level analysis
ibaq_file <- file.path(reporter$resultdir, paste0("IBAQ_", opt$workunit, ".xlsx"))
if (length(xd$lfqdata$config$hierarchy_keys_depth()) == 1) {
  ibaq <- compute_IBAQ_values(lfqdataIB, xd$protein_annotation)
  writexl::write_xlsx(
    ibaq$to_wide()$data,
    path = ibaq_file
  )
}

outdir$data_files$ibaq_file <- ibaq_file

# ---- Index HTML ----
prolfquapp::write_index_html(outdir, result_dir = reporter$ZIPDIR)

# ---- SummarizedExperiment ----
logger::log_info("Writing summarized experiment.")
SE <- reporter$make_SummarizedExperiment()
saveRDS(SE, file = file.path(reporter$resultdir, "SummarizedExperiment.rds"))

# ---- Archive inputs ----
logger::log_info("Creating directory with input files :", GRP2$get_input_dir())
dir.create(GRP2$get_input_dir())

prolfquapp::copy_DEA_Files(workdir = GRP2$get_input_dir())
prolfquapp::copy_shell_script(workdir = GRP2$get_input_dir())

file.copy(c(files$data, files$fasta, ymlfile, opt$dataset), GRP2$get_input_dir())

logger::log_info("Write yaml with parameters: ", file.path(GRP2$get_input_dir(), "minimal.yaml"))

GRP2$RES <- NULL
GRP2$pop <- NULL
yaml::write_yaml(prolfqua::R6_extract_values(GRP2), file = file.path(GRP2$get_input_dir(), "minimal.yaml"))
