if (!require("optparse", quietly = TRUE)) {
  install.packages("optparse", dependencies = TRUE)
}

option_list <- list(
  optparse::make_option(c("-i", "--indir"), type = "character", default = ".",
                        help = "folder containing fasta and diann-output files",
                        metavar = "path"),
  optparse::make_option(c("-d", "--dataset"), type = "character", default = "dataset.csv",
                        help = "file with annotation",
                        metavar = "character"),
  optparse::make_option(c("-y", "--yaml"), type = "character", default = "config.yaml",
                        help = "yaml configuration file",
                        metavar = "character"),
  optparse::make_option(c("-w","--workunit"), type = "character", default = NULL,
                        help = "yaml configuration file",
                        metavar = "character"),
  optparse::make_option(c("-s", "--software"), type = "character", default = NULL,
                        help = "possible options: DIANN, FP_TMT, FP_multisite, FP_combined_STY, MAXQUANT, MSSTATS, MSSTATS_FP_DIA",
                        metavar = "character"),
  optparse::make_option(c("-o", "--outdir"), type = "character", default = NULL,
                        help = "output directory",
                        metavar = "character"),
  optparse::make_option(c("--libPath"), type = "character", default = NULL,
                        help = " (optional) R library path",
                        metavar = "string")
)

parser <- optparse::OptionParser(usage = "%prog config.yaml --software DIANN --indir .", option_list = option_list)
arguments <- optparse::parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
ymlfile <- arguments$args

logger::log_appender(logger::appender_console)
logger::log_info("LIBRARY PATHS (.libPaths()):",paste(.libPaths(), collapse = "\n"))

prolfquapp::set_lib_path(opt$libPath);

library(prolfquapp)
logger::log_info("using : ", system.file(package = "prolfqua"))
logger::log_info("using : ", system.file(package = "prolfquapp"))

if (FALSE) {
  ymlfile <- "p35593_uniprot_paired/WholeProtUniprot.yaml"
  opt$dataset <- "p35593_uniprot_paired/dataset.xlsx"
  opt$indir <- "o35593_prot_ionquant/"
}
if (FALSE) {
  ymlfile <- "minimalPhosphoAnalysis.yaml"
  opt$dataset <- "p35540_WU313409_annotationTableXLSX.xlsx"
  opt$indir <- "Fragpipe_o35920_phospho/"
  opt$software <- "FP_multisite"
}
if (FALSE) {
  ymlfile <- "minimal.yaml"
  opt$indir <- "."
  opt$software <- "MSSTATS"
  opt$dataset <- "dataset.csv"
}
if (FALSE) {
  ymlfile <- "config.yaml"
  opt$indir <- "."
  opt$software <- "DIANN"
  opt$dataset <- "dataset.csv"
}
if (FALSE) {
  ymlfile <- "FragPipe_f20/msstats20.yaml"
  opt$indir <- "FragPipe_f20"
  opt$software <- "MSSTATS_FP_DIA_PEPTIDE"
  opt$dataset <- "FragPipe_f20/dataset_msstats20_parallel.xlsx"
  opt$workunit <- "testing_peptide"
}
if (FALSE) {
  ymlfile <- "FragPipe_all/msstatsall.yml"
  opt$indir <- "FragPipe_f20_diann"
  opt$software <- "DIANN_PEPTIDE"
  opt$dataset <- "FragPipe_all/dataset_msstats_all_interaction_no_Subject.xlsx"
  opt$workunit <- "f20_diann_PEPTIDE_with_subject"

}

ymlfile <- if ( length(ymlfile) == 0 ) { opt$yaml } else { ymlfile }

logger::log_info("YAML file read: ", ymlfile)
stopifnot(file.exists(ymlfile))

#undebug(get_config)
#undebug(list_to_R6_app_config)
GRP2 <- prolfquapp::get_config(ymlfile)

res <- prolfquapp::sync_opt_config(opt, GRP2)
opt <- res$opt
GRP2 <- res$config

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
logger::log_info( prolfquapp::capture_output( quote(lobstr::tree(R6_extract_values(GRP2)))))


annotation <- file.path(opt$dataset) |>
  prolfquapp::read_table_data() |> prolfquapp::read_annotation(prefix = GRP2$group)

logger::log_info("Contrasts: \n", paste(annotation$contrasts, collapse = "\n"))

logger::log_info("Factors : ",paste(annotation$atable$factor_keys_depth(), collapse = "\n"))

prolfquapp::copy_DEA_Files()
logger::log_info("Software: ", opt$software)

result <- tryCatch({
  # Attempt to run the function
  procsoft <- preprocess_software(
    opt$indir,
    annotation,
    prolfquapp::prolfq_preprocess_functions,
    pattern_contaminants = GRP2$processing_options$pattern_contaminants,
    pattern_decoys = GRP2$processing_options$pattern_decoys,
    software = opt$software
  )
  # Return the result if successful
  list(value = procsoft, error = NULL, stack_trace = NULL)
}, error = function(e) {
  # On error, capture the stack trace as text
  stack_trace <- capture.output(traceback())
  # Return the error message and stack trace
  list(
    value = NULL,
    error = conditionMessage(e),
    stack_trace = paste(stack_trace, collapse = "\n")
  )
})

if (!is.null(result$error)) {
  logger::log_error(result$error, "\n")
  logger::log_error("Stack trace:\n")
  logger::log_error(result$stack_trace, "\n")
  stop("error occured")
} else {
  xd <- result$value$xd
  files <- result$value$files
}


logger::log_info("Processing done:", opt$software)
logger::log_info(paste(c("Protein Annotation :\n",capture.output( print(xd$protein_annotation$get_summary()))),collapse = "\n"))
logger::log_info("AGGREGATING PEPTIDE DATA: {GRP2$processing_options$aggregate}.")
lfqdata <- prolfquapp::aggregate_data(xd$lfqdata, agg_method = GRP2$processing_options$aggregate)
logger::log_info("END OF PROTEIN AGGREGATION")
logger::log_info("RUN ANALYSIS")
grp <- prolfquapp::generate_DEA_reports2(
  lfqdata,
  GRP2,
  xd$protein_annotation,
  annotation$contrasts)
logger::log_info("Writing results to: " ,  GRP2$get_zipdir())

outdir <- prolfquapp::write_DEA_all(
  grp, name = "", boxplot = FALSE, markdown = "_Grp2Analysis_V2.Rmd")

lfqdataIB <- xd$lfqdata$get_subset(xd$protein_annotation$clean(
  contaminants = GRP2$processing_options$remove_cont,
  decoys = GRP2$processing_options$remove_decoys))

# do not write when peptide level analysis
if (length(xd$protein_annotation$pID) == 1) {
  ibaq <- compute_IBAQ_values(lfqdataIB, xd$protein_annotation)
  writexl::write_xlsx(
    ibaq$to_wide()$data,
    path = file.path(grp$get_result_dir(), paste0("IBAQ_",opt$workunit,".xlsx")))
}

logger::log_info("Writing summarized experiment.")
SE <- prolfquapp::make_SummarizedExperiment(grp)
saveRDS(SE, file = file.path( grp$get_result_dir(), "SummarizedExperiment.rds"))


logger::log_info("Creating directory with input files :", GRP2$get_input_dir())
dir.create(GRP2$get_input_dir())

prolfquapp::copy_DEA_Files(workdir = GRP2$get_input_dir())
prolfquapp::copy_shell_script(workdir = GRP2$get_input_dir())

file.copy(c(files$data, files$fasta, ymlfile, opt$dataset), GRP2$get_input_dir())

logger::log_info("Write yaml with parameters: ", file.path(GRP2$get_input_dir(), "minimal.yaml"))

GRP2$RES <- NULL
GRP2$pop <- NULL
yaml::write_yaml(prolfquapp::R6_extract_values(GRP2), file = file.path(GRP2$get_input_dir(), "minimal.yaml"))


