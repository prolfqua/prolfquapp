if (!require("optparse", quietly = TRUE)) {
  install.packages("optparse", dependencies = TRUE)
}

logger::log_info("LIBRARY PATHS (.libPaths()):",paste(.libPaths(), collapse = "\n"))


option_list <- list(
  make_option(c("-i", "--indir"), type = "character", default = NULL,
              help = "folder containing fasta and diann-output files",
              metavar = "path"),
  make_option(c("-d", "--dataset"), type = "character", default = "dataset.csv",
              help = "file with annotation",
              metavar = "character"),
  make_option(c("-y", "--yaml"), type = "character", default = "config.yaml",
              help = "yaml configuration file",
              metavar = "character"),
  make_option("--workunit", type = "character"),
  make_option(c("-s", "--software"), type = "character", default = "DIANN",
              help = "possible options DIANN, FP_TMT, MAXQUANT",
              metavar = "character"),
  make_option(c("-p", "--outdir"), type = "character", default = NULL,
              help = "yaml configuration file",
              metavar = "path"),
  make_option(c("--libPath"), type = "character", default = NULL,
              help = " (optional) R library path",
              metavar = "string")
)


library(prolfquapp)

parser <- OptionParser(usage = "%prog file [options] ", option_list = option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
lobstr::tree(arguments)

ymlfile <- arguments$args
ymlfile <- if ( length(ymlfile) == 0 ) { opt$yaml } else { ymlfile }

if (FALSE) {
  opt$software = "FP_TMT"
  #opt$indir = "2513659"
  #opt$dataset = "dataset_compare_concentrations.csv"
  #opt$yaml = "configuration_compare_concentrations.yaml"
}


yamlfile <- file.path(ymlfile)
logger::log_info("YAML file read: ", yamlfile)
stopifnot(file.exists(yamlfile))


GRP2 <- prolfquapp::get_config(yamlfile)
GRP2$set_zipdir_name()

if (!is.null(opt$workunit)) {
  logger::log_info("setting workunit to: " ,  opt$workunit)
  GRP2$project_spec$workunit_Id <- opt$workunit
  GRP2$set_zipdir_name()
}

logger::log_info("Writing results to: " ,  GRP2$zipdir)
lobstr::tree(R6_extract_values(GRP2))

annotation <- file.path(opt$dataset) |>
  readr::read_csv() |> prolfquapp::read_annotation(prefix = GRP2$group)
logger::log_info("Contrasts: \n", paste(annotation$contrasts, collapse = "\n"))

logger::log_info("Software :", opt$software)
if (opt$software == "DIANN") {
  prolfquapp::copy_DEA_DIANN(run_script = FALSE)
  files <- prolfquapp::get_DIANN_files(opt$indir)
  xd <- prolfquapp::preprocess_DIANN(
    quant_data = files$data, fasta_file = files$fasta,
    annotation = annotation, nrPeptides =  GRP2$processing_options$nr_peptides,
    q_value = 0.1,pattern_contaminants = GRP2$processing_options$pattern_contaminants,
    pattern_decoys = GRP2$processing_options$pattern_decoys)
} else if (opt$software == "FP_TMT") {
  prolfquapp::copy_DEA_FragPipe_TMT(run_script = FALSE)
  files <- prolfquapp::get_FP_PSM_files(opt$indir)
  psm <- prolfquapp::tidy_FragPipe_psm(files$data)

  xd <- prolfquapp::preprocess_FP_PSM(
    quant_data = files$data,
    fasta_file = files$fasta,
    annotation = annotation,
    purity_threshold = 0.5,
    PeptideProphetProb = 0.9
  )

} else if (opt$software == "MAXQUANT") {
  stop("support coming soon.")
} else {
  stop("no such software.")
}



logger::log_info("AGGREGATING PEPTIDE DATA: {GRP2$pop$aggregate}.")
lfqdata <- prolfquapp::aggregate_data(xd$lfqdata, agg_method = GRP2$processing_options$aggregate)
logger::log_info("END OF PROTEIN AGGREGATION")

logger::log_info("RUN ANALYSIS")



grp <- prolfquapp::generate_DEA_reports2(lfqdata, GRP2, xd$protein_annotation, annotation$contrasts)

outdir <- prolfquapp::write_DEA_all(
  grp, boxplot = FALSE, markdown = "_Grp2Analysis_V2.Rmd")
ibaq <- compute_IBAQ_values(xd$lfqdata, xd$protein_annotation)
writexl::write_xlsx(ibaq$to_wide()$data, path = file.path(outdir, "IBAQ.xlsx"))

logger::log_info("write results and summarized experiment")
SE <- prolfquapp::make_SummarizedExperiment(grp)
saveRDS(SE, file = file.path( outdir, "SummarizedExperiment.rds"))

inputs <- file.path(
  GRP2$zipdir,
  paste0("Inputs_DEA_WU", GRP2$project_spec$workunitID))
dir.create(inputs)

prolfquapp::copy_DEA_DIANN(workdir = inputs, run_script = TRUE)
file.copy(c(files$data, files$fasta, yamlfile, opt$dataset), inputs)

GRP2$RES <- NULL
GRP2$pop <- NULL
yaml::write_yaml(prolfquapp::R6_extract_values(GRP2), file = file.path(inputs, "minimal.yaml"))


