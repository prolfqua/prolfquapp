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
                        help = "possible options: DIANN, FP_TMT, FP_multisite, MAXQUANT, MSSTATS",
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
  ymlfile <- "p35593_uniprot_paired/MultisiteUniprot.yaml"
  opt$dataset <- "p35593_uniprot_paired/dataset.xlsx"
  opt$indir <- "o35593_phos_siteLocFiltered_ionquant/"
}
if (FALSE) {
  ymlfile <- "xenbaseAnalysis/WholeProtXenbase.yaml"
  opt$dataset <- "xenbaseAnalysis/dataset.xlsx"
  opt$indir <- "o35593_prot_ionquant_xenbase/"
}
if (FALSE) {
  ymlfile <- "xenbaseAnalysis/MultisiteXenbase.yaml"
  opt$dataset <- "xenbaseAnalysis/dataset.xlsx"
  opt$indir <- "o35593_phos_siteLocFiltered_ionquant_xenbase/"
}

if (FALSE) {
  ymlfile <- "config.yaml"
  opt$dataset <- "dataset_with_contrasts.xlsx"
  opt$indir <- "DIANN_19_all_18_50_50_MBR_v01//"
}
if (FALSE) {
  ymlfile <- "config.yaml"
  opt$dataset <- "dataset.csv"
  opt$software <- "MSSTATS"
  opt$indir <- "2543975/"
}
if (FALSE) {
  ymlfile <- "config.yaml"
  opt$dataset <- "dataset.csv"
  opt$indir <- "."
}



ymlfile <- if ( length(ymlfile) == 0 ) { opt$yaml } else { ymlfile }

logger::log_info("YAML file read: ", ymlfile)
stopifnot(file.exists(ymlfile))

GRP2 <- prolfquapp::get_config(ymlfile)
res <- prolfquapp::sync_opt_config(opt, GRP2)
opt <- res$opt
GRP2 <- res$config

appender_combined <- logger::appender_tee(file.path(opt$outdir, "prolfqua.log"))
logger::log_appender(appender_combined)
logger::log_info(prolfquapp::capture_output(quote(lobstr::tree(opt))))

logger::log_info("Writing to output directory : ", opt$outdir)
dir.create(opt$outdir)


logger::log_info("prolfquapp paramters : ")
logger::log_info( prolfquapp::capture_output( quote(lobstr::tree(R6_extract_values(GRP2)))))


annotation <- file.path(opt$dataset) |>
  prolfquapp::read_table_data() |> prolfquapp::read_annotation(prefix = GRP2$group)
logger::log_info("Contrasts: \n", paste(annotation$contrasts, collapse = "\n"))

logger::log_info("Factors : ",paste(annotation$atable$factor_keys_depth(), collapse = "\n"))

prolfquapp::copy_DEA_Files()
logger::log_info("Software: ", opt$software)

if (opt$software == "DIANN") {
  files <- prolfquapp::get_DIANN_files(opt$indir)
  logger::log_info("Files data: ", paste(files$data, collapse = "; "))
  logger::log_info("Files fasta: ", paste0(files$fasta, collapse = "; "))
  xd <- prolfquapp::preprocess_DIANN(
    quant_data = files$data,
    fasta_file = files$fasta,
    annotation = annotation,
    q_value = 0.1,
    pattern_contaminants = GRP2$processing_options$pattern_contaminants,
    pattern_decoys = GRP2$processing_options$pattern_decoys
  )
} else if (opt$software == "FP_TMT") {
  files <- prolfquapp::get_FP_PSM_files(opt$indir)
  logger::log_info("Files data: ", paste(files$data, collapse = "; "))
  logger::log_info("Files fasta: ", paste0(files$fasta, collapse = "; "))
  xd <- prolfquapp::preprocess_FP_PSM(
    quant_data = files$data,
    fasta_file = files$fasta,
    annotation = annotation,
    purity_threshold = 0.5,
    PeptideProphetProb = 0.9,
    pattern_contaminants = GRP2$processing_options$pattern_contaminants,
    pattern_decoys = GRP2$processing_options$pattern_decoys
  )
} else if (opt$software == "FP_multisite") {
  files <- prolfquapp::get_FP_multi_site_files(opt$indir)
  logger::log_info("Files data: ", paste(files$data, collapse = "; "))
  logger::log_info("Files fasta: ", paste0(files$fasta, collapse = "; "))
  xd <- prolfquapp::preprocess_FP_multi_site(
    files$data[1],
    files$fasta,
    annotation,
    pattern_contaminants = GRP2$processing_options$pattern_contaminants,
    pattern_decoys = GRP2$processing_options$pattern_decoys)
} else if (opt$software == "FP_combined_STY"){

} else if (opt$software == "MAXQUANT") {
  files <- prolfquapp::get_MQ_peptide_files(opt$indir)
  logger::log_info("Files data: ", paste(files$data, collapse = "; "))
  logger::log_info("Files fasta: ", paste0(files$fasta, collapse = "; "))

  xd <- prolfquapp::preprocess_MQ_peptide(
    quant_data = files$data,
    fasta_file = files$fasta,
    annotation = annotation,
    pattern_contaminants = GRP2$processing_options$pattern_contaminants,
    pattern_decoys = GRP2$processing_options$pattern_decoys
  )
} else if (opt$software == "MSSTATS") {
  files <- prolfquapp::get_MSstats_files(opt$indir)
  logger::log_info("Files data: ", paste(files$data, collapse = "; "))
  logger::log_info("Files fasta: ", paste0(files$fasta, collapse = "; "))

  xd <- prolfquapp::preprocess_MSstats(
    quant_data = files$data,
    fasta_file = files$fasta,
    annotation = annotation,
    pattern_contaminants = GRP2$processing_options$pattern_contaminants,
    pattern_decoys = GRP2$processing_options$pattern_decoys
  )
} else {
  logger::log_error("no such software :" , opt$software)
  stop("no such software.")
}

logger::log_info(paste(c("Protein Annotation :\n",capture.output( print(xd$protein_annotation$get_summary()))),collapse = "\n"))
logger::log_info("AGGREGATING PEPTIDE DATA: {GRP2$processing_options$aggregate}.")
lfqdata <- prolfquapp::aggregate_data(xd$lfqdata, agg_method = GRP2$processing_options$aggregate)
logger::log_info("END OF PROTEIN AGGREGATION")
logger::log_info("RUN ANALYSIS")
grp <- prolfquapp::generate_DEA_reports2(lfqdata, GRP2, xd$protein_annotation, annotation$contrasts)
logger::log_info("Writing results to: " ,  GRP2$get_zipdir())

outdir <- prolfquapp::write_DEA_all(
  grp, boxplot = FALSE, markdown = "_Grp2Analysis_V2.Rmd")

lfqdataIB <- xd$lfqdata$get_subset(xd$protein_annotation$clean(
  contaminants = GRP2$processing_options$remove_cont,
  decoys = GRP2$processing_options$remove_decoys))

# do not write when peptide level analysis.
if (length(xd$protein_annotation$pID) == 1) {
  ibaq <- compute_IBAQ_values(lfqdataIB, xd$protein_annotation)
  writexl::write_xlsx(ibaq$to_wide()$data, path = file.path(grp$get_result_dir(), "IBAQ.xlsx"))
}

logger::log_info("Writing summarized experiment.")
SE <- prolfquapp::make_SummarizedExperiment(grp)

saveRDS(SE, file = file.path( grp$get_result_dir(), "SummarizedExperiment.rds"))


logger::log_info("Creating directory with input files :", GRP2$get_input_dir())
dir.create(GRP2$get_input_dir())

prolfquapp::copy_DEA_Files(workdir = GRP2$get_input_dir())
prolfquapp::copy_shell_script(workdir = GRP2$get_input_dir())

file.copy(c(files$data, files$fasta, ymlfile, opt$dataset), GRP2$get_input_dir())

logger::log_info("Wirte yaml with parameters: ", file.path(GRP2$get_input_dir(), "minimal.yaml"))

GRP2$RES <- NULL
GRP2$pop <- NULL
yaml::write_yaml(prolfquapp::R6_extract_values(GRP2), file = file.path(GRP2$get_input_dir(), "minimal.yaml"))


