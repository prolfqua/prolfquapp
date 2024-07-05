if (!require("prolfqua", quietly = TRUE))
  remotes::install_github("prolfqua/prolfquapp", dependencies = TRUE)
if (!require("prolfqua", quietly = TRUE))
  remotes::install_github("prolfqua/prolfquapp", dependencies = TRUE)
if (!require("optparse", quietly = TRUE))
  install.packages("optparse", dependencies = TRUE)

logger::log_info("LIBRARY PATHS (.libPaths()):",paste(.libPaths(), collapse = "\n"))

library("optparse")
option_list <- list(
  make_option( c("-i", "--indir"), type = "character", default = ".",
               help = "folder containing fasta, diann-output.tsv and dataset.tsv file",
               metavar = "string"),
  make_option( c("-d", "--dataset"), type = "character", default = "dataset.csv",
               help = "name of annotation",
               metavar = "string"),
  make_option(c("-s", "--software"), type = "character", default = "DIANN",
              help = "possible options DIANN, FP_TMT, MAXQUANT",
              metavar = "character")
)


parser <- optparse::OptionParser(usage = "%prog config.yaml --software DIANN --indir .", option_list = option_list)
arguments <- optparse::parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
ymlfile <- arguments$args


if (opt$software == "DIANN") {
  files <- prolfquapp::get_DIANN_files(opt$indir)
  logger::log_info("Files data: ", files$data)
  logger::log_info("Files fasta: ", files$fasta)

  xd <- prolfquapp::preprocess_DIANN(
    quant_data = files$data,
    fasta_file = files$fasta,
    annotation = annotation,
    nrPeptides =  GRP2$processing_options$nr_peptides,
    q_value = 0.1,
    pattern_contaminants = GRP2$processing_options$pattern_contaminants,
    pattern_decoys = GRP2$processing_options$pattern_decoys
  )
} else if (opt$software == "FP_TMT") {
  files <- prolfquapp::get_FP_PSM_files(opt$indir)
  logger::log_info("Files data: ", files$data)
  logger::log_info("Files fasta: ", files$fasta)
  xd <- prolfquapp::preprocess_FP_PSM(
    quant_data = files$data,
    fasta_file = files$fasta,
    annotation = annotation,
    purity_threshold = 0.5,
    PeptideProphetProb = 0.9,
    pattern_contaminants = GRP2$processing_options$pattern_contaminants,
    pattern_decoys = GRP2$processing_options$pattern_decoys
  )
} else if (opt$software == "MAXQUANT") {
  files <- prolfquapp::get_MQ_peptide_files(opt$indir)
  logger::log_info("Files data: ", files$data)
  logger::log_info("Files fasta: ", files$fasta)

  xd <- prolfquapp::preprocess_MQ_peptide(
    quant_data = files$data,
    fasta_file = files$fasta,
    annotation = annotation,
    pattern_contaminants = GRP2$processing_options$pattern_contaminants,
    pattern_decoys = GRP2$processing_options$pattern_decoys
  )
} else if (opt$software == "FP_DDA") {

} else {
  logger::log_error("no such software :" , opt$software)
  stop("no such software.")
}
