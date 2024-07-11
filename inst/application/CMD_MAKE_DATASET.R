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
               help = "name of annotation file",
               metavar = "string"),
  make_option(c("-s", "--software"), type = "character", default = "DIANN",
              help = "possible options DIANN, FP_TMT, MAXQUANT",
              metavar = "character")
)

parser <- optparse::OptionParser(usage = "%prog config.yaml --software DIANN --indir .", option_list = option_list)
arguments <- optparse::parse_args(parser, positional_arguments = FALSE)
opt <- arguments
logger::log_info(prolfquapp::capture_output(quote(lobstr::tree(opt))))
if (FALSE) {
  opt$indir = "o35593_prot_ionquant"
  opt$dataset = "uniprotwhole/dataset.xlsx"
  opt$software = "FP_TMT"
}
logger::log_info(prolfquapp::capture_output(quote(lobstr::tree(opt))))
print(opt)



if (opt$software == "DIANN") {
  files <- prolfquapp::get_DIANN_files(opt$indir)
  logger::log_info("Files data: ", files$data)
  xx <- prolfquapp::diann_read_output(data, Q.Value = 0.01)
  datasetannot <- data.frame(raw.file = unique(xx$raw.file), Name = NA, Group = NA, Subject = NA, Control = NA)
  prolfquapp::write_annotation_file(datasetannot, opt$dataset)
} else if (opt$software == "FP_TMT") {
  files <- prolfquapp::get_FP_PSM_files(opt$indir)
  logger::log_info("Files data: ", files$data)
  x <- prolfquapp::tidy_FragPipe_psm(files$data)
  channel <- unique(x$data$channel)
  datasetannot <- data.frame(channel = channel, name = channel , group = NA, subject = NA, CONTROL = NA)
  prolfquapp::write_annotation_file(datasetannot, opt$dataset)
} else if (opt$software == "MAXQUANT") {
  #opt<- list()
  #opt$indir = "inst/application/MaxQuantDDA/C26109WU305220/DEA_20240703_PI_26109_OI_26109_WU_305220_robscale"
  files <- prolfquapp::get_MQ_peptide_files(opt$indir)
  logger::log_info("Files data: ", files$data)
  peptide <- prolfquapp::tidyMQ_Peptides( files$data, proteotypic_only = TRUE)
  head(peptide)
  datasetannot <- data.frame(raw.file = unique(peptide$raw.file), name = NA, group = NA, subject = NA, CONTROL = NA)
  prolfquapp::write_annotation_file(datasetannot, opt$dataset)
} else if (opt$software == "FP_DDA") {

} else {
  logger::log_error("no such software :" , opt$software)
  stop("no such software.")
}
