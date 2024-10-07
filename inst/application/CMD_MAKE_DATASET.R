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
               help = "folder containing fasta file and output of the quantification software.",
               metavar = "string"),
  make_option( c("-d", "--dataset"), type = "character", default = "dataset.csv",
               help = "name of annotation file",
               metavar = "string"),
  make_option(c("-s", "--software"), type = "character", default = "DIANN",
              help = "possible options DIANN, FP_TMT, MAXQUANT, MSSTATS, FP_multisite, FP_combined_STY",
              metavar = "character"),
  make_option(c("-f", "--function"), type = "character", default = "",
              help = "possible options: packagename::processing_function",
              metavar = "character")

)

parser <- optparse::OptionParser(usage = "%prog config.yaml --software DIANN --indir .", option_list = option_list)
arguments <- optparse::parse_args(parser, positional_arguments = FALSE)
opt <- arguments
logger::log_info(prolfquapp::capture_output(quote(lobstr::tree(opt))))

if (FALSE) {
  opt$indir = "DIANN_19_all_18_50_50_MBR_v01/"
  opt$dataset = "dataset.xlsx"
  opt$software = "DIANN"
}
if (FALSE) {
  opt$indir = "20240729_093759_o35116_SN19_outliersRemoved_MSStats_Report"
  opt$dataset = "dataset_2.xlsx"
  opt$software = "MSSTATS"
}
logger::log_info(prolfquapp::capture_output(quote(lobstr::tree(opt))))



if (opt$software == "DIANN") {
  files <- prolfquapp::get_DIANN_files(opt$indir)
  logger::log_info("Files data: ", files$data)
  logger::log_info("Files fasta: ", files$fasta)
  data <- readr::read_tsv(files$data)
  logger::log_info("Files: ", files$data, " loaded. Starting filtering.")
  xx <- prolfquapp::diann_read_output(data, Lib.PG.Q.Value = 0.01, PG.Q.Value = 0.01)
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
  files <- prolfquapp::get_MQ_peptide_files(opt$indir)
  logger::log_info("Files data: ", files$data)
  peptide <- prolfquapp::tidyMQ_Peptides( files$data, proteotypic_only = TRUE)
  head(peptide)
  datasetannot <- data.frame(raw.file = unique(peptide$raw.file), name = NA, group = NA, subject = NA, CONTROL = NA)
  prolfquapp::write_annotation_file(datasetannot, opt$dataset)
} else if (opt$software == "MSSTATS") {
  files <- prolfquapp::get_MSstats_files(opt$indir)
  logger::log_info("Files data: ", paste(files$data, collapse = "; "))
  logger::log_info("Files fasta: ", paste0(files$fasta, collapse = "; "))
  logger::log_info("Read data:", files$data)

  msstats_df <- prolfquapp::read_table_data(files$data)
  datasetannot <- msstats_df |> dplyr::select(raw.file = "Run", "Group" = "Condition", "Subject" = "BioReplicate") |> dplyr::distinct()
  datasetannot$Control <- ""
  datasetannot <- datasetannot |> tidyr::unite("Name", "Group", "Subject", sep = "_", remove = FALSE)
  prolfquapp::write_annotation_file(datasetannot, opt$dataset)
} else if (opt$software == "FP_multisite") {
} else if (opt$software == "FP_combined_STY") {
  logger::log_info("FP_combined_STY")
  fff <- prolfquapp::get_FP_combined_STY_files(opt$indir)
  res_data <- prolfquapp::read_combined_STY_file(fff$data)
  logger::log_info("using file: ", paste(fff, collapse = ";"))

  manifest <- readr::read_tsv(fff$fp.manifest, col_names = FALSE)
  colnames(manifest) <- c("raw.file", "Name", "Experiment", "Data_type")
  manifest$Experiment <- NULL
  manifest$Data_type <- NULL

  res_data <- dplyr::inner_join(manifest, res_data, by = c(Name = "SampleName" ))

  datasetannot <- res_data |> dplyr::select(tidyselect::all_of(c("raw.file", "Name"))) |> dplyr::distinct()
  datasetannot$Group <- ""
  datasetannot$Subject <- ""
  datasetannot$Control <- ""
  prolfquapp::write_annotation_file(datasetannot, opt$dataset)

}else {
  logger::log_error("no such software :" , opt$software)
  stop("no such software.")
}
