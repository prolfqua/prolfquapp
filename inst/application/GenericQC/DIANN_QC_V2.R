if (!require("prolfqua", quietly = TRUE))
  remotes::install_github("prolfqua/prolfquapp", dependencies = TRUE)
if (!require("prolfqua", quietly = TRUE))
  remotes::install_github("prolfqua/prolfquapp", dependencies = TRUE)
if (!require("optparse", quietly = TRUE))
  install.packages("optparse", dependencies = TRUE)


library("optparse")
parser <- OptionParser()
parser <- add_option(parser, c("-d", "--indir"), type = "character", default = ".",
                     help = "folder containing fasta, diann-output.tsv and dataset.tsv file",
                     metavar = "string")
parser <- add_option(parser, c("-p", "--projectId"), type = "character", default = "1234",
                     help = "your project identifier",
                     metavar = "string")
parser <- add_option(parser, c("-w", "--workunitId"), type = "character", default = "4321",
                     help = "workunit identifier",
                     metavar = "string")
parser <- add_option(parser, c("-w", "--dataset"), type = "character", default = "dataset.csv",
                     help = "file annotation",
                     metavar = "string")
parser <- add_option(parser, c("-o", "--output"), type = "character", default = "qc_dir",
                     help = "folder to write the results to.",
                     metavar = "string")
parser <- add_option(parser, c("--libPath"), type = "character", default = NULL,
                     help = " (optional) R library path",
                     metavar = "string")

opt <- parse_args(parser)


# set library path
libPath <- opt$libPath

if (!is.null(libPath) && dir.exists(libPath) ) {
  print(paste("Setting libPath:", libPath, collapse = " ;"))
  .libPaths(libPath)
  cat(.libPaths(), sep = "\n")
}


# this must be executed after the libPath is modified.
library(dplyr)
library(prolfquapp)
library(logger)


GRP2 <- prolfquapp::make_DEA_config_R6(ZIPDIR = opt$output, PROJECTID =  opt$projectId, WORKUNITID = opt$workunitId )
if (!dir.exists(GRP2$zipdir)) {
  dir.create(GRP2$zipdir)
}

output_dir <- GRP2$zipdir

source("utilities.R")
files <- find_default_files(opt$indir)
path <- opt$indir
annotation <- file.path(path,files$dataset) |>
  readr::read_csv() |> prolfquapp::read_annotation(QC = TRUE)

undebug(get_annot_from_fasta)
xd <- prolfquapp::preprocess_DIANN(
  quant_data = files$data,
  fasta_file = files$fasta,
  annotation = annotation,
  nrPeptides =  GRP2$processing_options$nr_peptides,
  q_value = 0.01)

head(xd$protein_annotation$row_annot)

xd$lfqdata$factors()


TABLES2WRITE <- list()
TABLES2WRITE$peptide_wide <- dplyr::left_join(xd$protein_annotation$row_annot,
                                       xd$lfqdata$to_wide()$data,
                                       multiple = "all")

TABLES2WRITE$annotation <- xd$lfqdata$factors()


lfqdataProt <- prolfquapp::aggregate_data(xd$lfqdata, agg_method = "medpolish")
TABLES2WRITE$proteins_wide <- dplyr::left_join(xd$protein_annotation$row_annot,
                                        lfqdataProt$to_wide()$data,
                                        multiple = "all")


lfqdataProtTotal <- prolfquapp::aggregate_data(xd$lfqdata, agg_method = "topN", N = 10000)
lfqdataProtTotal$data <- inner_join(lfqdataProtTotal$data , select(xd$protein_annotation$row_annot, protein_Id, protein_length, nr_tryptic_peptides))
lfqdataProtTotal$data <- lfqdataProtTotal$data |> mutate(IBAQV_proteinLength = srm_sum_N / protein_length)
lfqdataProtTotal$data <- lfqdataProtTotal$data |> mutate(IBAQV = srm_sum_N / nr_tryptic_peptides)

if (FALSE) {
  par(mfrow = c(2,3))
  plot( lfqdataProtTotal$data$srm_sum_N, lfqdataProtTotal$data$IBAQV, pch = ".", log = "xy")
  plot( lfqdataProtTotal$data$srm_sum_N, lfqdataProtTotal$data$IBAQV_proteinLength, pch = ".", log = "xy")
  plot( lfqdataProtTotal$data$IBAQV, lfqdataProtTotal$data$IBAQV_proteinLength, pch = ".", log = "xy")
  plot( lfqdataProtTotal$data$IBAQV, lfqdataProtTotal$data$IBAQV_proteinLength, pch = ".", log = "")
  plot(lfqdataProtTotal$data$protein_length, lfqdataProtTotal$data$nr_tryptic_peptides)
  xx <- with(lfqdataProtTotal$data, lm(nr_tryptic_peptides ~ protein_length))
  abline(xx, col = 2)
}

lfqdataProtTotal$config$table$set_response("IBAQV")

lfqdataProtTotal$hierarchy_counts()
summarizer <- lfqdataProtTotal$get_Summariser()
precabund <- summarizer$percentage_abundance()

precabund <- dplyr::inner_join(
  xd$protein_annotation$row_annot,
  precabund,
  multiple = "all",
  by = lfqdataProtTotal$config$table$hierarchy_keys_depth())

TABLES2WRITE$proteinAbundances <- precabund
TABLES2WRITE$IBAQ_abundances <-
  dplyr::left_join(xd$protein_annotation$row_annot,
                   lfqdataProtTotal$to_wide()$data,
                   multiple = "all")

writexl::write_xlsx(TABLES2WRITE, path = file.path(output_dir, "proteinAbundances.xlsx"))

file.copy(system.file("application/GenericQC/QC_ProteinAbundances.Rmd", package = "prolfquapp"),
          to = output_dir, overwrite = TRUE)


if (!is.null(lfqdataProt)) {
  rmarkdown::render(file.path(output_dir,"QC_ProteinAbundances.Rmd"),
                    params = list(lfqdataProt = lfqdataProt,
                                  precabund = precabund,
                                  project_info = GRP2$project_spec,
                                  factors = TRUE),
                    output_file = "proteinAbundances.html")

} else {
  str <- c("<!DOCTYPE html>",
           "<html>",
           "<head>","<title>",
           "There is a problem",
           "</title>","</head>",
           "<body>",
           "<h1>",
           paste0("the input file :" , files$data , " is empty"),
           "</h1>",
           "</body>",
           "</html>")
  cat(str, file = file.path(output_dir,"proteinAbundances.html"), sep = "\n")

}


if (nrow(lfqdataProt$factors()) > 1) {
  file.copy(system.file("application/GenericQC/QCandSSE.Rmd", package = "prolfquapp"),
            to = output_dir, overwrite = TRUE)
  rmarkdown::render(file.path(output_dir,"QCandSSE.Rmd"),
                    params = list(data = lfqdataProt$data,
                                  configuration = lfqdataProt$config,
                                  project_conf = GRP2$project_spec,
                                  pep = FALSE),
                    output_file = "QC_sampleSizeEstimation.html"
  )
} else{
  message("only a single sample: ", nrow(lfqdataProt$factors()))
}
