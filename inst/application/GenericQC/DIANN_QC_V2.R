if (!require("prolfqua", quietly = TRUE))
  remotes::install_github("prolfqua/prolfquapp", dependencies = TRUE)
if (!require("prolfqua", quietly = TRUE))
  remotes::install_github("prolfqua/prolfquapp", dependencies = TRUE)
if (!require("optparse", quietly = TRUE))
  install.packages("optparse", dependencies = TRUE)

logger::log_info("LIBRARY PATHS (.libPaths()):",paste(.libPaths(), collapse = "\n"))


library("optparse")
parser <- OptionParser()
parser <- add_option(parser, c("-i", "--indir"), type = "character", default = ".",
                     help = "folder containing fasta, diann-output.tsv and dataset.tsv file",
                     metavar = "string")
parser <- add_option(parser, c("-p", "--projectId"), type = "character", default = "1234",
                     help = "your project identifier",
                     metavar = "string")
parser <- add_option(parser, c("-w", "--workunitId"), type = "character", default = "4321",
                     help = "workunit identifier",
                     metavar = "string")
parser <- add_option(parser, c("-d", "--dataset"), type = "character", default = "dataset.csv",
                     help = "name of annotation",
                     metavar = "string")
parser <- add_option(parser, c("-o", "--output"), type = "character", default = "qc_dir",
                     help = "folder to write the results to.",
                     metavar = "string")
parser <- add_option(parser, c("--libPath"), type = "character", default = NULL,
                     help = " (optional) R library path",
                     metavar = "string")

opt <- parse_args(parser)
opt$indir <- "2494004"

# set library path
libPath <- opt$libPath

if (!is.null(libPath) && dir.exists(libPath) ) {
  logger::log_info(paste("Setting libPath:", libPath, collapse = " ;"))
  .libPaths(libPath)
  logger::log_info(.libPaths(), sep = "\n")
}


# this must be executed after the libPath is modified.
library(dplyr)
library(prolfquapp)
library(logger)


GRP2 <- prolfquapp::make_DEA_config_R6(ZIPDIR = opt$output, ORDERID = opt$projectId, PROJECTID =  opt$projectId, WORKUNITID = opt$workunitId )
if (!dir.exists(GRP2$zipdir)) {
  dir.create(GRP2$zipdir)
}

output_dir <- GRP2$zipdir


path <- opt$indir
files <- prolfquapp::get_DIANN_files(path)
xx <- get_annot_from_fasta(files$fasta)

annotfile <- file.path(path, opt$dataset)
if(!file.exists(annotfile)) {stop("No annotation file found",annotfile)}
annotation <- file.path(annotfile) |>
  readr::read_csv() |> prolfquapp::read_annotation(QC = TRUE)

xd <- prolfquapp::preprocess_DIANN(
  quant_data = files$data,
  fasta_file = files$fasta,
  annotation = annotation,
  nrPeptides =  GRP2$processing_options$nr_peptides,
  q_value = 0.01)

TABLES2WRITE <- list()
TABLES2WRITE$peptide_wide <- dplyr::left_join(xd$protein_annotation$row_annot,
                                       xd$lfqdata$to_wide()$data,
                                       multiple = "all")

TABLES2WRITE$annotation <- xd$lfqdata$factors()
lfqdataProt <- prolfquapp::aggregate_data(xd$lfqdata, agg_method = "medpolish")
TABLES2WRITE$proteins_wide <- dplyr::left_join(xd$protein_annotation$row_annot,
                                        lfqdataProt$to_wide()$data,
                                        multiple = "all")


lfqdataProtIBAQ <- compute_IBAQ_values(xd$lfqdata, xd$protein_annotation)

lfqdataProtIBAQ$hierarchy_counts()
summarizer <- lfqdataProtIBAQ$get_Summariser()
precabund <- summarizer$percentage_abundance()

precabund <- dplyr::inner_join(
  xd$protein_annotation$row_annot,
  precabund,
  multiple = "all",
  by = lfqdataProtIBAQ$config$table$hierarchy_keys_depth())

TABLES2WRITE$proteinAbundances <- precabund
TABLES2WRITE$IBAQ_abundances <-
  dplyr::left_join(xd$protein_annotation$row_annot,
                   lfqdataProtIBAQ$to_wide()$data,
                   multiple = "all")

writexl::write_xlsx(TABLES2WRITE, path = file.path(output_dir, "proteinAbundances.xlsx"))

file.copy(system.file("application/GenericQC/QC_ProteinAbundances.Rmd", package = "prolfquapp"),
          to = output_dir, overwrite = TRUE)


if (!is.null(lfqdataProtIBAQ)) {
  rmarkdown::render(file.path(output_dir,"QC_ProteinAbundances.Rmd"),
                    params = list(lfqdataProt = lfqdataProtIBAQ,
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
