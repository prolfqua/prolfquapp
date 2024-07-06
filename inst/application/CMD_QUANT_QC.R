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
  make_option( c("-p", "--projectId"), type = "character", default = "1234",
               help = "your project identifier",
               metavar = "string"),
  make_option( c("-w", "--workunitId"), type = "character", default = "4321",
               help = "workunit identifier",
               metavar = "string"),
  make_option( c("-d", "--dataset"), type = "character", default = "dataset.csv",
               help = "name of annotation",
               metavar = "string"),

  make_option( c("-o", "--outdir"), type = "character", default = "qc_dir",
               help = "folder to write the results to.",
               metavar = "string"),
  make_option(c("-s", "--software"), type = "character", default = "DIANN",
              help = "possible options DIANN, FP_TMT, MAXQUANT",
              metavar = "character"),
  make_option(c("--libPath"), type = "character", default = NULL,
              help = " (optional) R library path",
              metavar = "string")
)

parser <- OptionParser(usage = "%prog --indir . ", option_list = option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
lobstr::tree(arguments)

opt <- arguments$options

if (FALSE) {
  opt$indir <- "2532162/"
}

# set library path
prolfquapp::set_lib_path(opt$libPath);
library(prolfquapp)
library(logger)

logger::log_info("using : ", system.file(package = "prolfqua"))
logger::log_info("using : ", system.file(package = "prolfquapp"))

GRP2 <- prolfquapp::make_DEA_config_R6(
  PATH = opt$outdir,
  ORDERID = opt$projectId,
  PROJECTID =  opt$projectId,
  WORKUNITID = opt$workunitId,
  application = opt$software)

dir.create(GRP2$path)

logger::log_info(GRP2$get_zipdir())
if (!dir.exists(GRP2$get_zipdir())) {
  dir.create(GRP2$get_zipdir())
}

output_dir <- GRP2$get_zipdir()
path <- opt$indir

annotfile <- file.path(path, opt$dataset)
if (!file.exists(annotfile)) {stop("No annotation file found : ",annotfile)}
annotation <- file.path(annotfile) |>
  readr::read_csv() |> prolfquapp::read_annotation(QC = TRUE)
names(annotation)


if (opt$software == "DIANN") {
  files <- prolfquapp::get_DIANN_files(path)
  xd <- prolfquapp::preprocess_DIANN(
    quant_data = files$data,
    fasta_file = files$fasta,
    annotation = annotation,
    q_value = 0.01)
} else if (opt$software == "FP_TMT") {
  files <- prolfquapp::get_FP_PSM_files(path)
  xd <- prolfquapp::preprocess_FP_PSM(
    quant_data = files$data,
    fasta_file = files$fasta,
    annotation = annotation
  )
} else if (opt$software == "MAXQUANT") {

} else {
  stop("unknown software : ", opt$software)
}

pap <- ProteinAbundanceProcessor$new(xd$lfqdata, xd$protein_annotation, GRP2)

pap$write_xlsx()
pap$render_QC_protein_abundances()
pap$render_sample_size_QC()


