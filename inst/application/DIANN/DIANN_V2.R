if (!require("prolfqua", quietly = TRUE))
  remotes::install_github("prolfqua/prolfquapp", dependencies = TRUE)
if (!require("prolfqua", quietly = TRUE))
  remotes::install_github("prolfqua/prolfquapp", dependencies = TRUE)
if (!require("optparse", quietly = TRUE))
  install.packages("optparse", dependencies = TRUE)


library("optparse")
parser <- OptionParser()
parser <- add_option(parser, c("-d", "--dataset"), type = "character", default = "dataset.csv",
                     help = "file with annotation",
                     metavar = "character")
parser <- add_option(parser, c("-y", "--yaml"), type = "character", default = "config.yaml",
                     help = "yaml configuration file",
                     metavar = "character")
opt <- parse_args(parser)

library(prolfquapp)
prolfquapp::copy_DEA_DIANN()
path = "."

yamlfile <- file.path(path, opt$yaml)
GRP2 <- if (file.exists(yamlfile)) {
  yamlfile |> prolfquapp::read_BF_yamlR6(application = "DIANN")
} else {
  prolfquapp::make_DEA_config_R6(
    ZIPDIR = "DEA", PROJECTID = "1" ,ORDERID = "2", WORKUNITID = "HelloWorld" )
}

annotation <- file.path(path,opt$dataset) |>
  readr::read_csv() |> prolfquapp::read_annotation()

files <- prolfquapp::get_DIANN_files(path)

xd <- prolfquapp::preprocess_DIANN(
  quant_data = files$data, fasta_file = files$fasta,
  annotation = annotation, nrPeptides =  GRP2$processing_options$nr_peptides,
  q_value = 0.01)

logger::log_info("AGGREGATING PEPTIDE DATA: {GRP2$pop$aggregate}.")
lfqdata <- prolfquapp::aggregate_data(xd$lfqdata, agg_method = GRP2$processing_options$aggregate)
logger::log_info("END OF PROTEIN AGGREGATION")

logger::log_info("RUN ANALYSIS")
grp <- prolfquapp::generate_DEA_reports2(lfqdata, GRP2, xd$protein_annotation, annotation$contrasts)


outdir <- prolfquapp::write_DEA_all(
  grp, boxplot = FALSE, markdown = "_Grp2Analysis_V2.Rmd")



logger::log_info("write results and summarized experiment")
SE <- prolfquapp::make_SummarizedExperiment(grp)
saveRDS(SE, file = file.path( outdir, "SummarizedExperiment.rds"))

inputs <- file.path(
  GRP2$zipdir,
  paste0("Inputs_DEA_WU", GRP2$project_spec$workunitID))
dir.create(inputs)

prolfquapp::copy_DEA_DIANN(workdir = inputs, run_script = TRUE)
file.copy(c(files$data, files$fasta, yamlfile, "dataset.csv"), inputs)
yaml::write_yaml(GRP2, file = file.path(inputs, "minimal.yaml"))


