library(rlang)
library(prolfqua)
library(prolfquapp)


# function
# input params
# dataframe or path to annotation file
# fasta file path
# path to DIANN result
rm(list = ls())
source("utils.R")
source("utils_annotation.R")

prolfquapp::copy_DEA_DIANN()

path = "WU298535"
ymlfile <- file.path(path,"config.yaml")
GRP2 <- prolfquapp::read_BF_yamlR6(ymlfile, application = "DIANN")

dsf = "WU298535/dataset.csv"
dsf <- readr::read_csv(dsf)
annotation <- read_annotation(dsf)

files <- get_files(path)
xd <- preprocess_DIANN(quant_data = files$data,
                       fasta_file = files$fasta,
                       annotation = annotation)

xd$lfqdata$hierarchy_counts()

logger::log_info("AGGREGATING PEPTIDE DATA!")
lfqdata <- prolfquapp::aggregate_data(xd$lfqdata, agg_method = GRP2$processing_options$aggregate)
logger::log_info("data aggregated: {GRP2$pop$aggregate}.")
logger::log_info("END OF PROTEIN AGGREGATION")


grp <- generate_DEA_reports2(lfqdata, GRP2, xd$protein_annotation, annotation$contrasts)

prolfquapp::write_DEA_all(grp, names(grp), GRP2$zipdir , boxplot = FALSE)

SE <- prolfquapp::make_SummarizedExperiment(grp)
saveRDS(SE, file = file.path(GRP2$zipdir, paste0("DE_", names(grp)[i]) , paste0("SummarizedExperiment",".rds") ))

