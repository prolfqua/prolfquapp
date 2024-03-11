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

GRP2 <- prolfquapp::make_DEA_config_R6()
GRP2$processing_options$aggregate
GRP2$project_spec
GRP2$zipdir
GRP2$software

path = "WU298535"


files <- get_files(path)

dsf = "WU298535/dataset.csv"
dsf <- readr::read_csv(dsf)
annotation <- read_annotation(dsf)

xd <- preprocess_DIANN(quant_data = files$data,
                       fasta_file = files$fasta,
                       annotation = annotation)

xd$lfqdata$hierarchy_counts()

logger::log_info("AGGREGATING PEPTIDE DATA!")
lfqdata <- prolfquapp::aggregate_data(xd$lfqdata, agg_method = GRP2$processing_options$aggregate)
logger::log_info("data aggregated: {GRP2$pop$aggregate}.")
logger::log_info("END OF PROTEIN AGGREGATION")


grp <- generate_DEA_reports2(lfqdata, GRP2, xd$protein_annotation, annotation$contrasts)

prolfquapp::copy_DEA_DIANN()
for (i in seq_along(grp)) {
  prolfquapp::write_DEA_all(grp[[i]], names(grp)[i], GRP2$zipdir , boxplot = FALSE)
}

for (i in seq_along(grp)) {
  SE <- prolfquapp::make_SummarizedExperiment(grp[[i]])
  saveRDS(SE, file = file.path(GRP2$zipdir, paste0("DE_", names(grp)[i]) , paste0("SummarizedExperiment",".rds") ))
}

