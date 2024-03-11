#library(rlang)
#library(prolfqua)
#library(prolfquapp)


# function
# input params
# dataframe or path to annotation file
# fasta file path
# path to DIANN result

prolfquapp::copy_DEA_DIANN()

path = "WU298535"
ymlfile <- file.path(path,"config.yaml")
GRP2 <- prolfquapp::read_BF_yamlR6(ymlfile, application = "DIANN")

dsf = "WU298535/dataset.csv"
dsf <- readr::read_csv(dsf)
annotation <- read_annotation(dsf)
annotation


files <- get_DIANN_files(path)
xd <- prolfquapp::preprocess_DIANN(quant_data = files$data,
                       fasta_file = files$fasta,
                       annotation = annotation)

xd$lfqdata$hierarchy_counts()

logger::log_info("AGGREGATING PEPTIDE DATA!")
lfqdata <- prolfquapp::aggregate_data(xd$lfqdata, agg_method = GRP2$processing_options$aggregate)
logger::log_info("data aggregated: {GRP2$pop$aggregate}.")
logger::log_info("END OF PROTEIN AGGREGATION")


grp <- prolfquapp::generate_DEA_reports2(lfqdata, GRP2, xd$protein_annotation, annotation$contrasts)

prolfquapp::write_DEA_all(grp[[1]], names(grp)[1], GRP2$zipdir , boxplot = FALSE)

SE <- prolfquapp::make_SummarizedExperiment(grp[[1]])
saveRDS(SE, file = file.path(GRP2$zipdir, paste0("DE_", names(grp)[1]) , paste0("SummarizedExperiment",".rds") ))

