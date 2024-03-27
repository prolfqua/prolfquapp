prolfquapp::copy_DEA_DIANN()

path = "."
ymlfile <- file.path(path,"config.yaml")
GRP2 <- prolfquapp::read_BF_yamlR6(ymlfile, application = "DIANN")

dsf = file.path(path,"dataset.csv")
dsf <- readr::read_csv(dsf)
annotation <- read_annotation(dsf)

files <- get_DIANN_files(path)
xd <- prolfquapp::preprocess_FP_psm(quant_data = files$data,
                                   fasta_file = files$fasta,
                                   annotation = annotation)


logger::log_info("AGGREGATING PEPTIDE DATA!")
lfqdata <- prolfquapp::aggregate_data(xd$lfqdata, agg_method = GRP2$processing_options$aggregate)
logger::log_info("data aggregated: {GRP2$pop$aggregate}.")
logger::log_info("END OF PROTEIN AGGREGATION")

logger::log_info("run analysis")
grp <- prolfquapp::generate_DEA_reports2(lfqdata, GRP2, xd$protein_annotation, annotation$contrasts)

logger::log_info("write results and html reports")
prolfquapp::write_DEA_all(grp[[1]], names(grp)[1], GRP2$zipdir , boxplot = FALSE)

logger::log_info("write results and summarized experiment")
SE <- prolfquapp::make_SummarizedExperiment(grp[[1]])
saveRDS(SE, file = file.path(GRP2$zipdir, paste0("DE_", names(grp)[1]) , paste0("SummarizedExperiment",".rds") ))

