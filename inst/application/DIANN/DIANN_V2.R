library(prolfquapp)
prolfquapp::copy_DEA_DIANN()

path = "."
ymlfile <- file.path(path,"config.yaml")
GRP2 <- prolfquapp::read_BF_yamlR6(ymlfile, application = "DIANN")

dsf = file.path(path,"dataset.csv")
dsf <- readr::read_csv(dsf)
annotation <- prolfquapp::read_annotation(dsf)
files <- prolfquapp::get_DIANN_files(path)

xd <- prolfquapp::preprocess_DIANN(
  quant_data = files$data,
  fasta_file = files$fasta,
  annotation = annotation,
  nrPeptides = 1,
  q_value = 0.01)


logger::log_info("AGGREGATING PEPTIDE DATA!")
lfqdata <- prolfquapp::aggregate_data(xd$lfqdata, agg_method = GRP2$processing_options$aggregate)
logger::log_info("data aggregated: {GRP2$pop$aggregate}.")
logger::log_info("END OF PROTEIN AGGREGATION")

logger::log_info("run analysis")
grp <- prolfquapp::generate_DEA_reports2(lfqdata, GRP2, xd$protein_annotation, annotation$contrasts)


logger::log_info("write results and html reports")
prolfquapp::write_DEA_all(grp, "Groups_vs_Controls", GRP2$zipdir , boxplot = FALSE, markdown = "_Grp2Analysis_V2.Rmd")

logger::log_info("write results and summarized experiment")
SE <- prolfquapp::make_SummarizedExperiment(grp)

saveRDS(SE, file = file.path(GRP2$zipdir,
                             paste0("Results_DEA_WU", grp$project_spec$workunitID) ,
                             paste0("SummarizedExperiment",".rds") ))


### put all inputs into indir

inputs <- file.path(GRP2$zipdir,
                    paste0("Inputs_DEA_WU", GRP2$project_spec$workunitID))
dir.create(inputs)
prolfquapp::copy_DEA_DIANN(workdir = inputs, run_script = TRUE)
file.copy(files$data, inputs)
file.copy(files$fasta, inputs)


