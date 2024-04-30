if (!require("prolfquapp", quietly = TRUE))
  remotes::install_github("prolfqua/prolfquapp", dependencies = TRUE)

library(prolfquapp)
prolfquapp::copy_DEA_DIANN()

path = "."
GRP2 <- if (file.exists(file.path(path,"config.yaml"))) {
  file.path(path,"config.yaml") |> prolfquapp::read_BF_yamlR6(application = "DIANN")
} else {
  prolfquapp::make_DEA_config_R6(
    ZIPDIR = "DEA", PROJECTID = "1" ,ORDERID = "2", WORKUNITID = "HelloWorld" )
}

annotation <- file.path(path,"dataset.csv") |>
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
file.copy(c(files$data, files$fasta, ymlfile, "dataset.csv"), inputs)



