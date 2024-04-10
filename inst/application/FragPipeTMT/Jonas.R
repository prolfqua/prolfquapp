################################################################################
library(tidyverse)
library(prolfqua)
library(prolfquapp)

# params ideally taken from yaml
fgczProject <- "o34441"
OIDfgcz <- "o34441"
descri <- "FPguiTMTphospho_vsWTstim"
fracti <- "TotalProteome"
WUID <- "WU301538"

# data from
# # https://fgcz-bfabric.uzh.ch/bfabric/workunit/show.html?id=301538&tab=details


# v3
prolfquapp::copy_DEA_DIANN()

path = "."
ymlfile <- file.path(path,"config_FP_dea.yaml")
GRP2 <- prolfquapp::read_BF_yamlR6(ymlfile, application = "DIANN")

dsf = file.path(path,"o34441_FPguiTMTphospho__Dataset_TotalNEnriched_better.tsv")
# @jg: add raw.file to annotation
#dsf <- readr::read_csv(dsf) # in BF it is csv

dsf <- readr::read_tsv(dsf)

# now annotation contains 3 lists
# full annotation table
# all contrasts based on CONTROL column
# atable for config
annotation <- read_annotation(dsf)

path = "../o34441_FP_tsvFiles/philosopher_prot/"
files <- get_FP_PSM_files(path = path) # @WeW fix here in function to
specify path as argument as forseen but not working?

# quick hack
#mypsm <- dir(path = "../o34441_FP_tsvFiles/philosopher_prot/", pattern
  = "psm.tsv", recursive = TRUE, full.names = TRUE)
#files$data <- mypsm


# add option to pass to tidy_psm (last column before quant)
debug(preprocess_FP_PSM)

# hack (2024-04-09) -> in GUI the sample annotation in psm file is done
in "channel" with grouping.var
# in the annotation this is in Grouping.var
annotation$annot$channel <- annotation$annot$SampleName
# problem, the matching is done on raw.file that does not exist on psm!

xd <- prolfquapp::preprocess_FP_PSM(quant_data = files$data,
                                    fasta_file = files$fasta,
                                    annotation = annotation,
                                    column_before_quants = "Mapped
Proteins")

# tested till here!!



logger::log_info("AGGREGATING PEPTIDE DATA!")
lfqdata <- prolfquapp::aggregate_data(xd$lfqdata, agg_method =
                                        GRP2$processing_options$aggregate)
logger::log_info("data aggregated: {GRP2$pop$aggregate}.")
logger::log_info("END OF PROTEIN AGGREGATION")

logger::log_info("run analysis")
grp <- prolfquapp::generate_DEA_reports2(lfqdata, GRP2,
                                         xd$protein_annotation, annotation$contrasts)

logger::log_info("write results and html reports")
prolfquapp::write_DEA_all(grp[[1]], names(grp)[1], GRP2$zipdir , boxplot
                          = FALSE)

logger::log_info("write results and summarized experiment")
SE <- prolfquapp::make_SummarizedExperiment(grp[[1]])
saveRDS(SE, file = file.path(GRP2$zipdir, paste0("DE_", names(grp)[1]) ,
                             paste0("SummarizedExperiment",".rds") ))
