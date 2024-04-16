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

#GRP2 <- prolfquapp::R6_extract_values(ymlfile, application = "DIANN")
GRP2 <- prolfquapp::make_DEA_config_R6(ZIPDIR = "testresult",PROJECTID = fgczProject,
                                       ORDERID = OIDfgcz)

dsf = file.path(path,"o34441_FPguiTMTphospho__Dataset_TotalNEnriched_better.tsv")
# @jg: add raw.file to annotation
#dsf <- readr::read_csv(dsf) # in BF it is csv

dsf <- readr::read_tsv(dsf)
dsf$condition <- NULL
dsf$genotype <- NULL
dsf$treatment <- NULL
dsf$SampleName <- NULL
dsf$channel <- dsf$sample
# now annotation contains 3 lists
# full annotation table
# all contrasts based on CONTROL column
# atable for config

annotation <- read_annotation(dsf)
annotation$annot

path = "./philosopher_prot/"
dir(path)
files <- get_FP_PSM_files(path = path) # @WeW fix here in function to specify path as argument as forseen but not working?


debug(preprocess_FP_PSM)
xd <- prolfquapp::preprocess_FP_PSM(quant_data = files$data,
                                    fasta_file = files$fasta,
                                    annotation = annotation,
                                    column_before_quants = "Mapped Proteins")


# tested till here!!

lfqdata <- xd$lfqdata
lfqdata$hierarchy_counts()
lfqdata$config$table$hierarchyDepth <- 3
lfqdata$config$table$hierarchy_keys_depth()

lfqdata$config$table$ident_Score
lfqdata$config$table$ident_qValue

ag <- lfqdata$get_Aggregator()

ag$sum_topN(N = 10000)

lfqdata <- ag$lfq_agg
lfqdata$data
lfqdata$response()
lfqdata$hierarchy_counts()

lfqdata$config$table$hierarchyDepth <- 1

logger::log_info("AGGREGATING PEPTIDE DATA!")
lfqdata <- prolfquapp::aggregate_data(lfqdata, agg_method =
                                        GRP2$processing_options$aggregate)
logger::log_info("data aggregated: {GRP2$pop$aggregate}.")
logger::log_info("END OF PROTEIN AGGREGATION")

lfqdata$hierarchy_counts()

logger::log_info("run analysis")
grp <- prolfquapp::generate_DEA_reports2(lfqdata, GRP2,
                                         xd$protein_annotation, annotation$contrasts)

logger::log_info("write results and html reports")
prolfquapp::write_DEA_all(grp, "Groups_vs_Controls", GRP2$zipdir , boxplot
                          = FALSE, markdown = "_Grp2Analysis_V2.Rmd")

logger::log_info("write results and summarized experiment")
SE <- prolfquapp::make_SummarizedExperiment(grp)
saveRDS(SE, file = file.path(GRP2$zipdir,
                             paste0("Results_DEA_WU", grp$project_spec$workunitID) ,
                             paste0("SummarizedExperiment",".rds") ))
