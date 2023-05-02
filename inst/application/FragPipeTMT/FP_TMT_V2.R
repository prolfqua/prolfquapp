# This script starts from psm.tsv
# does basic filtering using purity.
#
library(tidyverse)
GRP2 <- prolfquapp::make_DEA_config()
###
dir.create(GRP2$zipdir)
###
# reading foreign data
REPEATED <- TRUE

# sanitize peptide csv.
psm_file <- dir(path = ".", pattern = "psm.tsv", recursive = TRUE, full.names = TRUE)
xx <- readr::read_tsv(psm_file)
length(grep("c",xx$`Modified Peptide`, value = TRUE))

xa <- xx$`Assigned Modifications`
tmp <- gsub(" ", "", unlist(str_split(xa, ",")))
tmp <- gsub("^[0-9]+","", tmp)
table(tmp)


scores <- xx |> select(all_of(c("Expectation","Hyperscore","Nextscore","PeptideProphet Probability")))
image(cor(scores))

psm <- prolfqua::tidy_FragPipe_psm(psm_file)


psm$qValue <- 1 - psm$PeptideProphet.Probability

ds_file <- "dataset_mock_H2O2.csv"
ds_file <- "dataset_wt_KO.csv"
annot <- read.csv(ds_file)
GRP2 <- prolfquapp::make_DEA_config(ZIPDIR = gsub(".csv", "", ds_file))


GRP2 <- prolfquapp::dataset_extract_contrasts(annot, GRP2)

channelCol  <- grep("channel", names(annot), ignore.case = TRUE, value = TRUE)
if (channelCol != "channel") {
  annot[["channel"]] <- annot[[channelCol]]
  annot[[channelCol]] <- NULL
}

nr <- sum(annot$channel %in% unique(psm$channel))
logger::log_info("nr : ", nr, " files annotated")
psm <- dplyr::inner_join(annot, psm, multiple = "all")

# Setup configuration

atable <- prolfqua::AnalysisTableAnnotation$new()

atable$ident_Score = "PeptideProphet.Probability"
atable$ident_qValue = "qValue"
atable$fileName = "channel"
atable$hierarchy[["protein_Id"]] <- c("Protein")
atable$hierarchy[["peptide_Id"]] <- c("Peptide")
atable$hierarchy[["mod_peptide_Id"]] <- c("Modified.Peptide","Assigned.Modifications")
atable$hierarchy[["Spectrum"]] <- c("Spectrum")

#
tmp <- prolfquapp::dataset_set_factors(atable, psm)
atable <- tmp$atable
atable$factors
psm <- tmp$msdata

head(psm)

# CREATE protein annotation.
prot_annot <- prolfquapp::dataset_protein_annot(
  psm,
  "Protein",
  protein_annot = "Protein.Description",
  more_columns = "nrPeptides")


atable$set_response("abundance")
# Preprocess data - aggregate proteins.

config <- prolfqua::AnalysisConfiguration$new(atable)
adata <- prolfqua::setup_analysis(psm, config)
colnames(adata)
lfqdata <- prolfqua::LFQData$new(adata, config)
lfqdata$hierarchy_counts()
lfqdata$remove_small_intensities(threshold = 1)
lfqdata$hierarchy_counts()
xx <- lfqdata$summarize_hierarchy()
xx$Spectrum_n |> max()

lfqdata$config$table$hierarchyDepth <- 3
lfqdata$config$table$hierarchy_keys_depth()

lfqdata$config$table$ident_Score
lfqdata$config$table$ident_qValue

ag <- lfqdata$get_Aggregator()

ag$sum_topN(N = 1000)

lfqdata <- ag$lfq_agg
lfqdata$data
lfqdata$response()
lfqdata$hierarchy_counts()


#logger::log_info("AGGREGATING PEPTIDE DATA!")
#lfqdata <- prolfquapp::aggregate_data(lfqdata, agg_method = GRP2$pop$aggregate)
#logger::log_info("data aggregated: {GRP2$pop$aggregate}.")
#lfqdata$factors()

logger::log_info("END OF DATA TRANSFORMATION.")
#debug(prolfquapp::generate_DEA_reports)
GRP2$pop$remove <- FALSE
grp <- prolfquapp::generate_DEA_reports(lfqdata, GRP2, prot_annot)

saveRDS(grp, file = paste0(GRP2$zipdir, ".rds"))
grp <- readRDS(paste0(GRP2$zipdir, ".rds"))

dir.create( GRP2$zipdir)
prolfquapp::copy_DEA_FragPipe_TMT()
for (i in seq_along(grp)) {
  prolfquapp::write_DEA_all(grp[[i]], names(grp)[i], GRP2$zipdir )
}


#prolfquapp::render_DEA(grp2, outpath = ".", htmlname = "check.html")
