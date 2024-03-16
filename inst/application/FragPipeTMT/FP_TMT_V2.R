# This script starts from psm.tsv
# does basic filtering using purity.
#
library(tidyverse)
library(prolfquapp)
# reading foreign data
REPEATED <- TRUE

# sanitize peptide csv.
psm_file <- dir(path = ".", pattern = "psm.tsv", recursive = TRUE, full.names = TRUE)
fasta.files <- grep("*.fasta$", dir(path = ".", recursive = TRUE), value = TRUE)


psm <- prolfqua::tidy_FragPipe_psm(psm_file)
nrowPSM <- nrow(psm)
fasta_annot <- get_annot_from_fasta(fasta.files)
psm <- left_join(psm, fasta_annot, by = c(Protein = "fasta.id"), multiple = "all")

stopifnot(nrow(psm) == nrowPSM)

prot_annot <- prolfquapp::dataset_protein_annot(
  psm,
  c("protein_Id" = "Protein"),
  protein_annot = "fasta.header",
  more_columns = "nrPeptides")

prot_annot |> filter(grepl("^zz",protein_Id)) |> head()

psm$qValue <- 1 - psm$PeptideProphet.Probability
ds_file <- "dataset.csv"
annot <- read.csv(ds_file)


### READ YAML
GRP2 <- prolfquapp::read_yaml("config.yaml")
GRP2 <- prolfquapp::dataset_extract_contrasts(annot, GRP2)

channelCol  <- grep("^channel", names(annot), ignore.case = TRUE, value = TRUE)

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
# CREATE protein annotation.



atable$set_response("abundance")
# Preprocess data - aggregate proteins.

config <- prolfqua::AnalysisConfiguration$new(atable)
adata <- prolfqua::setup_analysis(psm, config)
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
lfqdata$config$table$hkeysDepth()
lfqdata$config$table$hierarchyDepth <- 1
lfqdata <- prolfquapp::aggregate_data(lfqdata, agg_method = GRP2$pop$aggregate)
#logger::log_info("data aggregated: {GRP2$pop$aggregate}.")
#lfqdata$factors()

logger::log_info("END OF DATA TRANSFORMATION.")
#debug(prolfquapp::generate_DEA_reports)
protAnnot <- prolfqua::ProteinAnnotation$new(lfqdata , prot_annot)
grp <- prolfquapp::generate_DEA_reports(lfqdata, GRP2, protAnnot)

dir.create( GRP2$zipdir)
for (i in seq_along(grp)) {
  prolfquapp::write_DEA_all(grp[[i]], names(grp)[i], GRP2$zipdir, boxplot = FALSE )
}

for (i in seq_along(grp)) {
  SE <- prolfquapp::make_SummarizedExperiment(grp[[i]])
  saveRDS(SE, file = file.path(GRP2$zipdir, paste0("DE_", names(grp)[i]) , paste0("SummarizedExperiment",".rds") ))
}


