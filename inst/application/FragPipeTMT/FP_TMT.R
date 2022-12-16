library(tidyverse)

ymlfile <- "config.yaml"
GRP2 <- prolfquapp::read_yaml(ymlfile, application = "FragPipeTMT")
###
dir.create(GRP2$zipdir)
###
# reading foreign data
REPEATED <- TRUE
# sanitize peptide csv.
peptides <- read.csv(file = "peptide.tsv", sep = "\t")

peptides <- peptides |> dplyr::select(
  dplyr::all_of(c("Peptide", "Protein", "Protein.ID",   "Entry.Name",   "Gene", "Protein.Description")),
  dplyr::starts_with("channel_"))
peptides <- peptides |> tidyr::pivot_longer(dplyr::starts_with("channel"), names_to = "channel", values_to = "abundance")
peptides <- peptides |> dplyr::mutate(channel = gsub("channel_", "", .data$channel))


annot <- read.csv("dataset.csv")
GRP2 <- prolfquapp::dataset_extract_contrasts(annot, GRP2)
channelCol  <- grep("channel", names(annot), ignore.case = TRUE, value = TRUE)
if (channelCol != "channel") {
  annot[["channel"]] <- annot[[channelCol]]
  annot[[channelCol]] <- NULL
}

nr <- sum(annot$channel %in% unique(peptides$channel))
logger::log_info("nr : ", nr, " files annotated")
peptide <- dplyr::inner_join(annot, peptides)

# Setup configuration

atable <- prolfqua::AnalysisTableAnnotation$new()
atable$fileName = "channel"
atable$hierarchy[["protein_Id"]] <- c("Protein")
atable$hierarchy[["peptide_Id"]] <- c("Peptide")

#
tmp <- prolfquapp::dataset_set_factors(atable, peptide)
atable <- tmp$atable
peptide <- tmp$msdata

# CREATE protein annotation.
prot_annot <- prolfquapp::dataset_protein_annot(
  peptide,
  "Protein",
  protein_annot = "Protein.Description")

atable$set_response("abundance")
# Preprocess data - aggregate proteins.
config <- prolfqua::AnalysisConfiguration$new(atable)
adata <- prolfqua::setup_analysis(peptide, config)

lfqdata <- prolfqua::LFQData$new(adata, config)
lfqdata$remove_small_intensities()
lfqdata$filter_proteins_by_peptide_count()
GRP2$pop$nrPeptides <- 2


logger::log_info("AGGREGATING PEPTIDE DATA!")
lfqdata <- prolfquapp::aggregate_data(lfqdata, agg_method = GRP2$pop$aggregate)
logger::log_info("data aggregated: {GRP2$pop$aggregate}.")
lfqdata$factors()

logger::log_info("END OF DATA TRANSFORMATION.")
#debug(prolfquapp::generate_DEA_reports)
grp <- prolfquapp::generate_DEA_reports(lfqdata, GRP2, prot_annot)

for (i in seq_along(grp)) {
  prolfquapp::write_DEA_all(grp[[i]], names(grp)[i], GRP2$zipdir )
}
