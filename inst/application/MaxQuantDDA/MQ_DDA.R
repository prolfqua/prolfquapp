ymlfile <- "config.yaml"
GRP2 <- prolfquapp::read_yaml(ymlfile)
dir.create(GRP2$zipdir)

# internal <- data.frame(protein_Id = "sp|P60709|ACTB_HUMAN",
#                        "sp|P08670|VIME_HUMAN",
#                        "sp|P0C0S8|H2A1_HUMAN",
#                        "sp|P68431|H31_HUMAN",
#                        "sp|P62805|H4_HUMAN",
#                        "sp|P62807|H2B1C_HUMAN")
# GRP2$pop$internal = internal

peptidef <- "peptides.txt"
proteinf <- "proteinGroups.txt"
dsf <- "dataset.csv"
REPEATED <- TRUE
stopifnot(file.exists(peptidef), file.exists(proteinf), file.exists(dsf))
protein <- prolfqua::tidyMQ_ProteinGroups(proteinf)
fh <- unique(protein$fasta.headers)


peptide <- prolfqua::tidyMQ_Peptides(peptidef)
annot <- read.csv(dsf)
annot <- data.frame(lapply(annot, as.character))
annot <- annot |> dplyr::mutate(
  raw.file = gsub("^x|.d.zip$|.raw$","",
                  tolower(make.names(basename(annot$Relative.Path)))
  )
)

GRP2 <- prolfquapp::dataset_extract_contrasts(annot = annot, GRP2 = GRP2)

nr <- sum(annot$raw.file %in% unique(peptide$raw.file))
logger::log_info("nr : ", nr, " files annotated")
annot$Relative.Path <- NULL

proteinAnnot <- dplyr::select(protein, proteinID, fasta.headers ) |> dplyr::distinct()
peptide <- dplyr::inner_join(annot, peptide)
peptide <- dplyr::inner_join(proteinAnnot,
                             peptide,
                             by = c(proteinID = "leading.razor.protein"))
# Setup configuration
atable <- prolfqua::AnalysisTableAnnotation$new()
atable$fileName = "raw.file"
atable$hierarchy[["protein_Id"]] <- c("proteinID")
atable$hierarchy[["peptide_Id"]] <- c("sequence")
atable$set_response("peptide.intensity")
atable$hierarchyDepth <- 1

res <- prolfquapp::dataset_set_factors(atable, peptide)
atable <- res$atable
peptide <- res$msdata

# create protein annotation
prot_annot <- prolfquapp::dataset_protein_annot(peptide, "proteinID", protein_annot = "fasta.headers")

# Preprocess data - aggregate proteins.
config <- prolfqua::AnalysisConfiguration$new(atable)
adata <- prolfqua::setup_analysis(peptide, config)

lfqdata <- prolfqua::LFQData$new(adata, config)
lfqdata$remove_small_intensities()
lfqdata$filter_proteins_by_peptide_count()

GRP2$pop$nrPeptides <- 2
logger::log_info("AGGREGATING PEPTIDE DATA!")
lfqdata <- prolfquapp::aggregate_data(lfqdata, agg_method = GRP2$pop$aggregate)

grp <- prolfquapp::generate_DEA_reports(lfqdata, GRP2, prot_annot)
for (i in seq_along(grp)) {
  prolfquapp::write_DEA_all(grp[[i]], names(grp)[i], GRP2$zipdir )
}

