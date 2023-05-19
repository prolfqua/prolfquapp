ymlfile <- "config.yaml"
GRP2 <- prolfquapp::read_yaml(ymlfile)
dir.create(GRP2$zipdir)


peptidef <- "peptides.txt"
proteinf <- "proteinGroups.txt"
dsf <- "dataset.csv"
REPEATED <- TRUE
stopifnot(file.exists(peptidef), file.exists(proteinf), file.exists(dsf))

protein <- prolfqua::tidyMQ_ProteinGroups(proteinf)

prot_annot <- prolfquapp::dataset_protein_annot(
  protein,
  c(protein_Id = "proteinID"),
  protein_annot = "fasta.headers",
  more_columns = "nr.peptides")



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

peptide <- dplyr::inner_join(annot, peptide,multiple = "all")


# Setup configuration
atable <- prolfqua::AnalysisTableAnnotation$new()
atable$fileName = "raw.file"
atable$hierarchy[["protein_Id"]] <- c("leading.razor.protein")
atable$hierarchy[["peptide_Id"]] <- c("sequence")
atable$set_response("peptide.intensity")
atable$hierarchyDepth <- 1

res <- prolfquapp::dataset_set_factors(atable, peptide)
atable <- res$atable
peptide <- res$msdata
head(peptide)


# Preprocess data - aggregate proteins.
config <- prolfqua::AnalysisConfiguration$new(atable)

adata <- prolfqua::setup_analysis(peptide, config)

lfqdata <- prolfqua::LFQData$new(adata, config)
lfqdata$remove_small_intensities()
lfqdata$filter_proteins_by_peptide_count()

GRP2$pop$nrPeptides <- 2
logger::log_info("AGGREGATING PEPTIDE DATA!")
lfqdata <- prolfquapp::aggregate_data(lfqdata, agg_method = GRP2$pop$aggregate)

### XXX
protAnnot <- prolfqua::ProteinAnnotation$new(lfqdata , prot_annot)

grp <- prolfquapp::generate_DEA_reports(lfqdata, GRP2, protAnnot)

for (i in seq_along(grp)) {
  prolfquapp::write_DEA_all(grp[[i]], names(grp)[i], GRP2$zipdir, boxplot = TRUE )
}

