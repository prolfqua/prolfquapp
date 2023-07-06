library(rlang)
library(prolfqua)
library(prolfquapp)
library(R6)

ymlfile <- "WU289521/config.yaml"

GRP2 <- prolfquapp::read_yaml(ymlfile)


dsf <- "dataset.csv"
diann.path <- grep("diann-output.tsv", dir(recursive = TRUE), value = TRUE)
fasta.file <- grep("*.fasta", dir(recursive = TRUE), value = TRUE)



###
dir.create(GRP2$zipdir)

### read the fasta and DIANN ouptut
peptide <- read_DIANN_output(
  diann.path = diann.path[1],
  fasta.file = fasta.file[1],
  nrPeptides = 1,
  Q.Value = 0.1)

prot_annot <- prolfquapp::dataset_protein_annot(
  peptide,
  c("protein_Id" = "Protein.Group"),
  protein_annot = "fasta.header",
  more_columns = c("nrPeptides", "fasta.id", "Protein.Group.2")
)

atable <- prolfqua::AnalysisTableAnnotation$new()
atable$fileName = "raw.file"
atable$hierarchy[["protein_Id"]] <- c("Protein.Group")
atable$hierarchy[["peptide_Id"]] <- c("Stripped.Sequence")
atable$set_response("Peptide.Quantity")
atable$hierarchyDepth <- 1


### read annotation file

annot <- read.csv(dsf)
annot <- data.frame(lapply(annot, as.character))
annot <- annot |> dplyr::mutate(
  raw.file = gsub("^x|.d.zip$|.raw$","",
                  (basename(annot$Relative.Path))
  )
)
annot$Relative.Path <- NULL
GRP2 <- prolfquapp::dataset_extract_contrasts(annot, GRP2)

# Join peptide data and annotation

annot$raw.file[ !annot$raw.file %in% sort(unique(peptide$raw.file)) ]
nr <- sum(annot$raw.file %in% sort(unique(peptide$raw.file)))
logger::log_info("nr : ", nr, " files annotated")
peptide <- dplyr::inner_join(annot, peptide, multiple = "all")


# also massages some columns in the petpide table
res <- prolfquapp::dataset_set_factors(atable, peptide)
atable <- res$atable
peptide <- res$msdata

# Preprocess data - aggregate proteins.
config <- prolfqua::AnalysisConfiguration$new(atable)
adata <- prolfqua::setup_analysis(peptide, config)

lfqdata <- prolfqua::LFQData$new(adata, config)
lfqdata$remove_small_intensities()
GRP2$pop$nrPeptides <- 2
logger::log_info("AGGREGATING PEPTIDE DATA!")

lfqdata <- prolfquapp::aggregate_data(lfqdata, agg_method = GRP2$pop$aggregate)
logger::log_info("data aggregated: {GRP2$pop$aggregate}.")
logger::log_info("END OF DATA TRANSFORMATION.")

protAnnot <- prolfqua::ProteinAnnotation$new(lfqdata , prot_annot)
prolfquapp::copy_DEA_DIANN()

grp <- prolfquapp::generate_DEA_reports(lfqdata, GRP2, protAnnot)
saveRDS(grp, file = "DEAll.RDS")

for (i in seq_along(grp)) {
  prolfquapp::write_DEA_all(grp[[i]], names(grp)[i], GRP2$zipdir , boxplot = FALSE)
}



