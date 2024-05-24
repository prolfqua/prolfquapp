message("prolfqua Version :", packageVersion("prolfqua"), "\n")
message("prolfqua Version :", packageVersion("prolfquapp"), "\n")

ymlfile <- "config.yaml"
GRP2 <- prolfquapp::read_yaml(ymlfile, application = "FragPipe")

###
dir.create(GRP2$zipdir)
###

proteinf <- "combined_protein.tsv"
dsf <- "dataset.csv"
REPEATED <- TRUE
stopifnot( file.exists(proteinf), file.exists(dsf))
protein <- prolfqua::tidy_FragPipe_combined_protein("combined_protein.tsv")
# remove single hit wonders.
protein <- protein |> dplyr::filter(combined.total.peptides > 1 )
GRP2$pop$nrPeptides <- 2

annot <- read.csv(dsf)
annot <- data.frame(lapply(annot, as.character))



annot <- prolfquapp::sanitize_grouping_var(annot)
annot <- annot |> dplyr::mutate(
  raw.file = gsub(
    "^x|.d.zip$|.raw$",
    "",
    basename(annot$Relative.Path)
  )
)

GRP2 <- prolfquapp::dataset_set_factors_deprecated(annot = annot, GRP2 = GRP2)

nr <- sum(annot$raw.file %in% unique(protein$raw.file))
logger::log_info("nr : ", nr, " files annotated")
annot$Relative.Path <- NULL

protein <- dplyr::inner_join(annot, protein, multiple = "all")


################### annotations

#### Setup configuration ###

atable <- prolfqua::AnalysisTableAnnotation$new()
atable$fileName = "raw.file"
atable$hierarchy[["protein_Id"]] <- c("protein")
atable$hierarchyDepth <- 1
atable$set_response("razor.intensity")

res <- prolfquapp::dataset_set_factors(atable, protein)
atable <- res$atable
protein <- res$msdata


# Preprocess Data
config <- prolfqua::AnalysisConfiguration$new(atable)
adata <- prolfqua::setup_analysis(protein, config)
proteinID <- atable$hkeysDepth()

# create protein annotation
prot_annot <- prolfquapp::dataset_protein_annot(protein,
                                                idcol = c("protein_Id" = "protein"),
                                                protein_annot = "description",
                                                more_columns = c("combined.total.peptides","gene","coverage"))

lfqdata <- prolfqua::LFQData$new(adata, config)
lfqdata$remove_small_intensities()

prot_annot <- prolfquapp::ProteinAnnotation$new(lfqdata, prot_annot)

grp <- prolfquapp::generate_DEA_reports(lfqdata, GRP2, prot_annot)
for (i in seq_along(grp)) {
  prolfquapp::write_DEA_all(grp[[i]], names(grp)[i], GRP2$zipdir )
}

