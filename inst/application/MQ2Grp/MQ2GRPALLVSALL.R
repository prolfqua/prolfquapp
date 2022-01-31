library(prolfqua)

ymlfile <- "config.yaml"
yml = yaml::read_yaml(ymlfile)



WORKUNITID = yml$job_configuration$workunit_id
PROJECTID = yml$job_configuration$project_id
ORDERID = yml$job_configuration$order_id
INPUT_ID = yml$job_configuration$input[[1]][[1]]$resource_id
INPUT_URL = yml$job_configuration$input[[1]][[1]]$resource_url

peptidef <- "peptides.txt"
proteinf <- "proteinGroups.txt"
dsf <- "dataset.csv"
REPEATED <- TRUE


stopifnot(file.exists(peptidef), file.exists(proteinf), file.exists(dsf))

protein <- prolfqua::tidyMQ_ProteinGroups(proteinf)
peptide <- prolfqua::tidyMQ_Peptides(peptidef)
annot <- read.csv(dsf)

annot <- annot |> dplyr::mutate(
  raw.file = gsub("^x|.d.zip$|.raw$","",
                  tolower(make.names(basename(annot$Relative.Path)))
  )
)


annot$Relative.Path <- NULL

proteinAnnot <- dplyr::select(protein, proteinID, fasta.headers ) |> dplyr::distinct()
peptide <- dplyr::inner_join(annot, peptide)
peptide <- dplyr::inner_join(proteinAnnot, peptide, by = c(proteinID = "leading.razor.protein"))

################### annotations
GRP2 <- list()
GRP2$Bfabric <- list()
GRP2$Bfabric$projectID <- PROJECTID
GRP2$Bfabric$orderID <- ORDERID
GRP2$Bfabric$workunitID <- WORKUNITID

GRP2$Software <- "MaxQuant"

GRP2$Bfabric$inputID <- INPUT_ID
GRP2$Bfabric$inputURL <- INPUT_URL


GRP2$nrPeptides <- 2
GRP2$log2FCthreshold <- 1
GRP2$FDRthreshold <- 0.1

# Setup configuration

atable <- prolfqua::AnalysisTableAnnotation$new()
atable$fileName = "raw.file"
atable$hierarchy[["protein_Id"]] <- c("proteinID")
atable$hierarchy[["peptide_Id"]] <- c("sequence")

#
atable$hierarchyDepth <- 1
atable$factors[["Experiment_"]] = "Experiment"

if (!is.null(annot$Subject) & REPEATED) {
  atable$factors[["Subject"]] = "Subject"
}
atable$factorDepth <- 1
atable$setWorkIntensity("peptide.intensity")




# Compute all possible 2 Grps to avoid specifying reference.
levels <- annot$Experiment |> unique()
outdir <- "xyz"
dir.create(outdir)


if (TRUE) {
  i <- 1
  j <- 2
  cat(levels[i], levels[j], "\n")
  GRP2$Contrasts <- paste0("Experiment_",levels[i], " - ", "Experiment_",levels[j])
  names(GRP2$Contrasts) <- paste0("Experiment" , levels[i], "_vs_", levels[j])
  message(GRP2$Contrasts)
  fname <- paste0("Experiment_" , levels[i], "_vs_", levels[j])
  outpath <- file.path( outdir, fname)
  proteinF <- peptide |> dplyr::filter(.data$Experiment == levels[i] | .data$Experiment == levels[j])

  grp2 <- prolfqua::make2grpReport(proteinF,
                                   atable,
                                   GRP2,
                                   protein_annot = "fasta.headers",
                                   remove = TRUE)

  names(grp2)

  grp2$contrMore$get_contrasts()
  grp2$transformedlfqData$get_subset(grp2$contrMore)
  prolfqua::write_2GRP(grp2, outpath = outpath, xlsxname = fname)
  prolfqua::render_2GRP(grp2, outpath = outpath, htmlname = fname)


} else {
  for (i in 1:length(levels)) {
    for (j in 1:length(levels)) {
      if (i != j) {
        cat(levels[i], levels[j], "\n")
        GRP2$Contrasts <- paste0("Experiment_",levels[i], " - ", "Experiment_",levels[j])
        names(GRP2$Contrasts) <- paste0("Experiment" , levels[i], "_vs_", levels[j])
        message(GRP2$Contrasts)
        fname <- paste0("Experiment_" , levels[i], "_vs_", levels[j])
        outpath <- file.path( outdir, fname)
        proteinF <- peptide |> dplyr::filter(.data$Experiment == levels[i] | .data$Experiment == levels[j])

        grp2 <- prolfqua::make2grpReport(proteinF,
                                         atable,
                                         GRP2,
                                         protein_annot = "fasta.headers",
                                         remove = TRUE)

        prolfqua::write_2GRP(grp2, outpath = outpath, xlsxname = fname)
        prolfqua::render_2GRP(grp2, outpath = outpath, htmlname = fname)

      }
    }
  }

}
