ymlfile <- "config.yaml"
yml = yaml::read_yaml(ymlfile)

WORKUNITID = yml$job_configuration$workunit_id
PROJECTID = yml$job_configuration$project_id
ORDERID = yml$job_configuration$order_id

idxzip <- grep("*.zip",yml$application$input$MaxQuant_textfiles)

INPUT_ID = yml$job_configuration$input[[1]][[idxzip]]$resource_id
INPUT_URL = yml$job_configuration$input[[1]][[idxzip]]$resource_url
###

yml$application$parameters$`3|Normalization`

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

nr <- sum(annot$raw.file %in% unique(peptide$raw.file))
logger::log_info("nr", nr, " files annotated")
annot$Relative.Path <- NULL


proteinAnnot <- dplyr::select(protein, proteinID, fasta.headers ) |> dplyr::distinct()
peptide <- dplyr::inner_join(annot, peptide)
peptide <- dplyr::inner_join(proteinAnnot, peptide, by = c(proteinID = "leading.razor.protein"))


################### annotations
GRP2 <- list()
GRP2$Bfabric <- list()
GRP2$Bfabric$projectID <- PROJECTID


GRP2$Bfabric$projectName <- "" # workunit name in the future.
GRP2$Bfabric$orderID <- ORDERID

GRP2$Bfabric$workunitID <- WORKUNITID
GRP2$Bfabric$inputID <- INPUT_ID
GRP2$Bfabric$inputURL <- INPUT_URL

#at least 2 peptides per protein
GRP2$pop <- list()
GRP2$pop$transform <- yml$application$parameters$`3|Normalization`

GRP2$pop$aggregate <- yml$application$parameters$`2|Aggregation`
GRP2$pop$Diffthreshold <- as.numeric(yml$application$parameters$`4|Difference_threshold`)
GRP2$pop$FDRthreshold <- as.numeric(yml$application$parameters$`5|FDR_threshold`)
removeREV <- if(yml$application$parameters$`6|remConDec` == "true"){TRUE} else {FALSE}
revpattern <- yml$application$parameters$`7|REVpattern`
contpattern <- yml$application$parameters$`8|CONpattern`

GRP2$Software <- "MaxQuant"


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


for (i in 1:length(levels)) {
  for (j in 1:length(levels)) {
    if (i != j) {
      cat(levels[i], levels[j], "\n")
      GRP2$pop$Contrasts <- paste0("Experiment_",levels[i], " - ", "Experiment_",levels[j])
      names(GRP2$pop$Contrasts) <- paste0("Experiment" , levels[i], "_vs_", levels[j])
      message("CONTRAST", GRP2$pop$Contrasts)
      fname <- paste0("Experiment_" , levels[i], "_vs_", levels[j])
      outpath <- file.path( outdir, fname)
      proteinF <- peptide |> dplyr::filter(.data$Experiment == levels[i] | .data$Experiment == levels[j])
      grp2 <- prolfqua::make2grpReport(proteinF,
                                       atable,
                                       GRP2,
                                       protein_annot = "fasta.headers",
                                       revpattern = revpattern,
                                       contpattern = contpattern,
                                       remove = TRUE)

      prolfqua::write_2GRP(grp2, outpath = outpath, xlsxname = fname)
      prolfqua::render_2GRP(grp2, outpath = outpath, htmlname = fname)

    }
  }
}
