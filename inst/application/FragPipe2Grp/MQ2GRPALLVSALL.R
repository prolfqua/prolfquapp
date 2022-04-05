ymlfile <- "config.yaml"
yml = yaml::read_yaml(ymlfile)

WORKUNITID = yml$job_configuration$workunit_id
PROJECTID = yml$job_configuration$project_id
ORDERID = yml$job_configuration$order_id
ORDERID <- if(is.null(ORDERID)){PROJECTID}else{ORDERID}

ZIPDIR = paste0("C",ORDERID,"WU",WORKUNITID)

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

GRP2$Software <- "FragPipe"

###



dir.create(ZIPDIR)



idxzip <- grep("*.zip",yml$application$input[[1]])

INPUT_ID = yml$job_configuration$input[[1]][[idxzip]]$resource_id
INPUT_URL = yml$job_configuration$input[[1]][[idxzip]]$resource_url
###

NORMALIZATION <- yml$application$parameters$`3|Normalization`

proteinf <- "combined_protein.tsv"
dsf <- "dataset.csv"
REPEATED <- TRUE

stopifnot( file.exists(proteinf), file.exists(dsf))

protein <- prolfqua::tidy_MSFragger_combined_protein_V16("combined_protein.tsv")
# remove single hit wonders.
protein <- protein |> dplyr::filter(combined.total.peptides > 1 )

annot <- read.csv(dsf)
annot <- annot |> dplyr::mutate(
  raw.file = gsub("^x|.d.zip$|.raw$","",
                  basename(annot$Relative.Path)
  )
)

nr <- sum(annot$raw.file %in% unique(protein$raw.file))
logger::log_info("nr", nr, " files annotated")
annot$Relative.Path <- NULL

proteinAnnot <- dplyr::select(protein, protein, description ) |> dplyr::distinct()
protein <- dplyr::inner_join(annot, protein)


################### annotations

#### Setup configuration ###

atable <- prolfqua::AnalysisTableAnnotation$new()
atable$fileName = "raw.file"
atable$hierarchy[["protein_Id"]] <- c("protein")

if(sum(grepl("^name", colnames(annot), ignore.case = TRUE)) > 0){
  atable$sampleName <- grep("^name", colnames(annot) , value=TRUE, ignore.case = TRUE)
}

atable$hierarchyDepth <- 1

stopifnot(sum(grepl("^group", colnames(protein), ignore.case = TRUE)) == 1)
groupingVAR <- grep("^group", colnames(protein), value= TRUE, ignore.case = TRUE)
atable$factors[["Group_"]] = groupingVAR

if (!is.null(annot$Subject) & REPEATED) {
  atable$factors[["Subject"]] = grep("^subject", colnames(protein), value = TRUE, ignore.case = TRUE)
}

atable$factorDepth <- 1
atable$setWorkIntensity("razor.intensity")


# Compute all possible 2 Grps to avoid specifying reference.
levels <- protein[[groupingVAR]] |> unique()

logger::log_info("levels : ", paste(levels, collapse = " "))

if(! length(levels) > 1){
  logger::log_error("not enough group levels_ to make comparisons.")
}


GRP2$pop$nrPeptides <- 2


for (i in seq_along(levels)) {
  for (j in seq_along(levels)) {
    if (i != j) {
      cat(levels[i], levels[j], "\n")
      GRP2$pop$Contrasts <- paste0("Group_",levels[i], " - ", "Group_",levels[j])
      names(GRP2$pop$Contrasts) <- paste0("Group_" , levels[i], "_vs_", levels[j])
      message("CONTRAST : ", GRP2$pop$Contrasts)
      proteinF <- protein |> dplyr::filter(!!rlang::sym(groupingVAR) == levels[i] | !!rlang::sym(groupingVAR) == levels[j])
      grp2 <- prolfquapp::make2grpReport(proteinF,
                                         atable,
                                         GRP2,
                                         protein_annot = "description",
                                         revpattern = revpattern,
                                         contpattern = contpattern,
                                         remove = TRUE)

      fname <- paste0("Group_" , levels[i], "_vs_", levels[j])
      qcname <- paste0("QC_" , levels[i], "_vs_", levels[j])
      outpath <- file.path( ZIPDIR, fname)
      prolfquapp::write_2GRP(grp2, outpath = outpath, xlsxname = fname)
      prolfquapp::render_2GRP(grp2, outpath = outpath, htmlname = fname, stage=FALSE)
      prolfquapp::render_2GRP(grp2, outpath = outpath, htmlname = qcname, stage=FALSE,markdown = "_DiffExpQC.Rmd")

    }

  }
}
