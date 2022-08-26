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

idxzip <- grep("*.zip",yml$application$input[[1]])

GRP2$Bfabric$inputID <- yml$job_configuration$input[[1]][[idxzip]]$resource_id
GRP2$Bfabric$inputURL <- yml$job_configuration$input[[1]][[idxzip]]$resource_url

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

rm(yml)

###
dir.create(ZIPDIR)
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

annot <- annot |> dplyr::mutate(
  raw.file = gsub("^x|.d.zip$|.raw$","",
                  basename(annot$Relative.Path)
  )
)

if( all(c("ContrastName", "Contrast") %in% colnames(annot)) ) {
  contr <- annot |> dplyr::select(ContrastName, Contrast) |> dplyr::filter(nchar(Contrast) > 0)
  Contrasts <- contr$Contrast
  names(Contrasts) <- contr$ContrastName
  GRP2$pop$Contrasts <- Contrasts
}

nr <- sum(annot$raw.file %in% unique(protein$raw.file))
logger::log_info("nr : ", nr, " files annotated")
annot$Relative.Path <- NULL

protein <- dplyr::inner_join(annot, protein)


################### annotations

#### Setup configuration ###

atable <- prolfqua::AnalysisTableAnnotation$new()
atable$fileName = "raw.file"
atable$hierarchy[["protein_Id"]] <- c("protein")


atable$hierarchyDepth <- 1

if(sum(grepl("^name", colnames(annot), ignore.case = TRUE)) > 0){
  atable$sampleName <- grep("^name", colnames(annot) , value=TRUE, ignore.case = TRUE)
}

stopifnot(sum(grepl("^group", colnames(protein), ignore.case = TRUE)) == 1)
groupingVAR <- grep("^group", colnames(protein), value= TRUE, ignore.case = TRUE)
protein[[groupingVAR]]<- gsub("[[:space:]]", "", protein[[groupingVAR]])
protein[[groupingVAR]] <- gsub("[-\\+\\/\\*]", "_", protein[[groupingVAR]])

atable$factorDepth <- 1

atable$factors[["Group_"]] = groupingVAR

if (sum(grepl("^subject", colnames(protein), ignore.case = TRUE)) == 1 & REPEATED) {
  subvar <- grep("^subject", colnames(protein), value = TRUE, ignore.case = TRUE)
  atable$factors[["Subject"]] = subvar

  tmp <- data.frame(table(dplyr::distinct(protein[,c(groupingVAR,subvar)])) )
  if(all(tmp$Freq > 1)){
    atable$factorDepth <- 2
  }
}

if (sum(grepl("^control", colnames(protein), ignore.case = TRUE)) == 1) {
  atable$factors[["CONTROL"]] = grep("^control", colnames(protein), value = TRUE, ignore.case = TRUE)
}

atable$setWorkIntensity("razor.intensity")


# Preprocess Data
config <- prolfqua::AnalysisConfiguration$new(atable)
adata <- prolfqua::setup_analysis(protein, config)
proteinID <- atable$hkeysDepth()

# Create protein annotation.
protein_annot <- "description"
prot_annot <- dplyr::select(protein,
                            dplyr::all_of(c( atable$hierarchy[[proteinID]], protein_annot))) |>
  dplyr::distinct()
prot_annot <- dplyr::rename(prot_annot, description = !!rlang::sym(protein_annot))
prot_annot <- dplyr::rename(prot_annot, !!proteinID := (!!atable$hierarchy[[proteinID]]))


lfqdata <- prolfqua::LFQData$new(adata, config)
lfqdata$remove_small_intensities()

prolfquapp::generate_reports(lfqdata, GRP2, prot_annot, revpattern, contpattern, ZIPDIR)

