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

GRP2$pop$aggregate <- "medpolish"
GRP2$pop$Diffthreshold <- as.numeric(yml$application$parameters$`4|Difference_threshold`)
GRP2$pop$FDRthreshold <- as.numeric(yml$application$parameters$`5|FDR_threshold`)

GRP2$pop$remove <- if(yml$application$parameters$`6|remConDec` == "true"){TRUE} else {FALSE}
GRP2$pop$revpattern <- yml$application$parameters$`7|REVpattern`
GRP2$pop$contpattern <- yml$application$parameters$`8|CONpattern`

GRP2$Software <- "MaxQuant"

rm(yml)

###
dir.create(ZIPDIR)
###


peptidef <- "peptides.txt"
proteinf <- "proteinGroups.txt"
dsf <- "dataset.csv"
REPEATED <- TRUE

stopifnot(file.exists(peptidef), file.exists(proteinf), file.exists(dsf))

protein <- prolfqua::tidyMQ_ProteinGroups(proteinf)
peptide <- prolfqua::tidyMQ_Peptides(peptidef)
annot <- read.csv(dsf)
annot <- data.frame(lapply(annot, as.character))

annot <- annot |> dplyr::mutate(
  raw.file = gsub("^x|.d.zip$|.raw$","",
                  tolower(make.names(basename(annot$Relative.Path)))
  )
)

if( all(c("ContrastName", "Contrast") %in% colnames(annot)) ) {
  contr <- annot |> dplyr::select(ContrastName, Contrast) |> dplyr::filter(nchar(Contrast) > 0)
  Contrasts <- contr$Contrast
  names(Contrasts) <- contr$ContrastName
  GRP2$pop$Contrasts <- Contrasts
}


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
#
atable$hierarchyDepth <- 1

if(sum(grepl("^name", colnames(annot), ignore.case = TRUE)) > 0){
  atable$sampleName <- grep("^name", colnames(annot) , value=TRUE, ignore.case = TRUE)
}


stopifnot(sum(grepl("^group|^bait", colnames(peptide), ignore.case = TRUE)) == 1)

groupingVAR <- grep("^group|^bait", colnames(peptide), value= TRUE, ignore.case = TRUE)
peptide[[groupingVAR]]<- gsub("[[:space:]]", "", peptide[[groupingVAR]])
peptide[[groupingVAR]] <- gsub("[-\\+\\/\\*]", "_", peptide[[groupingVAR]])

atable$factors[["Group_"]] = groupingVAR
atable$factorDepth <- 1
if (sum(grepl("^subject", colnames(peptide), ignore.case = TRUE)) == 1 & REPEATED) {
  atable$factors[["Subject"]] = grep("^subject", colnames(peptide), value = TRUE, ignore.case = TRUE)
  atable$factorDepth <- 2
}

if (sum(grepl("^control", colnames(peptide), ignore.case = TRUE)) == 1) {
  atable$factors[["CONTROL"]] = grep("^control", colnames(peptide), value = TRUE, ignore.case = TRUE)
}


atable$setWorkIntensity("peptide.intensity")


# Preprocess data - aggregate proteins.
config <- prolfqua::AnalysisConfiguration$new(atable)
adata <- prolfqua::setup_analysis(peptide, config)
proteinID <- atable$hkeysDepth()

# CREATE protein annotation.
protein_annot <- "fasta.headers"
prot_annot <- dplyr::select(peptide ,
                            dplyr::all_of(c( atable$hierarchy[[proteinID]], protein_annot))) |>
  dplyr::distinct()
prot_annot <- dplyr::rename(prot_annot, description = !!rlang::sym(protein_annot))
prot_annot <- dplyr::rename(prot_annot, !!proteinID := (!!atable$hierarchy[[proteinID]]))


lfqdata <- prolfqua::LFQData$new(adata, config)
lfqdata$remove_small_intensities()
lfqdata$filter_proteins_by_peptide_count()
GRP2$pop$nrPeptides <- 2

logger::log_info("AGGREGATING PEPTIDE DATA!")

transformed <- lfqdata$get_Transformer()$log2()$lfq
aggregator <- transformed$get_Aggregator()

if (GRP2$pop$aggregate == "medpolish") {
  aggregator$medpolish()
} else if (GRP2$pop$aggregate == "topN") {
  aggregator$sum_topN()
} else if (GRP2$pop$aggregate == "lmrob") {
  aggregator$lmrob()
} else {
  logger::log_warn("no such aggregator {GRP2$pop$aggregate}.")
}

logger::log_info("data aggregated: {GRP2$pop$aggregate}.")
ag <- aggregator$lfq_agg
tr <- ag$get_Transformer()
tr <- tr$intensity_array(exp, force = TRUE)
lfqdata <- tr$lfq
lfqdata$is_transformed(FALSE)

logger::log_info("END OF DATA TRANSFORMATION.")

#debug(prolfquapp::generate_reports)
prolfquapp::generate_reports(lfqdata, GRP2, prot_annot, revpattern, contpattern, ZIPDIR)
