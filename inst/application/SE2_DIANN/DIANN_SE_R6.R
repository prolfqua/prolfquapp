# Author : Witold Wolski <wew@fgcz.ethz.ch>
# compatible with prolfqua 3.0.0 release available from https://github.com/wolski/prolfqua/releases/tag/v0.2.9


prolfquapp::copy_SAINT_express(run_script = FALSE)
# Read b-fabric related information
ymlfile <- "config.yaml"

GRP2 <- prolfquapp::make_DEA_config_R6()

GRP2$processing_options$other[["spc"]] <- FALSE

dataset <- dir(".", pattern = 'dataset.csv', full.names = TRUE, recursive = TRUE)[1]
annot <- read.csv(dataset)
annot <- data.frame(lapply(annot, as.character))
annot <- annot |> dplyr::mutate(
  raw.file = gsub("^x|.d.zip$|.raw$","",
                  (basename(annot$Relative.Path))
  )
)

annotation <- annot
colnames(annotation) <- tolower(make.names(colnames(annotation)))

# code to find input files
path = "."
diann.path <- grep("report\\.tsv$|diann-output\\.tsv", dir(path = path, recursive = TRUE), value = TRUE)
fasta.files <- grep("*\\.fasta$|*\\.fas$", dir(path = path, recursive = TRUE), value = TRUE)

if (any(grepl("database[0-9]*.fasta$", fasta.files))) {
  fasta.files <- grep("database[0-9]*.fasta$", fasta.files, value = TRUE)
}

if (length(fasta.files) == 0) {
  logger::log_error("No fasta file found!")
  stop()
}

peptide <- prolfquapp::read_DIANN_output(
  diann.path = diann.path,
  fasta.file = fasta.files,
  nrPeptides = 1,
  Q.Value = 0.1 )

prot_annot <- prolfquapp::dataset_protein_annot(
  peptide,
  c("protein_Id" = "Protein.Group"),
  protein_annot = "fasta.header",
  more_columns = c("nrPeptides", "fasta.id", "Protein.Group.2", "protein_length")
)

nr <- length(union(annot$raw.file , unique(peptide$raw.file)))
logger::log_info("nr : ", nr, " files annotated.")
nrNotData <- length(setdiff(annot$raw.file , unique(peptide$raw.file)))
logger::log_info("nr : ", nrNotData, " files in annotation but not in data. " )
nrNotInAnnot <- length(setdiff( unique(peptide$raw.file), annot$raw.file))
logger::log_info("nr : ", nrNotInAnnot, " files in data but not in annot. " )

annot$Relative.Path <- NULL
peptide <- dplyr::inner_join(annot, peptide, multiple = "all")

# Create config

#### this comes from DIANN
atable <- prolfqua::AnalysisTableAnnotation$new()
atable$fileName = "raw.file"
atable$hierarchy[["protein_Id"]] <- c("Protein.Group")
atable$hierarchy[["peptide_Id"]] <- c("Stripped.Sequence")
atable$set_response("Peptide.Quantity")
atable$hierarchyDepth <- 1


res <- prolfquapp::dataset_set_factors(atable, peptide)
atable <- res$atable
peptide <- res$msdata

# Preprocess data - aggregate proteins.
config <- prolfqua::AnalysisConfiguration$new(atable)
adata <- prolfqua::setup_analysis(peptide, config)
length(unique(adata$protein_Id))
lfqdata <- prolfqua::LFQData$new(adata, config)
lfqdata$remove_small_intensities()
logger::log_info("AGGREGATING PEPTIDE DATA!")

# aggregate data
lfqdataProt <- prolfquapp::aggregate_data(
  lfqdata,
  agg_method = "topN")


## Transformation on
lfqTrans <- prolfquapp::transform_lfqdata(lfqdataProt, method = GRP2$processing_options$transform)
lfqdataProt <- exp2(lfqTrans)
lfqdataProt$response()


lfqdata <- lfqdataProt
protAnnot <- prolfqua::ProteinAnnotation$new(lfqdata , prot_annot)



# CREATING RESULTS
RESULTS <- list() # RESULT is stored in excel table
RESULTS$annotation <- lfqdata$factors()

intdata <- dplyr::inner_join(protAnnot$row_annot, lfqdata$data, multiple = "all")
localSAINTinput <- prolfqua::protein_2localSaint(
  intdata,
  quantcolumn = lfqdata$config$table$get_response(),
  proteinID = "protein_Id",
  proteinLength = "protein_length",
  IP_name = "raw.file",
  baitCol = "Bait_",
  CorTCol = "CONTROL")


RESULTS <- c(RESULTS, localSAINTinput)

resSaint <- prolfqua::runSaint(
  localSAINTinput,
  spc = GRP2$processing_options$other[["spc"]] )

resSaint$list <- dplyr::inner_join(
  protAnnot$row_annot,
  resSaint$list,
  by = c(protein_Id = "Prey"),
  keep = TRUE , multiple = "all")

RESULTS <- c(RESULTS, resSaint)
# write analysis results
# Prepare result visualization and render report
cse <- prolfqua::ContrastsSAINTexpress$new(resSaint$list)


#### THIS NEEDS TO BE CLEANED UP #####


resContrasts <- cse$get_contrasts()



sig <- resContrasts |>
  dplyr::filter(.data$BFDR  <  REPORTDATA$FDRthreshold &
                  .data$log2_EFCs  >  log2(REPORTDATA$FCthreshold))

# Transform data for PCA visualization etc
tt <- lfqdata$get_Transformer()$log2()
lfqdata_transformed <- tt$lfq

REPORTDATA$pups <- prolfqua::UpSet_interaction_missing_stats(
  lfqdataProt$data,
  lfqdata$config,tr = 2)
RESULTS$InputData <- lfqdata$to_wide()$data

gs <- lfqdata$get_Summariser()
RESULTS$MissingInformation <- gs$interaction_missing_stats()$data
RESULTS$MissingInformation$isotopeLabel <- NULL
RESULTS$listFile <- NULL
writexl::write_xlsx(RESULTS, path = file.path(ZIPDIR,paste0(treat, "_data.xlsx")))

REPORTDATA$BFABRIC <- BFABRIC
REPORTDATA$lfqdata_transformed <- lfqdata_transformed
REPORTDATA$sig <- sig
REPORTDATA$resContrasts <- resContrasts
REPORTDATA$prot_annot <- dplyr::rename(prot_annot, protein = protein_Id)


tmp <- prolfqua::get_UniprotID_from_fasta_header(REPORTDATA$pups$data, "protein_Id")
if(is.null(tmp$UniprotID)) {tmp$UniprotID <- sapply(tmp$protein_Id, function(x){strsplit(x,";")[[1]][1]})}


write.table(data.frame(tmp$UniprotID), file = file.path(ZIPDIR,"ORA_background.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE )
sig |> dplyr::group_by(Bait) |> tidyr::nest() -> sigg

if (nrow(sigg) > 0) {
  for (i in 1:nrow(sigg)) {
    tmp <- prolfqua::get_UniprotID_from_fasta_header(sigg$data[[i]], "Prey")
    if(is.null(tmp$UniprotID)) {tmp$UniprotID <- sapply(tmp$Prey, function(x){strsplit(x,";")[[1]][1]})}
    filename <- paste0("ORA_Bait_", sigg$Bait[i] , ".txt")
    write.table(data.frame(tmp$UniprotID),
                file = file.path(ZIPDIR, filename),
                col.names = FALSE,
                row.names = FALSE,
                quote = FALSE )
  }
}

prolfqua::copy_SAINTe_doc(workdir = ZIPDIR)

SEP <- REPORTDATA

saveRDS(REPORTDATA, file = "REPORTDATA.rds")
rm(list = setdiff(ls(), c("REPORTDATA","ZIPDIR","treat","yml"))) # delete all variables not needed for rendering
SEP <- REPORTDATA

fragPipeDIA <- names(yml$application$input) == "FragPipeDIA-dataset"
text <-  c( "The LC-MS data was processed using the ",
            if(fragPipeDIA){"FragPipe-DIA proteomics pipeline [@yu2023analysis]."} else {"DIA-NN software [@demichev2020dia]."} ,
            "The quantification results were extracted from the DIANN main report, containing precursor (charged modified peptide) quantification results.",
            "The protein abundances were estimated by the sum of all the precursor abundances (Precursor.Quantity column) assigned to a protein. ")

text <- paste(text, collapse = "")

rmarkdown::render("SaintExpressReportMsFragger.Rmd",
                  params = list(sep = REPORTDATA, textpreprocessing = text),
                  output_format = bookdown::html_document2())


file.copy("SaintExpressReportMsFragger.html",
          file.path(ZIPDIR, paste0(treat, "SaintExpressReportMsFragger.html")),
          overwrite = TRUE)



