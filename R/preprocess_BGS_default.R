#' get BGS and fasta file location in folder
#' @return list with paths to data and fasta
#' @export
read_BGS <- function(file = "Experiment1_Report_BGS Factory Report (Normal).tsv"){
  bgs <- readr::read_tsv(file)
  colnames(bgs)
  colnames(bgs) <- colnames(bgs) |>
    stringr::str_replace_all( "[[:space:]\\(\\)\\-]", "_") |>
    stringr::str_replace_all( "_+", "_") |>
    stringr::str_replace("_$", "")

  ctoselect <- c(
    "R.FileName",
    "PG.ProteinGroups",
    "PG.ProteinAccessions",
    "PG.Qvalue",
    "PG.QValue_Run_Wise",
    "PG.Quantity",
    "PEP.GroupingKey",
    "PEP.IsProteotypic",
    "PEP.RunEvidenceCount",
    "PEP.NrOfMissedCleavages",
    "EG.ModifiedSequence",
    "EG.Qvalue",
    #"EG.TotalQuantity_Settings",
    "FG.Qvalue",
    "FG.Charge",
    #"FG.LabeledSequence",
    "FG.Quantity")
  bgsf <- bgs[,ctoselect]
  return(bgsf)
}


#' get BGS and fasta file location in folder
#' @return list with paths to data and fasta
#' @export
#' @examples
#' \dontrun{
#' x <- get_DIANN_files("inst/application/DIANN/2517219/")
#' }
get_BGS_files <- function(path, bgs_pattern = "*BGS Factory Report \\(Normal\\).tsv|_Report.tsv"){
  diann.path <- grep(bgs_pattern, dir(path = path, recursive = TRUE, full.names = TRUE), value = TRUE)
  fasta.files <- grep("*.fasta$|*.fas$", dir(path = path, recursive = TRUE, full.names = TRUE), value = TRUE)
  if (any(grepl("database[0-9]*.fasta$", fasta.files))) {
    fasta.files <- grep("database[0-9]*.fasta$", fasta.files, value = TRUE)
  }
  fasta.files <- fasta.files[!grepl("first-pass",fasta.files)]

  if (length(fasta.files) == 0) {
    logger::log_error("No fasta file found!")
    stop()
  }
  return(list(data = diann.path, fasta = fasta.files))
}

#' create templte dataset for BGS data
#' @return data.frame
#' @export
dataset_template_BGS <- function(files){
  xt <- readr::read_tsv(files$data)
  ds <- xt |> dplyr::select(raw.file = "R.FileName",
                            Group = "R.Condition",
                            name = "R.Replicate")
  ds <- ds |>
    tidyr::unite("Name", c("Group", "name"), remove = FALSE) |>
    dplyr::distinct()
  ds$name <- NULL
  return(ds)
}

#' preprocess DIANN ouput, filter by q_value and nr_peptides
#' @return list with lfqdata and protein annotation
#' @export
#' @examples
#' \dontrun{
#' x <- get_BGS_files("DefaultParsing")
#' bgs <- read_BGS(x$data)
#' annot <- data.frame(raw.file = bgs$R.FileName |> unique(),
#'  Name = paste(c(rep("A",3),rep("B",3)),1:6, sep="_"),
#' group = c(rep("A",3),rep("B",3)))
#' annotation <- annot |> prolfquapp::read_annotation(QC = TRUE)
#' #debug(preprocess_BGS)
#' xd <- preprocess_BGS(x$data, x$fasta, annotation)
#' }
preprocess_BGS <- function(quant_data,
                             fasta_file,
                             annotation,
                             pattern_contaminants = "^zz|^CON|Cont_",
                             pattern_decoys = "^REV_|^rev",
                             q_value = 0.01,
                             hierarchy_depth = 2){

  annot <- annotation$annot
  atable <- annotation$atable$clone(deep = FALSE)
  annot <- annot |> dplyr::mutate(
    raw.file = gsub("^x|.d.zip$|.raw$","",
                    (basename(annot[[atable$fileName]]))
    ))
  report2 <- read_BGS(quant_data)

  nrPEP <- report2 |> dplyr::select("PG.ProteinGroups", "PEP.GroupingKey") |> dplyr::distinct() |>
    dplyr::group_by(dplyr::across("PG.ProteinGroups")) |> dplyr::summarize(nrPeptides = n())

  nrPEP$Protein.Group.2 <- sapply(nrPEP$PG.ProteinGroups, function(x){ unlist(strsplit(x, "[ ;]"))[1]} )

  nr <- sum(annot$raw.file %in% sort(unique(report2$R.FileName)))
  logger::log_info("nr : ", nr, " files annotated out of ", length(unique(report2$R.FileName)))
  if (nr == 0) { stop("No files are annotated. The annotation file is not compatible withe quant data.") }

  atable$fileName = "raw.file"
  #atable$nr_children = "nr_children"
  atable$ident_qValue = "FG.Qvalue"
  atable$hierarchy[["protein_Id"]] <- c("PG.ProteinGroups")
  atable$hierarchy[["peptide_Id"]] <- c("PEP.GroupingKey")
  atable$hierarchy[["elution_group"]] <- c("EG.ModifiedSequence", "FG.Charge")
  atable$set_response("FG.Quantity")
  atable$hierarchyDepth <- hierarchy_depth

  report2 <- dplyr::inner_join(annot, report2, multiple = "all", by=c("raw.file" = "R.FileName"))
  config <- prolfqua::AnalysisConfiguration$new(atable)
  adata <- prolfqua::setup_analysis(report2, config)
  lfqdata <- prolfqua::LFQData$new(adata, config)
  lfqdata$remove_small_intensities()

  # build protein annotation
  logger::log_info("start reading fasta.")
  fasta_annot <- get_annot_from_fasta(fasta_file, pattern_decoys = pattern_decoys, isUniprot = TRUE)
  logger::log_info("reading fasta done, creating protein annotation.")

  prot_annot <- dplyr::left_join(nrPEP, fasta_annot, by = c(Protein.Group.2 = "proteinname"))
  prot_annot <- dplyr::rename(prot_annot, IDcolumn = "Protein.Group.2",
                              description = "fasta.header",
                              protein_Id = "PG.ProteinGroups" )

  protAnnot <- prolfquapp::ProteinAnnotation$new(
    lfqdata , prot_annot, description = "description",
    cleaned_ids = "IDcolumn",
    full_id = "fasta.id",
    exp_nr_children = "nrPeptides",
    pattern_contaminants = pattern_contaminants,
    pattern_decoys = pattern_decoys
  )
  logger::log_info("protein annotation done.")
  return(list(lfqdata = lfqdata , protein_annotation = protAnnot))
}

