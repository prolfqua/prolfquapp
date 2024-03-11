#' get report.tsv and fasta file location in folder
#' @return list with paths to data and fasta
#' @export
get_DIANN_files <- function(path){
  diann.path <- grep("report\\.tsv$|diann-output\\.tsv", dir(path = path, recursive = TRUE, full.names = TRUE), value = TRUE)
  fasta.files <- grep("*.fasta$", dir(path = path, recursive = TRUE,full.names = TRUE), value = TRUE)
  if (any(grepl("database[0-9]*.fasta$", fasta.files))) {
    fasta.files <- grep("database[0-9]*.fasta$", fasta.files, value = TRUE)
  }
  if (length(fasta.files) == 0) {
    logger::log_error("No fasta file found!")
    stop()
  }
  return(list(data = diann.path, fasta = fasta.files))
}

#' preprocess DIANN ouput, filter by q_value and nr_peptides
#' @return list with lfqdata and protein annotation
#' @export
preprocess_DIANN <- function(quant_data,
                             fasta_file,
                             annotation,
                             q_value = 0.01,
                             nrPeptides = 1){

  annot <- annotation$annot
  atable <- annotation$atable
  annot <- annot |> dplyr::mutate(
    raw.file = gsub("^x|.d.zip$|.raw$","",
                    (basename(annot[[atable$fileName]]))
    ))


  peptide <- prolfquapp::read_DIANN_output(
    diann.path = quant_data,
    fasta.file = fasta_file,
    nrPeptides = 1,
    Q.Value = q_value)

  prot_annot <- prolfquapp::dataset_protein_annot(
    peptide,
    c("protein_Id" = "Protein.Group"),
    protein_annot = "fasta.header",
    more_columns = c("nrPeptides", "fasta.id", "Protein.Group.2")
  )

  nr <- sum(annot$raw.file %in% sort(unique(peptide$raw.file)))
  logger::log_info("nr : ", nr, " files annotated out of ", length(unique(peptide$raw.file)))

  atable$fileName = "raw.file"
  atable$hierarchy[["protein_Id"]] <- c("Protein.Group")
  atable$hierarchy[["peptide_Id"]] <- c("Stripped.Sequence")
  atable$set_response("Peptide.Quantity")
  atable$hierarchyDepth <- 1

  peptide <- dplyr::inner_join(annot, peptide, multiple = "all")
  config <- prolfqua::AnalysisConfiguration$new(atable)
  adata <- prolfqua::setup_analysis(peptide, config)
  lfqdata <- prolfqua::LFQData$new(adata, config)
  lfqdata$remove_small_intensities()
  protAnnot <- prolfqua::ProteinAnnotation$new(lfqdata , prot_annot)
  return(list(lfqdata = lfqdata , protein_annotation = protAnnot))
}


