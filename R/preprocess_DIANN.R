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
  atable <- annotation$atable$clone(deep = FALSE)
  annot <- annot |> dplyr::mutate(
    raw.file = gsub("^x|.d.zip$|.raw$","",
                    (basename(annot[[atable$fileName]]))
    ))

  peptide <- prolfquapp::read_DIANN_output(
    diann.path = quant_data,
    fasta.file = fasta_file,
    nrPeptides = nrPeptides,
    Q.Value = q_value)

  nr <- sum(annot$raw.file %in% sort(unique(peptide$raw.file)))
  logger::log_info("nr : ", nr, " files annotated out of ", length(unique(peptide$raw.file)))
  if (nr == 0) { stop("No files are annotated. The annotation file is not compatible withe quant data.") }

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

  protAnnot <- build_protein_annot(
    lfqdata,
    peptide,
    c("protein_Id" = "Protein.Group"),
    cleaned_protein_id = "Protein.Group.2",
    protein_description = "fasta.header",
    nr_children = "nrPeptides",
    more_columns = c("nrPeptides", "fasta.id"))

  return(list(lfqdata = lfqdata , protein_annotation = protAnnot))
}


