#' get report.tsv and fasta file location in folder
#' @return list with paths to data and fasta
#' @export
get_FP_PSM_files <- function(path){

  psm_file <- dir(path = ".", pattern = "psm.tsv", recursive = TRUE, full.names = TRUE)
  fasta.files <- grep("*.fasta$", dir(path = path, recursive = TRUE,full.names = TRUE), value = TRUE)
  if (any(grepl("database[0-9]*.fasta$", fasta.files))) {
    fasta.files <- grep("database[0-9]*.fasta$", fasta.files, value = TRUE)
  }
  if (length(fasta.files) == 0) {
    logger::log_error("No fasta file found!")
    stop()
  }
  return(list(data = psm_file, fasta = fasta.files))
}

#' preprocess DIANN ouput, filter by q_value and nr_peptides
#' @return list with lfqdata and protein annotation
#' @export
preprocess_FP_PSM <- function(quant_data,
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

  psm <- prolfqua::tidy_FragPipe_psm(quant_data)
  psm$qValue <- 1 - psm$PeptideProphet.Probability

  nrowPSM <- nrow(psm)
  fasta_annot <- get_annot_from_fasta(fasta.files)
  psm <- dplyr::left_join(psm, fasta_annot, by = c(Protein = "fasta.id"), multiple = "all")
  stopifnot(nrow(psm) == nrowPSM)

  prot_annot <- prolfquapp::dataset_protein_annot(
    psm,
    c("protein_Id" = "Protein"),
    protein_annot = "fasta.header",
    more_columns = "nrPeptides")



  nr <- sum(annot$raw.file %in% sort(unique(peptide$raw.file)))
  logger::log_info("nr : ", nr, " files annotated out of ", length(unique(peptide$raw.file)))

  atable$ident_Score = "PeptideProphet.Probability"
  atable$ident_qValue = "qValue"
  atable$fileName = "channel"
  atable$hierarchy[["protein_Id"]] <- c("Protein")
  atable$hierarchy[["peptide_Id"]] <- c("Peptide")
  atable$hierarchy[["mod_peptide_Id"]] <- c("Modified.Peptide","Assigned.Modifications")
  atable$hierarchy[["Spectrum"]] <- c("Spectrum")

  peptide <- dplyr::inner_join(annot, peptide, multiple = "all")
  config <- prolfqua::AnalysisConfiguration$new(atable)
  adata <- prolfqua::setup_analysis(peptide, config)
  lfqdata <- prolfqua::LFQData$new(adata, config)
  lfqdata$remove_small_intensities()
  protAnnot <- prolfqua::ProteinAnnotation$new(lfqdata , prot_annot)
  return(list(lfqdata = lfqdata , protein_annotation = protAnnot))
}


