#' get report.tsv and fasta file location in folder
#' @return list with paths to data and fasta
#' @export
get_FP_PSM_files <- function(path){

  psm_file <- dir(path = path, pattern = "psm.tsv", recursive = TRUE, full.names = TRUE)
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
                             purity_threshold = 0.5,
                             PeptideProphetProb = 0.9,
                             column_before_quants = "Quan Usage",
                             pattern_contaminants = "^zz|^CON",
                             pattern_decoys = "REV_"){

  annot <- annotation$annot
  atable <- annotation$atable
  annot <- annot |> dplyr::mutate(
    raw.file = gsub("^x|.d.zip$|.raw$","",
                    (basename(annot[[atable$fileName]]))
    ))

  psm <- prolfquapp::tidy_FragPipe_psm(quant_data, column_before_quants = column_before_quants)
  nrPeptides_exp <- psm$nrPeptides
  psm <- psm$data
  psm$qValue <- 1 - psm$PeptideProphet.Probability

  nr <- sum(annot[[annotation$atable$fileName]] %in% sort(unique(psm$channel)))
  logger::log_info("nr : ", nr, " files annotated out of ", length(unique(psm$channel)))
  stopifnot(nr > 0)
  logger::log_info("channels in annotation which are not in psm.tsv file : ", paste(setdiff(annot[[annotation$atable$fileName]],sort(unique(psm$channel))), collapse = " ; ") )
  logger::log_info("channels in psm.tsv which are not in annotation file : ", paste(setdiff(sort(unique(psm$channel)),annot[[annotation$atable$fileName]]), collapse = " ; ") )

  atable$ident_Score = "PeptideProphet.Probability"
  atable$ident_qValue = "qValue"
  atable$hierarchy[["protein_Id"]] <- c("Protein")
  atable$hierarchy[["peptide_Id"]] <- c("Peptide")
  atable$hierarchy[["mod_peptide_Id"]] <- c("Modified.Peptide","Assigned.Modifications")
  # atable$hierarchy[["Spectrum"]] <- c("Spectrum")
  atable$set_response("abundance")

  bycol <- c("channel")
  names(bycol) <- atable$fileName
  psma <- dplyr::inner_join(annot, psm, multiple = "all", by = bycol)
  config <- prolfqua::AnalysisConfiguration$new(atable)
  adata <- prolfqua::setup_analysis(psma, config)
  lfqdata <- prolfqua::LFQData$new(adata, config)
  lfqdata$data

  # build rowAnnotation.
  fasta_annot <- get_annot_from_fasta(fasta_file)
  #psm <- dplyr::left_join(psm, fasta_annot, by = c(Protein = "fasta.id"), multiple = "all")
  #stopifnot(nrow(psm) == nrowPSM)
  fasta_annot <- dplyr::left_join(nrPeptides_exp, fasta_annot, by = c("Protein" = "fasta.id"))

  fasta_annot <- fasta_annot |> dplyr::rename(!!lfqdata$config$table$hierarchy_keys_depth()[1] := !!sym("Protein"))
  fasta_annot <- fasta_annot |> dplyr::rename(description = fasta.header)

  prot_annot <- prolfquapp::ProteinAnnotation$new(
    lfqdata ,
    fasta_annot,
    description = "description",
    cleaned_ids = "proteinname",
    full_id = "protein_Id",
    exp_nr_children = "nrPeptides",
    pattern_contaminants = pattern_contaminants,
    pattern_decoys = pattern_decoys
  )

  lfqdata$remove_small_intensities()
  return(list(lfqdata = lfqdata , protein_annotation = prot_annot))
}


