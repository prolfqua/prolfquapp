#' get petpide.txt and fasta file location in folder
#' @return list with paths to data and fasta
#' @export
get_MSstats_files <- function(path){
  msstats.path <- grep("msstats.*\\.(csv|tsv)$", dir(path = path, recursive = TRUE, full.names = TRUE),
                       value = TRUE,
                       ignore.case = TRUE)
  fasta.files <- grep("*\\.fasta$|*\\.fas$", dir(path = path, recursive = TRUE, full.names = TRUE),
                      ignore.case = TRUE, value = TRUE)

  if (any(grepl("database[0-9]*.fasta$", fasta.files))) {
    fasta.files <- grep("database[0-9]*.fasta$", fasta.files, value = TRUE)
  }
  if (length(fasta.files) == 0) {
    logger::log_error("No fasta file found!")
    stop()
  }
  if (length(msstats.path) > 1) {
    logger::log_warn("more then 1 msstats.tsv file found :", length(msstats.path), ". Returning first.")
  }
  return(list(data = msstats.path[1], fasta = fasta.files))
}


#' read MSstats.csv files and rollup to ProteinSequence level.
#' @export
read_msstats <- function(file){
  msstats <- readr::read_csv(file)
  msstats <- msstats |> dplyr::select(-all_of(c("Condition","BioReplicate")))

  peptideLevelInt <- msstats |>
    dplyr::group_by(dplyr::across(c("ProteinName", "PeptideSequence", "IsotopeLabelType", "Run" ))) |>
    dplyr::summarise(nr_children = dplyr::n(),
                     Intensity = sum(Intensity, na.rm = TRUE),
                     .groups = "drop")
  peptideLevelInt <- peptideLevelInt |> dplyr::mutate(Intensity = ifelse(Intensity < 1e-10, NA, Intensity))
  return(peptideLevelInt)
}

#' preprocess MSstats fragpipe
#' @export
#'
preprocess_MSstats_FPDIA <- function(quant_data,
                                  fasta_file,
                                  annotation,
                                  pattern_contaminants = "",
                                  pattern_decoys = ""){

  annot <- annotation$annot
  atable <- annotation$atable
  annot <- annot |> dplyr::mutate(
    !!annotation$atable$fileName := (gsub("^x|.d.zip$|.raw$","",
                                          (basename(annot[[atable$fileName]])))
    ))

  peptide <- read_msstats(quant_data)
  peptide$nr_peptides <- 1

  nrPeptides_exp <- peptide |> dplyr::select(all_of(c("ProteinName", "PeptideSequence"))) |>
    dplyr::distinct() |>
    dplyr::group_by(dplyr::across("ProteinName")) |>
    dplyr::summarize(nrPeptides = dplyr::n())


  nr <- sum(annot[[annotation$atable$fileName]] %in% sort(unique(peptide$Run)))
  logger::log_info("nr : ", nr, " files annotated out of ", length(unique(peptide$Run)))
  stopifnot(nr > 0)
  logger::log_info("channels in annotation which are not in peptide.txt file : ",
                   paste(setdiff(annot[[annotation$atable$fileName]],sort(unique(peptide$Run))), collapse = " ; ") )
  logger::log_info("channels in peptide.txt which are not in annotation file : ",
                   paste(setdiff(sort(unique(peptide$Run)),annot[[annotation$atable$fileName]]), collapse = " ; ") )


  peptide$qValue <- 0
  atable$ident_qValue = "qValue"
  atable$hierarchy[["protein_Id"]] <- c("ProteinName")
  atable$hierarchy[["peptide_Id"]] <- c("PeptideSequence")
  atable$nr_children = "nr_peptides"
  atable$set_response("Intensity")
  atable$hierarchyDepth <- 1

  bycol <- c("Run")
  names(bycol) <- atable$fileName
  apeptide <- dplyr::inner_join(annot, peptide, multiple = "all", by = bycol)

  config <- prolfqua::AnalysisConfiguration$new(atable)
  adata <- prolfqua::setup_analysis(apeptide, config)
  lfqdata <- prolfqua::LFQData$new(adata, config)
  logger::log_info("Start reading fasta: ", fasta_file )
  fasta_annot <- get_annot_from_fasta(fasta_file, pattern_decoys = pattern_decoys)
  logger::log_info("Finished reading fasta: ", fasta_file )

  fasta_annot <- dplyr::left_join(nrPeptides_exp, fasta_annot, by = c("ProteinName" = "proteinname"))

  fasta_annot <- fasta_annot |> dplyr::rename(!!lfqdata$config$table$hierarchy_keys_depth()[1] := !!rlang::sym("ProteinName"))
  fasta_annot <- fasta_annot |> dplyr::rename(description = fasta.header)

  prot_annot <- prolfquapp::ProteinAnnotation$new(
    lfqdata ,
    fasta_annot,
    description = "description",
    cleaned_ids = "protein_Id",
    full_id = "fasta.id",
    exp_nr_children = "nrPeptides",
    pattern_contaminants = pattern_contaminants,
    pattern_decoys = pattern_decoys
  )
  lfqdata$remove_small_intensities()
  return(list(lfqdata = lfqdata , protein_annotation = prot_annot))
}




#' preprocess MQ peptide
#' @export
#'
preprocess_MSstats <- function(quant_data,
                                  fasta_file,
                                  annotation,
                                  pattern_contaminants = "^zz|^CON|Cont_",
                                  pattern_decoys = "^REV_|^rev_"){

  annot <- annotation$annot
  atable <- annotation$atable
  annot <- annot |> dplyr::mutate(
    !!annotation$atable$fileName := (gsub("^x|.d.zip$|.raw$","",
                                          (basename(annot[[atable$fileName]])))
    ))

  peptide <- read_msstats(quant_data)
  peptide$nr_peptides <- 1

  nrPeptides_exp <- peptide |> dplyr::select(all_of(c("ProteinName", "PeptideSequence"))) |>
    dplyr::distinct() |>
    dplyr::group_by(dplyr::across("ProteinName")) |>
    dplyr::summarize(nrPeptides = dplyr::n())


  nr <- sum(annot[[annotation$atable$fileName]] %in% sort(unique(peptide$Run)))
  logger::log_info("nr : ", nr, " files annotated out of ", length(unique(peptide$Run)))
  stopifnot(nr > 0)
  logger::log_info("channels in annotation which are not in peptide.txt file : ",
                   paste(setdiff(annot[[annotation$atable$fileName]],sort(unique(peptide$Run))), collapse = " ; ") )
  logger::log_info("channels in peptide.txt which are not in annotation file : ",
                   paste(setdiff(sort(unique(peptide$Run)),annot[[annotation$atable$fileName]]), collapse = " ; ") )


  peptide$qValue <- 0
  atable$ident_qValue = "qValue"
  atable$hierarchy[["protein_Id"]] <- c("ProteinName")
  atable$hierarchy[["peptide_Id"]] <- c("PeptideSequence")
  atable$nr_children = "nr_peptides"
  atable$set_response("Intensity")
  atable$hierarchyDepth <- 1

  bycol <- c("Run")
  names(bycol) <- atable$fileName
  apeptide <- dplyr::inner_join(annot, peptide, multiple = "all", by = bycol)

  config <- prolfqua::AnalysisConfiguration$new(atable)
  adata <- prolfqua::setup_analysis(apeptide, config)
  lfqdata <- prolfqua::LFQData$new(adata, config)
  logger::log_info("Start reading fasta: ", fasta_file )
  fasta_annot <- get_annot_from_fasta(fasta_file, pattern_decoys = pattern_decoys)
  logger::log_info("Finished reading fasta: ", fasta_file )

  fasta_annot <- dplyr::left_join(nrPeptides_exp, fasta_annot, by = c("ProteinName" = "fasta.id"))

  fasta_annot <- fasta_annot |> dplyr::rename(!!lfqdata$config$table$hierarchy_keys_depth()[1] := !!rlang::sym("ProteinName"))
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


