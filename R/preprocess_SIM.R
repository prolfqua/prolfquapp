#' Discover simulated data files (no real files needed)
#'
#' Returns placeholder paths. Used by the SIM preprocessor for testing
#' CMD scripts without real quantification data.
#'
#' @param path ignored
#' @return list with \code{data} and \code{fasta} placeholder strings
#' @export
get_SIM_files <- function(path) {
  list(data = "simulated", fasta = "simulated")
}

#' Preprocess simulated data
#'
#' Returns a simulated LFQData + ProteinAnnotation, bypassing all file I/O.
#' Designed for integration testing of CMD scripts via \code{--software SIM}.
#'
#' The simulated data is reconfigured to use the annotation's factor prefix
#' so that contrasts derived from the annotation match the model terms.
#'
#' @param quant_data ignored (placeholder)
#' @param fasta_file ignored (placeholder)
#' @param annotation annotation list from \code{\link{read_annotation}}
#' @param pattern_contaminants regex for contaminant proteins
#' @param pattern_decoys regex for decoy proteins
#' @param hierarchy_depth 1 = protein level, 2 = peptide level
#' @return list with \code{lfqdata} (LFQData) and
#'   \code{protein_annotation} (ProteinAnnotation)
#' @export
preprocess_SIM <- function(
  quant_data,
  fasta_file,
  annotation,
  pattern_contaminants = "^zz|^CON|Cont_",
  pattern_decoys = "^REV_|^rev_",
  hierarchy_depth = 1
) {
  sim <- prolfqua::sim_lfq_data_peptide_config(Nprot = 50)

  # Clone annotation atable and configure for simulated data columns
  atable <- annotation$atable$clone(deep = TRUE)
  atable$fileName <- "sample"
  atable$sampleName <- "sampleName"
  atable$hierarchy[["protein_Id"]] <- "protein_Id"
  atable$hierarchy[["peptide_Id"]] <- "peptide_Id"
  atable$set_response("abundance")
  atable$hierarchyDepth <- hierarchy_depth
  atable$nr_children <- "nr_children"
  atable$ident_qValue <- "qValue"

  # Map annotation factor columns onto simulated data
  # atable$factors maps e.g. G_ -> "group", meaning:
  # "create factor column G_ from source column group"
  # The sim data has group_ but not group, so add the source column
  raw <- sim$data
  for (fkey in names(atable$factors)) {
    src_col <- atable$factors[[fkey]]
    if (!src_col %in% colnames(raw)) {
      raw[[src_col]] <- raw[["group_"]]
    }
  }

  config <- prolfqua::AnalysisConfiguration$new(atable)
  adata <- prolfqua::setup_analysis(raw, config)
  lfqdata <- prolfqua::LFQData$new(adata, config)

  # Add contaminant/decoy prefixes (same as sim_data_protAnnot)
  lfqdata$data$protein_Id <- prolfquapp::add_RevCon(
    lfqdata$data$protein_Id
  )

  pids <- grep(
    "^zz|^REV",
    unique(lfqdata$data$protein_Id),
    value = TRUE,
    invert = TRUE
  )
  addannot <- data.frame(
    protein_Id = pids,
    description = stringi::stri_rand_strings(length(pids), 13)
  )
  addannot <- addannot |>
    tidyr::separate(
      protein_Id, c("cleanID", NA), remove = FALSE
    )
  pannot <- prolfquapp::ProteinAnnotation$new(
    lfqdata,
    addannot,
    description = "description",
    cleaned_ids = "cleanID",
    pattern_contaminants = pattern_contaminants,
    pattern_decoys = pattern_decoys
  )
  pannot$row_annot$nr_tryptic_peptides <- pannot$row_annot$nr_peptides * 2
  pannot$row_annot$protein_length <- pannot$row_annot$nr_peptides * 10
  pannot$row_annot$nrPeptides <- pannot$row_annot$nr_peptides

  list(
    lfqdata = lfqdata,
    protein_annotation = pannot
  )
}
