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
#' @param nr_peptides minimum number of distinct (stripped) peptides per protein (>= 1, default 1)
#' @return list with \code{lfqdata} (LFQData) and
#'   \code{protein_annotation} (ProteinAnnotation)
#' @export
preprocess_SIM <- function(
  quant_data,
  fasta_file,
  annotation,
  pattern_contaminants = "^zz|^CON|Cont_",
  pattern_decoys = "^REV_|^rev_",
  hierarchy_depth = 1,
  nr_peptides = 1
) {
  sim <- prolfqua::sim_lfq_data_peptide_config(Nprot = 50)

  config <- annotation$atable$clone(deep = TRUE)
  config$file_name <- "sample"
  config$sample_name <- "sampleName"
  config$hierarchy[["protein_Id"]] <- "protein_Id"
  config$hierarchy[["peptide_Id"]] <- "peptide_Id"
  config$set_response("abundance")
  config$hierarchy_depth <- hierarchy_depth
  config$nr_children <- "nr_children"
  config$ident_q_value <- "qValue"

  # config$factors maps a factor key to its annotation source column (e.g. G_ -> "group");
  # the sim data only has `group_`, so add the missing source columns.
  raw <- sim$data
  annot <- annotation$annot
  group_src <- config$factors[[setdiff(names(config$factors), c("CONTROL", "Subject_"))[1]]]
  for (fkey in names(config$factors)) {
    src_col <- config$factors[[fkey]]
    if (identical(fkey, "CONTROL") && all(c(src_col, group_src) %in% colnames(annot))) {
      control_map <- dplyr::distinct(annot[, c(group_src, src_col), drop = FALSE])
      raw[[src_col]] <- control_map[[src_col]][match(raw[["group_"]], control_map[[group_src]])]
    } else if (!src_col %in% colnames(raw)) {
      raw[[src_col]] <- raw[["group_"]]
    }
  }

  raw <- prolfquapp::filter_by_peptide_count(raw, "protein_Id", "peptide_Id", nr_peptides)
  adata <- prolfqua::setup_analysis(raw, config)
  lfqdata <- prolfqua::LFQData$new(adata, config)

  # Add contaminant/decoy prefixes (same as sim_data_protAnnot)
  tmp_data <- lfqdata$data_long()
  tmp_data$protein_Id <- prolfquapp::add_RevCon(tmp_data$protein_Id)
  lfqdata$set_data(tmp_data)

  pids <- grep("^zz|^REV", unique(lfqdata$data_long()$protein_Id), value = TRUE, invert = TRUE)
  addannot <- data.frame(protein_Id = pids, description = stringi::stri_rand_strings(length(pids), 13)) |>
    tidyr::separate(protein_Id, c("cleanID", NA), remove = FALSE)
  pannot <- prolfquapp::ProteinAnnotation$new(
    lfqdata,
    addannot,
    description = "description",
    cleaned_ids = "cleanID",
    pattern_contaminants = pattern_contaminants,
    pattern_decoys = pattern_decoys
  )
  pannot$row_annot$nr_tryptic_peptides <- pannot$row_annot$nrPeptides * 2
  pannot$row_annot$protein_length <- pannot$row_annot$nrPeptides * 10
  list(lfqdata = lfqdata, protein_annotation = pannot)
}
