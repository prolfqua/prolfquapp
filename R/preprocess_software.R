#' collects preprocess methods for various software
#' @export
#'
preprocess_software <- function(indir,
                                annotation,
                                pattern_contaminants ="^zz|^CON|Cont_" , pattern_decoys = "^rev_|^REV_",
                                software = c("DIANN", "FP_TMT", "FP_multisite", "MAXQUANT", "MSSTATS", "MSSTATS_FP_DIA")) {
  software <- match.arg(software)
  if (software == "DIANN") {
    files <- prolfquapp::get_DIANN_files(indir)
    logger::log_info("Files data: ", paste(files$data, collapse = "; "))
    logger::log_info("Files fasta: ", paste0(files$fasta, collapse = "; "))
    xd <- prolfquapp::preprocess_DIANN(
      quant_data = files$data,
      fasta_file = files$fasta,
      annotation = annotation,
      q_value = 0.1,
      pattern_contaminants = pattern_contaminants,
      pattern_decoys = pattern_decoys
    )
  } else if (software == "FP_TMT") {
    files <- prolfquapp::get_FP_PSM_files(indir)
    logger::log_info("Files data: ", paste(files$data, collapse = "; "))
    logger::log_info("Files fasta: ", paste0(files$fasta, collapse = "; "))
    xd <- prolfquapp::preprocess_FP_PSM(
      quant_data = files$data,
      fasta_file = files$fasta,
      annotation = annotation,
      purity_threshold = 0.5,
      PeptideProphetProb = 0.9,
      pattern_contaminants = pattern_contaminants,
      pattern_decoys = pattern_decoys
    )
  } else if (software == "FP_multisite") {
    files <- prolfquapp::get_FP_multi_site_files(indir)
    logger::log_info("Files data: ", paste(files$data, collapse = "; "))
    logger::log_info("Files fasta: ", paste0(files$fasta, collapse = "; "))
    xd <- prolfquapp::preprocess_FP_multi_site(
      files$data[1],
      files$fasta,
      annotation,
      pattern_contaminants = pattern_contaminants,
      pattern_decoys = pattern_decoys)
  } else if (software == "FP_combined_STY") {
    return(NULL)
  } else if (software == "MAXQUANT") {
    files <- prolfquapp::get_MQ_peptide_files(indir)
    logger::log_info("Files data: ", paste(files$data, collapse = "; "))
    logger::log_info("Files fasta: ", paste0(files$fasta, collapse = "; "))

    xd <- prolfquapp::preprocess_MQ_peptide(
      quant_data = files$data,
      fasta_file = files$fasta,
      annotation = annotation,
      pattern_contaminants = pattern_contaminants,
      pattern_decoys = pattern_decoys
    )
  } else if (software == "MSSTATS") {
    files <- prolfquapp::get_MSstats_files(indir)
    logger::log_info("Files data: ", paste(files$data, collapse = "; "))
    logger::log_info("Files fasta: ", paste0(files$fasta, collapse = "; "))
    xd <- prolfquapp::preprocess_MSstats(
      quant_data = files$data,
      fasta_file = files$fasta,
      annotation = annotation,
      pattern_contaminants = pattern_contaminants,
      pattern_decoys = pattern_decoys
    )
  } else if (software == "MSSTATS_FP_DIA") {
    files <- prolfquapp::get_MSstats_files(indir)
    logger::log_info("Files data: ", paste(files$data, collapse = "; "))
    logger::log_info("Files fasta: ", paste0(files$fasta, collapse = "; "))
    xd <- prolfquapp::preprocess_MSstats_FPDIA(
      quant_data = files$data,
      fasta_file = files$fasta,
      annotation = annotation,
      pattern_contaminants = pattern_contaminants,
      pattern_decoys = pattern_decoys
    )
  } else {
    logger::log_error("no such software :" , software)
    stop("no such software.")
    return(NULL)

  }
  return(list(xd = xd, files = files))
}
