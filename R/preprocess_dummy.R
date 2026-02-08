#' get_dummy_files
#' @param path path to data directory
#' @export
get_dummy_files <- function(path) {
  return(list(data = "data.path", fasta = "fasta.files.path"))
}

#' preprocess_dummy
#' @param quant_data path to quantification data file
#' @param fasta_file path to fasta file(s)
#' @param annotation annotation list from read_annotation
#' @param pattern_contaminants regex pattern for contaminants
#' @param pattern_decoys regex pattern for decoys
#' @export
preprocess_dummy <- function(quant_data,
                             fasta_file,
                             annotation,
                             pattern_contaminants = "^zz|^CON|Cont_",
                             pattern_decoys = "^REV_|^rev") {
  return(list(lfqdata = "prolfqua::LFQData$new()", protein_annotation = "prolfquapp::ProteinAnnotation$new()"))
}
