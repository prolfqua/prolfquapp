#' get_dummy_files
#' @export
get_dummy_files <- function(path) {
  return(list(data = "data.path", fasta = "fasta.files.path"))
}

#' preprocess_dummy
#' @export
preprocess_dummy <- function(quant_data,
                             fasta_file,
                             annotation,
                             pattern_contaminants = "^zz|^CON|Cont_",
                             pattern_decoys = "^REV_|^rev"
){
  return(list(lfqdata = "prolfqua::LFQData$new()" , protein_annotation = "prolfquapp::ProteinAnnotation$new()" ))
}
