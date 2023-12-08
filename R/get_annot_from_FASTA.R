#' get_annot_from_fasta
#'
#' @export
#' @examples
#'
#' #prolfquapp::get_annot_from_fasta(fasta.files)
#'
#'
get_annot_from_fasta <- function(
    fasta.files,
    rev= "REV_",
    isUniprot = TRUE) {
  fasta <- list()

  if ("connection" %in% class(fasta.files) ) {
    fasta <- seqinr::read.fasta(file = fasta.files, as.string = TRUE, seqtype = "AA")
  } else {
    for (fasta.file in fasta.files) {
      x <- seqinr::read.fasta(file = fasta.file, as.string = TRUE, seqtype = "AA")
      fasta <- c(fasta, x)
    }
  }

  fasta_annot <- prolfqua::matrix_to_tibble(
    data.frame(annot = sapply(fasta, seqinr::getAnnot)), preserve_row_names = NULL
  )
  fasta_annot$protein_length <- vapply(fasta, nchar, 0)

  fasta_annot <- fasta_annot |> tidyr::separate(.data$annot,
                                                c("fasta.id","fasta.header"),
                                                sep = "\\s", extra = "merge")
  fasta_annot <- fasta_annot |> dplyr::mutate(fasta.id = gsub(">","", .data$fasta.id) )
  fasta_annot <- fasta_annot[!grepl(rev,fasta_annot$fasta.id), ]


  if (isUniprot) {
    fasta_annot <- fasta_annot |> dplyr::mutate(proteinname = gsub(".+\\|(.+)\\|.*","\\1", .data$fasta.id) )
  } else {
    fasta_annot <- fasta_annot |> dplyr::mutate(proteinname = .data$fasta.id )
  }
  # remove duplicated id's
  fasta_annot <- fasta_annot[!duplicated(fasta_annot$proteinname),]
  return(fasta_annot)
}
