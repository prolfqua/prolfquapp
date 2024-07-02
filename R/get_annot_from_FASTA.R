
.find_cleavage_sites <- function(sequence, pattern = "(K|R)(?!P|$)") {
  # Convert sequence to uppercase
  sequence <- toupper(sequence)
  # Find positions of tryptic sites
  positions <- gregexpr(pattern, sequence, perl = TRUE)[[1]]
  cleavage_sites <- unlist(positions)
  return(cleavage_sites)
}

.compute_peptide_lengths <- function(sequence, cleavage_sites) {
  # Compute the lengths of the peptides
  start_positions <- c(1, cleavage_sites + 1)
  end_positions <- c(cleavage_sites, nchar(sequence))
  peptide_lengths <- end_positions - start_positions + 1
  return(peptide_lengths)
}


#' Compute number of tryptic peptides
#' @export
#' @examples
#' # example code
#'
#' sequence <- "MKGLPRAKSHGSTGWGKRKRNKPK"
#' nr_tryptic_peptides(sequence, min_length=5)
#'
nr_tryptic_peptides <- function(sequence, min_length = 6, max_length = 30){
  peptide_lengths <- .compute_peptide_lengths(sequence, .find_cleavage_sites(sequence))
  res <- sum(peptide_lengths >= min_length & peptide_lengths < max_length)
  return(res)
}


#' extract gene names from uniprot 1sp fasta.headers
#' @param fasta.headers vector with
#' @export
extract_GN <- function(fasta.headers){
  res <- vector(mode = "character", length = length(fasta.headers))
  for (i in seq_along(fasta.headers)) {
    header <- fasta.headers[i]
    res[i] <- if (grepl(".+ GN=(.+) PE=.+", header)) {
      gsub(".+ GN=(.+) PE=.*","\\1",header)} else {""}
  }
  return(res)
}

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
    isUniprot = TRUE,
    min_length = 7,
    max_length = 30) {
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
  fasta_annot$nr_tryptic_peptides <- vapply(fasta, nr_tryptic_peptides, 0, min_length = min_length, max_length = max_length)

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

  if (mean(grepl(".+ GN=(.+) PE=.+",fasta_annot$fasta.header)) > 0.5) {
    fasta_annot <- fasta_annot |> dplyr::mutate(gene_name = extract_GN(.data$fasta.header))
  }

  # remove duplicated id's
  fasta_annot <- fasta_annot[!duplicated(fasta_annot$proteinname),]
  return(fasta_annot)
}


