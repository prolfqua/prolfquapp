#' Compute number of tryptic peptides
#' @param sequence amino acid sequence
#' @param min_length minimum peptide length
#' @param max_length maximum peptide length
#' @export
#' @examples
#' # example code
#'
#' sequence <- "MKGLPRAKSHGSTGWGKRKRNKPK"
#' nr_tryptic_peptides(sequence, min_length=5)
#'
nr_tryptic_peptides <- function(sequence, min_length = 6, max_length = 30) {
  cleavage_sites <- stringr::str_locate_all(toupper(sequence), "(K|R)(?!P|$)")[[1]][, "end"]
  peptide_lengths <- c(cleavage_sites, nchar(sequence)) - c(0, cleavage_sites)
  sum(peptide_lengths >= min_length & peptide_lengths < max_length)
}


#' extract gene names from uniprot 1sp fasta.headers
#' @param fasta.headers vector with
#' @export
extract_GN <- function(fasta.headers) {
  res <- character(length(fasta.headers))
  hit <- grepl(".+ GN=(.+) PE=.+", fasta.headers)
  res[hit] <- gsub(".+ GN=(.+) PE=.*", "\\1", fasta.headers[hit])
  res
}


#' get_annot_from_fasta
#'
#' @param fasta.files path to fasta file(s) or connection
#' @param pattern_decoys regex for decoy sequence IDs
#' @param isUniprot if TRUE parse UniProt-style headers
#' @param min_length minimum tryptic peptide length
#' @param max_length maximum tryptic peptide length
#' @param include_seq if TRUE include protein sequences
#' @export
#' @examples
#' fasta_text <- c(
#'   ">sp|P00001|TEST Protein OS=Human GN=TEST PE=1 SV=1",
#'   "MKRISTTITTT",
#'   ">REV_sp|P00002|DECOY Protein OS=Human GN=DECOY PE=1 SV=1",
#'   "MPEPTIDER"
#' )
#' fasta_conn <- textConnection(fasta_text)
#' get_annot_from_fasta(fasta_conn, pattern_decoys = "^REV_")
#' close(fasta_conn)
#'
get_annot_from_fasta <- function(
  fasta.files,
  pattern_decoys = "^REV_|^rev_",
  isUniprot = TRUE,
  min_length = 7,
  max_length = 30,
  include_seq = FALSE
) {
  read_fasta <- function(file) seqinr::read.fasta(file = file, as.string = TRUE, seqtype = "AA")
  fasta <- if (inherits(fasta.files, "connection")) {
    read_fasta(fasta.files)
  } else {
    do.call(
      c,
      lapply(fasta.files, function(fasta.file) {
        logger::log_info("get_annot : ", fasta.file)
        read_fasta(fasta.file)
      })
    )
  }
  logger::log_info("get_annot : finished reading")

  fasta_annot <- data.frame(annot = vapply(fasta, seqinr::getAnnot, ""), sequence = as.character(fasta)) |>
    tidyr::separate(.data$annot, c("fasta.id", "fasta.header"), sep = "\\s", extra = "merge") |>
    dplyr::mutate(fasta.id = gsub("^>|;", "", .data$fasta.id))
  logger::log_info("get_annot : all seq : ", nrow(fasta_annot))
  # Pure parser: decoy removal and protein-ID uniqueness are resolved downstream by
  # ProteinAnnotation. `pattern_decoys` is accepted for backward compatibility only.
  logger::log_info("get_annot : isUniprot : ", isUniprot)
  fasta_annot$proteinname <- if (isUniprot) {
    gsub(".+\\|(.+)\\|.*", "\\1", fasta_annot$fasta.id)
  } else {
    fasta_annot$fasta.id
  }
  if (sum(grepl(".+ GN=(.+) PE=.+", fasta_annot$fasta.header)) > 1) {
    fasta_annot$gene_name <- extract_GN(fasta_annot$fasta.header)
    logger::log_info("get_annot : extracted gene names")
  }
  fasta_annot$protein_length <- vapply(fasta_annot$sequence, nchar, 0)
  fasta_annot$nr_tryptic_peptides <- vapply(
    fasta_annot$sequence,
    nr_tryptic_peptides,
    0,
    min_length = min_length,
    max_length = max_length
  )
  logger::log_info("get_annot : nr of tryptic peptides per protein computed.")
  if (!include_seq) {
    fasta_annot$sequence <- NULL
  }
  fasta_annot
}
