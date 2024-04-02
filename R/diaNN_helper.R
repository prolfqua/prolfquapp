#' DIANN helper functions
#' @param diann.path path to diann-output.tsv
#' @param fasta.file path to fasta file
#' @param nrPeptides peptide threshold
#' @param Q.Value q value threshold
#' @import data.table
#' @export
#'
read_DIANN_output <- function(diann.path,
                              fasta.file,
                              nrPeptides = 2,
                              Q.Value = 0.01,
                              isUniprot = TRUE,
                              rev = "REV_"
                              ) {
  report2 <- prolfquapp::diann_read_output(diann.path, nrPeptides = nrPeptides, Q.Value = Q.Value)
  if (nrow(report2) == 0) {
    return(NULL)
  }
  report2$raw.file <- gsub("^x|.d.zip$|.d$|.raw$|.mzML$","",basename(gsub("\\\\","/",report2$File.Name)))
  report2$Protein.Group <- sub("zz\\|(.+)\\|.+", "\\1", report2$Protein.Group )

  peptide <- prolfquapp::diann_output_to_peptide(report2)
  peptide$Protein.Group.2 <- sapply(peptide$Protein.Group, function(x){ unlist(strsplit(x, "[ ;]"))[1]} )
  # we need to add the fasta.header information.
  fasta_annot <- get_annot_from_fasta(fasta.file)
  message("Percent of Proteins with description:" ,mean(peptide$Protein.Group.2 %in% fasta_annot$proteinname) * 100)
  # add fasta headers.
  if (nrow(peptide) == 0) {
    return(NULL)
  }
  peptide <- dplyr::left_join(dtplyr::lazy_dt(peptide), dtplyr::lazy_dt(fasta_annot),
                              by = c( Protein.Group.2 = "proteinname")) |>
    dplyr::as_tibble()
  return(peptide)
}
