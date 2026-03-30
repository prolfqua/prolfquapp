#' Create an example DEAnalyse object from simulated data
#'
#' Builds a complete DEAnalyse R6 object using simulated peptide data.
#' Useful for vignette defaults, examples, and testing.
#'
#' @param Nprot number of simulated proteins (default 100)
#' @return a \code{DEAnalyse} R6 object with contrasts computed and annotated
#' @export
#' @examples
#' dea <- example_deanalyse(Nprot = 10)
#' dea$contrast_results[[dea$default_model]]$get_contrasts()
#'
example_deanalyse <- function(Nprot = 100) {
  pep <- prolfqua::sim_lfq_data_peptide_config(Nprot = Nprot)
  pep <- prolfqua::LFQData$new(pep$data, pep$config)

  pA <- data.frame(protein_Id = unique(pep$data$protein_Id))
  pA$description <- paste0(pA$protein_Id, "_description")
  pA <- ProteinAnnotation$new(pep, row_annot = pA, description = "description")

  GRP2 <- make_DEA_config_R6()
  GRP2$processing_options$diff_threshold <- 0.2
  GRP2$processing_options$transform <- "robscale"

  contrasts <- c(
    "AVsC" = "group_A - group_Ctrl",
    "BVsC" = "group_B - group_Ctrl"
  )

  data_prep <- ProteinDataPrep$new(pep, pA, GRP2)
  data_prep$cont_decoy_summary()
  data_prep$remove_cont_decoy()
  data_prep$aggregate()
  data_prep$transform_data()

  deanalyse <- data_prep$build_deanalyse(contrasts)
  deanalyse$build_default()
  deanalyse$get_annotated_contrasts()
  deanalyse$filter_contrasts()

  deanalyse
}

#' Create an example QC_generator object from simulated data
#'
#' Builds a QC_generator R6 object with simulated peptide data and
#' fake FASTA-derived columns (protein_length, nr_tryptic_peptides)
#' so that IBAQ computation works. Useful for vignette defaults and testing.
#'
#' @param Nprot number of simulated proteins (default 100)
#' @return a \code{QC_generator} R6 object
#' @export
#' @examples
#' pap <- example_qc_generator(Nprot = 10)
#' pap$get_prot_data()
#'
example_qc_generator <- function(Nprot = 100) {
  res <- sim_data_protAnnot(Nprot = Nprot)
  GRP2 <- make_DEA_config_R6()
  GRP2$set_zipdir_name()
  dir.create(GRP2$get_zipdir(), showWarnings = FALSE, recursive = TRUE)
  # Alias nr_peptides -> nrPeptides (expected by QC report templates)
  res$pannot$row_annot$nrPeptides <- res$pannot$row_annot$nr_peptides
  pap <- QC_generator$new(res$lfqdata, res$pannot, GRP2)
  pap$get_prot_data()
  pap
}
