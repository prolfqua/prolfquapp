#' extect contrasts from dataset
#' @export
#' @examples
#'
#' file <- system.file("application/dataset_csv/dataset 25.csv", package = "prolfquapp")
#' res <- readr::read_csv(file)
#' GRP2 <- make_DEA_config()
#' GRP2 <- dataset_extract_contrasts(res,GRP2)
#' stopifnot(length(GRP2$pop$Contrasts) == 0)
#'
#' file <- system.file("application/dataset_csv/dataset 26.csv", package = "prolfquapp")
#' res <- readr::read_csv(file)
#' GRP2 <- dataset_extract_contrasts(res,GRP2)
#' stopifnot(length(GRP2$pop$Contrasts) == 3)
#'
dataset_extract_contrasts <- function(annot, GRP2) {
  warning("DEPRECATED")
  if ( all(c("ContrastName", "Contrast") %in% colnames(annot)) ) {
    contr <- annot |>
      dplyr::select(all_of(c("ContrastName", "Contrast"))) |>
      dplyr::filter(nchar(!!rlang::sym("Contrast")) > 0)
    Contrasts <- contr$Contrast
    names(Contrasts) <- contr$ContrastName
    GRP2$pop$Contrasts <- Contrasts
  }
  return(GRP2)
}

#' Sanitize grouping variable in annotation file
#' @export
sanitize_grouping_var <- function(annot){
  stopifnot(sum(grepl("^group|^bait|^Experiment", colnames(annot), ignore.case = TRUE)) >= 1)
  groupingVAR <- grep("^group|^bait|^Experiment", colnames(annot), value = TRUE, ignore.case = TRUE)
  if (any(grepl("^bait", groupingVAR, ignore.case = TRUE))) {
    groupingVAR <- grep("^bait", groupingVAR, value = TRUE, ignore.case = TRUE)[1]
  } else {
    groupingVAR <- groupingVAR[1]
  }
  annot[[groupingVAR]] <- gsub("\\.","_",gsub("\\.\\.","\\.",make.names(annot[[groupingVAR]]))) # sanitize group variable
  return(annot)
}


#' Dataset protein annot
#'
#' @export
#' @param msdata data frame
#' @param idcolName name of column with ID's
#' @param protein_annot fasta haeder column
#' @param more_columns more columns to include

build_protein_annot <- function(
    lfqdata,
    msdata,
    idcol = c("protein_Id" = "Protein.Group"),
    cleaned_protein_id = "Protein.Group.2",
    protein_description = "fasta.header",
    nr_children = "nrPeptides",
    more_columns = c("fasta.id")
) {
  proteinID_column = names(idcol)[1]
  msdata <- dplyr::mutate(msdata, !!proteinID_column := !!rlang::sym(idcol) )
  length_protIDs <- length(unique(msdata[[proteinID_column]]))
  prot_annot <- dplyr::select(
    msdata ,
    dplyr::all_of(c( proteinID_column, protein_description, cleaned_protein_id, nr_children, more_columns))) |>
    dplyr::distinct()
  stopifnot( length_protIDs == nrow(prot_annot) )
  prot_annot <- dplyr::rename(prot_annot, description = !!rlang::sym(protein_description))
  prot_annot <- dplyr::rename(prot_annot, IDcolumn = !!rlang::sym(cleaned_protein_id))
  protAnnot <- prolfqua::ProteinAnnotation$new(
    lfqdata , prot_annot,description = "description",
    ids = "IDcolumn",
    nr_children = nr_children)
  return(protAnnot)
}

#' Dataset protein annot
#'
#' @export
#' @param msdata data frame
#' @param idcolName name of column with ID's
#' @param protein_annot fasta haeder column
#' @param more_columns more columns to include
dataset_protein_annot <- function(
    msdata,
    idcol = c("protein_Id" = "Protein.Group"),
    protein_annot = "fasta.header",
    more_columns = c("nrPeptides", "fasta.id")) {
  warning("deprecated! use build_protein_annot")
  proteinID_column = names(idcol)[1]
  msdata <- dplyr::rename(msdata, !!proteinID_column := !!rlang::sym(idcol) )
  prot_annot <- dplyr::select(
    msdata ,
    dplyr::all_of(c( proteinID_column, protein_annot, more_columns))) |>
    dplyr::distinct()
  prot_annot <- dplyr::rename(prot_annot, description = !!rlang::sym(protein_annot))
  # figure out if this is an uniprot database.

  UNIPROT <- mean(grepl("^sp\\||^tr\\|", prot_annot[[proteinID_column]])) > 0.8
  message("uniprot database : ", UNIPROT)

  if (UNIPROT) {
    prot_annot <- prolfqua::get_UniprotID_from_fasta_header(prot_annot, idcolumn = proteinID_column)
    prot_annot <- prot_annot |> dplyr::rename(!!"IDcolumn" := !!rlang::sym("UniprotID"))
  } else {
    prot_annot$IDcolumn <- prot_annot[[proteinID_column]]
  }
  return(prot_annot)
}
#' dataset transform data
#'
#' @export
#'
aggregate_data <- function(lfqdata,
                           agg_method = c("medpolish", "lmrob", "topN"), N = 3) {
  agg_method <- match.arg(agg_method)

  if (agg_method == "topN") {
    aggregator <- lfqdata$get_Aggregator()
    aggregator$sum_topN(N = N)
    lfqdata <- aggregator$lfq_agg
  } else if (agg_method == "lmrob" || agg_method == "medpolish") {
    transformed <- lfqdata$get_Transformer()$intensity_array(log)$lfq
    aggregator <- transformed$get_Aggregator()
    if (agg_method == "lmrob" ) {
      aggregator$lmrob()
    } else if (agg_method == "medpolish") {
      aggregator$medpolish()
    }
    ag <- aggregator$lfq_agg
    tr <- ag$get_Transformer()
    tr <- tr$intensity_array(exp, force = TRUE)
    lfqdata <- tr$lfq
    lfqdata$is_transformed(FALSE)
  } else {
    logger::log_warn("no such aggregator {agg_method}.")
  }
  return(lfqdata)
}

#' compute ibaq values
#' @export
compute_IBAQ_values <- function(lfqdata, protein_annotation) {
  stopifnot(all(c("protein_length", c("nr_tryptic_peptides") %in% colnames(protein_annotation$row_annot))))
  rel_annot <- dplyr::select(protein_annotation$row_annot, protein_Id, protein_length, nr_tryptic_peptides)
  lfqdata$config$table$hierarchyDepth <- 1 # you want to roll up to portein
  lfqdataProtTotal <- prolfquapp::aggregate_data(lfqdata, agg_method = "topN", N = 10000)
  lfqdataProtTotal$data <- dplyr::inner_join(lfqdataProtTotal$data , rel_annot, by = protein_annotation$pID)
  lfqdataProtTotal$data <- lfqdataProtTotal$data |>
    dplyr::mutate(IBAQValue_proteinLength = .data$srm_sum_N / .data$protein_length)
  lfqdataProtTotal$data <- lfqdataProtTotal$data |>
    dplyr::mutate(IBAQValue = .data$srm_sum_N / .data$nr_tryptic_peptides)
  lfqdataProtTotal$config$table$set_response("IBAQValue")
  return(lfqdataProtTotal)
}

