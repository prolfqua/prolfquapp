# ProteinAnnotation ----
#' Decorates LFQData with a row annotation and some protein specific functions.
#'
#' @export
#' @family LFQData
#' @examples
#'
#' istar <-prolfqua::sim_lfq_data_peptide_config()
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' pannot <- ProteinAnnotation$new( lfqdata )
#' pannot$annotate_decoys()
#' pannot$annotate_contaminants()
#' dd <- pannot$clean()
#' tmp <- lfqdata$get_subset(dd)
#' pannot$row_annot
#'
ProteinAnnotation <-
  R6::R6Class("ProteinAnnotation",
              public = list(
                #' @field row_annot data.frame containing further information
                row_annot = NULL,
                #' @field pID column with protein ids
                pID = character(),
                #' @field description name of column containing descriptions
                description = "description",
                #' @field cleaned_ids vector with columns containing addition IDs
                cleaned_ids = character(),
                #' @field nr_children name of columns with the number of peptides
                nr_children = character(),
                contaminant_pattern = character(),
                decoy_pattern = character(),
                #' @description initialize
                #' @param lfqdata data frame from \code{\link{setup_analysis}}
                #' @param row_annot data frame with row annotation. Must have columns matching \code{config$table$hierarchy_keys_depth()}
                #' @param description name of column with description
                #' @param ids names of columns with cleaned Ids
                #' @param nr_peptides additional peptides
                #' @param nr_children column with the number of children
                initialize = function(lfqdata,
                                      row_annot = NULL,
                                      description = NULL,
                                      cleaned_ids = NULL,
                                      nr_children = "nr_peptides",
                                      contaminant_pattern = "^zz|^CON",
                                      decoy_pattern = "^REV_"){
                  self$pID = lfqdata$config$table$hierarchy_keys_depth()[1]
                  self$nr_children = nr_children
                  self$contaminant_pattern = contaminant_pattern
                  self$decoy_pattern = decoy_pattern
                  if ( !is.null(cleaned_ids)) {self$cleaned_ids = cleaned_ids} else {self$cleaned_ids = self$pID}
                  if ( !is.null(description)) {self$description = description} else {self$description = self$pID}
                  if ( !is.null(row_annot)) {
                    stopifnot(self$pID %in% colnames(row_annot))
                    row_annot <- dplyr::filter(row_annot, !!sym(self$pID) %in% lfqdata$data[[self$pID]] )
                    self$row_annot <- row_annot
                  } else {
                    self$row_annot <- distinct(select(lfqdata$data, self$pID))
                  }

                  stopifnot(self$cleaned_ids %in% colnames(self$row_annot))
                  stopifnot(self$description %in% colnames(self$row_annot))

                  if (!self$nr_children %in% colnames(row_annot) ) {
                    warning("no nr_children column specified, computing using nr_obs_experiment function")
                    self$row_annot <- inner_join(
                      self$row_annot,
                      nr_obs_experiment(lfqdata$data, lfqdata$config, name_nr_child = self$nr_children),
                      by = self$pID)
                  }

                },
                #' @description
                #' annotate rev sequences
                #' @param pattern default "REV_"
                annotate_decoys = function() {
                  self$row_annot <- self$row_annot |> mutate(
                    REV = case_when(grepl(self$decoy_pattern, !!sym(self$pID), ignore.case = TRUE) ~ TRUE,
                                    TRUE ~ FALSE))

                  return(sum(self$row_annot$REV))
                },
                #' @description
                #' annotate contaminants
                #' @param pattern default "^zz|^CON"
                annotate_contaminants = function() {
                  self$row_annot <- self$row_annot |> mutate(
                    CON = case_when(grepl(self$contaminant_pattern, !!sym(self$pID), ignore.case = TRUE) ~ TRUE,
                                    TRUE ~ FALSE))
                  return(sum(self$row_annot$CON))
                },
                #' @description
                #' get summary
                #' @param
                get_summary = function() {
                  allProt <- nrow(self$row_annot)
                  contdecoySummary <- data.frame(
                    totalNrOfProteins = allProt,
                    percentOfContaminants = round(self$annotate_contaminants() / allProt * 100, digits = 2),
                    percentOfFalsePositives = round(self$annotate_decoys() / allProt * 100, digits = 2),
                    NrOfProteinsNoDecoys = self$nr_clean()
                  )
                  return(contdecoySummary)
                },
                #' @description get number of neither contaminants nor decoys
                #' @param contaminants remove contaminants
                #' @param decoys remove decoys
                #' return number of cleans
                nr_clean = function(contaminants = TRUE, decoys = TRUE){

                  if (decoys && !("REV" %in% colnames(self$row_annot)) ) { stop("annotate REV") }
                  if (contaminants & !("CON" %in% colnames(self$row_annot)) ) { stop("annotate CON") }

                  res <- if (decoys && contaminants) {
                    sum(!self$row_annot$REV & !self$row_annot$CON)
                  } else if (contaminants) {
                    sum(!self$row_annot$CON)
                  } else if (decoys) {
                    sum(!self$row_annot$REV)
                  } else {
                    nrow(self$row_annot)
                  }
                  return(res)
                },
                #' @description remove REV and CON sequences
                #' @param contaminants remove contaminants
                #' @param decoys remove decoys
                #'
                clean = function(contaminants = TRUE, decoys = TRUE){
                  if (contaminants && !("REV" %in% colnames(self$row_annot)) ) { stop("annotate REV") }
                  if (decoys && !("CON" %in% colnames(self$row_annot)) ) { stop("annotate CON") }
                  res <- if (decoys && contaminants) {
                    filter(self$row_annot , !self$row_annot$REV & !self$row_annot$CON )
                  } else if (contaminants) {
                    filter(self$row_annot , !self$row_annot$CON)
                  } else if (decoys) {
                    filter(self$row_annot , !self$row_annot$REV )
                  } else {
                    self$row_annot
                  }
                  return(res)
                }

              )
  )


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
    more_columns = c("fasta.id"),
    contaminant_pattern = "^zz|^CON",
    decoy_pattern="REV_"
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
  protAnnot <- prolfquapp::ProteinAnnotation$new(
    lfqdata , prot_annot, description = "description",
    cleaned_ids = "IDcolumn",
    nr_children = nr_children,
    contaminant_pattern = contaminant_pattern,
    decoy_pattern = decoy_pattern
    )
  return(protAnnot)
}


