#' add REV and zz entries - used for testing
#' @export
#'
add_RevCon <- function(stringsAll, pattern_decoys= "REV_", pattern_contaminants = "zz"){
  # Set seed for reproducibility
  set.seed(123)

  dd <- data.frame(idx = unique(stringsAll), tomod = unique(stringsAll))
  # Determine the number of elements to prefix
  n <- length(dd$idx)
  num_rev <- ceiling(0.10 * n) # 10 percent with "REV_"
  num_zz <- ceiling(0.05 * n)  # 5 percent with "ZZ"

  # Randomly select indices for "REV_" prefix
  indices_rev <- sample(1:n, num_rev, replace = FALSE)

  # Apply the "REV_" prefix
  dd$tomod[indices_rev] <- paste(pattern_decoys, dd$tomod[indices_rev], sep = "")

  # Exclude already modified strings and select indices for "ZZ" prefix
  available_indices <- setdiff(1:n, indices_rev)
  indices_zz <- sample(available_indices, num_zz, replace = FALSE)

  # Apply the "ZZ" prefix
  dd$tomod[indices_zz] <- paste(pattern_contaminants, dd$tomod[indices_zz], sep = "")

  res <- merge(data.frame(idx = stringsAll), dd, by = "idx")

  return(res$tomod)
}


# ProteinAnnotation ----
#' Decorates LFQData with a row annotation and some protein specific functions.
#'
#' @export
#' @family LFQData
#' @examples
#'
#' istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 100)
#' lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
#'
#' lfqdata$data$protein_Id <- add_RevCon(lfqdata$data$protein_Id)
#' pids <- grep("^zz|^REV",unique(lfqdata$data$protein_Id),value=TRUE, invert=TRUE)
#' addannot <- data.frame(protein_Id = pids,
#' description = stringi::stri_rand_strings(length(pids), 13))
#' addannot <- addannot |> tidyr::separate(protein_Id, c("cleanID",NA), remove=FALSE)
#' pannot <- ProteinAnnotation$new( lfqdata, addannot, description = "description", cleaned_ids ="cleanID",
#' pattern_contaminants = "^zz", pattern_decoys="^REV" )
#'
#' stopifnot(pannot$annotate_decoys() == 10)
#' stopifnot(pannot$annotate_contaminants() == 5)
#' dd <- pannot$clean()
#' pannot$nr_clean()
#' pannot$get_summary()
#' stopifnot(nrow(dd) == 85)
#' tmp <- lfqdata$get_subset(dd)
#' dx <- pannot$clean(contaminants = TRUE, decoys = FALSE)
#' stopifnot(nrow(dx) == 95)
#' dx <- pannot$clean(contaminants = FALSE, decoys = TRUE)
#' stopifnot(nrow(dx) == 90)
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
                #' @field pattern_contaminants pattern_contaminants
                pattern_contaminants = character(),
                #' @field pattern_decoys pattern_decoys
                pattern_decoys = character(),
                #' @description initialize
                #' @param lfqdata data frame from \code{\link{setup_analysis}}
                #' @param row_annot data frame with row annotation. Must have columns matching \code{config$table$hierarchy_keys_depth()}
                #' @param description name of column with description
                #' @param cleaned_ids names of columns with cleaned Ids
                #' @param nr_children column with the number of children
                #' @param pattern_contaminants pattern_contaminants
                #' @param pattern_decoys pattern_decoys
                initialize = function(lfqdata,
                                      row_annot = NULL,
                                      description = NULL,
                                      cleaned_ids = NULL,
                                      nr_children = "nr_peptides",
                                      pattern_contaminants = "^zz|^CON",
                                      pattern_decoys = "^REV_"){
                  self$pID = lfqdata$config$table$hierarchy_keys_depth()[1]
                  self$nr_children = nr_children
                  self$pattern_contaminants = pattern_contaminants
                  self$pattern_decoys = pattern_decoys
                  if ( !is.null(cleaned_ids)) {self$cleaned_ids = cleaned_ids} else {self$cleaned_ids = self$pID}
                  if ( !is.null(description)) {self$description = description} else {self$description = self$pID}

                  self$row_annot <- dplyr::distinct(dplyr::select(lfqdata$data, self$pID))
                  if ( !is.null(row_annot)) {
                    stopifnot(self$pID %in% colnames(row_annot))
                    #row_annot <- dplyr::filter(row_annot, !!sym(self$pID) %in% lfqdata$data[[self$pID]] )
                    self$row_annot <- dplyr::left_join(self$row_annot, row_annot, by = self$pID)
                  }
                  stopifnot(self$cleaned_ids %in% colnames(self$row_annot))
                  stopifnot(self$description %in% colnames(self$row_annot))
                  if (!self$nr_children %in% colnames(row_annot) ) {
                    warning("no nr_children column specified, computing using nr_obs_experiment function")
                    cf <- lfqdata$config$clone(deep=TRUE)
                    cf$table$hierarchyDepth <- 1
                    self$row_annot <- dplyr::inner_join(
                      self$row_annot,
                      prolfqua::nr_obs_experiment(lfqdata$data, cf, name_nr_child = self$nr_children),
                      by = self$pID)
                  }
                },
                #' @description
                #' annotate rev sequences
                #' @param pattern default "REV_"
                annotate_decoys = function() {
                  self$row_annot <- self$row_annot |> dplyr::mutate(
                    REV = dplyr::case_when(grepl(self$pattern_decoys, !!sym(self$pID), ignore.case = TRUE) ~ TRUE,
                                           TRUE ~ FALSE))

                  return(sum(self$row_annot$REV))
                },
                #' @description
                #' annotate contaminants
                #' @param pattern default "^zz|^CON"
                annotate_contaminants = function() {
                  self$row_annot <- self$row_annot |> dplyr::mutate(
                    CON = dplyr::case_when(grepl(self$pattern_contaminants, !!sym(self$pID), ignore.case = TRUE) ~ TRUE,
                                           TRUE ~ FALSE))
                  return(sum(self$row_annot$CON))
                },
                #' @description
                #' get summary
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
                    dplyr::filter(self$row_annot , !self$row_annot$REV & !self$row_annot$CON )
                  } else if (contaminants) {
                    dplyr::filter(self$row_annot , !self$row_annot$CON)
                  } else if (decoys) {
                    dplyr::filter(self$row_annot , !self$row_annot$REV )
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
#' @examples
#' # example code
#'
build_protein_annot <- function(
    lfqdata,
    msdata,
    idcol = c("protein_Id" = "Protein.Group"),
    cleaned_protein_id = "Protein.Group.2",
    protein_description = "fasta.header",
    nr_children = "nrPeptides",
    more_columns = c("fasta.id"),
    pattern_contaminants = "^zz|^CON",
    pattern_decoys="REV_"
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
    pattern_contaminants = pattern_contaminants,
    pattern_decoys = pattern_decoys
  )
  return(protAnnot)
}


