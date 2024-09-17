#' Generate differential expression analysis reports
#' @param prot_annot ProteinAnnotation$new
#' @export
#'
generate_DEA_reports <- function(lfqdata, GRP2, prot_annot) {
  warning("DEPRECATED generate_DEA_reports -> use generate_DEA_reports2")
  # Compute all possible 2 Grps to avoid specifying reference.
  res <- list()
  levels  <- lfqdata$factors() |>
    dplyr::select(
      Group_ = Group_,
      control = dplyr::starts_with("control", ignore.case = TRUE)) |>
    dplyr::distinct()
  logger::log_info("levels : ", paste(levels, collapse = " "))
  if (!length(levels$Group_) > 1) {
    logger::log_error("not enough group levels_ to make comparisons.")
  }

  # Contrasts are already defined
  if (!is.null(GRP2$pop$Contrasts)) {
    logger::log_info("CONTRAST : ", paste( GRP2$pop$Contrasts, collapse = "\n"))
    lfqwork <- lfqdata$get_copy()
    lfqwork$data <- lfqdata$data |> dplyr::filter(.data$Group_ %in% levels$Group_)

    grp2 <- prolfquapp::make_DEA_report(
      lfqwork,
      prot_annot,
      GRP2)
    res[["Groups_vs_Controls"]] <- grp2
    return(res)
  }

  ## Generate contrasts from dataset
  if (!is.null(levels$control)) {
    Contrasts <- character()
    Names <- character()
    for (i in 1:nrow(levels)) {
      for (j in 1:nrow(levels)) {
        if (i != j) {
          if (levels$control[j] == "C") {
            cat(levels$Group_[i], levels$Group_[j], "\n")
            Contrasts <- c(Contrasts, paste0("Group_",levels$Group_[i], " - ", "Group_",levels$Group_[j]))
            Names <- c(Names, paste0(levels$Group_[i], "_vs_", levels$Group_[j]))
          }
        }
      }
    }

    names(Contrasts) <- Names
    GRP2$pop$Contrasts <- Contrasts
    logger::log_info("CONTRAST : ", paste( GRP2$pop$Contrasts, collapse = " "))
    lfqwork <- lfqdata$get_copy()
    lfqwork$data <- lfqdata$data |> dplyr::filter(.data$Group_ %in% levels$Group_)

    grp2 <- prolfquapp::make_DEA_report(
      lfqwork,
      prot_annot,
      GRP2)
    res[["Groups_vs_Controls"]] <- grp2
    return(res)
  } else {
    # create all possible 2 grp comparisons
    for (i in seq_along(levels$Group_)) {
      for (j in seq_along(levels$Group_)) {
        if (i != j) {
          cat("COMPARING : ", levels$Group_[i], " vs " , levels$Group_[j], "\n")
          GRP2$pop$Contrasts <- paste0("Group_",levels$Group_[i], " - ", "Group_",levels$Group_[j])
          names(GRP2$pop$Contrasts) <- paste0(levels$Group_[i], "_vs_", levels$Group_[j])
          logger::log_info("CONTRAST : ", GRP2$pop$Contrasts, collapse = " ")
          lfqwork <- lfqdata$get_copy()
          lfqwork$data <- lfqdata$data |> dplyr::filter(.data$Group_ == levels$Group_[i] | .data$Group_ == levels$Group_[j])
          grp2 <- prolfquapp::make_DEA_report(lfqwork,
                                              prot_annot,
                                              GRP2)

          name <- paste0(levels$Group_[i], "_vs_", levels$Group_[j])
          res[[name]] <- grp2
        }# end for 1
      }# end for 2
    }
    return(res)
  }


}




#' Create DEA report in html and write data to xlsx table
#'
#' For use examples see run_scripts directory
#' @rdname make_DEA_report
#' @param startdata table in long format
#' @param atable AnalysisTableAnnotation annotate startdata table
#' @param GRP2 list with named arguments i.e. Contrasts, projectID, projectName, workunitID, nrPeptides, Diffthreshold, FDRthreshold
#' @param protein_annot column with protein description e.g. (fasta header)
#' @export
#' @family workflow
#'
make_DEA_report <- function(lfqdata,
                            protAnnot,
                            GRP2
) {
  ### Do some type of data normalization (or do not)
  warning("DEPRECATED make_DEA_report -> use make_DEA_report2")
  transformed <- transform_lfqdata(
    lfqdata,
    method = GRP2$pop$transform,
    internal = GRP2$pop$internal
  )


  allProt <- nrow( protAnnot$row_annot )
  GRP2$RES <- list()
  GRP2$RES$Summary <- data.frame(
    totalNrOfProteins = allProt,
    percentOfContaminants = round(protAnnot$annotate_contaminants()/allProt * 100 , digits = 2),
    percentOfFalsePositives  = round(protAnnot$annotate_decoys()/allProt * 100 , digits = 2),
    NrOfProteinsNoDecoys = protAnnot$nr_clean()
  )
  GRP2$RES$rowAnnot <- protAnnot

  if (GRP2$pop$removeCon || GRP2$pop$removeDecoys) {
    lfqdata <- lfqdata$get_subset(protAnnot$clean(contaminants = GRP2$pop$removeCon, decoys = GRP2$pop$removeDecoys))
    transformed <- transformed$get_subset(protAnnot$clean())
    logger::log_info(
      paste0("removing contaminants and reverse sequences with patterns:",
             GRP2$pop$contpattern,
             GRP2$pop$revpattern ))
  }

  lfqdata$rename_response("protein_abundance")
  transformed$rename_response("normalized_protein_abundance")

  GRP2$RES$lfqData <- lfqdata
  GRP2$RES$transformedlfqData <- transformed

  ################## Run Modelling ###############
  factors <- transformed$config$table$factor_keys()[!grepl("^control", transformed$config$table$factor_keys() , ignore.case = TRUE)]

  if (is.null(GRP2$pop$interaction) || !GRP2$pop$interaction ) {
    formula <- paste0(transformed$config$table$get_response(), " ~ ",
                      paste(factors, collapse = " + "))
  } else {
    formula <- paste0(transformed$config$table$get_response(), " ~ ",
                      paste(factors, collapse = " * "))
  }

  message("FORMULA :",  formula)
  GRP2$RES$formula <- formula
  formula_Condition <-  prolfqua::strategy_lm(formula)
  # specify model definition
  modelName  <- "Model"

  mod <- prolfqua::build_model(
    transformed,
    formula_Condition,
    subject_Id = transformed$config$table$hierarchy_keys() )

  logger::log_info("fitted model with formula : {formula}")
  GRP2$RES$models <- mod

  contr <- prolfqua::Contrasts$new(mod, GRP2$pop$Contrasts,
                                   modelName = "Linear_Model")

  conrM <- prolfqua::ContrastsModerated$new(
    contr)

  if (is.null(GRP2$pop$missing) || GRP2$pop$missing ) {
    mC <- prolfqua::ContrastsMissing$new(
      lfqdata = transformed,
      contrasts = GRP2$pop$Contrasts,
      modelName = "Imputed_Mean")

    conMI <- prolfqua::ContrastsModerated$new(
      mC)

    res <- prolfqua::merge_contrasts_results(conrM, conMI)
    GRP2$RES$contrMerged <- res$merged
    GRP2$RES$contrMore <- res$more

  } else {
    GRP2$RES$contrMerged <- conrM
    GRP2$RES$contrMore <- NULL
  }


  datax <- GRP2$RES$contrMerged$get_contrasts()
  datax <- dplyr::inner_join(GRP2$RES$rowAnnot$row_annot, datax, multiple = "all")
  GRP2$RES$contrastsData  <- datax

  datax <- datax |>
    dplyr::filter(.data$FDR < GRP2$pop$FDRthreshold & abs(.data$diff) > GRP2$pop$Diffthreshold )
  GRP2$RES$contrastsData_signif <- datax
  return(GRP2)
}




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


#' read yaml file
#' @export
#' @return list with applications parameters
read_yaml <- function(ymlfile, application = "FragPipeTMT" ) {
  yml = yaml::read_yaml(ymlfile)

  WORKUNITID = yml$job_configuration$workunit_id
  PROJECTID = yml$job_configuration$project_id
  ORDERID = yml$job_configuration$order_id
  ORDERID <- if (is.null(ORDERID)) { PROJECTID }else{ ORDERID }
  ZIPDIR = paste0("C",ORDERID,"WU",WORKUNITID)

  GRP2 <- list()
  GRP2$project_spec <- list()
  GRP2$project_spec$project_Id <- PROJECTID


  GRP2$project_spec$project_name <- "" # workunit name in the future.
  GRP2$project_spec$order_Id <- ORDERID

  GRP2$project_spec$workunit_Id <- WORKUNITID

  idxzip <- grep("[0-9]{7,7}.zip",yml$application$input[[1]])

  GRP2$project_spec$input_Id <- yml$job_configuration$input[[1]][[idxzip]]$resource_id
  GRP2$project_spec$input_URL <- yml$job_configuration$input[[1]][[idxzip]]$resource_url

  #at least 2 peptides per protein
  GRP2$pop <- list()
  GRP2$pop$transform <- yml$application$parameters$`3|Normalization`

  GRP2$pop$aggregate <- "medpolish"
  GRP2$pop$Diffthreshold <- as.numeric(yml$application$parameters$`4|Difference_threshold`)
  GRP2$pop$FDRthreshold <- as.numeric(yml$application$parameters$`5|FDR_threshold`)

  GRP2$pop$removeCon <- if (yml$application$parameters$`6|remConDec` == "true") { TRUE } else { FALSE }
  GRP2$pop$removeDecoys <- if (yml$application$parameters$`6|remConDec` == "true") { TRUE } else { FALSE }

  GRP2$pop$revpattern <- yml$application$parameters$`7|REVpattern`
  GRP2$pop$contpattern <- yml$application$parameters$`8|CONpattern`

  GRP2$Software <- application
  GRP2$zipdir <- ZIPDIR
  return(GRP2)
}


#' create GRP2 configuration.
#' Use this function if there is no Yaml Input.
#' @param patternDecoys default "^REV_"
#' @param patternContaminants default "^zz_"
#' @export
#' @examples
#' DEAconfig <- make_DEA_config()
#' DEAconfig
make_DEA_config <- function(
    ZIPDIR  = ".",
    PROJECTID = "",
    ORDERID ="",
    WORKUNITID ="",
    Normalization = c("none","vsn", "quantile", "robscale"),
    aggregation = c("medpolish" , "top3", "lmrob"),
    Diffthreshold = 1,
    FDRthreshold = 0.1,
    removeContaminants = FALSE,
    removeDecoys = FALSE,
    patternDecoys = "^REV_",
    patternContaminants = "^zz",
    application = "FragPipeTMT",
    nrPeptides = 2){
  warning("DEPRECATED")
  Normalization <- match.arg(Normalization)
  aggregation <- match.arg(aggregation)


  GRP2 <- list()
  GRP2$project_spec <- list()
  GRP2$project_spec$project_Id <- PROJECTID
  GRP2$project_spec$project_name <- "" # workunit name in the future.
  GRP2$project_spec$order_Id <- ORDERID
  GRP2$project_spec$workunit_Id <- WORKUNITID
  GRP2$project_spec$input_URL <- NULL

  #at least 2 peptides per protein
  GRP2$pop <- list()
  GRP2$pop$transform <- Normalization

  GRP2$pop$aggregate <- aggregation
  GRP2$pop$Diffthreshold <- Diffthreshold
  GRP2$pop$FDRthreshold <- FDRthreshold

  GRP2$pop$removeCon <- removeContaminants
  GRP2$pop$removeDecoys <- removeDecoys

  GRP2$pop$revpattern <- patternDecoys
  GRP2$pop$contpattern <- patternContaminants
  GRP2$pop$nr_peptdes <- nrPeptides
  GRP2$Software <- application
  GRP2$zipdir <- ZIPDIR

  return(GRP2)
}
