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


#' Generate differential expression analysis reports
#'
#' Writes results of DEA see \code{\link{generate_DEA_reports}}
#' @export
#'
write_DEA_all <- function(
    grp2,
    ZIPDIR = grp2$zipdir,
    name = "Groups_vs_Controls",
    boxplot = TRUE ,
    render = TRUE,
    markdown ="_Grp2Analysis.Rmd"){
  dir.create(GRP2$zipdir)
  fname <- paste0("DE_",  name, "_WU", grp2$project_spec$workunit_Id )
  qcname <- paste0("QC_", name, "_WU", grp2$project_spec$workunit_Id )
  outpath <- file.path( ZIPDIR, paste0("Results_DEA_WU", grp2$project_spec$workunit_Id))
  logger::log_info("writing into : ", outpath, " <<<<")

  prolfquapp::write_DEA(grp2, outpath = outpath, xlsxname = fname)

  if (render) {
    prolfquapp::render_DEA(grp2, outpath = outpath, htmlname = fname, markdown = markdown)
    prolfquapp::render_DEA(grp2, outpath = outpath, htmlname = qcname, markdown = "_DiffExpQC.Rmd")
  }

  bb <- grp2$RES$transformedlfqData
  grsizes <- bb$factors() |>
    dplyr::group_by(dplyr::across(bb$config$table$factor_keys_depth())) |>
    dplyr::summarize(n = n()) |>
    dplyr::pull(n)
  if (boxplot) {
    if (sum(!grepl("^control",bb$config$table$factor_keys(), ignore.case = TRUE))  > 1 &
        all(grsizes == 1)
    ) {
      prolfquapp::writeLinesPaired(bb, outpath)
    } else {
      pl <- bb$get_Plotter()
      pl$write_boxplots(outpath)
    }
  }
  return(outpath)
}


#' Render DEA analysis report
#' @rdname make_DEA_report
#' @param GRP2 return value of \code{\link{make_DEA_report}}
#' @param outpath path to place output
#' @param htmlname name for html file
#' @param word default FALSE, if true create word document.s
#' @param markdown which file to render
#' @export
#' @family workflow
render_DEA <- function(GRP2,
                       outpath,
                       htmlname="Result2Grp",
                       word = FALSE,
                       markdown = "_Grp2Analysis.Rmd"){
  dir.create(outpath)

  rmarkdown::render(
    markdown,
    params = list(grp = GRP2) ,
    output_format = if (word) {
      bookdown::word_document2(toc = TRUE, toc_float = TRUE) } else {
        bookdown::html_document2(toc = TRUE, toc_float = TRUE)
      }
  )
  fname <- paste0(tools::file_path_sans_ext(markdown), if (word) {".docx"} else {".html"})
  if (file.copy(fname, file.path(outpath, paste0(htmlname,if (word) {".docx"} else {".html"})), overwrite = TRUE)) {
    file.remove(fname)
  }
}


