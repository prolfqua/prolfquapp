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
#' @examples
#'
#' istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 100)
#' data <- istar$data
#' data$Description <-"AAAAA"
#' protein_annot <- data |> dplyr::select(protein_Id, description = Description) |> dplyr::distinct()
#' GRP2 <- list()
#' GRP2$Bfabric <- list()
#' GRP2$Bfabric$projectID <- "3765"
#' GRP2$Bfabric$projectName <- "Order_26863"
#' GRP2$Bfabric$orderID <- "3765"
#'
#' GRP2$Bfabric$workunitID <- "2057368.zip"
#' GRP2$Bfabric$inputID <- "2057368.zip"
#' GRP2$Bfabric$inputURL <- "https://www.fgcz.ch"
#'
#' #at least 2 peptides per protein
#' GRP2$pop <- list()
#' GRP2$pop$transform <- "vsn"
#' GRP2$pop$aggregate <- "medpolish"
#' GRP2$pop$Diffthreshold <- 0.5
#' GRP2$pop$FDRthreshold <- 0.25
#' GRP2$pop$Contrasts <- c(b_vs_a = "group_A - group_Ctrl")
#' GRP2$Software <- "MaxQuant"
#'
#'
#' lfqdata <- prolfqua::LFQData$new(data, istar$config)
#' lfqdata$remove_small_intensities()
#' aggregator <- lfqdata$get_Aggregator()
#'
#' aggregator$sum_topN()
#' lfqdata <- aggregator$lfq_agg
#'
#' GRP2$pop$nrPeptides <- 2
#'
#' GRP2$pop$revpattern = "^REV_"
#' GRP2$pop$contpattern = "^zz|^CON__"
#'
#' GRP2$pop$removeCon = TRUE
#' GRP2$pop$removeDecoys = TRUE
#'
#'
#' protAnnot <- prolfqua::ProteinAnnotation$new(
#'     lfqdata,
#'     row_annot = protein_annot)
#'
#' grp <- make_DEA_report(lfqdata, protAnnot, GRP2)
#' st <- grp$RES$transformedlfqData$get_Stats()
#' bb <- st$stats()
#' sr <- grp$RES$transformedlfqData$get_Summariser()
#' int <- sr$interaction_missing_stats()
#' res <- write_DEA(grp,".", write=FALSE)
#' se <- make_SummarizedExperiment(grp)
#' \dontrun{
#' render_DEA(grp, ".")
#' render_DEA(grp, "." ,word = TRUE)
#'
#' }
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
    percentOfContaminants = round(protAnnot$annotate_contaminants(GRP2$pop$contpattern)/allProt * 100 , digits = 2),
    percentOfFalsePositives  = round(protAnnot$annotate_decoys(GRP2$pop$revpattern)/allProt * 100 , digits = 2),
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


