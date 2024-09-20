# simplified version of

#' make DEA
#' @export
#'
#' @examples
#' # example code
#'
#' pep <- prolfqua::sim_lfq_data_protein_config()
#' pep <- prolfqua::LFQData$new(pep$data, pep$config)
#' pA <- data.frame(protein_Id = unique(pep$data$protein_Id))
#' pA <- pA |> dplyr::mutate(fasta.annot = paste0(pA$protein_Id, "_description"))
#' pA <- prolfquapp::ProteinAnnotation$new(pep,row_annot = pA ,description = "fasta.annot")
#' GRP2 <- prolfquapp::make_DEA_config_R6()
#'
#' pep$factors()
#' GRP2$pop <- list(Contrasts = c("AVsC" = "group_A - group_Ctrl", BVsC = "group_B - group_Ctrl"))
#' grp <- make_DEA_report2(pep, pA, GRP2)
#'
make_DEA_report2 <- function(lfqdata,
                             protAnnot,
                             GRP2
) {
  ### Do some type of data normalization (or do not)

  allProt <- nrow( protAnnot$row_annot )
  GRP2$RES <- list()
  GRP2$RES$Summary <- data.frame(
    totalNrOfProteins = allProt,
    percentOfContaminants = round(protAnnot$annotate_contaminants()/allProt * 100 , digits = 2),
    percentOfFalsePositives  = round(protAnnot$annotate_decoys()/allProt * 100 , digits = 2),
    NrOfProteinsNoDecoys = protAnnot$nr_clean()
  )
  GRP2$RES$rowAnnot <- protAnnot

  if (GRP2$processing_options$remove_cont || GRP2$processing_options$remove_decoys) {
    lfqdata <- lfqdata$get_subset(protAnnot$clean(
      contaminants = GRP2$processing_options$remove_cont,
      decoys = GRP2$processing_options$remove_decoys))
    logger::log_info(
      paste0("removing contaminants and reverse sequences with patterns:",
             GRP2$processing_options$pattern_contaminants,
             GRP2$processing_options$pattern_decoys ))
  }

  transformed <- prolfquapp::transform_lfqdata(
    lfqdata,
    method = GRP2$processing_options$transform,
    internal = GRP2$pop$internal
  )

  lfqdata$rename_response("protein_abundance")
  transformed$rename_response("normalized_protein_abundance")

  GRP2$RES$lfqData <- lfqdata
  GRP2$RES$transformedlfqData <- transformed

  ################## Run Modelling ###############
  # remove control column from factors.
  factors <- transformed$config$table$factor_keys_depth()[
    !grepl("^control", transformed$config$table$factor_keys_depth() , ignore.case = TRUE)
  ]

  # model with or without interactions
  if (is.null(GRP2$processing_options$interaction) || !GRP2$processing_options$interaction ) {
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

  if (is.null(GRP2$processing_options$missing) || GRP2$processing_options$missing ) {
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
    dplyr::filter(.data$FDR < GRP2$processing_options$FDR_threshold &
                    abs(.data$diff) > GRP2$processing_options$diff_threshold )
  GRP2$RES$contrastsData_signif <- datax
  return(GRP2)
}


#' will replace generate_DEA_reports
#' @export
generate_DEA_reports2 <- function(lfqdata, GRP2, prot_annot, Contrasts) {
  # Compute all possible 2 Grps to avoid specifying reference.
  GRP2$pop$Contrasts <- Contrasts
  logger::log_info("CONTRAST : ", paste( GRP2$pop$Contrasts, collapse = " "))
  lfqwork <- lfqdata$get_copy()
  # lfqwork$data <- lfqdata$data |> dplyr::filter(.data$Group_ %in% levels$Group_)
  grp2 <- make_DEA_report2(
    lfqwork,
    prot_annot,
    GRP2)
  return(grp2)
}



prep_result_list <- function(GRP2){
  rd <- GRP2$RES$lfqData
  tr <- GRP2$RES$transformedlfqData
  ra <- GRP2$RES$rowAnnot
  contrasts <- data.frame(
    contrast_name = names(GRP2$pop$Contrasts),
    contrast = GRP2$pop$Contrasts
  )

  wideraw <- dplyr::inner_join(ra$row_annot, rd$to_wide()$data, multiple = "all")
  widetr <- dplyr::inner_join(ra$row_annot , tr$to_wide()$data, multiple = "all")

  ctr <- dplyr::inner_join(ra$row_annot , GRP2$RES$contrMerged$get_contrasts(), multiple = "all")
  ctr_wide <- dplyr::inner_join(ra$row_annot , GRP2$RES$contrMerged$to_wide(), multiple = "all")

  resultList <- list()

  resultList$annotation <- dplyr::inner_join(
    rd$factors(),
    rd$get_Summariser()$hierarchy_counts_sample(),
    by = rd$config$table$sampleName,
    multiple = "all")

  resultList$normalized_abundances = dplyr::inner_join(ra$row_annot, tr$data, multiple = "all")
  resultList$raw_abundances_matrix = wideraw
  resultList$normalized_abundances_matrix = widetr
  resultList$diff_exp_analysis = ctr
  resultList$diff_exp_analysis_wide = ctr_wide
  resultList$formula = data.frame(formula = GRP2$RES$formula)
  resultList$summary = GRP2$RES$Summary
  resultList$missing_information = prolfqua::UpSet_interaction_missing_stats(rd$data, rd$config, tr = 1)$data
  resultList$contrasts <- contrasts

  # add protein statistics
  st <- GRP2$RES$transformedlfqData$get_Stats()
  resultList$stats_normalized <- st$stats()
  resultList$stats_normalized_wide <- st$stats_wide()

  st <- GRP2$RES$lfqData$get_Stats()
  resultList$stats_raw <- st$stats()
  resultList$stats_raw_wide <- st$stats_wide()
  return(resultList)
}


.write_ORA <- function(fg, outpath, workunit_id, id_column = "IDcolumn") {
  fg <- fg |> dplyr::mutate(updown = paste0(contrast, ifelse(diff > 0 , "_up", "_down")))
  ora_sig <- split(fg[[id_column]], fg$updown)

  for (i in names(ora_sig)) {
    ff <- file.path(outpath, paste0("ORA_",i,"_WU",workunit_id,".txt" ))
    logger::log_info("Writing File ", ff)
    write.table(ora_sig[[i]],file = ff, col.names = FALSE,
                row.names = FALSE, quote = FALSE)
  }
}


write_result_list <- function(outpath,
                              GRP2,
                              resultList,
                              xlsxname,
                              GSEA = TRUE,
                              ORA = TRUE) {
  workunit_id <- GRP2$project_spec$workunit_Id
  dir.create(outpath)
  id_column = GRP2$RES$rowAnnot$cleaned_ids

  if (ORA) {
    bkg <- GRP2$RES$lfqData$hierarchy()
    ff <- file.path(outpath ,paste0("ORA_background_WU",workunit_id,".txt"))
    write.table(bkg,file = ff, col.names = FALSE,
                row.names = FALSE, quote = FALSE)
    .write_ORA(GRP2$RES$contrastsData_signif, outpath, workunit_id,id_column = id_column)
  }
  if (GSEA) {
    fg <- GRP2$RES$contrastsData
    gsea <- dplyr::select(fg , c("contrast",id_column , "statistic")) |>
      dplyr::arrange( .data$statistic )
    gsea <- split(dplyr::select( gsea, c(id_column, "statistic")), gsea$contrast)

    for (i in names(gsea)) {
      ff <- file.path(outpath, paste0("GSEA_",i,"_WU",workunit_id,".rnk" ))
      logger::log_info("Writing File ", ff)
      write.table(na.omit(gsea[[i]]),file = ff, col.names = FALSE,
                  row.names = FALSE, quote = FALSE, sep = "\t")
    }
  }
  if (nrow(resultList$normalized_abundances) > 1048575) {
    resultList$normalized_abundances <- NULL
  }
  writexl::write_xlsx(resultList, path = file.path(outpath, paste0(xlsxname, ".xlsx")))
}

#' Write differential expression analysis results
#'
#' @param GRP2 return value of \code{\link{make_DEA_report}}
#' @param outpath path to place output
#' @param xlsxname file name for xlsx
#' @export
#' @family workflow
#'
write_DEA <- function(GRP2, outpath, ORA = TRUE, GSEA = TRUE, xlsxname = "AnalysisResults", write = TRUE){
  resultList <- prep_result_list(GRP2)
  if (write) {
    write_result_list(outpath, GRP2, resultList, xlsxname = xlsxname, ORA = ORA, GSEA = GSEA)
  }
  return(resultList)
}




#' Generate differential expression analysis reports
#'
#' Writes results of DEA see \code{\link{generate_DEA_reports}}
#' @export
#'
write_DEA_all <- function(
    grp2,
    name = "",
    boxplot = TRUE ,
    render = TRUE,
    ORA = TRUE,
    GSEA = TRUE,
    markdown ="_Grp2Analysis.Rmd"){
  dir.create(grp2$get_zipdir())
  name <- if (name != "") {paste0(name,"_")}
  fname <- paste0("DE_",  name, "WU", grp2$project_spec$workunit_Id )
  qcname <- paste0("QC_", name, "WU", grp2$project_spec$workunit_Id )
  outpath <- grp2$get_result_dir()

  prolfquapp::write_DEA(grp2, outpath = outpath, xlsxname = fname, ORA = ORA, GSEA = GSEA)

  if (render) {
    prolfquapp::render_DEA(grp2, outpath = outpath, htmlname = fname, markdown = markdown)
    prolfquapp::render_DEA(grp2, outpath = outpath, htmlname = qcname, markdown = "_DiffExpQC.Rmd")
  }

  bb <- grp2$RES$transformedlfqData
  grsizes <- bb$factors() |>
    dplyr::group_by(dplyr::across(bb$config$table$factor_keys_depth())) |>
    dplyr::summarize(n = n()) |>
    dplyr::pull( n )
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

#' convert tibble to data.frame with rownames
#' @param .data a tibble or data.frame
#' @param var name of the column with new row.names
#' @return a data.frame with rownames
#' @export
#' @examples
#' ind <- tibble::tibble(a = 1:3, rowname = letters[1:3])
#' column_to_rownames(ind)
column_to_rownames <- function(.data, var = "rowname", sep = "~lfq~"){
  #.data[,var] |> tidyr::unite("id", everything, sep = sep) |> pull("id")
  res <- as.data.frame(.data)
  rownames(res) <- .data[,var] |> tidyr::unite("id", tidyselect::everything(), sep = sep) |> dplyr::pull("id")
  return(res)
}

strip_rownames <- function(.data, strip="~lfq~light$"){
  newrnames <- gsub(strip, "", rownames(.data))
  rownames(.data) <- newrnames
  return(.data)
}


#' build bfabric urls
#' @export
#' @examples
#'
#'
#' ps <- ProjectSpec$new()
#' ps$project_Id <- 32258
#' ps$order_Id <- 34628
#' ps$workunit_Id <- 302212
#' bfabric_url_builder(ps)
#'
#' ps <- ProjectSpec$new()
#' ps$order_Id <- 34628
#' ps$workunit_Id <- 302212
#' bfabric_url_builder(ps)
#'
bfabric_url_builder <- function(project_spec){

  orderURL <- NULL
  workunitURL <- NULL
  projectURL <- NULL
  orderID <- as.numeric(project_spec$order_Id)
  if ((length(orderID) > 0) && !is.na(orderID)) {
    orderURL <- paste0("https://fgcz-bfabric.uzh.ch/bfabric/order/show.html?id=", orderID, "&tab=details")
  }
  workunitID <- as.numeric(project_spec$workunit_Id)
  if ((length(workunitID) > 0) && !is.na(workunitID)) {
    workunitURL <- paste0("https://fgcz-bfabric.uzh.ch/bfabric/workunit/show.html?id=", orderID, "&tab=details")
  }
  projectID <- as.numeric(project_spec$project_Id)
  if ((length(projectID) > 0) && !is.na(projectID)) {
    projectURL <- paste0("https://fgcz-bfabric.uzh.ch/bfabric/project/show.html?id=", orderID, "&tab=details")
  }

  return(list(orderURL = orderURL, projectURL = projectURL, workunitURL = workunitURL))
}

#' Convert prolfqua differential expression analysis results to SummarizedExperiment
#'
#' @param GRP2 return value of \code{\link{make_DEA_report}}
#' @return SummarizedExperiment
#' @export
#' @family workflow
#'
#'
#'
make_SummarizedExperiment <- function(GRP2,
                                      colname = NULL,
                                      rowname = NULL,
                                      strip="~lfq~light",
                                      .url_builder = bfabric_url_builder){
  if (is.null(colname)) {
    colname <- GRP2$RES$lfqData$config$table$sampleName
  }
  if (is.null(rowname)) {
    rowname <- GRP2$RES$lfqData$config$table$hierarchyKeys()
  }
  resTables <- prep_result_list(GRP2)
  matTr <- GRP2$RES$transformedlfqData$to_wide(as.matrix = TRUE)
  matRaw <- GRP2$RES$transformedlfqData$to_wide(as.matrix = TRUE)

  mat.raw <- strip_rownames(matRaw$data, strip)
  mat.trans <- strip_rownames(matTr$data, strip)
  col.data <- column_to_rownames(matRaw$annotation, var = colname)
  col.data <- col.data[colnames(mat.raw),]
  x <- SummarizedExperiment::SummarizedExperiment(
    assays = list(rawData = mat.raw, transformedData = mat.trans),
    colData = col.data, metadata = list(bfabric_urls = .url_builder(GRP2$project_spec), contrasts = resTables$contrasts, formula = resTables$formula )
  )

  diffbyContrast <- split(resTables$diff_exp_analysis, resTables$diff_exp_analysis$contrast)
  for (i in names(diffbyContrast)) {
    row.data <- column_to_rownames(diffbyContrast[[i]], var = rowname)
    row.data <- row.data[rownames(mat.raw),]
    if ("nr_compounds" %in% colnames(row.data)) {
      row.data <- dplyr::rename(nr_peptides = "nr_compounds")
    }
    SummarizedExperiment::rowData(x)[[paste0("constrast_",i)]] <- row.data
  }

  SummarizedExperiment::rowData(x)[["stats_normalized_wide"]] <- column_to_rownames(resTables$stats_normalized_wide, var = rowname)[rownames(mat.raw),]
  SummarizedExperiment::rowData(x)[["stats_raw_wide"]] <- column_to_rownames(resTables$stats_raw_wide, var = rowname)[rownames(mat.raw),]
  return(x)
}










