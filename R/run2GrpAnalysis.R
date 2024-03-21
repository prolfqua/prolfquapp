#' transform lfq data using robscale, vsn or non,
#'
#' It will also run internal but then robscale must be used.
#' @param lfqdata \code{\link{LFQData}}
#' @param method normalization method to use
#' @param internal a data.frame with protein ids to be used for internal calibration, column name must be the same as
#' @export
#' @examples
#' istar <- prolfqua_data('data_ionstar')$filtered()
#' config <- prolfqua:::old2new(istar$config)
#' tmp <- prolfqua::LFQData$new(istar$data, config)
#' d <- istar$d
#' internal <- dplyr::filter(d, protein_Id %in% sample(unique(d$protein_Id), 3 )) |>
#'   dplyr::select(all_of(tmp$config$table$hierarchy_keys()[1])) |> dplyr::distinct()
#' tmp2 <- transform_lfqdata(tmp, internal = internal)
#' tmp2 <- transform_lfqdata(tmp)
#'
transform_lfqdata <- function(lfqdata, method = c("robscale", "vsn", "none"), internal = NULL) {
  method <- match.arg(method)
  lt <- lfqdata$get_Transformer()
  if (method == "robscale") {
    transformed <- lt$log2()$robscale()$lfq
    if (!is.null(internal)) {
      logger::log_info("Attempt of internal calibration.")
      subset <- lfqdata$get_subset(internal)$get_Transformer()$log2()$lfq
      transformed <- lfqdata$get_Transformer()$log2()$robscale_subset(subset)$lfq
    }
  } else if (method == "vsn") {
    transformed <- lt$intensity_matrix( .func = vsn::justvsn)$lfq
  } else if (method == "none") {
    transformed <- lt$log2()$lfq
  } else {
    logger::log_warn("no such transformaton : {method}")
    return(NULL)
  }
  logger::log_info("data transformed : {method}.")
  return(transformed)
}

#' transform lfq data with x^2 - apply if non log data is needed
#' @param lfqTrans transformed LFQData
#' @export
exp2 <- function( lfqTrans ){
  if(!lfqTrans$config$table$is_response_transformed) {
    warning("Data not transformed.")
  }
  tr <- lfqTrans$get_Transformer()
  .exp2 <- function(x){
    2^x
  }
  tr$intensity_array(.exp2, force = TRUE)
  tr$lfq$config$table$is_response_transformed <- FALSE
  lfqdataProt <- tr$lfq
  return(lfqdataProt)
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
#' @examples
#'
#' library(prolfquapp)
#' library(prolfqua)
#' istar <- prolfqua_data('data_ionstar')$filtered()
#' data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' data$Description <-"AAAAA"
#' protein_annot <- data |> dplyr::select(protein_Id, description = Description) |> dplyr::distinct()
#' protein_annot2 <- prolfqua::get_UniprotID_from_fasta_header(protein_annot, idcolumn = "protein_Id")
#' protein_annot <- protein_annot2 |> dplyr::rename(IDcolumn = UniprotID)
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
#' GRP2$pop$Contrasts <- c(b_vs_a = "dilution.b - dilution.a")
#'
#' GRP2$Software <- "MaxQuant"
#'
#' data <- dplyr::filter(data, dilution. == "a" |  dilution. == "b")
#' atab <- AnalysisTableAnnotation$new()
#' atab$ident_qValue = "pep"
#' atab$fileName = "raw.file"
#' atab$hierarchy["protein_Id"] = "protein_Id"
#' atab$hierarchy["peptide_Id"] = "peptide_Id"
#' atab$factors["dilution."] = "dilution."
#' atab$set_response("peptide.intensity")
#' atab$isotopeLabel = "isotope"
#' config <- prolfqua::AnalysisConfiguration$new(atab)
#'
#' lfqdata <- prolfqua::LFQData$new(data, config)
#' lfqdata$remove_small_intensities()
#' aggregator <- lfqdata$get_Aggregator()
#'
#' aggregator$sum_topN()
#' lfqdata <- aggregator$lfq_agg
#'
#' GRP2$pop$nrPeptides <- 2
#'
#' GRP2$pop$Contrasts <- c("avsb" = "dilution.a - dilution.b")
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
  warning("DEPRECATED")
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

#' convert tibble to data.frame with rownames
#' @param .data a tibble or data.frame
#' @param var name of the column with new row.names
#' @return a data.frame with rownames
#' @export
#' @examples
#' ind <- tibble::tibble(a = 1:3, rowname = letters[1:3])
#' column_to_rownames(ind)
column_to_rownames <- function(.data, var = "rowname"){
  res <- as.data.frame(.data)
  rownames(res) <- .data[[var]]
  return(res)
}

strip_rownames <- function(.data, strip="~lfq~light$"){
  newrnames <- gsub(strip, "", rownames(.data))
  rownames(.data) <- newrnames
  return(.data)
}


#' Convert prolfqua differential expression analysis results to SummarizedExperiment
#'
#' @rdname make_DEA_report
#' @param GRP2 return value of \code{\link{make_DEA_report}}
#' @return SummarizedExperiment
#' @export
#' @family workflow
#'
make_SummarizedExperiment <- function(GRP2, colname = NULL, rowname = NULL, strip="~lfq~light"){
  if (is.null(colname)) {
    colname <- GRP2$RES$lfqData$config$table$sampleName
  }
  if (is.null(rowname)) {
    rowname <- GRP2$RES$lfqData$config$table$hierarchyKeys()
  }
  resTables <- write_DEA(GRP2,".", write = FALSE)
  matTr <- GRP2$RES$transformedlfqData$to_wide(as.matrix = TRUE)
  matRaw <- GRP2$RES$transformedlfqData$to_wide(as.matrix = TRUE)

  mat.raw <- strip_rownames(matRaw$data, strip)
  mat.trans <- strip_rownames(matTr$data, strip)
  col.data <- column_to_rownames(matRaw$annotation, var = colname)
  col.data <- col.data[colnames(mat.raw),]
  x <- SummarizedExperiment::SummarizedExperiment(
    assays = list(rawData = mat.raw, transformedData = mat.trans),
    colData = col.data, metadata = as.list(resTables$formula)
  )

  diffbyContrast <- split(resTables$diff_exp_analysis, resTables$diff_exp_analysis$contrast)
  for (i in names(diffbyContrast)) {
    row.data <- column_to_rownames(diffbyContrast[[i]], var = rowname)
    row.data <- row.data[rownames(mat.raw),]

    SummarizedExperiment::rowData(x)[[paste0("constrast_",i)]] <- row.data
  }

  SummarizedExperiment::rowData(x)[["stats_normalized_wide"]] <- column_to_rownames(resTables$stats_normalized_wide, var = rowname)[rownames(mat.raw),]
  SummarizedExperiment::rowData(x)[["stats_raw_wide"]] <- column_to_rownames(resTables$stats_raw_wide, var = rowname)[rownames(mat.raw),]
  return(x)
}


#' Write differential expression analysis results
#'
#' @rdname make_DEA_report
#' @param GRP2 return value of \code{\link{make_DEA_report}}
#' @param outpath path to place output
#' @param xlsxname file name for xlsx
#' @export
#' @family workflow
#'
write_DEA <- function(GRP2, outpath, xlsxname = "AnalysisResults", write = TRUE){
  if (write) {dir.create(outpath)}
  rd <- GRP2$RES$lfqData
  tr <- GRP2$RES$transformedlfqData
  ra <- GRP2$RES$rowAnnot
  formula <- data.frame(
    formula = GRP2$RES$formula,
    contrast_name = names(GRP2$pop$Contrasts),
    contrast = GRP2$pop$Contrasts)

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

  resultList$normalized_abundances = dplyr::inner_join(ra$row_annot, tr$data,multiple = "all")
  resultList$raw_abundances_matrix = wideraw
  resultList$normalized_abundances_matrix = widetr
  resultList$diff_exp_analysis = ctr
  resultList$diff_exp_analysis_wide = ctr_wide
  resultList$formula = formula
  resultList$summary = GRP2$RES$Summary
  resultList$missing_information = prolfqua::UpSet_interaction_missing_stats(rd$data, rd$config, tr = 1)$data

  # add protein statistics
  st <- GRP2$RES$transformedlfqData$get_Stats()
  resultList$stats_normalized <- st$stats()
  resultList$stats_normalized_wide <- st$stats_wide()

  st <- GRP2$RES$lfqData$get_Stats()
  resultList$stats_raw <- st$stats()
  resultList$stats_raw_wide <- st$stats_wide()

  if (write) {

    bkg <- GRP2$RES$rowAnnot$row_annot$IDcolumn
    ff <- file.path(outpath ,"ORA_background.txt")
    write.table(bkg,file = ff, col.names = FALSE,
                row.names = FALSE, quote = FALSE)

    fg <- GRP2$RES$contrastsData_signif
    ora_sig <- split(fg$IDcolumn, fg$contrast)

    for (i in names(ora_sig)) {
      ff <- file.path(outpath, paste0("Ora_",i,".txt" ))
      logger::log_info("Writing File ", ff)
      write.table(ora_sig[[i]],file = ff, col.names = FALSE,
                  row.names = FALSE, quote = FALSE)
    }

    fg <- GRP2$RES$contrastsData
    gsea <- fg |> dplyr::select( contrast, IDcolumn, statistic) |> dplyr::arrange( statistic )
    gsea <- split(dplyr::select( gsea, IDcolumn, statistic ), gsea$contrast)

    for (i in names(gsea)) {
      ff <- file.path(outpath, paste0("GSEA_",i,".rnk" ))
      logger::log_info("Writing File ", ff)
      write.table(na.omit(gsea[[i]]),file = ff, col.names = FALSE,
                  row.names = FALSE, quote = FALSE, sep = "\t")
    }
    if (nrow(resultList$normalized_abundances) > 1048575) {
      resultList$normalized_abundances <- NULL
    }
    writexl::write_xlsx(resultList, path = file.path(outpath, paste0(xlsxname, ".xlsx")))
  }
  return(resultList)
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


#' Generate differential expression analysis reports
#'
#' Writes results of DEA see \code{\link{generate_DEA_reports}}
#' @export
#'
write_DEA_all <- function(grp2, name, ZIPDIR, boxplot = TRUE , render = TRUE){
  dir.create(GRP2$zipdir)
  fname <- paste0("DE_", name)
  qcname <- paste0("QC_", name)
  outpath <- file.path( ZIPDIR, fname)
  logger::log_info("writing into : ", outpath, " <<<<")
  prolfquapp::write_DEA(grp2, outpath = outpath, xlsxname = fname)

  if (render) {
    prolfquapp::render_DEA(grp2, outpath = outpath, htmlname = fname)
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
}


