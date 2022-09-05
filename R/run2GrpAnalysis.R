#' Create 2 grp report in html and write data to xlsx table
#'
#' For use examples see run_scripts directory
#' @rdname make2grpReport
#' @param startdata table in long format
#' @param atable AnalysisTableAnnotation annotate startdata table
#' @param GRP2 list with named arguments i.e. Contrasts, projectID, projectName, workunitID, nrPeptides, Diffthreshold, FDRthreshold
#' @param protein_annot column with protein description e.g. (fasta header)
#' @export
#' @family workflow
#' @examples
#' library(prolfquapp)
#' library(prolfqua)
#' istar <- prolfqua_data('data_ionstar')$filtered()
#' data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' data$Description <-"AAAAA"
#' protein_annot <- data |> dplyr::select(protein_Id, Description) |> dplyr::distinct()
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
#' atab$setWorkIntensity("peptide.intensity")
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
#' GRP2$pop$remove = TRUE
#' grp <- make2grpReport(lfqdata, protein_annot, GRP2)
#' st <- grp$RES$transformedlfqData$get_Stats()
#' bb <- st$stats()
#' sr <- grp$RES$transformedlfqData$get_Summariser()
#' int <- sr$interaction_missing_stats()
#'
#' \dontrun{
#'
#' render_2GRP(grp, ".")
#' render_2GRP(grp, "." ,word = TRUE)
#' write_2GRP(grp,".")
#' }
make2grpReport <- function(lfqdata,
                           prot_annot,
                           GRP2
) {



  ### Do some type of data normalization (or do not)
  lt <- lfqdata$get_Transformer()
  if (GRP2$pop$transform == "robscale") {
    transformed <- lt$log2()$robscale()$lfq
  } else if (GRP2$pop$transform == "vsn") {
    transformed <- lt$intensity_matrix( .func = vsn::justvsn)$lfq
  } else if (GRP2$pop$transform == "none") {
    transformed <- lt$log2()$lfq
  } else {
    logger::log_warn("no such transformaton : {GRP2$pop$transform}")
  }
  logger::log_info("data transformed : {GRP2$pop$transform}.")

  ## count contaminants.
  protAnnot <- prolfqua::RowAnnotProtein$new(
    transformed,
    row_annot = prot_annot)

  allProt <- nrow( protAnnot$row_annot )
  GRP2$RES <- list()
  GRP2$RES$Summary <- data.frame(
    totalNrOfProteins = allProt,
    percentOfContaminants = round(protAnnot$annotateCON(GRP2$pop$contpattern)/allProt * 100 , digits = 2),
    percentOfFalsePositives  = round(protAnnot$annotateREV(GRP2$pop$revpattern)/allProt * 100 , digits = 2),
    NrOfProteinsNoDecoys = protAnnot$nr_clean()
  )
  GRP2$RES$rowAnnot <- protAnnot

  if (GRP2$pop$remove) {
    lfqdata <- lfqdata$get_subset(protAnnot$clean())
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

  formula <- paste0(transformed$config$table$getWorkIntensity(), " ~ ",
                    paste(transformed$config$table$factorKeys()[!grepl("^control", transformed$config$table$factorKeys() , ignore.case=TRUE)], collapse = " + "))
  message("FORMULA :",  formula)
  GRP2$RES$formula <- formula
  formula_Condition <-  prolfqua::strategy_lm(formula)
  # specify model definition
  modelName  <- "Model"

  mod <- prolfqua::build_model(
    transformed,
    formula_Condition,
    subject_Id = transformed$config$table$hierarchyKeys() )

  logger::log_info("fitted model with formula : {formula}")
  GRP2$RES$models <- mod

  contr <- prolfqua::Contrasts$new(mod, GRP2$pop$Contrasts)
  conrM <- prolfqua::ContrastsModerated$new(
    contr,
    modelName = "Linear_Model_Moderated")
  #saveRDS(list(transformed = transformed ,
  #             Contrasts = GRP2$pop$Contrasts),
  #        file = "totest.Rds")
  mC <- prolfqua::ContrastsSimpleImpute$new(
    lfqdata = transformed,
    contrasts = GRP2$pop$Contrasts)

  conMI <- prolfqua::ContrastsModerated$new(
    mC,
    modelName = "Imputed_Mean")

  res <- prolfqua::addContrastResults(conrM, conMI)
  GRP2$RES$contrMerged <- res$merged
  GRP2$RES$contrMore <- res$more

  datax <- GRP2$RES$contrMerged$get_contrasts()
  datax <- dplyr::inner_join(GRP2$RES$rowAnnot$row_annot, datax)
  GRP2$RES$contrastsData  <- datax

  datax <- datax |>
    dplyr::filter(.data$FDR < GRP2$pop$FDRthreshold & abs(.data$diff) > GRP2$pop$Diffthreshold )
  GRP2$RES$contrastsData_signif <- datax

  return(GRP2)
}


#' write 2 grp results
#' @rdname make2grpReport
#' @param GRP2 return value of \code{\link{make2grpReport}}
#' @param outpath path to place output
#' @param xlsxname file name for xlsx
#' @export
#' @family workflow
#'
write_2GRP <- function(GRP2, outpath, xlsxname = "AnalysisResults"){
  dir.create(outpath)
  rd <- GRP2$RES$lfqData
  tr <- GRP2$RES$transformedlfqData
  ra <- GRP2$RES$rowAnnot
  formula <- data.frame(formula = GRP2$RES$formula, contrast_name = names(GRP2$pop$Contrasts), contrast = GRP2$pop$Contrasts)
  wideraw <- dplyr::inner_join(ra$row_annot, rd$to_wide()$data)
  widetr <- dplyr::inner_join(ra$row_annot , tr$to_wide()$data )
  ctr <- dplyr::inner_join(ra$row_annot , GRP2$RES$contrMerged$get_contrasts())
  resultList <- list()
  resultList$annotation = tr$to_wide()$annot
  resultList$normalized_abundances = dplyr::inner_join(ra$row_annot, rd$data)
  resultList$raw_abundances_matrix = wideraw
  resultList$normalized_abundances_matrix = widetr
  resultList$diff_exp_analysis = ctr
  resultList$formula = formula
  resultList$summary = GRP2$RES$Summary
  resultList$missing_information = prolfqua::UpSet_interaction_missing_stats(rd$data, rd$config, tr=1)

  # add protein statistics
  st <- GRP2$RES$transformedlfqData$get_Stats()
  resultList$protein_variances <- st$stats()

  bkg <- prolfqua::get_UniprotID_from_fasta_header(
    GRP2$RES$rowAnnot$row_annot,
    idcolumn = "protein_Id")$UniprotID
  ff <- file.path(outpath ,"ORA_background.txt")
  write.table(bkg,file = ff, col.names = FALSE,
              row.names = FALSE, quote=FALSE)

  fg <- prolfqua::get_UniprotID_from_fasta_header(
    GRP2$RES$contrastsData_signif,
  idcolumn = "protein_Id")
  ora_sig <- split(fg$UniprotID, fg$contrast)
  for(i in names(ora_sig)){
    ff <- file.path(outpath, paste0("Ora_",i,".txt" ))
    logger::log_info("Writing File ", ff)
    write.table(ora_sig[[i]],file = ff, col.names = FALSE,
               row.names = FALSE, quote=FALSE)
  }
  fg <- prolfqua::get_UniprotID_from_fasta_header(
    GRP2$RES$contrastsData,
    idcolumn = "protein_Id")
  gsea <- fg |> dplyr::select(contrast, UniprotID, statistic) |> dplyr::arrange(statistic)
  gsea <- split(dplyr::select(gsea,UniprotID, statistic ), gsea$contrast)
  for(i in names(gsea)){
    ff <- file.path(outpath, paste0("GSEA_",i,".rnk" ))
    logger::log_info("Writing File ", ff)
    write.table(na.omit(gsea[[i]]),file = ff, col.names = FALSE,
                row.names = FALSE, quote=FALSE, sep = "\t")
  }
  writexl::write_xlsx(resultList, path = file.path(outpath, paste0(xlsxname, ".xlsx")))
}

#' render 2GRP analysis report
#' @rdname make2grpReport
#' @param GRP2 return value of \code{\link{make2grpReport}}
#' @param outpath path to place output
#' @param htmlname name for html file
#' @param word default FALSE, if true create word document.
#' @param markdown which file to render
#' @export
#' @family workflow
render_2GRP <- function(GRP2,
                        outpath,
                        htmlname="Result2Grp",
                        word = FALSE,
                        markdown = "_Grp2Analysis.Rmd"){
  dir.create(outpath)

  rmarkdown::render(
    markdown,
    params = list(grp = GRP2) ,
    output_format = if(word){
      bookdown::word_document2(toc = TRUE, toc_float = TRUE) } else {
        bookdown::html_document2(toc = TRUE, toc_float = TRUE)
      }
  )
  fname <- paste0(tools::file_path_sans_ext(markdown), if(word) {".docx"} else {".html"})
  if (file.copy(fname, file.path(outpath, paste0(htmlname,if(word) {".docx"} else {".html"})), overwrite = TRUE)) {
    file.remove(fname)
  }
}

