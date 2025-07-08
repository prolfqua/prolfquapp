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
#' pA <- prolfquapp::ProteinAnnotation$new(pep, row_annot = pA, description = "fasta.annot")
#' GRP2 <- prolfquapp::make_DEA_config_R6()
#'
#' pep$factors()
#' GRP2$pop <- list(Contrasts = c("AVsC" = "group_A - group_Ctrl", BVsC = "group_B - group_Ctrl"))
#' grp <- make_DEA_report2(pep, pA, GRP2)
#'
make_DEA_report2 <- function(lfqdata,
                             protAnnot,
                             GRP2) {
  ### Do some type of data normalization (or do not)

  allProt <- nrow(protAnnot$row_annot)
  GRP2$RES <- list()
  GRP2$RES$Summary <- data.frame(
    totalNrOfProteins = allProt,
    percentOfContaminants = round(protAnnot$annotate_contaminants() / allProt * 100, digits = 2),
    percentOfFalsePositives = round(protAnnot$annotate_decoys() / allProt * 100, digits = 2),
    NrOfProteinsNoDecoys = protAnnot$nr_clean()
  )
  GRP2$RES$rowAnnot <- protAnnot

  if (GRP2$processing_options$remove_cont || GRP2$processing_options$remove_decoys) {
    lfqdata <- lfqdata$get_subset(protAnnot$clean(
      contaminants = GRP2$processing_options$remove_cont,
      decoys = GRP2$processing_options$remove_decoys
    ))
    logger::log_info(
      paste0(
        "removing contaminants and reverse sequences with patterns:",
        GRP2$processing_options$pattern_contaminants,
        GRP2$processing_options$pattern_decoys
      )
    )
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

  formula <- get_formula(transformed$config$table, interaction = GRP2$processing_options$interaction)
  message("FORMULA :", formula)
  GRP2$RES$formula <- formula
  formula_Condition <- prolfqua::strategy_lm(formula)
  # specify model definition
  mod <- prolfqua::build_model(
    transformed,
    formula_Condition,
    subject_Id = transformed$config$table$hierarchy_keys()
  )

  logger::log_info("fitted model with formula : {formula}")
  GRP2$RES$models <- mod

  contr <- prolfqua::Contrasts$new(mod, GRP2$pop$Contrasts,
                                   modelName = "Linear_Model"
  )
  conrM <- prolfqua::ContrastsModerated$new(
    contr
  )

  if (is.null(GRP2$processing_options$model_missing) || GRP2$processing_options$model_missing) {
    mC <- prolfqua::ContrastsMissing$new(
      lfqdata = transformed,
      contrasts = GRP2$pop$Contrasts,
      modelName = "Imputed_Mean"
    )

    conMI <- prolfqua::ContrastsModerated$new(
      mC
    )

    res <- prolfqua::merge_contrasts_results(conrM, conMI)
    GRP2$RES$contrMerged <- res$merged
    GRP2$RES$contrMore <- res$more
  } else {
    GRP2$RES$contrMerged <- conrM
    GRP2$RES$contrMore <- NULL
  }


  datax <- GRP2$RES$contrMerged$get_contrasts()
  datax <- dplyr::inner_join(GRP2$RES$rowAnnot$row_annot, datax, multiple = "all")
  GRP2$RES$contrastsData <- datax

  datax <- datax |>
    dplyr::filter(.data$FDR < GRP2$processing_options$FDR_threshold &
                    abs(.data$diff) > GRP2$processing_options$diff_threshold)
  GRP2$RES$contrastsData_signif <- datax
  return(GRP2)
}


#' will generates formula
#' @export
get_formula <- function(table, interaction = FALSE) {
  ################## Run Modelling ###############
  # remove control column from factors.
  factors <- table$factor_keys_depth()[
    !grepl("^control", table$factor_keys_depth(), ignore.case = TRUE)
  ]

  # model with or without interactions
  if (is.null(interaction) || !interaction) {
    formula <- paste0(
      table$get_response(), " ~ ",
      paste(factors, collapse = " + ")
    )
  } else {
    formula <- paste0(
      table$get_response(), " ~ ",
      paste(factors, collapse = " * ")
    )
  }
}


#' will replace generate_DEA_reports
#' @export
generate_DEA_reports2 <- function(lfqdata,
                                  GRP2,
                                  prot_annot,
                                  Contrasts) {
  # Compute all possible 2 Grps to avoid specifying reference.
  GRP2$pop$Contrasts <- Contrasts
  logger::log_info("CONTRAST : ", paste(GRP2$pop$Contrasts, collapse = " "))
  lfqwork <- lfqdata$get_copy()
  # lfqwork$data <- lfqdata$data |> dplyr::filter(.data$Group_ %in% levels$Group_)
  grp2 <- make_DEA_report2(
    lfqwork,
    prot_annot,
    GRP2
  )
  return(grp2)
}



prep_result_list <- function(GRP2) {
  rd <- GRP2$RES$lfqData
  tr <- GRP2$RES$transformedlfqData
  ra <- GRP2$RES$rowAnnot
  contrasts <- data.frame(
    contrast_name = names(GRP2$pop$Contrasts),
    contrast = GRP2$pop$Contrasts
  )

  wideraw <- dplyr::inner_join(ra$row_annot, rd$to_wide()$data, multiple = "all")
  widetr <- dplyr::inner_join(ra$row_annot, tr$to_wide()$data, multiple = "all")

  ctr <- dplyr::inner_join(ra$row_annot, GRP2$RES$contrMerged$get_contrasts(), multiple = "all")
  ctr_wide <- dplyr::inner_join(ra$row_annot, GRP2$RES$contrMerged$to_wide(), multiple = "all")

  resultList <- list()

  resultList$annotation <- dplyr::inner_join(
    rd$factors(),
    rd$get_Summariser()$hierarchy_counts_sample(),
    by = rd$config$table$sampleName,
    multiple = "all"
  )

  resultList$normalized_abundances <- dplyr::inner_join(ra$row_annot, tr$data, multiple = "all")
  resultList$raw_abundances_matrix <- wideraw
  resultList$normalized_abundances_matrix <- widetr
  resultList$diff_exp_analysis <- ctr
  resultList$diff_exp_analysis_wide <- ctr_wide
  resultList$formula <- data.frame(formula = GRP2$RES$formula)
  resultList$summary <- GRP2$RES$Summary
  resultList$missing_information <- prolfqua::UpSet_interaction_missing_stats(rd$data, rd$config, tr = 1)$data
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
  fg <- fg |> dplyr::mutate(updown = paste0(contrast, ifelse(diff > 0, "_up", "_down")))
  ora_sig <- split(fg[[id_column]], fg$updown)
  ora_files <- list()
  for (i in names(ora_sig)) {
    filename <- paste0("ORA_", i, "_WU", workunit_id, ".txt")
    ff <- file.path(outpath, filename)
    ora_files[[filename]] <- ff
    logger::log_info("Writing File ", ff)
    write.table(unique(ora_sig[[i]]),
                file = ff, col.names = FALSE,
                row.names = FALSE, quote = FALSE
    )
  }
  return(ora_files)
}

.write_GSEA <- function(fg, outpath, workunit_id, id_column) {
  gsea <- dplyr::select(fg, c("contrast", id_column, "statistic")) |>
    dplyr::arrange(.data$statistic)
  gsea <- gsea |>
    dplyr::group_by(dplyr::across(c("contrast", id_column))) |>
    dplyr::summarize(statistic = mean(.data$statistic)) |>
    dplyr::ungroup()
  gsea <- split(dplyr::select(gsea, c(id_column, "statistic")), gsea$contrast)
  gsea_files <- list()
  for (i in names(gsea)) {
    filernk <- paste0("GSEA_", i, "_WU", workunit_id, ".rnk")
    ff <- file.path(outpath, filernk)
    gsea_files[[filernk]] <- ff
    logger::log_info("Writing File ", ff)
    write.table(na.omit(gsea[[i]]),
                file = ff, col.names = FALSE,
                row.names = FALSE, quote = FALSE, sep = "\t"
    )
  }
  return(gsea_files)
}

write_result_list <- function(outpath,
                              GRP2,
                              resultList,
                              xlsxname,
                              GSEA = TRUE,
                              ORA = TRUE) {
  workunit_id <- GRP2$project_spec$workunit_Id
  dir.create(outpath)
  id_column <- GRP2$RES$rowAnnot$cleaned_ids

  if (ORA) {
    ff <- file.path(outpath, paste0("ORA_background_WU", workunit_id, ".txt"))
    write.table(GRP2$RES$rowAnnot$row_annot[[id_column]],
                file = ff, col.names = FALSE,
                row.names = FALSE, quote = FALSE
    )
    ora_files <- .write_ORA(GRP2$RES$contrastsData_signif, outpath, workunit_id, id_column = id_column)
  }
  if (GSEA) {
    gsea_files <- .write_GSEA(GRP2$RES$contrastsData, outpath, workunit_id, id_column)
  }
  if (nrow(resultList$normalized_abundances) > 1048575) {
    resultList$normalized_abundances <- NULL
  }
  xlsx_file <- file.path(outpath, paste0(xlsxname, ".xlsx"))
  writexl::write_xlsx(resultList, path = xlsx_file)
  return(list(xlsx_file = xlsx_file, ora_files = ora_files, gsea_files = gsea_files))
}

#' Write differential expression analysis results
#'
#' @param GRP2 return value of \code{\link{make_DEA_report}}
#' @param outpath path to place output
#' @param xlsxname file name for xlsx
#' @export
#' @family workflow
#'
write_DEA <- function(GRP2,
                      outpath,
                      ORA = TRUE,
                      GSEA = TRUE,
                      xlsxname = "AnalysisResults",
                      write = TRUE) {
  resultList <- prep_result_list(GRP2)
  if (write) {
    outfile <- write_result_list(outpath, GRP2, resultList, xlsxname = xlsxname, ORA = ORA, GSEA = GSEA)
  }
  return(list(resultList = resultList, outfile = outfile))
}


#' Generate differential expression analysis reports
#'
#' Writes results of DEA see \code{\link{generate_DEA_reports}}
#' @export
#'
write_DEA_all <- function(
    grp2,
    name = "",
    boxplot = TRUE,
    render = TRUE,
    ORA = TRUE,
    GSEA = TRUE,
    markdown = "_Grp2Analysis.Rmd") {
  dir.create(grp2$get_zipdir())
  name <- if (name != "") {
    paste0(name, "_")
  }
  fname <- paste0("DE_", name, "WU", grp2$project_spec$workunit_Id)
  qcname <- paste0("QC_", name, "WU", grp2$project_spec$workunit_Id)
  outpath <- grp2$get_result_dir()

  data_files <- prolfquapp::write_DEA(grp2, outpath = outpath, xlsxname = fname, ORA = ORA, GSEA = GSEA)$outfile

  if (render) {
    dea_file <- prolfquapp::render_DEA(grp2, outpath = outpath, htmlname = fname, markdown = markdown)
    qc_file <- prolfquapp::render_DEA(grp2, outpath = outpath, htmlname = qcname, markdown = "_DiffExpQC.Rmd")
  }

  bb <- grp2$RES$transformedlfqData
  grsizes <- bb$factors() |>
    dplyr::group_by(dplyr::across(bb$config$table$factor_keys_depth())) |>
    dplyr::summarize(n = n()) |>
    dplyr::pull(n)
  if (boxplot) {
    nr_controls <- sum(!grepl("^control", bb$config$table$factor_keys(), ignore.case = TRUE))
    if (nr_controls > 1 && all(grsizes == 1)
    ) {
      prolfquapp::writeLinesPaired(bb, outpath)
    } else {
      pl <- bb$get_Plotter()
      pl$write_boxplots(outpath)
    }
  }
  return(list(dea_file = dea_file, qc_file = qc_file, data_files = data_files))
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
                       htmlname = "Result2Grp",
                       word = FALSE,
                       markdown = "_Grp2Analysis.Rmd") {
  dir.create(outpath)

  rmarkdown::render(
    markdown,
    params = list(grp = GRP2),
    output_format = if (word) {
      bookdown::word_document2(toc = TRUE, toc_float = TRUE)
    } else {
      bookdown::html_document2(toc = TRUE, toc_float = TRUE)
    }
  )
  extension <- if (word) {
    ".docx"
  } else {
    ".html"
  }
  fname <- paste0(tools::file_path_sans_ext(markdown), extension)
  outfile <- file.path(outpath, paste0(htmlname, extension))
  if (file.copy(fname, outfile, overwrite = TRUE)) {
    file.remove(fname)
  }
  return(outfile)
}

#' convert tibble to data.frame with rownames
#' @param .data a tibble or data.frame
#' @param var name of the column with new row.names
#' @return a data.frame with rownames
#' @export
#' @examples
#' ind <- tibble::tibble(a = 1:3, rowname = letters[1:3])
#' column_to_rownames(ind)
column_to_rownames <- function(.data, var = "rowname", sep = "~lfq~") {
  # .data[,var] |> tidyr::unite("id", everything, sep = sep) |> pull("id")
  res <- as.data.frame(.data)
  rownames(res) <- .data[, var] |>
    tidyr::unite("id", tidyselect::everything(), sep = sep) |>
    dplyr::pull("id")
  return(res)
}

strip_rownames <- function(.data, strip = "~lfq~light$") {
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
bfabric_url_builder <- function(project_spec) {
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
                                      strip = "~lfq~light",
                                      .url_builder = bfabric_url_builder) {
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
  col.data <- col.data[colnames(mat.raw), ]
  x <- SummarizedExperiment::SummarizedExperiment(
    assays = list(rawData = mat.raw, transformedData = mat.trans),
    colData = col.data, metadata = list(bfabric_urls = .url_builder(GRP2$project_spec), contrasts = resTables$contrasts, formula = resTables$formula)
  )

  diffbyContrast <- split(resTables$diff_exp_analysis, resTables$diff_exp_analysis$contrast)
  for (i in names(diffbyContrast)) {
    row.data <- column_to_rownames(diffbyContrast[[i]], var = rowname)
    row.data <- row.data[rownames(mat.raw), ]
    if ("nr_compounds" %in% colnames(row.data)) {
      row.data <- dplyr::rename(nr_peptides = "nr_compounds")
    }
    SummarizedExperiment::rowData(x)[[paste0("constrast_", i)]] <- row.data
  }

  SummarizedExperiment::rowData(x)[["stats_normalized_wide"]] <- column_to_rownames(resTables$stats_normalized_wide, var = rowname)[rownames(mat.raw), ]
  SummarizedExperiment::rowData(x)[["stats_raw_wide"]] <- column_to_rownames(resTables$stats_raw_wide, var = rowname)[rownames(mat.raw), ]
  return(x)
}


.test_links <- list(dea_file = "./DEA_20250704_PI35298_O38953_WUtotal_proteome_none/Results_WU_total_proteome/DE_WUtotal_proteome.html",
                    qc_file = "./DEA_20250704_PI35298_O38953_WUtotal_proteome_none/Results_WU_total_proteome/QC_WUtotal_proteome.html",
                    data_files = list(xlsx_file = "./DEA_20250704_PI35298_O38953_WUtotal_proteome_none/Results_WU_total_proteome/DE_WUtotal_proteome.xlsx",
                                      ora_files = list(ORA_Treated_vs_Control_down_WUtotal_proteome.txt = "./DEA_20250704_PI35298_O38953_WUtotal_proteome_none/Results_WU_total_proteome/ORA_Treated_vs_Control_down_WUtotal_proteome.txt",
                                                       ORA_Treated_vs_Control_up_WUtotal_proteome.txt = "./DEA_20250704_PI35298_O38953_WUtotal_proteome_none/Results_WU_total_proteome/ORA_Treated_vs_Control_up_WUtotal_proteome.txt"),
                                      gsea_files = list(`GSEA_Treated_vs_Control_WUtotal_proteome.rnk` = "./DEA_20250704_PI35298_O38953_WUtotal_proteome_none/Results_WU_total_proteome/GSEA_Treated_vs_Control_WUtotal_proteome.rnk"),
                                      ibaq_file = "./DEA_20250704_PI35298_O38953_WUtotal_proteome_none/Results_WU_total_proteome/IBAQ_total_proteome.xlsx"))

#' write index.html file with links to all relevant files:
#' @export
#' @examples
#' .resdir <- "."
#' write_index_html(prolfquapp:::.test_links,tempdir())
write_index_html <- function(file_path_list, result_dir) {
  # Determine top‐level directory and name
  dea_path    <- file_path_list$dea_file
  topdir_path <- dirname(dea_path)
  topdir_name <- basename(topdir_path)

  # Helper to make a section
  make_section <- function(title, paths) {
    lines <- c(
      sprintf("  <h2>%s</h2>", title),
      "  <ul>"
    )
    for (i in seq_along(paths)) {
      p    <- paths[i]
      name <- names(paths)[i]
      text <- if (!is.null(name) && nzchar(name)) name else basename(p)
      rel  <- gsub(result_dir, "", p)
      rel <- paste0("./",rel)
      lines <- c(
        lines,
        sprintf('    <li><a href="%s">%s</a></li>', rel, text)
      )
    }
    c(lines, "  </ul>")
  }

  style_block <- R"(
  <style>
    .footer {
      text-align: center;
      margin-top: 40px;
      padding-top: 20px;
      border-top: 1px solid #ecf0f1;
      color: #7f8c8d;
    }
    h2 {
      text-align: left;
      margin-top: 20px;
      padding-bottom: 10px;
      border-bottom: 1px solid #ecf0f1;
      color: #2c3e50;
    }

    h3 {
      text-align: center;
      margin-top: 10px;
      color: #7f8c8d;
    }
    body {
      font-family: Arial, sans-serif;
    }

    h1 {
      text-align: center;
    }
  </style>
  )"


  # Build html
  title <- paste0("Differential Expression Analysis (DEA)")
  html_lines <- c(
    "<!DOCTYPE html>",
    "<html>",
    "<head>",
    "  <meta charset='UTF-8'>",
    sprintf("  <title>%s results</title>", title),

    style_block,
    "</head>",
    "<body>",
    sprintf("  <h1>%s results</h1>", title),
    sprintf(" <h3>for: %s</h3>", topdir_name),
    "<h3>Functional Genomics Center Zurich</h3>"
  )

  # 1. Top‐level files: dea + qc
  html_lines <- c(
    html_lines,
    make_section("Overview files", c(
      "<b>DEA Report (read first)</b>" = file_path_list$dea_file,
      "QC report" = file_path_list$qc_file
    )),
    make_section("Spreadsheet exports", c(
      "DEA results (XLSX)" = file_path_list$data_files$xlsx_file,
      "Intensity Based Absolute Quantitation (XLSX)" = file_path_list$data_files$ibaq_file
    )),
    make_section("ORA inputs:", unlist(file_path_list$data_files$ora_files)),
    make_section("GSEA ranklists:", unlist(file_path_list$data_files$gsea_files))
  )


  analysis_date <- format(Sys.time(), "%B %d, %Y")
  version       <- as.character(packageVersion("prolfquapp"))

  # Close out
  html_lines <- c(
    html_lines,
    "  <div class='footer'>",
    sprintf("    <p><strong>Analysis Date:</strong> %s</p>",    analysis_date),
    "    <p><strong>Project:</strong> Differential Expression analysis</p>",
    sprintf("    <p><strong>Generated by:</strong> prolfquapp package v%s</p>", version),
    '    <p><strong>Publication:</strong> <a href="https://doi.org/10.1021/acs.jproteome.4c00911" target="_blank">Wolski et al., J Proteome Res. 2025;24(2):955-965</a></p>',
    "  </div>",
    "</body>",
    "</html>"
  )
  # Write it
  index_file <- file.path(result_dir, "index.html")
  writeLines(html_lines, con = index_file)
  message("Wrote HTML index to: ", index_file)
  invisible(html_lines)
}


