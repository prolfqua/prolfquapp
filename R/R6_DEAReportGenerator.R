#' DEAReportGenerator
#'
#' Generates all output files for a differential expression analysis.
#' Uses a DEAnalyse object as data source instead of the legacy GRP2$RES list.
#'
#' @export
#'
DEAReportGenerator <- R6::R6Class(
  "DEAReportGenerator",
  public = list(
    #' @field deanalyse DEAnalyse object containing all analysis results
    deanalyse = NULL,
    #' @field GRP2 ProlfquAppConfig object containing analysis configuration
    GRP2 = NULL,
    #' @field fname filename prefix for DEA results
    fname = "",
    #' @field qcname filename prefix for QC results
    qcname = "",
    #' @field resultdir directory for storing results
    resultdir = "",
    #' @field ZIPDIR zip directory path
    ZIPDIR = "",

    #' @description
    #' Initialize DEAReportGenerator
    #' @param deanalyse DEAnalyse R6 object with completed analysis
    #' @param GRP2 ProlfquAppConfig R6 object
    #' @param name optional name prefix for output files
    initialize = function(deanalyse, GRP2, name = "") {
      stopifnot("DEAnalyse" %in% class(deanalyse))
      stopifnot("ProlfquAppConfig" %in% class(GRP2))
      self$deanalyse <- deanalyse
      self$GRP2 <- GRP2
      self$ZIPDIR <- GRP2$get_zipdir()
      name_prefix <- if (nchar(name) > 0) paste0(name, "_") else ""
      self$fname <- paste0("DE_", name_prefix, "WU", GRP2$project_spec$workunit_Id)
      self$qcname <- paste0("QC_", name_prefix, "WU", GRP2$project_spec$workunit_Id)
      self$resultdir <- GRP2$get_result_dir()
      logger::log_info("writing into : ", self$resultdir, " <<<<")
      dir.create(self$ZIPDIR, showWarnings = FALSE, recursive = TRUE)
      dir.create(self$resultdir, showWarnings = FALSE, recursive = TRUE)
    },

    #' @description
    #' Prepare result list with all analysis outputs for XLSX
    #' @return list containing all analysis results (14 sheets)
    prep_result_list = function() {
      dea <- self$deanalyse
      rd <- dea$lfq_data
      tr <- dea$lfq_data_transformed
      ra <- dea$rowAnnot

      contrasts_df <- data.frame(
        contrast_name = names(dea$contrasts),
        contrast = dea$contrasts
      )

      wideraw <- dplyr::inner_join(ra$row_annot, rd$to_wide()$data, multiple = "all")
      widetr <- dplyr::inner_join(ra$row_annot, tr$to_wide()$data, multiple = "all")

      contr_obj <- dea$contrast_results[[dea$default_model]]
      ctr <- dplyr::inner_join(ra$row_annot, contr_obj$get_contrasts(), multiple = "all")
      ctr_wide <- dplyr::inner_join(ra$row_annot, contr_obj$to_wide(), multiple = "all")

      resultList <- list()

      resultList$annotation <- dplyr::inner_join(
        rd$factors(),
        rd$get_Summariser()$hierarchy_counts_sample(),
        by = rd$config$sampleName,
        multiple = "all"
      )

      resultList$normalized_abundances <- dplyr::inner_join(ra$row_annot, tr$data, multiple = "all")
      resultList$raw_abundances_matrix <- wideraw
      resultList$normalized_abundances_matrix <- widetr
      resultList$diff_exp_analysis <- ctr
      resultList$diff_exp_analysis_wide <- ctr_wide
      resultList$formula <- data.frame(formula = dea$formula)
      resultList$summary <- dea$summary
      resultList$missing_information <- prolfqua::UpSet_interaction_missing_stats(rd$data, rd$config, tr = 1)$data
      resultList$contrasts <- contrasts_df

      # add protein statistics
      st <- tr$get_Stats()
      resultList$stats_normalized <- st$stats()
      resultList$stats_normalized_wide <- st$stats_wide()

      st <- rd$get_Stats()
      resultList$stats_raw <- st$stats()
      resultList$stats_raw_wide <- st$stats_wide()
      return(resultList)
    },

    #' @description
    #' Write DEA results (XLSX, ORA, GSEA files)
    #' @param ORA if TRUE write ORA gene lists
    #' @param GSEA if TRUE write GSEA rank files
    #' @return list with xlsx_file, ora_files, gsea_files paths
    write_DEA = function(ORA = TRUE, GSEA = TRUE) {
      resultList <- self$prep_result_list()
      dea <- self$deanalyse
      outpath <- self$resultdir
      workunit_id <- self$GRP2$project_spec$workunit_Id
      id_column <- dea$rowAnnot$cleaned_ids

      dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

      ora_files <- list()
      if (ORA) {
        ff <- file.path(outpath, paste0("ORA_background_WU", workunit_id, ".txt"))
        write.table(dea$rowAnnot$row_annot[[id_column]],
          file = ff, col.names = FALSE,
          row.names = FALSE, quote = FALSE
        )
        ora_files <- .write_ORA(dea$annotated_contrasts_signif, outpath, workunit_id, id_column = id_column)
      }

      gsea_files <- list()
      if (GSEA) {
        gsea_files <- .write_GSEA(dea$annotated_contrasts, outpath, workunit_id, id_column)
      }

      if (nrow(resultList$normalized_abundances) > 1048575) {
        resultList$normalized_abundances <- NULL
      }

      xlsx_file <- file.path(outpath, paste0(self$fname, ".xlsx"))
      writexl::write_xlsx(resultList, path = xlsx_file)
      return(list(xlsx_file = xlsx_file, ora_files = ora_files, gsea_files = gsea_files))
    },

    #' @description
    #' Render DEA report using R Markdown
    #' @param htmlname name for the output HTML file
    #' @param markdown path to the R Markdown template file
    #' @param word logical, if TRUE output Word document, otherwise HTML
    #' @param toc logical, if TRUE include table of contents
    #' @return path to the output file
    render_DEA = function(htmlname, markdown = "_Grp2Analysis_V2_R6.Rmd", word = FALSE, toc = TRUE) {
      dir.create(self$resultdir, showWarnings = FALSE, recursive = TRUE)

      rmarkdown::render(
        markdown,
        params = list(deanalyse = self$deanalyse),
        output_format = if (word) {
          bookdown::word_document2(toc = toc, toc_float = toc)
        } else {
          bookdown::html_document2(toc = toc, toc_float = toc)
        }
      )

      extension <- if (word) ".docx" else ".html"
      tmpname <- paste0(tools::file_path_sans_ext(markdown), extension)
      outfile <- file.path(self$resultdir, paste0(htmlname, extension))
      if (file.copy(tmpname, outfile, overwrite = TRUE)) {
        file.remove(tmpname)
      }
      return(outfile)
    },

    #' @description
    #' Generate boxplots for quality control
    #' @param boxplot logical, if TRUE write boxplots
    make_boxplots = function(boxplot = TRUE) {
      if (!boxplot) return(invisible(NULL))
      bb <- self$deanalyse$lfq_data_transformed
      grsizes <- bb$factors() |>
        dplyr::group_by(dplyr::across(bb$config$factor_keys_depth())) |>
        dplyr::summarize(n = dplyr::n(), .groups = "drop") |>
        dplyr::pull(n)
      nr_controls <- sum(!grepl("^control", bb$config$factor_keys(), ignore.case = TRUE))
      if (nr_controls > 1 && all(grsizes == 1)) {
        prolfquapp::writeLinesPaired(bb, self$resultdir)
      } else {
        pl <- bb$get_Plotter()
        pl$write_boxplots(self$resultdir)
      }
    },

    #' @description
    #' Write all DEA results: XLSX, ORA, GSEA, HTML reports, boxplots
    #' @param boxplot if TRUE generate boxplots
    #' @param render if TRUE render HTML reports
    #' @param ORA if TRUE write ORA gene lists
    #' @param GSEA if TRUE write GSEA rank files
    #' @param markdown Rmd template for main DEA report
    #' @param markdown_qc Rmd template for QC report
    #' @param toc if TRUE include table of contents
    #' @return list with dea_file, qc_file, data_files paths
    write_DEA_all = function(
        boxplot = TRUE,
        render = TRUE,
        ORA = TRUE,
        GSEA = TRUE,
        markdown = "_Grp2Analysis_V2_R6.Rmd",
        markdown_qc = "_DiffExpQC_R6.Rmd",
        toc = TRUE) {
      data_files <- self$write_DEA(ORA = ORA, GSEA = GSEA)

      dea_file <- NULL
      qc_file <- NULL
      if (render) {
        dea_file <- self$render_DEA(
          htmlname = self$fname,
          markdown = markdown,
          toc = toc
        )
        qc_file <- self$render_DEA(
          htmlname = self$qcname,
          markdown = markdown_qc,
          toc = toc
        )
      }

      self$make_boxplots(boxplot = boxplot)

      return(list(dea_file = dea_file, qc_file = qc_file, data_files = data_files))
    },

    #' @description
    #' Create SummarizedExperiment object from analysis results
    #' @param strip pattern to strip from rownames
    #' @param .url_builder function to build URLs for bfabric
    #' @return SummarizedExperiment object
    make_SummarizedExperiment = function(strip = "~lfq~light",
                                         .url_builder = prolfquapp::bfabric_url_builder) {
      dea <- self$deanalyse
      colname <- dea$lfq_data$config$sampleName
      rowname <- dea$lfq_data$config$hierarchyKeys()
      resTables <- self$prep_result_list()

      matTr <- dea$lfq_data_transformed$to_wide(as.matrix = TRUE)
      matRaw <- dea$lfq_data$to_wide(as.matrix = TRUE)

      mat.raw <- prolfquapp::strip_rownames(matRaw$data, strip)
      mat.trans <- prolfquapp::strip_rownames(matTr$data, strip)
      col.data <- prolfquapp::column_to_rownames(matRaw$annotation, var = colname)
      col.data <- col.data[colnames(mat.raw), ]
      x <- SummarizedExperiment::SummarizedExperiment(
        assays = list(rawData = mat.raw, transformedData = mat.trans),
        colData = col.data,
        metadata = list(
          bfabric_urls = .url_builder(self$GRP2$project_spec),
          contrasts = resTables$contrasts,
          formula = resTables$formula
        )
      )

      diffbyContrast <- split(resTables$diff_exp_analysis, resTables$diff_exp_analysis$contrast)
      for (i in names(diffbyContrast)) {
        row.data <- prolfquapp::column_to_rownames(diffbyContrast[[i]], var = rowname)
        row.data <- row.data[rownames(mat.raw), ]
        SummarizedExperiment::rowData(x)[[paste0("constrast_", i)]] <- row.data
      }

      SummarizedExperiment::rowData(x)[["stats_normalized_wide"]] <-
        prolfquapp::column_to_rownames(resTables$stats_normalized_wide, var = rowname)[rownames(mat.raw), ]
      SummarizedExperiment::rowData(x)[["stats_raw_wide"]] <-
        prolfquapp::column_to_rownames(resTables$stats_raw_wide, var = rowname)[rownames(mat.raw), ]
      return(x)
    }
  )
)
