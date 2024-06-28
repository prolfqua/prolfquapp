#'

DEAReportGenerator <- R6::R6Class(
  "DEAReportGenerator",
  public = list(
    lfqdata = NULL,
    GRP2 = NULL,
    prot_annot = NULL,
    Contrasts = NULL,
    fname = "",
    qcname = "",
    resultdir = "",
    initialize = function(lfqdata,
                          GRP2,
                          prot_annot,
                          Contrasts
                          ) {
      self$lfqdata = lfqdata
      self$GRP2 = GRP2
      self$prot_annot = prot_annot
      self$Contrasts = Contrasts
      self$ZIPDIR = GRP2$zipdir
      self$fname = paste0("DE_",  name, "_WU", grp2$project_spec$workunit_Id )
      self$qcname = paste0("QC_", name, "_WU", grp2$project_spec$workunit_Id )
      self$resultdir =  file.path( self$ZIPDIR, paste0("Results_DEA_WU", grp2$project_spec$workunit_Id))
      logger::log_info("writing into : ", self$resultdir, " <<<<")
      dir.create( self$ZIPDIR )
      dir.create( self$resultdir )
    },

    write_DEA_all = function(){

    },

    render_DEA = function(htmlname ,markdown = "_Grp2Analysis.Rmd", word = FALSE){

      rmarkdown::render(
        markdown,
        params = list(grp = self$GRP2) ,
        output_format = if (word) {
          bookdown::word_document2(toc = TRUE, toc_float = TRUE) } else {
            bookdown::html_document2(toc = TRUE, toc_float = TRUE)
          }
      )
      tmpname <- paste0(tools::file_path_sans_ext(markdown), if (word) {".docx"} else {".html"})
      if (file.copy(tmpname, file.path(self$resultdir, paste0(htmlname , if (word) {".docx"} else {".html"})), overwrite = TRUE)) {
        file.remove(tmpname)
      }
    },

    make_boxplots = function(){
      bb <- self$GRP2$RES$transformedlfqData
      grsizes <- bb$factors() |>
        dplyr::group_by(dplyr::across(bb$config$table$factor_keys_depth())) |>
        dplyr::summarize(n = n()) |>
        dplyr::pull( n )
      if (boxplot) {
        if (sum(!grepl("^control",bb$config$table$factor_keys(), ignore.case = TRUE))  > 1 &
            all(grsizes == 1)
        ) {
          prolfquapp::writeLinesPaired(bb, self$resultdir)
        } else {
          pl <- bb$get_Plotter()
          pl$write_boxplots( self$resultdir)
        }
      }
    },

    prep_result_list = function(){
      rd <- self$GRP2$RES$lfqData
      tr <- self$GRP2$RES$transformedlfqData
      ra <- self$GRP2$RES$rowAnnot

      contrasts <- data.frame(
        contrast_name = names(self$GRP2$pop$Contrasts),
        contrast = self$GRP2$pop$Contrasts
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
    },

    make_SummarizedExperiment = function(strip="~lfq~light",
                                          .url_builder = bfabric_url_builder){

      colname <- self$GRP2$RES$lfqData$config$table$sampleName
      rowname <- self$GRP2$RES$lfqData$config$table$hierarchyKeys()
      resTables <- self$prep_result_list(GRP2)

      matTr <- self$GRP2$RES$transformedlfqData$to_wide(as.matrix = TRUE)
      matRaw <- self$GRP2$RES$transformedlfqData$to_wide(as.matrix = TRUE)

      mat.raw <- strip_rownames(matRaw$data, strip)
      mat.trans <- strip_rownames(matTr$data, strip)
      col.data <- column_to_rownames(matRaw$annotation, var = colname)
      col.data <- col.data[colnames(mat.raw),]
      x <- SummarizedExperiment::SummarizedExperiment(
        assays = list(rawData = mat.raw, transformedData = mat.trans),
        colData = col.data, metadata = list(bfabric_urls = .url_builder(self$GRP2$project_spec), contrasts = resTables$contrasts, formula = resTables$formula )
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






  ))
