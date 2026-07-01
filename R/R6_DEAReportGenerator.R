.safe_enrichment_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

.map_enrichment_ids <- function(data, row_annot, subject_id, id_column) {
  if (id_column %in% colnames(data)) {
    return(data)
  }
  if (identical(subject_id, id_column)) {
    data[[id_column]] <- data[[subject_id]]
    return(data)
  }
  id_map <- row_annot |>
    dplyr::select(dplyr::all_of(c(subject_id, id_column))) |>
    dplyr::distinct()
  dplyr::left_join(data, id_map, by = subject_id, multiple = "all")
}

.as_enrichment_contrasts <- function(contrast_obj, subject_id) {
  has_rank <- is.function(contrast_obj$get_rank)
  has_ora <- is.function(contrast_obj$get_ora)
  if (has_rank && has_ora) {
    return(contrast_obj)
  }
  prolfqua::ContrastsTable$new(
    contrast_obj$get_contrasts(),
    subject_id = subject_id
  )
}

.write_ORA <- function(
  contrast_obj,
  row_annot,
  outpath,
  workunit_id,
  id_column = "IDcolumn",
  FDR_threshold = 0.05,
  diff_threshold = 1
) {
  cfg <- if (is.function(contrast_obj$get_config)) contrast_obj$get_config() else NULL
  subject_id <- if (!is.null(cfg) && length(cfg$subject_id) > 0) {
    cfg$subject_id
  } else {
    contrast_obj$subject_id
  }
  # Directional backends (SAINT) only output the "up" ORA list and use
  # a contrast-column-named filename prefix; symmetric backends (LM)
  # output both up/down lists.
  saint <- isTRUE(cfg$significance_directional)
  ora_up <- contrast_obj$get_ora(
    up = TRUE,
    FDR_threshold = FDR_threshold,
    diff_threshold = diff_threshold
  )
  ora_up <- .map_enrichment_ids(ora_up, row_annot, subject_id, id_column)
  ora_up <- ora_up |>
    dplyr::filter(!is.na(.data[[id_column]]))
  if (saint) {
    ora_sig <- split(ora_up[[id_column]], ora_up$contrast)
  } else {
    ora_down <- contrast_obj$get_ora(
      up = FALSE,
      FDR_threshold = FDR_threshold,
      diff_threshold = diff_threshold
    )
    ora_down <- .map_enrichment_ids(ora_down, row_annot, subject_id, id_column)
    ora_down <- ora_down |>
      dplyr::filter(!is.na(.data[[id_column]]))
    ora_up$updown <- paste0(ora_up$contrast, "_up")
    ora_down$updown <- paste0(ora_down$contrast, "_down")
    ora_all <- dplyr::bind_rows(ora_up, ora_down)
    ora_sig <- split(ora_all[[id_column]], ora_all$updown)
  }
  ora_files <- list()
  for (i in names(ora_sig)) {
    contrast_name <- .safe_enrichment_name(i)
    filename <- if (saint) {
      paste0("ORA_Bait_", contrast_name, "_WU", workunit_id, ".txt")
    } else {
      paste0("ORA_", contrast_name, "_WU", workunit_id, ".txt")
    }
    ff <- file.path(outpath, filename)
    ora_files[[filename]] <- ff
    logger::log_info("Writing File ", ff)
    write.table(
      unique(ora_sig[[i]]),
      file = ff,
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE
    )
  }
  return(ora_files)
}

.write_GSEA <- function(
  contrast_obj,
  row_annot,
  outpath,
  workunit_id,
  id_column
) {
  cfg <- if (is.function(contrast_obj$get_config)) contrast_obj$get_config() else NULL
  subject_id <- if (!is.null(cfg) && length(cfg$subject_id) > 0) {
    cfg$subject_id
  } else {
    contrast_obj$subject_id
  }
  saint <- isTRUE(cfg$significance_directional)
  gsea <- contrast_obj$get_rank()
  gsea <- .map_enrichment_ids(gsea, row_annot, subject_id, id_column)
  gsea <- gsea |>
    dplyr::filter(!is.na(.data[[id_column]]))
  gsea <- dplyr::select(
    gsea,
    dplyr::all_of(c("contrast", id_column, "score"))
  ) |>
    dplyr::arrange(.data$score)
  gsea <- gsea |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c("contrast", id_column)))) |>
    dplyr::summarize(score = mean(.data$score), .groups = "drop") |>
    dplyr::ungroup()
  gsea <- split(
    dplyr::select(gsea, dplyr::all_of(c(id_column, "score"))),
    gsea$contrast
  )
  gsea_files <- list()
  for (i in names(gsea)) {
    contrast_name <- .safe_enrichment_name(i)
    filernk <- if (saint) {
      paste0("Bait_", contrast_name, ".rnk")
    } else {
      paste0("GSEA_", contrast_name, "_WU", workunit_id, ".rnk")
    }
    ff <- file.path(outpath, filernk)
    gsea_files[[filernk]] <- ff
    logger::log_info("Writing File ", ff)
    write.table(
      na.omit(gsea[[i]]),
      file = ff,
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )
  }
  return(gsea_files)
}

custom_round <- function(arr) {
  cr <- function(x) {
    if (x == 0) {
      return(0)
    } else if (abs(x) >= 1) {
      return(round(x, 2))
    } else {
      return(signif(x, 2))
    }
  }
  return(vapply(arr, cr, numeric(1)))
}

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
      self$fname <- paste0(
        "DE_",
        name_prefix,
        "WU",
        GRP2$project_spec$workunit_Id
      )
      self$qcname <- paste0(
        "QC_",
        name_prefix,
        "WU",
        GRP2$project_spec$workunit_Id
      )
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
      rd <- dea$lfq_data_raw
      tr <- dea$lfq_data
      ra <- dea$rowAnnot

      contrasts_df <- data.frame(
        contrast_name = names(dea$contrasts),
        contrast = dea$contrasts
      )

      wideraw <- .join_annotation(ra$row_annot, rd$data_wide()$data, ra$pID)
      widetr <- .join_annotation(ra$row_annot, tr$data_wide()$data, ra$pID)

      contr_obj <- dea$contrast_results[[dea$default_model]]
      ctr <- .join_annotation(ra$row_annot, contr_obj$get_contrasts(), ra$pID)
      ctr_wide <- .join_annotation(ra$row_annot, contr_obj$to_wide(), ra$pID)

      resultList <- list()

      resultList$annotation <- dplyr::inner_join(
        rd$factors(),
        rd$get_Summariser()$hierarchy_counts_sample(),
        by = rd$sample_name(),
        multiple = "all"
      )

      resultList$normalized_abundances <- .join_annotation(
        ra$row_annot, tr$data_long(), ra$pID
      )
      resultList$raw_abundances_matrix <- wideraw
      resultList$normalized_abundances_matrix <- widetr
      resultList$diff_exp_analysis <- ctr
      resultList$diff_exp_analysis_wide <- ctr_wide
      resultList$formula <- data.frame(formula = dea$formula)
      resultList$summary <- dea$summary
      resultList$missing_information <- prolfqua::upset_interaction_missing_stats(
        rd,
        tr = 1
      )$data
      resultList$contrasts <- contrasts_df
      # Backend-specific extras (SAINT input tables, etc.) are surfaced
      # through ContrastsInterface$extra_artifacts(); the default
      # returns an empty list, so LM-style backends add nothing here.
      extras <- contr_obj$extra_artifacts()
      if (length(extras) > 0) {
        resultList <- c(resultList, extras)
      }

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
      contrast_obj <- dea$contrast_results[[dea$default_model]]
      contrast_obj <- .as_enrichment_contrasts(
        contrast_obj,
        subject_id = dea$lfq_data$subject_id()
      )

      dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

      ora_files <- list()
      if (ORA) {
        ff <- file.path(
          outpath,
          paste0("ORA_background_WU", workunit_id, ".txt")
        )
        write.table(
          dea$rowAnnot$row_annot[[id_column]],
          file = ff,
          col.names = FALSE,
          row.names = FALSE,
          quote = FALSE
        )
        ora_files <- .write_ORA(
          contrast_obj,
          dea$rowAnnot$row_annot,
          outpath,
          workunit_id,
          id_column = id_column,
          FDR_threshold = dea$FDR_threshold,
          diff_threshold = dea$diff_threshold
        )
      }

      gsea_files <- list()
      if (GSEA) {
        gsea_files <- .write_GSEA(
          contrast_obj,
          dea$rowAnnot$row_annot,
          outpath,
          workunit_id,
          id_column
        )
      }

      if (nrow(resultList$normalized_abundances) > 1048575) {
        resultList$normalized_abundances <- NULL
      }

      xlsx_file <- file.path(outpath, paste0(self$fname, ".xlsx"))
      writexl::write_xlsx(resultList, path = xlsx_file)
      return(list(
        xlsx_file = xlsx_file,
        ora_files = ora_files,
        gsea_files = gsea_files
      ))
    },

    #' @description
    #' Render DEA report using R Markdown
    #' @param htmlname name for the output HTML file
    #' @param markdown path to the R Markdown template file
    #' @param word logical, if TRUE output Word document, otherwise HTML
    #' @param toc logical, if TRUE include table of contents
    #' @return path to the output file
    render_DEA = function(
      htmlname,
      markdown = "Grp2Analysis_V2_R6.Rmd",
      word = FALSE,
      toc = TRUE
    ) {
      dir.create(self$resultdir, showWarnings = FALSE, recursive = TRUE)

      rmarkdown::render(
        markdown,
        params = list(deanalyse = self$deanalyse),
        output_format = if (word) {
          bookdown::word_document2(toc = toc, toc_float = toc)
        } else {
          bookdown::html_document2(
            toc = toc,
            toc_float = toc,
            self_contained = TRUE,
            includes = rmarkdown::includes(
              in_header = system.file(
                "templates/fgcz_header.html",
                package = "prolfquapp"
              )
            ),
            css = system.file("templates/fgcz.css", package = "prolfquapp")
          )
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
    #' Generate sample-level boxplots for quality control
    #' @param boxplot logical, if TRUE write boxplots
    make_boxplots = function(boxplot = TRUE) {
      if (!boxplot) {
        return(invisible(NULL))
      }
      bb <- self$deanalyse$lfq_data
      # Paired layout (connect samples across the second factor) applies only
      # when there is a pairing factor (factor_keys()[2], used by
      # writeLinesPaired) and every factor-combination cell holds a single
      # sample. Group sizes and the factor count must use the SAME accessor.
      grsizes <- bb$factors() |>
        dplyr::group_by(dplyr::across(bb$factor_keys())) |>
        dplyr::summarize(n = dplyr::n(), .groups = "drop") |>
        dplyr::pull(n)
      nr_factors <- length(bb$factor_keys())
      if (nr_factors > 1 && all(grsizes == 1)) {
        prolfquapp::writeLinesPaired(bb, self$resultdir)
      } else {
        pl <- bb$get_Plotter()
        pl$write_boxplots(self$resultdir)
      }
    },

    #' @description
    #' Get subset of transformed data for significant proteins
    filter_data = function() {
      dea <- self$deanalyse
      dx <- dea$filter_contrasts()
      invisible(dea$lfq_data$get_subset(dx))
    },

    #' @description
    #' Get per-protein boxplots for significant proteins
    get_protein_boxplots = function() {
      self$filter_data()$get_Plotter()$boxplots()
    },

    #' @description
    #' Convert significant contrast results to table grobs. Column
    #' selection and rounding are driven by the contrast object's
    #' \code{ContrastConfiguration} so SAINT and LM backends both
    #' produce grobs with canonical \code{contrast}/\code{effect}/
    #' \code{score}/\code{fdr} columns without backend-specific code.
    contrasts_to_Grob = function() {
      dea <- self$deanalyse
      contrast_obj <- dea$contrast_results[[dea$default_model]]
      cfg <- contrast_obj$get_config()
      datax <- dea$filter_contrasts()
      hkeys <- dea$lfq_data$relevant_hierarchy_keys()

      canonical <- dplyr::transmute(
        datax,
        !!!rlang::syms(hkeys),
        contrast = .data[[cfg$contrast_col]],
        effect = custom_round(.data[[cfg$effect_col]]),
        score = custom_round(.data[[cfg$score_col]]),
        fdr = custom_round(.data[[cfg$fdr_col]])
      )
      xdn <- canonical |> dplyr::nest_by(!!!rlang::syms(hkeys))
      grobs <- vector(mode = "list", length = nrow(xdn))
      pb <- progress::progress_bar$new(total = nrow(xdn))
      for (i in seq_len(nrow(xdn))) {
        pb$tick()
        grobs[[i]] <- gridExtra::tableGrob(xdn$data[[i]])
      }
      xdn$grobs <- grobs
      return(xdn)
    },

    #' @description
    #' Get per-protein boxplots combined with contrast summary tables
    get_protein_boxplots_contrasts = function() {
      ctrG <- self$contrasts_to_Grob()
      bp <- self$get_protein_boxplots()
      stopifnot(nrow(ctrG) == nrow(bp))
      res <- vector(mode = "list", length = nrow(ctrG))
      pb <- progress::progress_bar$new(total = nrow(ctrG))
      for (i in seq_len(nrow(ctrG))) {
        res[[i]] <- gridExtra::arrangeGrob(
          bp$boxplot[[i]],
          ctrG$grobs[[i]],
          nrow = 2,
          heights = c(2 / 3, 1 / 3)
        )
        pb$tick()
      }
      ctrG$bxpl_grobs <- res
      return(ctrG)
    },

    #' @description
    #' Write per-protein boxplots with contrast tables to PDF
    #' @param filename base filename (without extension)
    write_protein_boxplots = function(filename = "boxplots") {
      dea <- self$deanalyse
      ctrG <- self$get_protein_boxplots_contrasts()
      filename <- paste0(
        filename,
        "_FDR_",
        dea$FDR_threshold,
        "_diff_",
        dea$diff_threshold,
        ".pdf"
      )
      logger::log_info("start writing boxplots into file : ", filename)
      pdf(file = file.path(self$ZIPDIR, filename))
      pb <- progress::progress_bar$new(total = length(ctrG$bxpl_grobs))
      for (i in seq_along(ctrG$bxpl_grobs)) {
        pb$tick()
        grid::grid.newpage()
        grid::grid.draw(ctrG$bxpl_grobs[[i]])
      }
      dev.off()
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
      markdown = "Grp2Analysis_V2_R6.Rmd",
      markdown_qc = "DiffExpQC_R6.Rmd",
      toc = TRUE
    ) {
      data_files <- self$write_DEA(ORA = ORA, GSEA = GSEA)

      dea_file <- NULL
      qc_file <- NULL
      if (render) {
        dea_file <- self$render_DEA(
          htmlname = self$fname,
          markdown = markdown,
          toc = toc
        )
        contrast_obj <- self$deanalyse$contrast_results[[
          self$deanalyse$default_model
        ]]
        if (isTRUE(contrast_obj$get_config()$supports_dea_qc)) {
          qc_file <- self$render_DEA(
            htmlname = self$qcname,
            markdown = markdown_qc,
            toc = toc
          )
        }
      }

      self$make_boxplots(boxplot = boxplot)

      return(list(
        dea_file = dea_file,
        qc_file = qc_file,
        data_files = data_files
      ))
    },

    #' @description
    #' Create SummarizedExperiment object from analysis results
    #' @param strip pattern to strip from rownames
    #' @param .url_builder function to build URLs for bfabric
    #' @return SummarizedExperiment object
    make_SummarizedExperiment = function(
      strip = "~lfq~light",
      .url_builder = prolfquapp::bfabric_url_builder
    ) {
      dea <- self$deanalyse
      colname <- dea$lfq_data_raw$sample_name()
      rowname <- dea$lfq_data_raw$hierarchy_keys()
      resTables <- self$prep_result_list()

      matTr <- dea$lfq_data$data_wide(as.matrix = TRUE)
      matRaw <- dea$lfq_data_raw$data_wide(as.matrix = TRUE)

      mat.raw <- prolfquapp::strip_rownames(matRaw$data, strip)
      mat.trans <- prolfquapp::strip_rownames(matTr$data, strip)
      assays <- list(rawData = mat.raw, transformedData = mat.trans)

      nr_children_col <- dea$lfq_data_raw$nr_children_col()
      if (
        length(nr_children_col) == 1 &&
          nzchar(nr_children_col) &&
          nr_children_col %in% colnames(dea$lfq_data_raw$data_long())
      ) {
        mat_children <- dea$lfq_data_raw$data_wide(
          as.matrix = TRUE,
          value = nr_children_col
        )
        assays[["nr_children"]] <- prolfquapp::strip_rownames(
          mat_children$data,
          strip
        )[rownames(mat.raw), colnames(mat.raw), drop = FALSE]
      }

      col.data <- prolfquapp::column_to_rownames(
        matRaw$annotation,
        var = colname
      )
      col.data <- col.data[colnames(mat.raw), ]
      x <- SummarizedExperiment::SummarizedExperiment(
        assays = assays,
        colData = col.data,
        metadata = list(
          bfabric_urls = .url_builder(self$GRP2$project_spec),
          contrasts = resTables$contrasts,
          formula = resTables$formula,
          default_model = dea$default_model,
          analysis_configuration_raw = prolfqua::R6_extract_values(dea$lfq_data_raw$get_config()),
          analysis_configuration_transformed = prolfqua::R6_extract_values(dea$lfq_data$get_config())
        )
      )

      contrast_obj <- dea$contrast_results[[dea$default_model]]
      contrast_column <- contrast_obj$get_config()$contrast_col
      diffbyContrast <- split(
        resTables$diff_exp_analysis,
        resTables$diff_exp_analysis[[contrast_column]]
      )
      for (i in names(diffbyContrast)) {
        row.data <- prolfquapp::column_to_rownames(
          diffbyContrast[[i]],
          var = rowname
        )
        row.data <- row.data[rownames(mat.raw), ]
        SummarizedExperiment::rowData(x)[[paste0("constrast_", i)]] <- row.data
      }

      SummarizedExperiment::rowData(x)[["stats_normalized_wide"]] <-
        prolfquapp::column_to_rownames(
          resTables$stats_normalized_wide,
          var = rowname
        )[rownames(mat.raw), ]
      SummarizedExperiment::rowData(x)[["stats_raw_wide"]] <-
        prolfquapp::column_to_rownames(resTables$stats_raw_wide, var = rowname)[
          rownames(mat.raw),
        ]
      return(x)
    }
  )
)
