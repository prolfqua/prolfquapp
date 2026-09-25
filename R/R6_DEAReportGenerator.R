# Enrichment rows carrying `id_column` (joined from the row annotation when the
# backend does not provide it), without rows lacking an id.
.map_enrichment_ids <- function(data, row_annot, subject_id, id_column) {
  if (!id_column %in% colnames(data)) {
    id_map <- dplyr::distinct(dplyr::select(row_annot, dplyr::all_of(c(subject_id, id_column))))
    data <- dplyr::left_join(data, id_map, by = subject_id, multiple = "all")
  }
  dplyr::filter(data, !is.na(.data[[id_column]]))
}

# Writes one file per element of `sets`, named by `filename(<element name>)`.
.write_enrichment_files <- function(sets, outpath, filename, ...) {
  files <- list()
  for (i in names(sets)) {
    fname <- filename(gsub("[^A-Za-z0-9_.-]+", "_", i))
    files[[fname]] <- file.path(outpath, fname)
    logger::log_info("Writing File ", files[[fname]])
    write.table(sets[[i]], file = files[[fname]], col.names = FALSE, row.names = FALSE, quote = FALSE, ...)
  }
  files
}

custom_round <- function(arr) {
  cr <- function(x) {
    if (x == 0) {
      0
    } else if (abs(x) >= 1) {
      round(x, 2)
    } else {
      signif(x, 2)
    }
  }
  vapply(arr, cr, numeric(1))
}

#' DEAReportGenerator
#'
#' Generates all output files for a differential expression analysis from a DEAnalyse object.
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
      stopifnot("DEAnalyse" %in% class(deanalyse), "ProlfquAppConfig" %in% class(GRP2))
      self$deanalyse <- deanalyse
      self$GRP2 <- GRP2
      self$ZIPDIR <- GRP2$get_zipdir()
      suffix <- paste0(if (nchar(name) > 0) paste0(name, "_"), "WU", GRP2$project_spec$workunit_Id)
      self$fname <- paste0("DE_", suffix)
      self$qcname <- paste0("QC_", suffix)
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
      join <- function(x) .join_annotation(dea$rowAnnot$row_annot, x, dea$rowAnnot$pID)
      contr_obj <- dea$contrast_results[[dea$default_model]]

      resultList <- list(
        annotation = dplyr::inner_join(
          rd$factors(),
          rd$get_Summariser()$hierarchy_counts_sample(),
          by = rd$sample_name(),
          multiple = "all"
        ),
        normalized_abundances = join(tr$data_long()),
        raw_abundances_matrix = join(rd$data_wide()$data),
        normalized_abundances_matrix = join(tr$data_wide()$data),
        diff_exp_analysis = join(contr_obj$get_contrasts()),
        diff_exp_analysis_wide = join(contr_obj$to_wide()),
        formula = data.frame(formula = dea$formula)
      )
      resultList$summary <- dea$summary
      resultList$missing_information <- prolfqua::upset_interaction_missing_stats(rd, tr = 1)$data
      resultList$contrasts <- data.frame(contrast_name = names(dea$contrasts), contrast = dea$contrasts)
      # Backend-specific extras (e.g. SAINT input tables); empty for LM-style backends.
      resultList <- c(resultList, contr_obj$extra_artifacts())

      st <- tr$get_Stats()
      resultList$stats_normalized <- st$stats()
      resultList$stats_normalized_wide <- st$stats_wide()
      st <- rd$get_Stats()
      resultList$stats_raw <- st$stats()
      resultList$stats_raw_wide <- st$stats_wide()
      resultList
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
      wu <- self$GRP2$project_spec$workunit_Id
      id_column <- dea$rowAnnot$cleaned_ids
      row_annot <- dea$rowAnnot$row_annot
      contrast_obj <- dea$contrast_results[[dea$default_model]]
      cfg <- contrast_obj$get_config()
      subject_id <- if (length(cfg$subject_id) > 0) cfg$subject_id else contrast_obj$subject_id
      # Directional backends (SAINT) write only the "up" ORA list and bait-named files.
      saint <- isTRUE(cfg$significance_directional)
      dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

      ora_files <- list()
      if (ORA) {
        ff <- file.path(outpath, paste0("ORA_background_WU", wu, ".txt"))
        write.table(row_annot[[id_column]], file = ff, col.names = FALSE, row.names = FALSE, quote = FALSE)
        ora <- function(up) {
          res <- contrast_obj$get_ora(up = up, FDR_threshold = dea$FDR_threshold, diff_threshold = dea$diff_threshold)
          res <- .map_enrichment_ids(res, row_annot, subject_id, id_column)
          if (!saint) {
            res$contrast <- paste0(res$contrast, if (up) "_up" else "_down")
          }
          res
        }
        ora_all <- if (saint) ora(TRUE) else dplyr::bind_rows(ora(TRUE), ora(FALSE))
        ora_files <- .write_enrichment_files(
          lapply(split(ora_all[[id_column]], ora_all$contrast), unique),
          outpath,
          function(x) paste0(if (saint) "ORA_Bait_" else "ORA_", x, "_WU", wu, ".txt")
        )
      }

      gsea_files <- list()
      if (GSEA) {
        gsea <- .map_enrichment_ids(
          contrast_obj$get_rank(score = .gsea_rank_column(cfg)),
          row_annot,
          subject_id,
          id_column
        )
        gsea <- dplyr::arrange(gsea, .data$score) |>
          dplyr::group_by(dplyr::across(dplyr::all_of(c("contrast", id_column)))) |>
          dplyr::summarize(score = mean(.data$score), .groups = "drop")
        gsea_files <- .write_enrichment_files(
          lapply(split(gsea[, c(id_column, "score")], gsea$contrast), na.omit),
          outpath,
          function(x) if (saint) paste0("Bait_", x, ".rnk") else paste0("GSEA_", x, "_WU", wu, ".rnk"),
          sep = "\t"
        )
      }

      if (nrow(resultList$normalized_abundances) > 1048575) {
        resultList$normalized_abundances <- NULL
      }
      xlsx_file <- file.path(outpath, paste0(self$fname, ".xlsx"))
      writexl::write_xlsx(resultList, path = xlsx_file)
      list(xlsx_file = xlsx_file, ora_files = ora_files, gsea_files = gsea_files)
    },

    #' @description
    #' Generate sample-level boxplots for quality control
    #' @param boxplot logical, if TRUE write boxplots
    make_boxplots = function(boxplot = TRUE) {
      if (!boxplot) {
        return(invisible(NULL))
      }
      bb <- self$deanalyse$lfq_data
      # Paired layout only with a pairing factor (factor_keys()[2], used by
      # writeLinesPaired) and a single sample per factor-combination cell.
      grsizes <- bb$factors() |>
        dplyr::group_by(dplyr::across(bb$factor_keys())) |>
        dplyr::summarize(n = dplyr::n(), .groups = "drop") |>
        dplyr::pull(n)
      if (length(bb$factor_keys()) > 1 && all(grsizes == 1)) {
        prolfquapp::writeLinesPaired(bb, self$resultdir)
      } else {
        bb$get_Plotter()$write_boxplots(self$resultdir)
      }
    },

    #' @description
    #' Get subset of transformed data for significant proteins
    filter_data = function() {
      invisible(self$deanalyse$lfq_data$get_subset(self$deanalyse$filter_contrasts()))
    },

    #' @description
    #' Get per-protein boxplots for significant proteins
    get_protein_boxplots = function() {
      self$filter_data()$get_Plotter()$boxplots()
    },

    #' @description
    #' Convert significant contrast results to table grobs with the canonical
    #' \code{contrast}/\code{effect}/\code{score}/\code{fdr} columns, selected
    #' through the contrast object's \code{ContrastConfiguration}.
    contrasts_to_Grob = function() {
      dea <- self$deanalyse
      cfg <- dea$contrast_results[[dea$default_model]]$get_config()
      hkeys <- rlang::syms(dea$lfq_data$relevant_hierarchy_keys())
      xdn <- dplyr::transmute(
        dea$filter_contrasts(),
        !!!hkeys,
        contrast = .data[[cfg$contrast_col]],
        effect = custom_round(.data[[cfg$effect_col]]),
        score = custom_round(.data[[cfg$score_col]]),
        fdr = custom_round(.data[[cfg$fdr_col]])
      ) |>
        dplyr::nest_by(!!!hkeys)
      pb <- progress::progress_bar$new(total = nrow(xdn))
      xdn$grobs <- lapply(xdn$data, function(d) {
        pb$tick()
        gridExtra::tableGrob(d)
      })
      xdn
    },

    #' @description
    #' Get per-protein boxplots combined with contrast summary tables
    get_protein_boxplots_contrasts = function() {
      ctrG <- self$contrasts_to_Grob()
      bp <- self$get_protein_boxplots()
      stopifnot(nrow(ctrG) == nrow(bp))
      pb <- progress::progress_bar$new(total = nrow(ctrG))
      ctrG$bxpl_grobs <- lapply(seq_len(nrow(ctrG)), function(i) {
        pb$tick()
        gridExtra::arrangeGrob(bp$boxplot[[i]], ctrG$grobs[[i]], nrow = 2, heights = c(2 / 3, 1 / 3))
      })
      ctrG
    },

    #' @description
    #' Write per-protein boxplots with contrast tables to PDF
    #' @param filename base filename (without extension)
    write_protein_boxplots = function(filename = "boxplots") {
      dea <- self$deanalyse
      ctrG <- self$get_protein_boxplots_contrasts()
      filename <- paste0(filename, "_FDR_", dea$FDR_threshold, "_diff_", dea$diff_threshold, ".pdf")
      logger::log_info("start writing boxplots into file : ", filename)
      pdf(file = file.path(self$ZIPDIR, filename))
      pb <- progress::progress_bar$new(total = length(ctrG$bxpl_grobs))
      for (grob in ctrG$bxpl_grobs) {
        pb$tick()
        grid::grid.newpage()
        grid::grid.draw(grob)
      }
      dev.off()
    },

    #' @description
    #' Write DEA data outputs: XLSX, ORA gene lists, GSEA rank files, and
    #' boxplots. HTML reports are rendered by `render_dea_reports()`.
    #' @param boxplot if TRUE generate boxplots
    #' @param ORA if TRUE write ORA gene lists
    #' @param GSEA if TRUE write GSEA rank files
    #' @return list with `data_files` paths; `dea_file` / `qc_file` are NULL as
    #'   the Quarto reports are produced by `render_dea_reports()`
    write_DEA_all = function(boxplot = TRUE, ORA = TRUE, GSEA = TRUE) {
      data_files <- self$write_DEA(ORA = ORA, GSEA = GSEA)
      self$make_boxplots(boxplot = boxplot)
      list(dea_file = NULL, qc_file = NULL, data_files = data_files)
    },

    #' @description
    #' Create SummarizedExperiment object from analysis results
    #'
    #' For the \code{lm_impute} default model the object also carries the
    #' assay \code{imputedData}, \code{transformedData} with every missing
    #' cell filled by \code{prolfqua::impute_from_model()}, and the rowData
    #' block \code{imputation} with \code{n_observed}, \code{n_imputed} and
    #' \code{route} per feature. Decoys are not modelled and stay NA in both.
    #'
    #' Given \code{ibaq}, the object also carries the assay \code{ibaq};
    #' features without an IBAQ value are NA.
    #' @param strip pattern to strip from rownames
    #' @param .url_builder function to build URLs for bfabric
    #' @param ibaq optional protein-level \code{LFQData} with IBAQ values, as
    #'   returned by \code{compute_IBAQ_values()}
    #' @return SummarizedExperiment object
    make_SummarizedExperiment = function(
      strip = "~lfq~light",
      .url_builder = prolfquapp::bfabric_url_builder,
      ibaq = NULL
    ) {
      dea <- self$deanalyse
      raw <- dea$lfq_data_raw
      colname <- raw$sample_name()
      rowname <- raw$hierarchy_keys()
      resTables <- self$prep_result_list()
      wide <- function(lfq, ...) prolfquapp::strip_rownames(lfq$data_wide(as.matrix = TRUE, ...)$data, strip)

      matRaw <- raw$data_wide(as.matrix = TRUE)
      mat.raw <- prolfquapp::strip_rownames(matRaw$data, strip)
      features <- rownames(mat.raw)
      assays <- list(rawData = mat.raw, transformedData = wide(dea$lfq_data))

      nr_children_col <- raw$nr_children_col()
      if (length(nr_children_col) == 1 && nr_children_col %in% colnames(raw$data_long())) {
        assays[["nr_children"]] <- wide(raw, value = nr_children_col)[features, colnames(mat.raw), drop = FALSE]
      }
      contrast_obj <- dea$contrast_results[[dea$default_model]]
      imputation <- NULL
      if (identical(dea$default_model, "lm_impute")) {
        imputation <- .imputed_assay(contrast_obj, mat.raw, rowname, strip)
        assays[["imputedData"]] <- imputation$matrix
      }
      if (!is.null(ibaq)) {
        assays[["ibaq"]] <- .align_to_features(wide(ibaq), mat.raw)
      }

      col.data <- prolfquapp::column_to_rownames(matRaw$annotation, var = colname)[colnames(mat.raw), ]
      ps <- self$GRP2$project_spec
      x <- SummarizedExperiment::SummarizedExperiment(
        assays = assays,
        colData = col.data,
        metadata = list(
          artifact_type = "dea_results",
          schema_version = "2.1.0",
          source_software = as.character(self$GRP2$software),
          feature_keys = rowname,
          sample_key = colname,
          # The annotation column carrying the identifier enrichment tools are
          # given (STRING, ORA), so a consumer need not guess at column names.
          identifier_key = dea$rowAnnot$cleaned_ids,
          bfabric_urls = .url_builder(ps),
          provenance = list(
            project_Id = ps$project_Id,
            project_name = ps$project_name,
            order_Id = ps$order_Id,
            workunit_Id = ps$workunit_Id,
            input_URL = ps$input_URL,
            software = self$GRP2$software,
            model = dea$default_model
          ),
          contrasts = resTables$contrasts,
          formula = resTables$formula,
          default_model = dea$default_model,
          analysis_configuration_raw = prolfqua::R6_extract_values(raw$get_config()),
          analysis_configuration = prolfqua::R6_extract_values(dea$lfq_data$get_config()),
          contrast_configuration = prolfqua::R6_extract_values(contrast_obj$get_config()),
          processing_options = prolfqua::R6_extract_values(self$GRP2$processing_options)
        )
      )

      # Feature annotation is stored once, as its own rowData frame; the
      # contrast frames carry the feature keys and the results only.
      SummarizedExperiment::rowData(x)[["annotation"]] <- .feature_rows(
        dea$rowAnnot$row_annot,
        dea$rowAnnot$pID,
        features
      )
      contrasts <- contrast_obj$get_contrasts()
      diffbyContrast <- split(contrasts, contrasts[[contrast_obj$get_config()$contrast_col]])
      for (i in names(diffbyContrast)) {
        SummarizedExperiment::rowData(x)[[paste0("constrast_", i)]] <- .feature_rows(
          diffbyContrast[[i]],
          rowname,
          features
        )
      }
      for (stats in c("stats_normalized_wide", "stats_raw_wide")) {
        stats_rows <- prolfquapp::column_to_rownames(resTables[[stats]], var = rowname)
        SummarizedExperiment::rowData(x)[[stats]] <- stats_rows[features, ]
      }
      if (!is.null(imputation)) {
        SummarizedExperiment::rowData(x)[["imputation"]] <- imputation$summary
      }
      x
    }
  )
)

# Rows of `df` keyed by `var`, in the order of `features`; features it lacks are NA rows.
.feature_rows <- function(df, var, features) {
  df <- prolfquapp::column_to_rownames(df, var = var)[features, , drop = FALSE]
  rownames(df) <- features
  df
}

# `values` on the rows and columns of `mat.raw`; features it lacks are NA.
.align_to_features <- function(values, mat.raw) {
  aligned <- matrix(NA_real_, nrow(mat.raw), ncol(mat.raw), dimnames = dimnames(mat.raw))
  shared <- intersect(rownames(mat.raw), rownames(values))
  aligned[shared, ] <- values[shared, colnames(mat.raw), drop = FALSE]
  aligned
}

# The imputedData assay and its per-feature imputation summary, aligned to the
# raw matrix. Decoys are removed before the fit, so their rows stay NA.
.imputed_assay <- function(facade, mat.raw, rowname, strip) {
  filled <- prolfqua::impute_from_model(facade$model, facade$.lfqdata)
  mat_filled <- prolfquapp::strip_rownames(filled$lfqdata$data_wide(as.matrix = TRUE)$data, strip)
  imputed <- .align_to_features(mat_filled, mat.raw)
  modelled <- intersect(rownames(mat.raw), rownames(mat_filled))
  n_unfilled <- sum(rowSums(is.na(imputed[modelled, , drop = FALSE])) > 0)
  if (n_unfilled > 0) {
    stop("lm_impute left missing values for ", n_unfilled, " feature(s); imputedData must be complete.", call. = FALSE)
  }
  list(matrix = imputed, summary = .feature_rows(filled$summary, rowname, rownames(mat.raw)))
}
