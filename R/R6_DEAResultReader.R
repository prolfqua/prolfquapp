# DEAResultReader ----
#
# Reads back what DEAReportGenerator wrote. The artifact is the single source
# of truth: the abundance layers rebuild the LFQData objects, the nested
# rowData frames rebuild the contrast table, and the serialized
# ContrastConfiguration says which column plays which role -- so nothing here,
# and nothing in the reports, needs to know which modelling backend ran.
#
# Every table is keyed by the artifact's feature keys. The rowData frames carry
# them as columns; an abundance matrix names its rows by feature id only, so
# its long form takes the keys from the annotation frame, joined on that id.

#' Read prolfquapp DEA results into prolfqua objects
#'
#' @description
#' Rebuilds the standard prolfqua structures from a DEA result artifact
#' written by [DEAReportGenerator]: the raw, transformed and (for
#' \code{lm_impute}) imputed \code{prolfqua::LFQData}, the feature and sample
#' annotation, the per-feature imputation summary, the contrast table with the
#' column names the modelling backend produced, and a
#' \code{prolfqua::ContrastsTable} carrying the backend's
#' \code{prolfqua::ContrastConfiguration}. Every table is keyed by the
#' artifact's \code{feature_keys}. Accepts a \code{SummarizedExperiment}, an
#' \code{anndataR} AnnData, or a path to an \code{.rds} or \code{.h5ad} file.
#'
#' Because the column roles travel with the artifact, consumers filter, rank
#' and plot through the configuration instead of hard-coding column names.
#'
#' @export
#' @examples
#' dea <- prolfquapp::example_deanalyse(Nprot = 12)
#' se <- prolfquapp::DEAReportGenerator$new(dea, dea$prolfq_app_config)$make_SummarizedExperiment()
#' reader <- prolfquapp::DEAResultReader$new(se)
#' reader$contrast_config$effect_col
#' nrow(reader$contrast_table)
#' nrow(reader$significant(FDR_threshold = 0.25, diff_threshold = 0.5))
DEAResultReader <- R6::R6Class(
  "DEAResultReader",
  public = list(
    #' @field se the \code{SummarizedExperiment} the results were read from
    se = NULL,
    #' @field metadata the artifact metadata (\code{S4Vectors::metadata})
    metadata = NULL,
    #' @field contrast_config \code{prolfqua::ContrastConfiguration} of the
    #'   backend that produced the results
    contrast_config = NULL,
    #' @field lfq_raw \code{prolfqua::LFQData} with the raw abundances
    lfq_raw = NULL,
    #' @field lfq_transformed \code{prolfqua::LFQData} with the transformed
    #'   abundances
    lfq_transformed = NULL,
    #' @field lfq_imputed \code{prolfqua::LFQData} with the transformed
    #'   abundances and every missing value filled by the feature's model, or
    #'   \code{NULL} when the artifact has no \code{imputedData}
    lfq_imputed = NULL,
    #' @field annotation feature annotation, one row per feature
    annotation = NULL,
    #' @field samples sample annotation, one row per sample
    samples = NULL,
    #' @field imputation per feature, how many values were observed and
    #'   imputed and by which route, or \code{NULL} when the artifact has no
    #'   imputation block
    imputation = NULL,
    #' @field contrast_table all contrasts stacked into one tibble
    contrast_table = NULL,
    #' @field contrasts \code{prolfqua::ContrastsTable} over
    #'   \code{contrast_table}, or \code{NULL} when the artifact holds no
    #'   usable contrast results
    contrasts = NULL,
    #' @field subject_id the feature-key columns identifying a contrast subject
    subject_id = NULL,
    #' @description
    #' Read a DEA result artifact.
    #' @param source a \code{SummarizedExperiment}, an AnnData, or a path to
    #'   an \code{.rds} or \code{.h5ad} file
    initialize = function(source) {
      self$se <- .dea_result_source_to_se(source)
      assay_names <- SummarizedExperiment::assayNames(self$se)
      missing <- setdiff(c("rawData", "transformedData"), assay_names)
      if (length(missing) > 0) {
        stop("SummarizedExperiment is missing required assay(s): ", paste(missing, collapse = ", "), call. = FALSE)
      }
      self$metadata <- S4Vectors::metadata(self$se)
      row_data <- SummarizedExperiment::rowData(self$se)
      self$subject_id <- as.character(unlist(self$metadata$feature_keys, use.names = FALSE))
      annotation <- as.data.frame(row_data[["annotation"]])
      self$annotation <- tibble::as_tibble(annotation)
      self$samples <- tibble::as_tibble(as.data.frame(SummarizedExperiment::colData(self$se)))
      feature_ids <- tibble::as_tibble(annotation[self$subject_id], rownames = ".feature")

      lfq <- function(assay_name, response, config, prefix) {
        .se_report_lfqdata_from_assay(self$se, feature_ids, assay_name, response, config, prefix)
      }
      self$lfq_raw <- lfq("rawData", "abundance", self$metadata$analysis_configuration_raw, "raw_")
      self$lfq_transformed <- lfq(
        "transformedData",
        "transformedIntensity",
        self$metadata$analysis_configuration,
        "transformed_"
      )
      if ("imputedData" %in% assay_names) {
        self$lfq_imputed <- lfq("imputedData", "transformedIntensity", self$metadata$analysis_configuration, "imputed_")
      }
      if (!is.null(row_data[["imputation"]])) {
        self$imputation <- .se_report_feature_rows(row_data[["imputation"]], self$subject_id)
      }

      self$contrast_config <- prolfqua::list_to_ContrastConfiguration(self$metadata$contrast_configuration %||% list())
      self$contrast_config$subject_id <- self$subject_id
      self$contrast_table <- .se_report_contrast_table(row_data, self$contrast_config)
      if (
        self$has_roles(self$contrast_config$contrast_col, self$contrast_config$effect_col, self$contrast_config$fdr_col)
      ) {
        self$contrasts <- prolfqua::ContrastsTable$new(
          self$contrast_table,
          subject_id = self$subject_id,
          model_name = self$metadata$default_model %||% "SummarizedExperiment"
        )
        self$contrasts$config <- self$contrast_config
      }
    },
    #' @description
    #' Does the contrast table carry these columns?
    #' @param ... column names, or vectors of column names
    has_roles = function(...) {
      needed <- unlist(list(...), use.names = FALSE)
      nrow(self$contrast_table) > 0 &&
        all(needed %in% colnames(self$contrast_table))
    },
    #' @description
    #' Contrast rows passing both thresholds. One-sided on the effect column
    #' for backends whose negative effects are not interpretable (SAINT).
    #' @param FDR_threshold false discovery rate threshold
    #' @param diff_threshold effect-size threshold
    significant = function(FDR_threshold = 0.05, diff_threshold = 1) {
      if (is.null(self$contrasts)) {
        return(self$contrast_table[0, , drop = FALSE])
      }
      tibble::as_tibble(
        self$contrasts$filter_significant(FDR_threshold, diff_threshold)
      )
    },
    #' @description
    #' \code{prolfqua::ContrastsPlotter} over the contrast table, with every
    #' score panel resolved through the column roles. Returns \code{NULL}
    #' when the artifact holds no usable contrast results.
    #' @param fc_threshold fold-change threshold
    #' @param fdr_threshold false discovery rate threshold
    get_Plotter = function(fc_threshold = 1, fdr_threshold = 0.1) {
      if (is.null(self$contrasts)) {
        return(NULL)
      }
      cfg <- self$contrast_config
      volcano <- list(list(score = cfg$fdr_col, thresh = fdr_threshold))
      histogram <- list(list(score = cfg$fdr_col, xlim = c(0, 1, 0.05)))
      score <- list()
      if (cfg$has_pvalue() && self$has_roles(cfg$pvalue_col)) {
        volcano <- c(list(list(score = cfg$pvalue_col)), volcano)
        histogram <- c(
          list(list(score = cfg$pvalue_col, xlim = c(0, 1, 0.05))),
          histogram
        )
      }
      if (self$has_roles(cfg$score_col)) {
        if (cfg$has_pvalue()) {
          score <- list(list(score = cfg$score_col, thresh = NULL))
        } else {
          # A backend without a p-value reports a bounded probability score,
          # so it gets its own histogram and a score panel of its own.
          histogram <- c(
            histogram,
            list(list(score = cfg$score_col, xlim = c(0, 1, 0.05)))
          )
          score <- list(list(score = cfg$score_col, thresh = 0.75))
        }
      }
      prolfqua::ContrastsPlotter$new(
        self$contrast_table,
        subject_id = self$subject_id,
        fcthresh = fc_threshold,
        volcano = volcano,
        histogram = histogram,
        score = score,
        modelName = cfg$model_name_col,
        diff = cfg$effect_col,
        contrast = cfg$contrast_col,
        avg.abundance = cfg$avg_abundance_col
      )
    }
  )
)

.dea_result_source_to_se <- function(source) {
  if (is.character(source) && length(source) == 1) {
    if (grepl("[.]h5ad$", source, ignore.case = TRUE)) {
      return(anndata_to_summarized_experiment(anndataR::read_h5ad(source)))
    }
    return(.dea_result_source_to_se(readRDS(source)))
  }
  if (inherits(source, "AbstractAnnData")) {
    return(anndata_to_summarized_experiment(source))
  }
  if (!inherits(source, "SummarizedExperiment")) {
    stop("Expected a SummarizedExperiment, an AnnData, or a path to an .rds or .h5ad file.", call. = FALSE)
  }
  source
}

# One assay in long form, keyed by the feature keys, as prolfqua::LFQData.
# `config` is the artifact's serialized AnalysisConfiguration; its hierarchy is
# replaced by the feature keys, since the assay holds one row per feature.
.se_report_lfqdata_from_assay <- function(se, feature_ids, assay_name, response, config, prefix) {
  config <- prolfqua::list_to_AnalysisConfiguration(config)
  feature_keys <- setdiff(colnames(feature_ids), ".feature")
  config$hierarchy <- stats::setNames(as.list(feature_keys), feature_keys)
  config$hierarchy_depth <- length(feature_keys)
  config$work_intensity <- character()
  config$set_response(response)
  config$is_response_transformed <- !identical(assay_name, "rawData")
  sample_col <- config$sample_name

  long <- function(name, values_to) {
    tidyr::pivot_longer(
      tibble::as_tibble(SummarizedExperiment::assay(se, name), rownames = ".feature"),
      cols = -".feature",
      names_to = sample_col,
      values_to = values_to
    )
  }
  long_data <- long(assay_name, response)
  if ("nr_children" %in% SummarizedExperiment::assayNames(se)) {
    long_data <- dplyr::left_join(long_data, long("nr_children", config$nr_children), by = c(".feature", sample_col))
  }
  long_data <- dplyr::left_join(long_data, feature_ids, by = ".feature")
  long_data$.feature <- NULL
  col_data <- as.data.frame(SummarizedExperiment::colData(se))
  col_data[[sample_col]] <- as.character(col_data[[sample_col]])
  long_data <- dplyr::left_join(long_data, col_data, by = sample_col)

  if (is.null(config$file_name) || !nzchar(config$file_name)) {
    config$file_name <- "raw.file"
  }
  long_data[[config$file_name]] <- long_data[[config$file_name]] %||% long_data[[sample_col]]
  long_data[[config$isotope_label]] <- long_data[[config$isotope_label]] %||% "light"
  long_data[[config$nr_children]] <- dplyr::coalesce(long_data[[config$nr_children]] %||% NA_integer_, 1L)
  for (factor_col in setdiff(config$factor_keys(), colnames(long_data))) {
    long_data[[factor_col]] <- "all"
  }
  prolfqua::LFQData$new(tibble::as_tibble(long_data), config, prefix = prefix)
}

.se_report_contrast_table <- function(row_data, contrast_config) {
  contrast_names <- grep("^constrast_", colnames(row_data), value = TRUE)
  res <- lapply(contrast_names, function(name) {
    df <- .se_report_feature_rows(row_data[[name]], contrast_config$subject_id)
    # Rows that predate estimate_type are observed estimates.
    defaults <- stats::setNames(
      list(sub("^constrast_", "", name), "SummarizedExperiment", "observed"),
      c(contrast_config$contrast_col, contrast_config$model_name_col, "estimate_type")
    )
    for (col in setdiff(names(defaults), colnames(df))) {
      df[[col]] <- defaults[[col]]
    }
    dplyr::filter(df, !is.na(.data[[contrast_config$contrast_col]]))
  })
  if (length(res) == 0) tibble::tibble() else dplyr::bind_rows(res)
}

# rowData frames are padded to the feature axis: a feature without a result
# has NA in every column, its ids included.
.se_report_feature_rows <- function(frame, feature_keys) {
  dplyr::filter(
    tibble::as_tibble(as.data.frame(frame)),
    dplyr::if_all(dplyr::all_of(feature_keys), ~ !is.na(.x))
  )
}
