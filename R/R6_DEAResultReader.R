# DEAResultReader ----
#
# Reads back what DEAReportGenerator wrote. The artifact is the single source
# of truth: the abundance layers rebuild the LFQData objects, the nested
# rowData frames rebuild the contrast table, and the serialized
# ContrastConfiguration says which column plays which role -- so nothing here,
# and nothing in the reports, needs to know which modelling backend ran.

#' Read prolfquapp DEA results into prolfqua objects
#'
#' @description
#' Rebuilds the standard prolfqua structures from a DEA result artifact
#' written by [DEAReportGenerator]: the raw and transformed
#' \code{prolfqua::LFQData}, the contrast table with the column names the
#' modelling backend produced, and a \code{prolfqua::ContrastsTable} carrying
#' the backend's \code{prolfqua::ContrastConfiguration}. Accepts a
#' \code{SummarizedExperiment}, an \code{anndataR} AnnData, or a path to an
#' \code{.rds} or \code{.h5ad} file.
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
    #' @field contrast_table all contrasts stacked into one tibble
    contrast_table = NULL,
    #' @field contrasts \code{prolfqua::ContrastsTable} over
    #'   \code{contrast_table}, or \code{NULL} when the artifact holds no
    #'   usable contrast results
    contrasts = NULL,
    #' @field subject_id feature-axis column identifying a contrast subject
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
        stop(
          "SummarizedExperiment is missing required assay(s): ",
          paste(missing, collapse = ", "),
          call. = FALSE
        )
      }

      self$metadata <- S4Vectors::metadata(self$se)
      col_data <- .se_report_col_data(self$se)

      self$lfq_raw <- .se_report_lfqdata_from_assay(
        self$se,
        col_data = col_data,
        assay_name = "rawData",
        response = "abundance",
        metadata_config = self$metadata$analysis_configuration_raw,
        is_transformed = FALSE,
        prefix = "raw_"
      )
      self$lfq_transformed <- .se_report_lfqdata_from_assay(
        self$se,
        col_data = col_data,
        assay_name = "transformedData",
        response = "transformedIntensity",
        metadata_config = self$metadata$analysis_configuration,
        is_transformed = TRUE,
        prefix = "transformed_"
      )

      # The artifact flattens the feature hierarchy into the row names, so the
      # subject of a contrast is that single key, whatever the model used.
      self$subject_id <- self$lfq_transformed$relevant_hierarchy_keys()[[1]]
      self$contrast_config <- prolfqua::list_to_ContrastConfiguration(
        self$metadata$contrast_configuration %||% list()
      )
      self$contrast_config$subject_id <- self$subject_id

      self$contrast_table <- .se_report_contrast_table(
        SummarizedExperiment::rowData(self$se),
        rownames(self$se),
        self$contrast_config
      )
      if (
        self$has_roles(
          self$contrast_config$contrast_col,
          self$contrast_config$effect_col,
          self$contrast_config$fdr_col
        )
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
    stop(
      "Expected a SummarizedExperiment, an AnnData, or a path to an ",
      ".rds or .h5ad file.",
      call. = FALSE
    )
  }
  source
}

.se_report_col_data <- function(se) {
  col_data <- as.data.frame(SummarizedExperiment::colData(se))
  if (!"sampleName" %in% colnames(col_data)) {
    sample_names <- rownames(col_data)
    if (is.null(sample_names) || any(!nzchar(sample_names))) {
      sample_names <- colnames(SummarizedExperiment::assay(se, "rawData"))
    }
    col_data$sampleName <- sample_names
  }
  col_data
}

.se_report_lfqdata_from_assay <- function(
  se,
  col_data,
  assay_name,
  response,
  metadata_config = NULL,
  is_transformed = FALSE,
  prefix = ""
) {
  mat <- SummarizedExperiment::assay(se, assay_name)
  config <- .se_report_config(
    col_data,
    response,
    metadata_config,
    is_transformed
  )
  sample_col <- config$sample_name
  hierarchy_col <- config$hierarchy_keys_depth()[[1]]

  long_data <- as.data.frame(mat, check.names = FALSE)
  long_data[[hierarchy_col]] <- rownames(mat)
  long_data <- tidyr::pivot_longer(
    long_data,
    cols = -dplyr::all_of(hierarchy_col),
    names_to = sample_col,
    values_to = response
  )

  long_data <- .se_report_add_nr_children(
    se,
    long_data,
    config = config,
    hierarchy_col = hierarchy_col,
    sample_col = sample_col
  )

  col_join <- col_data
  col_join[[sample_col]] <- as.character(col_join[[sample_col]])
  long_data[[sample_col]] <- as.character(long_data[[sample_col]])
  long_data <- dplyr::left_join(long_data, col_join, by = sample_col)

  if (is.null(config$file_name) || !nzchar(config$file_name)) {
    config$file_name <- "raw.file"
  }
  if (!config$file_name %in% colnames(long_data)) {
    long_data[[config$file_name]] <- long_data[[sample_col]]
  }
  if (!config$isotope_label %in% colnames(long_data)) {
    long_data[[config$isotope_label]] <- "light"
  }
  if (!config$nr_children %in% colnames(long_data)) {
    long_data[[config$nr_children]] <- 1L
  } else {
    long_data[[config$nr_children]][is.na(long_data[[
      config$nr_children
    ]])] <- 1L
  }
  for (factor_col in config$factor_keys()) {
    if (!factor_col %in% colnames(long_data)) {
      long_data[[factor_col]] <- "all"
    }
  }

  prolfqua::LFQData$new(tibble::as_tibble(long_data), config, prefix = prefix)
}

.se_report_add_nr_children <- function(
  se,
  long_data,
  config,
  hierarchy_col,
  sample_col
) {
  if (!"nr_children" %in% SummarizedExperiment::assayNames(se)) {
    return(long_data)
  }

  nr_children <- SummarizedExperiment::assay(se, "nr_children")
  child_data <- as.data.frame(nr_children, check.names = FALSE)
  child_data[[hierarchy_col]] <- rownames(nr_children)
  child_data <- tidyr::pivot_longer(
    child_data,
    cols = -dplyr::all_of(hierarchy_col),
    names_to = sample_col,
    values_to = config$nr_children
  )
  child_data[[sample_col]] <- as.character(child_data[[sample_col]])
  dplyr::left_join(long_data, child_data, by = c(hierarchy_col, sample_col))
}

.se_report_config <- function(
  col_data,
  response,
  metadata_config = NULL,
  is_transformed = FALSE
) {
  if (!is.null(metadata_config)) {
    config <- prolfqua::list_to_AnalysisConfiguration(metadata_config)
    config$hierarchy <- list("protein_Id" = "protein_Id")
    config$hierarchy_depth <- 1
    config$work_intensity <- character()
    config$set_response(response)
    config$is_response_transformed <- is_transformed
    return(config)
  }

  factor_cols <- .se_report_factor_cols(col_data)
  config <- prolfqua::AnalysisConfiguration$new()
  config$file_name <- .se_report_first_existing(
    c("raw.file", "fileName", "File.Name", "filename", "sampleName"),
    colnames(col_data),
    default = "raw.file"
  )
  config$sample_name <- "sampleName"
  config$isotope_label <- "isotopeLabel"
  config$hierarchy[["protein_Id"]] <- "protein_Id"
  config$hierarchy_depth <- 1
  config$nr_children <- "nr_children"
  config$set_response(response)
  config$is_response_transformed <- is_transformed
  config$factors <- stats::setNames(as.list(factor_cols), factor_cols)
  config$factor_depth <- max(
    1L,
    length(.se_report_primary_factor_cols(factor_cols))
  )
  config
}

.se_report_factor_cols <- function(col_data) {
  exclude <- c(
    "sampleName",
    "raw.file",
    "fileName",
    "File.Name",
    "filename",
    "Name",
    "CONTROL",
    "control",
    "isotopeLabel",
    "nr_children",
    "n_proteins"
  )
  candidates <- setdiff(colnames(col_data), exclude)
  candidates <- candidates[vapply(
    col_data[candidates],
    .se_report_is_factor_like,
    logical(1)
  )]
  if (length(candidates) == 0) {
    col_data$group_ <- "all"
    return("group_")
  }
  primary <- .se_report_primary_factor_cols(candidates)
  unique(c(primary, setdiff(candidates, primary)))
}

.se_report_primary_factor_cols <- function(cols) {
  primary <- grep(
    "group|condition|treatment|background|genotype",
    cols,
    ignore.case = TRUE,
    value = TRUE
  )
  if (length(primary) > 0) {
    return(primary)
  }
  head(cols, 1)
}

.se_report_is_factor_like <- function(x) {
  n <- length(x)
  distinct <- dplyr::n_distinct(x, na.rm = TRUE)
  if (distinct <= 1 || distinct >= n) {
    return(FALSE)
  }
  is.character(x) ||
    is.factor(x) ||
    is.logical(x) ||
    distinct <= min(10, ceiling(n / 2))
}

.se_report_first_existing <- function(candidates, values, default) {
  hit <- candidates[candidates %in% values]
  if (length(hit) > 0) {
    return(hit[[1]])
  }
  default
}

.se_report_contrast_table <- function(row_data, feature_ids, contrast_config) {
  contrast_names <- grep("^constrast_", colnames(row_data), value = TRUE)
  if (length(contrast_names) == 0) {
    return(tibble::tibble())
  }
  subject_col <- contrast_config$subject_id
  contrast_col <- contrast_config$contrast_col

  res <- lapply(contrast_names, function(name) {
    df <- .se_report_row_data_frame(row_data[[name]], feature_ids, subject_col)
    if (!subject_col %in% colnames(df)) {
      df[[subject_col]] <- feature_ids
    }
    if (!contrast_col %in% colnames(df)) {
      df[[contrast_col]] <- sub("^constrast_", "", name)
    }
    if (!contrast_config$model_name_col %in% colnames(df)) {
      df[[contrast_config$model_name_col]] <- "SummarizedExperiment"
    }
    # Preserve the rescue-state column when round-tripping; legacy rows that
    # predate estimate_type are treated as observed.
    if (!"estimate_type" %in% colnames(df)) {
      df$estimate_type <- "observed"
    }
    # Contrast frames are padded to the feature axis, so rows without an
    # estimate carry no subject and no contrast label.
    dplyr::filter(
      df,
      !is.na(.data[[subject_col]]),
      !is.na(.data[[contrast_col]])
    )
  })
  dplyr::bind_rows(res)
}

.se_report_row_data_frame <- function(value, feature_ids, subject_col) {
  if (is.null(value)) {
    return(tibble::tibble())
  }
  df <- as.data.frame(value)
  if (nrow(df) == length(feature_ids) && !subject_col %in% colnames(df)) {
    df[[subject_col]] <- feature_ids
  }
  tibble::as_tibble(df)
}
