# Internal helpers for SummarizedExperiment-backed Quarto reports.

se_report_lfqdata <- function(se) {
  if (is.character(se) && length(se) == 1) {
    se <- readRDS(se)
  }
  if (!inherits(se, "SummarizedExperiment")) {
    stop("Expected a SummarizedExperiment object or path to an .rds file.", call. = FALSE)
  }

  assay_names <- SummarizedExperiment::assayNames(se)
  required <- c("rawData", "transformedData")
  missing <- setdiff(required, assay_names)
  if (length(missing) > 0) {
    stop(
      "SummarizedExperiment is missing required assay(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  meta <- S4Vectors::metadata(se)
  col_data <- .se_report_col_data(se)
  row_data <- SummarizedExperiment::rowData(se)
  contrast_table <- .se_report_contrast_table(row_data, rownames(se))

  lfq_raw <- .se_report_lfqdata_from_assay(
    se,
    col_data = col_data,
    assay_name = "rawData",
    response = "abundance",
    metadata_config = meta$analysis_configuration_raw,
    is_transformed = FALSE,
    prefix = "raw_"
  )
  lfq_transformed <- .se_report_lfqdata_from_assay(
    se,
    col_data = col_data,
    assay_name = "transformedData",
    response = "transformedIntensity",
    metadata_config = meta$analysis_configuration_transformed,
    is_transformed = TRUE,
    prefix = "transformed_"
  )

  contrast_object <- NULL
  if (nrow(contrast_table) > 0) {
    contrast_object <- prolfqua::ContrastsTable$new(
      contrast_table,
      subject_id = lfq_transformed$relevant_hierarchy_keys(),
      model_name = "SummarizedExperiment"
    )
  }

  list(
    se = se,
    metadata = meta,
    col_data = col_data,
    row_data = row_data,
    lfq_raw = lfq_raw,
    lfq_transformed = lfq_transformed,
    contrast_table = contrast_table,
    contrast_object = contrast_object,
    stats_raw = .se_report_row_data_frame(row_data[["stats_raw_wide"]], rownames(se)),
    stats_normalized = .se_report_row_data_frame(row_data[["stats_normalized_wide"]], rownames(se))
  )
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
  config <- .se_report_config(col_data, response, metadata_config, is_transformed)
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
  }
  for (factor_col in config$factor_keys()) {
    if (!factor_col %in% colnames(long_data)) {
      long_data[[factor_col]] <- "all"
    }
  }

  prolfqua::LFQData$new(tibble::as_tibble(long_data), config, prefix = prefix)
}

.se_report_config <- function(col_data, response, metadata_config = NULL, is_transformed = FALSE) {
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
  config$factor_depth <- max(1L, length(.se_report_primary_factor_cols(factor_cols)))
  config
}

.se_report_factor_cols <- function(col_data) {
  exclude <- c(
    "sampleName", "raw.file", "fileName", "File.Name", "filename", "Name",
    "CONTROL", "control", "isotopeLabel", "nr_children", "n_proteins"
  )
  candidates <- setdiff(colnames(col_data), exclude)
  candidates <- candidates[vapply(col_data[candidates], .se_report_is_factor_like, logical(1))]
  if (length(candidates) == 0) {
    col_data$group_ <- "all"
    return("group_")
  }
  primary <- .se_report_primary_factor_cols(candidates)
  unique(c(primary, setdiff(candidates, primary)))
}

.se_report_primary_factor_cols <- function(cols) {
  primary <- grep("group|condition|treatment|background|genotype", cols, ignore.case = TRUE, value = TRUE)
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
  is.character(x) || is.factor(x) || is.logical(x) || distinct <= min(10, ceiling(n / 2))
}

.se_report_first_existing <- function(candidates, values, default) {
  hit <- candidates[candidates %in% values]
  if (length(hit) > 0) {
    return(hit[[1]])
  }
  default
}

.se_report_contrast_table <- function(row_data, feature_ids) {
  contrast_names <- grep("^constrast_", colnames(row_data), value = TRUE)
  if (length(contrast_names) == 0) {
    return(tibble::tibble())
  }

  res <- lapply(contrast_names, function(name) {
    df <- .se_report_row_data_frame(row_data[[name]], feature_ids)
    if (!"protein_Id" %in% colnames(df)) {
      df$protein_Id <- feature_ids
    }
    if (!"contrast" %in% colnames(df)) {
      df$contrast <- sub("^constrast_", "", name)
    }
    if (!"modelName" %in% colnames(df)) {
      df$modelName <- "SummarizedExperiment"
    }
    df
  })
  dplyr::bind_rows(res)
}

.se_report_row_data_frame <- function(value, feature_ids) {
  if (is.null(value)) {
    return(tibble::tibble())
  }
  df <- as.data.frame(value)
  if (nrow(df) == length(feature_ids) && !"protein_Id" %in% colnames(df)) {
    df$protein_Id <- feature_ids
  }
  tibble::as_tibble(df)
}
