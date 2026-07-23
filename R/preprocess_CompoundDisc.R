#' massage CD output compound table.
#' @param in_file path to Compound Discoverer file or data frame
#' @param EXCEL if TRUE read Excel format
#' @export
massage_CD <- function(in_file, EXCEL = TRUE) {
  xd <- if (is.character(in_file) && file.exists(in_file)) {
    readxl::read_excel(in_file)
  } else if (is.data.frame(in_file)) {
    in_file
  } else {
    stopifnot("expecting data frame or path got : ", class(in_file))
  }
  xd$my_C_ID <- seq_len(nrow(xd))

  if (EXCEL) {
    annot <- xd |>
      dplyr::select(
        "my_C_ID",
        "Checked",
        "Tags",
        "description" = "Name",
        "Formula",
        "Annot. Source: Predicted Compositions",
        "Annot. Source: mzCloud Search",
        "Annot. Source: mzVault Search",
        "Annot. Source: ChemSpider Search",
        "Annot. Source: MassList Search",
        "Calc. MW",
        mz = "m/z",
        RT_min = "RT [min]"
      )
    columns <- c("Area:", "Gap Status:", "Gap Fill Status:", "Peak Rating:")
    deselect <- NULL
    npatt <- "(.*)\\: (.*)(\\s\\(F\\d+\\))"
  } else {
    annot <- xd |>
      dplyr::select(dplyr::any_of(c(
        "my_C_ID",
        "Checked",
        "Tags",
        "description" = "Name",
        "Formula",
        "Annot Source Predicted Compositions",
        "Annot Source mzCloud Search",
        "Annot Source mzVault Search",
        "Annot Source ChemSpider Search",
        "Annot Source MassList Search",
        "Calc MW",
        "mz",
        RT_min = "RT in min"
      )))
    columns <- c("Area", "Gap Status", "Gap Fill Status", "Peak Rating")
    deselect <- c("Area Max", "Area SD", "Area CV in Percent")
    npatt <- "(.*)\\s(.*)\\s(F\\d+)"
  }
  colnames(annot) <- gsub("[[:space:].:/]+", "_", colnames(annot))
  colnames(annot) <- gsub("\\[|\\]", "", colnames(annot))
  annot <- annot |>
    dplyr::mutate(FormulaB = stringr::str_replace_all(Formula, " ", ""))

  annot <- annot |>
    tidyr::unite(
      "metabolite_feature_Id",
      c("my_C_ID", "FormulaB", "mz", "RT_min"),
      sep = "_",
      remove = FALSE
    )

  tolong <- xd |> dplyr::select("my_C_ID", tidyselect::starts_with(columns))
  if (!is.null(deselect)) {
    tolong <- dplyr::select(tolong, -all_of(deselect))
  }
  sum(grepl(columns[1], colnames(tolong)))
  sum(grepl(columns[2], colnames(tolong)))
  sum(grepl(columns[3], colnames(tolong)))
  sum(grepl(columns[4], colnames(tolong)))

  xdl <- tolong |>
    tidyr::pivot_longer(
      cols = tidyselect::starts_with(columns),
      names_to = c(".value", "filename", "file_id"),
      names_pattern = npatt
    )
  xdl$`Gap Fill Status` |> table()

  colnames(xdl) <- gsub("# ", "", colnames(xdl))
  colnames(xdl) <- gsub("[[:space:]]", "_", colnames(xdl))
  xdl <- xdl |>
    dplyr::mutate(file_id = gsub(" ", "", gsub("\\(|\\)", "", file_id)))

  # use nr_children to encode gap status.
  xdl <- xdl |>
    dplyr::mutate(
      nr_children = dplyr::case_when(
        Gap_Status == "Full gap" ~ 0,
        Gap_Status == "Missing ions" ~ 1,
        Gap_Status == "No gap" ~ 2,
        TRUE ~ 3
      )
    )

  xdl <- dplyr::inner_join(annot, xdl, by = "my_C_ID")
  return(xdl)
}

#' load compound discoverer (CD) files
#' @param xdl data frame from massage_CD
#' @param annotation list returned by `read_annotation` function
#' @return A list containing the prepared \code{LFQData} and
#'   \code{ProteinAnnotation} objects.
#' @export
preprocess_CD <- function(
  xdl,
  annotation
) {
  annot <- annotation$annot
  nr <- sum(annot$File.Name %in% sort(unique(xdl$filename)))
  logger::log_info(
    "nr : ",
    nr,
    " files annotated out of ",
    length(unique(xdl$filename))
  )
  stopifnot(nr > 0)

  config <- annotation$atable$clone(deep = TRUE)
  config$hierarchy[["metabolite_feature_Id"]] <- c("metabolite_feature_Id")
  config$set_response("Area")
  config$opt_rt <- "RT_min"
  config$opt_mz <- "mz"
  byv <- c("filename")
  names(byv) <- config$file_name
  byv <- c(byv, intersect(colnames(annot), colnames(xdl)))

  peptide <- dplyr::inner_join(annot, xdl, by = byv, multiple = "all")
  adata <- prolfqua::setup_analysis(peptide, config)
  lfqdata <- prolfqua::LFQData$new(adata, config)
  lfqdata$remove_small_intensities()

  m_annot <- xdl |>
    dplyr::select(
      "metabolite_feature_Id",
      "description",
      "Formula",
      starts_with("Annot_")
    ) |>
    dplyr::distinct()

  m_annot <- m_annot |> dplyr::mutate(IDcolumn = metabolite_feature_Id)
  prot_annot <- prolfquapp::ProteinAnnotation$new(
    lfqdata,
    m_annot,
    description = "description",
    cleaned_ids = "IDcolumn",
    full_id = "IDcolumn",
    exp_nr_children = "nr_compounds",
    pattern_contaminants = NULL,
    pattern_decoys = NULL
  )
  return(list(lfqdata = lfqdata, protein_annotation = prot_annot))
}

#' Locate files inside a CompoundDiscoverer prolfqua ZIP export
#' @param input ZIP file path or directory containing one ZIP export
#' @return list with zip, data, samples, and tempdir entries
#' @export
get_CD_export_files <- function(input) {
  if (!file.exists(input)) {
    stop("CompoundDiscoverer input not found: ", input, call. = FALSE)
  }

  zipfile <- input
  if (dir.exists(input)) {
    zips <- list.files(
      input,
      pattern = "[.]zip$",
      full.names = TRUE,
      ignore.case = TRUE
    )
    if (length(zips) != 1) {
      stop(
        "Expected exactly one CompoundDiscoverer ZIP in directory '",
        input,
        "', found ",
        length(zips),
        ".",
        call. = FALSE
      )
    }
    zipfile <- zips[[1]]
  }

  tempdir_export <- tempfile("compound_discoverer_")
  dir.create(tempdir_export, recursive = TRUE)
  utils::unzip(zipfile, exdir = tempdir_export)

  long_file <- list.files(
    tempdir_export,
    pattern = "_long[.]csv$",
    full.names = TRUE,
    recursive = TRUE
  )
  sample_file <- list.files(
    tempdir_export,
    pattern = "_prolfqua_samples[.]csv$",
    full.names = TRUE,
    recursive = TRUE
  )

  if (length(long_file) != 1 || length(sample_file) != 1) {
    unlink(tempdir_export, recursive = TRUE)
    stop(
      "CompoundDiscoverer ZIP must contain exactly one *_long.csv and one ",
      "*_prolfqua_samples.csv. Found long=",
      length(long_file),
      ", samples=",
      length(sample_file),
      ".",
      call. = FALSE
    )
  }

  list(
    zip = zipfile,
    data = long_file[[1]],
    samples = sample_file[[1]],
    tempdir = tempdir_export
  )
}

#' Sanitize a CompoundDiscoverer subset name for file paths
#' @param x subset name
#' @return filesystem-safe subset name
#' @export
sanitize_CD_subset_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nchar(x) == 0, "subset", x)
}

#' Detect subset columns in a CompoundDiscoverer long export
#' @param long_file path to *_long.csv
#' @return character vector of columns after Group
#' @export
get_CD_subset_columns <- function(long_file) {
  header <- readr::read_csv(long_file, n_max = 0, show_col_types = FALSE)
  columns <- colnames(header)
  group_pos <- match("Group", columns)
  if (is.na(group_pos) || group_pos >= length(columns)) {
    return(character())
  }
  columns[(group_pos + 1):length(columns)]
}

.is_cd_subset_selected <- function(x) {
  x <- trimws(as.character(x))
  !is.na(x) &
    x != "" &
    !tolower(x) %in% c("false", "f", "no", "n", "0", "na", "nan")
}

#' Make duplicated CompoundDiscoverer feature/sample rows unique
#' @param data long CompoundDiscoverer data frame
#' @param feature_col feature identifier column
#' @param sample_col sample identifier column
#' @return data frame with unique feature identifiers for affected rows
#' @export
make_CD_duplicate_features_unique <- function(
  data,
  feature_col = "Feature_ID",
  sample_col = "Sample"
) {
  stopifnot(feature_col %in% colnames(data))
  stopifnot(sample_col %in% colnames(data))

  data$.cd_row_order <- seq_len(nrow(data))
  duplicate_features <- data |>
    dplyr::count(
      dplyr::across(dplyr::all_of(c(feature_col, sample_col))),
      name = ".cd_n"
    ) |>
    dplyr::filter(.data$.cd_n > 1) |>
    dplyr::distinct(dplyr::across(dplyr::all_of(feature_col))) |>
    dplyr::pull(dplyr::all_of(feature_col))

  if (length(duplicate_features) == 0) {
    data$.cd_row_order <- NULL
    return(data)
  }

  logger::log_info(
    "Making duplicated CompoundDiscoverer Feature_ID x Sample rows unique for ",
    length(duplicate_features),
    " feature(s)."
  )

  data <- data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(feature_col, sample_col)))) |>
    dplyr::arrange(.data$.cd_row_order, .by_group = TRUE) |>
    dplyr::mutate(.cd_duplicate_rank = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::arrange(.data$.cd_row_order)

  original_col <- paste0(feature_col, "_original")
  data[[original_col]] <- data[[feature_col]]
  is_duplicated_feature <- data[[feature_col]] %in% duplicate_features
  data[[feature_col]][is_duplicated_feature] <- paste0(
    data[[feature_col]][is_duplicated_feature],
    " [duplicate ",
    data$.cd_duplicate_rank[is_duplicated_feature],
    "]"
  )

  data$.cd_row_order <- NULL
  data
}

.read_cd_annotation <- function(sample_file, config) {
  samples <- prolfquapp::read_table_data(sample_file)
  if (!"Sample" %in% colnames(samples)) {
    stop(
      "CD sample file must contain a Sample column.",
      call. = FALSE
    )
  }
  samples |>
    dplyr::rename(file = "Sample") |>
    prolfquapp::read_annotation(prefix = config$group)
}

.read_cd_long <- function(long_file) {
  data <- readr::read_csv(long_file, show_col_types = FALSE)
  required <- c("Feature_ID", "Sample", "Intensity")
  missing <- setdiff(required, colnames(data))
  if (length(missing) > 0) {
    stop(
      "CD long file is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  data
}

.subset_cd_long <- function(data, subset_column) {
  if (is.null(subset_column) || nchar(subset_column) == 0) {
    return(data)
  }
  if (!subset_column %in% colnames(data)) {
    stop("CD subset column not found: ", subset_column, call. = FALSE)
  }

  before <- nrow(data)
  data <- data[
    .is_cd_subset_selected(data[[subset_column]]),
    ,
    drop = FALSE
  ]
  logger::log_info(
    "CD subset '",
    subset_column,
    "': keeping ",
    nrow(data),
    " of ",
    before,
    " long-format rows."
  )
  if (nrow(data) == 0) {
    stop(
      "CD subset column has no selected rows: ",
      subset_column,
      call. = FALSE
    )
  }
  data
}

.prepare_cd_features <- function(data) {
  data <- prolfquapp::make_CD_duplicate_features_unique(data)
  data <- dplyr::rename(
    data,
    metabolite_feature_Id = "Feature_ID"
  )
  feature_source <- if ("Feature_ID_original" %in% colnames(data)) {
    data$Feature_ID_original
  } else {
    data$metabolite_feature_Id
  }
  parsed_mz_rt <- stringr::str_match(
    feature_source,
    "_mz([0-9.]+)_rt([0-9.]+)"
  )
  if (!"mz" %in% colnames(data)) {
    data$mz <- as.numeric(parsed_mz_rt[, 2])
  }
  if (!"RT_min" %in% colnames(data)) {
    data$RT_min <- as.numeric(parsed_mz_rt[, 3])
  }
  data
}

.configure_cd_analysis <- function(annotation) {
  config <- annotation$atable$clone(deep = TRUE)
  config$hierarchy[["metabolite_feature_Id"]] <-
    "metabolite_feature_Id"
  config$set_response("Intensity")
  config$hierarchy_depth <- 1
  config$opt_rt <- "RT_min"
  config$opt_mz <- "mz"
  config
}

.join_cd_features <- function(annotation, data, config) {
  by <- "Sample"
  names(by) <- config$file_name
  by <- c(
    by,
    intersect(colnames(annotation$annot), colnames(data))
  )
  feature_data <- dplyr::inner_join(
    annotation$annot,
    data,
    by = by,
    multiple = "all"
  )
  if (nrow(feature_data) == 0) {
    stop(
      "No CD feature rows matched the sample annotation.",
      call. = FALSE
    )
  }
  if (!"isotopeLabel" %in% colnames(feature_data)) {
    feature_data$isotopeLabel <- "light"
  }
  if (!"qValue" %in% colnames(feature_data)) {
    feature_data$qValue <- 0
  }
  if (!"nr_children" %in% colnames(feature_data)) {
    feature_data$nr_children <- 1
  }
  feature_data
}

.cd_feature_annotation <- function(data) {
  data |>
    dplyr::select(
      "metabolite_feature_Id",
      dplyr::any_of(
        c("Feature_ID_original", ".cd_duplicate_rank")
      )
    ) |>
    dplyr::distinct() |>
    dplyr::mutate(
      description = .data$metabolite_feature_Id,
      IDcolumn = .data$metabolite_feature_Id,
      exp_children = 1,
      nrPeptides = 1,
      protein_length = 1,
      nr_tryptic_peptides = 1
    )
}

#' Preprocess a CompoundDiscoverer prolfqua ZIP export
#' @param long_file path to *_long.csv
#' @param sample_file path to *_prolfqua_samples.csv
#' @param config ProlfquAppConfig object
#' @param subset_column optional column in long_file marking features to keep
#' @return list with lfqdata, protein_annotation, and annotation
#' @export
preprocess_CD_export <- function(long_file, sample_file, config, subset_column = NULL) {
  annotation <- .read_cd_annotation(sample_file, config)
  xdl <- .read_cd_long(long_file) |>
    .subset_cd_long(subset_column) |>
    .prepare_cd_features()
  analysis_config <- .configure_cd_analysis(annotation)
  feature_data <- .join_cd_features(
    annotation,
    xdl,
    analysis_config
  )
  adata <- prolfqua::setup_analysis(feature_data, analysis_config)
  lfqdata <- prolfqua::LFQData$new(adata, analysis_config)
  lfqdata$remove_small_intensities()
  prot_annot <- prolfquapp::ProteinAnnotation$new(
    lfqdata,
    .cd_feature_annotation(xdl),
    description = "description",
    cleaned_ids = "IDcolumn",
    full_id = "metabolite_feature_Id",
    exp_nr_children = "exp_children",
    pattern_contaminants = NULL,
    pattern_decoys = NULL
  )

  list(
    lfqdata = lfqdata,
    protein_annotation = prot_annot,
    annotation = annotation
  )
}
