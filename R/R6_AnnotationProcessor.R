#' read dataset file in csv, tsv or xlsx format
#' @param file_path path to csv, tsv, or xlsx file
#' @return A data frame read from \code{file_path}.
#' @export
read_table_data <- function(file_path) {
  reader <- list(
    csv = readr::read_csv,
    tsv = readr::read_tsv,
    xlsx = readxl::read_xlsx
  )[[tools::file_ext(file_path)]]
  if (is.null(reader)) {
    stop("Unsupported file extension")
  }
  reader(file_path)
}

#' Write dataset to file in csv, tsv, or xlsx format
#' @param data data frame to write
#' @param file_path output file path (csv, tsv, or xlsx)
#' @export
#' @examples
#'
#' ds <- data.frame(channel = c("A","B","C"), Name = NA, Subject = NA, Group = NA, Control = NA)
#' write_annotation_file(ds, file_path = file.path(tempdir(),"test.xlsx"))
#'
write_annotation_file <- function(data, file_path) {
  writer <- list(
    csv = readr::write_csv,
    tsv = readr::write_tsv,
    xlsx = writexl::write_xlsx
  )[[tools::file_ext(file_path)]]
  if (is.null(writer)) {
    stop("Unsupported file extension")
  }
  writer(data, file_path)
}

# AnnotationProcessor  -----
#' AnnotationProcessor
#' @export
#' @examples
#'
#' # AnnotationProcessor$debug("read_annotation")
#' ap <- AnnotationProcessor$new(prefix = "G_")
#'
#' annot <- data.frame(
#' file = c("a1.raw","a2.raw","a3.raw","a4.raw"),
#' group = c("a","a","b","b"),
#' CONTROL = c("C","C","T","T"),
#' Subject = c("X","Y","X","Y"))
#' ap$check_annotation(annot)
#' af <- annot
#' af$file <- NULL
#' testthat::expect_error(ap$check_annotation(af), "column starting with :")
#' af <- annot
#' af$group <- NULL
#' testthat::expect_error(ap$check_annotation(af),"column starting with :")
#' aa <- ap$read_annotation(annot)
#' stopifnot(length(aa$atable$factor_keys_depth()) == 2)
#' stopifnot(all(c("atable", "annot", "contrasts") %in% names(aa)))
#' stopifnot(aa$contrasts == "G_b - G_a")
#' af <- annot
#' af$CONTROL <- NULL
#' testthat::expect_error(ap$check_annotation(af),"you must specify a CONTROL column")
#' af <- annot
#' af$Subject <- NULL
#' testthat::expect_warning(ap$check_annotation(af),"column starting with")
#'
#'
#' # should not throw exception since QC does not require group or subject
#' ap <- AnnotationProcessor$new(QC = TRUE)
#' af <- annot
#' af$group <- NULL
#' af$CONTROL <- NULL
#' af$Subject <- NULL
#' ap$check_annotation(af)
#' aa <- ap$read_annotation(af)
#'
#' stopifnot(aa$atable$factor_keys() == "G_")
#' stopifnot(aa$atable$factors == "group")
#' aa <- ap$read_annotation(annot)
#' aa$atable$file_name
#' aa$atable$sample_name
#' stopifnot(is.null(aa$annotation))
#'
#' annot <- data.frame(
#' file = c("a1.raw","a2.raw","a3.raw","a4.raw"),
#' Name = c("a1.raw","a2.raw","a3.raw","a4.raw"),
#' "Grouping Var" = c("a","a","b","b"),
#' CONTROL = c("C","C","T","T"),
#' Subject = c("X","Y","X","Y"))
#' ax <- ap$read_annotation(annot)

AnnotationProcessor <- R6::R6Class(
  "AnnotationProcessor",

  public = list(
    #' @field QC is it a QC run
    QC = FALSE,
    #' @field prefix name for one factor designs
    prefix = "G_",
    #' @field repeated is it a repeated measurement
    repeated = TRUE,
    #' @field SAINT is it a AP MS experiment, then use Bait_ as prefix
    SAINT = FALSE,
    #' @field file_pattern colnames for file
    file_pattern = "^channel|^Relative|^raw|^file|^run",
    #' @field grouping_pattern colnames grouping variable
    grouping_pattern = "^group|^bait|^Experiment",
    #' @field subject_pattern colnames for pairing variable
    subject_pattern = "^subject|^BioReplicate",
    #' @field control_pattern contrast specification columns
    control_pattern = "ContrastName|Contrast|control",
    #' @field control_col_pattern columns which contains C or T.
    control_col_pattern = "^control",
    #' @field sample_name_pattern sample name column
    sample_name_pattern = "^name",
    #' @field sample_name_suffix_length maximum suffix length for display sample names
    sample_name_suffix_length = 14L,
    #' @field sample_name_display_column preferred derived display sample-name column
    sample_name_display_column = "sampleName",
    #' @field shorten_sample_names derive short display names for long sample names
    shorten_sample_names = TRUE,
    #' @field norm_value_pattern normalization value column (e.g., Creatinine)
    norm_value_pattern = "^creatinine|^normvalue",
    #' @field strict should name check be strict
    strict = FALSE,

    #' @description initialize
    #' @param QC default FALSE
    #' @param prefix default "G_"
    #' @param repeated default TRUE
    #' @param SAINT default FALSE
    #' @param shorten_sample_names derive short display sample names from long names
    #' @param sample_name_suffix_length suffix length used for derived sample names
    #' @param sample_name_display_column preferred derived display sample-name column
    initialize = function(
      QC = FALSE,
      prefix = "G_",
      repeated = TRUE,
      SAINT = FALSE,
      shorten_sample_names = TRUE,
      sample_name_suffix_length = 14L,
      sample_name_display_column = "sampleName"
    ) {
      self$QC <- QC
      self$prefix <- prefix
      self$repeated <- repeated
      self$SAINT <- SAINT
      self$shorten_sample_names <- shorten_sample_names
      self$sample_name_suffix_length <- sample_name_suffix_length
      self$sample_name_display_column <- sample_name_display_column
    },
    #' @description
    #' check annotation
    #' @param annot annotation
    check_annotation = function(annot) {
      warn_multiple <- function(cols) {
        if (length(cols) > 1) {
          warning(
            "there are more than one column for sample: ",
            paste(cols, collapse = ", ")
          )
        }
      }
      filename <- private$find_cols(annot, self$file_pattern)
      if (length(filename) < 1) {
        stop("column starting with :", self$file_pattern, " is missing.")
      }
      warn_multiple(filename)

      samples <- private$find_cols(annot, self$sample_name_pattern)
      if (length(samples) < 1) {
        warning("column starting with :", self$sample_name_pattern, " is missing.")
      }
      warn_multiple(samples)

      grouping <- private$grouping_cols(annot)
      # QC does not require a grouping column: set_grouping_var() injects a
      # single dummy group. For DEA a grouping column is mandatory.
      if (length(grouping) < 1 && self$QC) {
        warning(
          "no grouping column (",
          self$grouping_pattern,
          ") found; QC will use a single group."
        )
      } else if (length(grouping) < 1) {
        stop("column starting with :", self$grouping_pattern, " is missing.")
      }
      warn_multiple(grouping)

      if (!self$QC && length(private$find_cols(annot, self$control_pattern)) < 1) {
        stop("you must specify a CONTROL column.")
      }
      if ("CONTROL" %in% colnames(annot)) {
        stopifnot(all(c("C", "T") %in% annot[["CONTROL"]]))
      }
    },
    #' @description
    #' read annotation
    #' @param dsf either dataframe or file path.
    read_annotation = function(dsf) {
      annot <- if (inherits(dsf, "data.frame")) dsf else read_table_data(dsf)
      annot <- data.frame(lapply(annot, as.character), check.names = FALSE)
      self$check_annotation(annot)
      res <- private$dataset_set_factors(annot)
      if (!self$QC) {
        res$contrasts <- self$extract_contrasts(
          res$annot,
          group = res$atable$factors[[private$primary_factor_key()]]
        )
      }
      res
    },
    #' @description
    #' check annotation
    #' @param annot annotation
    #' @param group group column e.g. group
    extract_contrasts = function(annot, group) {
      factor_key <- private$primary_factor_key()
      levels <- annot |>
        dplyr::select(
          !!factor_key := starts_with(group, ignore.case = TRUE),
          control = starts_with("control", ignore.case = TRUE)
        ) |>
        dplyr::distinct()
      logger::log_info("levels: ", paste(levels, collapse = " "))
      if (length(levels[[factor_key]]) <= 1) {
        logger::log_error("not enough group levels to make comparisons.")
      }
      if (all(c("ContrastName", "Contrast") %in% colnames(annot))) {
        contr <- dplyr::filter(annot, nchar(!!rlang::sym("Contrast")) > 0)
        Contrasts <- contr$Contrast
        names(Contrasts) <- contr$ContrastName
        if (!any(grepl(paste0("\\b", factor_key), Contrasts))) {
          stop(
            "Group prefix should be: ",
            factor_key,
            "; but contrasts look like this: ",
            paste(Contrasts, collapse = "\n")
          )
        }
        return(Contrasts)
      }
      if (ncol(levels) != 2) {
        stop(
          "either column ",
          group,
          " or column control are missing. We found only column: ",
          paste(colnames(levels), collapse = " ")
        )
      }
      lv <- levels[[factor_key]]
      Contrasts <- character()
      Names <- character()
      for (i in seq_along(lv)) {
        for (j in setdiff(which(levels$control == "C"), i)) {
          logger::log_info("contrast: {lv[i]} vs {lv[j]}")
          Contrasts <- c(Contrasts, paste0(factor_key, lv[i], " - ", factor_key, lv[j]))
          Names <- c(Names, paste0(lv[i], "_vs_", lv[j]))
        }
      }
      if (!is.null(levels$control)) {
        names(Contrasts) <- Names
      }
      Contrasts
    },
    #' @description
    #' add vector of contrasts to annot table
    #' @param annot annotation
    #' @param Contrasts vector with contrasts
    add_contrasts_vec = function(annot, Contrasts) {
      if (length(Contrasts) > nrow(annot)) {
        warning("There are more Contrasts than samples.")
        return(annot)
      }
      pad <- rep(NA, nrow(annot) - length(Contrasts))
      annot$CONTROL <- NULL
      annot$ContrastName <- c(names(Contrasts), pad)
      annot$Contrast <- c(Contrasts, pad)
      annot
    }
  ),

  private = list(
    primary_factor_key = function() if (self$SAINT) "Bait_" else self$prefix,

    find_cols = function(annot, pattern) {
      grep(pattern, colnames(annot), value = TRUE, ignore.case = TRUE)
    },

    # Grouping candidates, dropping columns that carry no information (all NA /
    # blank): datasets often ship an empty "Bait ID" beside a populated
    # "Grouping Var", and the bait preference would otherwise pick the empty one.
    grouping_cols = function(annot) {
      cols <- private$find_cols(annot, self$grouping_pattern)
      non_empty <- vapply(
        cols,
        function(col) any(!is.na(annot[[col]]) & trimws(annot[[col]]) != ""),
        logical(1)
      )
      if (any(non_empty)) cols[non_empty] else cols
    },

    dataset_set_factors = function(annot) {
      atable <- prolfqua::AnalysisConfiguration$new()
      annot <- private$set_sample_name(annot, atable)
      atable$file_name <- private$find_cols(annot, self$file_pattern)[1]
      if (any(duplicated(annot[[atable$file_name]]))) {
        stop("file Names must be unique.")
      }
      annot <- private$set_grouping_var(annot, atable)
      private$process_subject_var(annot, atable)
      ctrl <- private$find_cols(annot, self$control_col_pattern)
      if (length(ctrl) == 1) {
        atable$factors[["CONTROL"]] <- ctrl
        stopifnot(all(annot[[ctrl]] %in% c("C", "T")))
      }
      norm_col <- private$find_cols(annot, self$norm_value_pattern)
      if (length(norm_col) >= 1) {
        atable$norm_value <- norm_col[1]
      }
      if (length(norm_col) > 1) {
        warning(
          "Multiple normalization value columns found: ",
          paste(norm_col, collapse = ", "),
          ". Using: ",
          norm_col[1]
        )
      }
      list(atable = atable, annot = annot)
    },

    set_sample_name = function(annot, atable) {
      source_sample_name <- private$find_cols(annot, self$sample_name_pattern)[1]
      if (is.na(source_sample_name)) {
        return(annot)
      }
      atable$sample_name <- source_sample_name
      display_names <- annot[[source_sample_name]]
      n_chars <- nchar(display_names, type = "chars")
      needs_shortening <- self$shorten_sample_names &&
        any(n_chars > self$sample_name_suffix_length, na.rm = TRUE)
      needs_unique_display <- any(duplicated(display_names))
      if (self$strict && needs_unique_display) {
        stop("sample Names must be unique.")
      }
      if (!needs_shortening && !needs_unique_display) {
        return(annot)
      }

      display_col <- self$sample_name_display_column
      candidate <- display_col
      index <- 0L
      while (candidate %in% colnames(annot) && candidate != source_sample_name) {
        index <- index + 1L
        candidate <- paste0(display_col, "_", index)
      }
      if (needs_shortening) {
        starts <- pmax(1L, n_chars - self$sample_name_suffix_length + 1L)
        display_names <- substring(display_names, starts, n_chars)
      }
      display_names[is.na(display_names) | !nzchar(display_names)] <- "NA"
      annot[[candidate]] <- make.unique(display_names, sep = "_")
      atable$sample_name <- candidate
      logger::log_info("Using derived sample display names in column '{candidate}'.")
      annot
    },

    set_grouping_var = function(annot, atable) {
      groupingVAR <- private$grouping_cols(annot)
      # QC datasets may carry no grouping column at all: synthesize one, which
      # the NA coercion below turns into a single "NA" group.
      if (length(groupingVAR) < 1) {
        annot[["group"]] <- NA_character_
        groupingVAR <- "group"
      }
      groupingVAR <- c(grep("^bait", groupingVAR, value = TRUE, ignore.case = TRUE), groupingVAR)[1]
      # Missing / blank entries become the literal "NA" group so the grouping
      # factor is never all-NA (which crashes the missingness heatmap).
      vals <- as.character(annot[[groupingVAR]])
      vals[is.na(vals) | trimws(vals) == ""] <- "NA"
      vals <- gsub("[[:space:]]", "", vals)
      annot[[groupingVAR]] <- gsub("[-\\+\\/\\*\\(\\)]", "_", vals)
      atable$factors[[private$primary_factor_key()]] <- groupingVAR
      atable$factor_depth <- 1
      annot
    },

    process_subject_var = function(annot, atable) {
      subvar <- private$find_cols(annot, self$subject_pattern)
      if (length(subvar) == 1 && self$repeated) {
        atable$factors[["Subject_"]] <- subvar
        group <- atable$factors[[private$primary_factor_key()]]
        fct <- dplyr::distinct(annot[, c(atable$file_name, group, subvar)])
        if (all(table(fct[, c(group, subvar)]) >= 1)) {
          atable$factor_depth <- 2
        }
      }
    }
  )
)

# read_annotation -----
#' read annotation files
#' @return list with annot (annotation table), atable (analtysis table configuration), contrasts list with contrasts.
#' @param dsf annotation table
#' @param repeated is this a repeated measurement
#' @param SAINT is this a SAINTexpress analysis
#' @param prefix prefix for group levels
#' @param QC if TRUE, read as QC annotation
#' @param shorten_sample_names derive short display sample names from long names
#' @param sample_name_suffix_length suffix length used for derived sample names
#' @param sample_name_display_column preferred derived display sample-name column
#' @export
#' @examples
#' annot <- data.frame(
#' file = c("a1.raw","a2.raw","a3.raw","a4.raw"),
#' name = c("aa","ba","aa","ba"),
#' group = c("a","a","b","b"))
#' read_annotation(annot, QC = TRUE)
#'
read_annotation <- function(
  dsf,
  repeated = TRUE,
  SAINT = FALSE,
  prefix = "G_",
  QC = FALSE,
  shorten_sample_names = TRUE,
  sample_name_suffix_length = 14L,
  sample_name_display_column = "sampleName"
) {
  AnnotationProcessor$new(
    repeated = repeated,
    SAINT = SAINT,
    prefix = prefix,
    QC = QC,
    shorten_sample_names = shorten_sample_names,
    sample_name_suffix_length = sample_name_suffix_length,
    sample_name_display_column = sample_name_display_column
  )$read_annotation(dsf)
}

#' extract contrast from annotation file
#' @param annot annotation data frame
#' @param prefix prefix for group levels
#' @param group name of the group column
#' @export
#' @examples
#'
#' annot <- data.frame(names = c("a1","b1"), group= c("a","b"), ddd = c("T","C"))
#' testthat::expect_error(extract_contrasts(annot))
#' annot$control <- annot$ddd
#' contrast <- extract_contrasts(annot)
#' stopifnot(contrast == "G_a - G_b")
#'
#' annot$Contrast <- c("G_a - G_b","G_b - G_a")
#' annot$ContrastName <- c("a_vs_b","b_vs_a")
#' annot$control <- NULL
#' ct <- extract_contrasts(annot)
#' stopifnot(length(ct) == 2)
extract_contrasts <- function(annot, prefix = "G_", group = "group") {
  AnnotationProcessor$new(prefix = prefix)$extract_contrasts(
    annot,
    group = group
  )
}

#' add vector of contrasts to annotation data frame
#' @param xx annotation data frame
#' @param Contrasts character vector of contrasts
#' @export
#' @examples
#' annot <- data.frame(Group = rep(c("A","B","C"), each = 3))
#' annot$Name
add_contrasts_vec <- function(xx, Contrasts) {
  AnnotationProcessor$new()$add_contrasts_vec(xx, Contrasts)
}
