# Create a named list of functions
#' read dataset file in csv, tsv or xlsx format
#' @param file_path path to csv, tsv, or xlsx file
#' @export
read_table_data <- function(file_path) {
  read_functions <- list(
    csv = readr::read_csv,
    tsv = readr::read_tsv,
    xlsx = readxl::read_xlsx
  )

  # Get the file extension
  file_extension <- tools::file_ext(file_path)

  # Check if the file extension is supported
  if (!file_extension %in% names(read_functions)) {
    stop("Unsupported file extension")
  }

  # Call the appropriate reading function
  data <- read_functions[[file_extension]](file_path)

  return(data)
}

# Create a named list of functions
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
  write_functions <- list(
    csv = readr::write_csv,
    tsv = readr::write_tsv,
    xlsx = writexl::write_xlsx
  )

  # Get the file extension
  file_extension <- tools::file_ext(file_path)

  # Check if the file extension is supported
  if (!file_extension %in% names(write_functions)) {
    stop("Unsupported file extension")
  }

  # Call the appropriate writing function
  write_functions[[file_extension]](data, file_path)
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
#' as <- annot
#' as$sample <- c("s1","s2","s3","s4")
#' aa <- ap$read_annotation(annot)
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
      filename <- grep(
        self$file_pattern,
        colnames(annot),
        ignore.case = TRUE,
        value = TRUE
      )
      if (length(filename) < 1) {
        stop("column starting with :", self$file_pattern, " is missing.")
      }
      if (length(filename) > 1) {
        warning(
          "there are more than one column for sample: ",
          paste(filename, collapse = ", ")
        )
      }

      samples <- grep(
        self$sample_name_pattern,
        colnames(annot),
        ignore.case = TRUE,
        value = TRUE
      )
      if (length(samples) < 1) {
        warning(
          "column starting with :",
          self$sample_name_pattern,
          " is missing."
        )
      }
      if (length(samples) > 1) {
        warning(
          "there are more than one column for sample: ",
          paste(samples, collapse = ", ")
        )
      }

      grouping <- grep(
        self$grouping_pattern,
        colnames(annot),
        ignore.case = TRUE,
        value = TRUE
      )
      non_empty_grouping <- vapply(
        grouping,
        function(col) {
          any(!is.na(annot[[col]]) & trimws(as.character(annot[[col]])) != "")
        },
        logical(1)
      )
      if (any(non_empty_grouping)) {
        grouping <- grouping[non_empty_grouping]
      }
      if (length(grouping) < 1) {
        # QC does not require a grouping variable: a single dummy group is
        # injected later in set_grouping_var(). For DEA a grouping column is
        # mandatory. (Message names grouping_pattern, not sample_name_pattern.)
        if (self$QC) {
          warning(
            "no grouping column (",
            self$grouping_pattern,
            ") found; QC will use a single group."
          )
        } else {
          stop("column starting with :", self$grouping_pattern, " is missing.")
        }
      }
      if (length(grouping) > 1) {
        warning(
          "there are more than one column for sample: ",
          paste(grouping, collapse = ", ")
        )
      }

      if (!self$QC) {
        contrast <- grep(
          self$control_pattern,
          colnames(annot),
          ignore.case = TRUE,
          value = TRUE
        )
        if (length(contrast) < 1) {
          stop(paste0("you must specify a CONTROL column."))
        }
      }

      if ("CONTROL" %in% colnames(annot)) {
        stopifnot(all(c("C", "T") %in% annot[["CONTROL"]]))
      }
    },
    #' @description
    #' read annotation
    #' @param dsf either dataframe or file path.
    read_annotation = function(dsf) {
      if ("data.frame" %in% class(dsf)) {
        annot <- dsf
      } else {
        annot <- prolfquapp::read_table_data(dsf)
      }
      annot <- data.frame(lapply(annot, as.character), check.names = FALSE)
      self$check_annotation(annot)
      res <- private$dataset_set_factors(annot)
      if (!self$QC) {
        factor_key <- private$primary_factor_key()
        contrasts <- self$extract_contrasts(
          res$annot,
          group = res$atable$factors[[factor_key]]
        )
        res[["contrasts"]] <- contrasts
      }
      return(res)
    },
    #' @description
    #' check annotation
    #' @param annot annotation
    #' @param group group column e.g. group
    extract_contrasts = function(annot, group) {
      levels <- private$get_levels(annot, group)
      logger::log_info("levels: ", paste(levels, collapse = " "))
      factor_key <- private$primary_factor_key()
      if (!(length(levels[[factor_key]]) > 1)) {
        logger::log_error("not enough group levels to make comparisons.")
      }
      if (all(c("ContrastName", "Contrast") %in% colnames(annot))) {
        return(private$get_defined_contrasts(annot))
      } else {
        return(private$generate_contrasts(annot, levels, group))
      }
    },
    #' @description
    #' add vector of contrasts to annot table
    #' @param annot annotation
    #' @param Contrasts vector with contrasts
    add_contrasts_vec = function(annot, Contrasts) {
      if (length(Contrasts) <= nrow(annot)) {
        annot$CONTROL <- NULL
        annot$ContrastName <- c(
          names(Contrasts),
          rep(NA, nrow(annot) - length(Contrasts))
        )
        annot$Contrast <- c(Contrasts, rep(NA, nrow(annot) - length(Contrasts)))
      } else {
        warning("There are more Contrasts than samples.")
      }
      return(annot)
    }
  ),

  private = list(
    primary_factor_key = function() {
      if (self$SAINT) {
        return("Bait_")
      }
      self$prefix
    },

    dataset_set_factors = function(annot) {
      atable <- prolfqua::AnalysisConfiguration$new()
      annot <- private$set_sample_name(annot, atable)
      private$set_file_name(annot, atable)
      annot <- private$set_grouping_var(annot, atable)
      private$process_subject_var(annot, atable)
      private$set_control_var(annot, atable)
      private$set_norm_value(annot, atable)
      return(list(atable = atable, annot = annot))
    },

    set_sample_name = function(annot, atable) {
      sample_cols <- grep(
        self$sample_name_pattern,
        colnames(annot),
        value = TRUE,
        ignore.case = TRUE
      )
      if (length(sample_cols) == 0) {
        return(annot)
      }

      source_sample_name <- sample_cols[1]
      atable$sample_name <- source_sample_name
      sample_names <- annot[[source_sample_name]]
      needs_shortening <- self$shorten_sample_names &&
        any(
          nchar(sample_names, type = "chars") > self$sample_name_suffix_length,
          na.rm = TRUE
        )
      needs_unique_display <- any(duplicated(sample_names))

      if (self$strict && needs_unique_display) {
        stop("sample Names must be unique.")
      }

      if (!needs_shortening && !needs_unique_display) {
        return(annot)
      }

      display_col <- private$available_sample_name_column(
        annot,
        source_sample_name
      )
      display_names <- sample_names
      if (needs_shortening) {
        display_names <- private$suffix_sample_names(display_names)
      }
      display_names[is.na(display_names) | !nzchar(display_names)] <- "NA"
      if (any(duplicated(display_names))) {
        display_names <- make.unique(display_names, sep = "_")
      }

      annot[[display_col]] <- display_names
      atable$sample_name <- display_col
      logger::log_info(
        "Using derived sample display names in column '{display_col}'."
      )
      return(annot)
    },

    available_sample_name_column = function(annot, source_sample_name) {
      display_col <- self$sample_name_display_column
      if (
        !display_col %in% colnames(annot) || display_col == source_sample_name
      ) {
        return(display_col)
      }

      candidate <- display_col
      index <- 1L
      while (
        candidate %in% colnames(annot) && candidate != source_sample_name
      ) {
        candidate <- paste0(display_col, "_", index)
        index <- index + 1L
      }
      candidate
    },

    suffix_sample_names = function(sample_names) {
      n_chars <- nchar(sample_names, type = "chars")
      starts <- pmax(1L, n_chars - self$sample_name_suffix_length + 1L)
      substring(sample_names, starts, n_chars)
    },

    set_file_name = function(annot, atable) {
      fileName <- grep(
        self$file_pattern,
        colnames(annot),
        value = TRUE,
        ignore.case = TRUE
      )[1]
      atable$file_name <- fileName
      if (any(duplicated(annot[[atable$file_name]]))) {
        stop("file Names must be unique.")
      }
    },

    set_grouping_var = function(annot, atable) {
      groupingVAR <- grep(
        self$grouping_pattern,
        colnames(annot),
        value = TRUE,
        ignore.case = TRUE
      )
      # Drop candidate columns that carry no information (all NA / all blank).
      # Datasets frequently ship an empty "Bait ID" column alongside a populated
      # "Grouping Var"; without this filter the bait-preference below would pick
      # the empty column, yielding an all-NA grouping factor that crashes the
      # missingness heatmap (pheatmap: "'gpar' element 'fill' must not be length 0").
      non_empty <- vapply(
        groupingVAR,
        function(col) {
          any(!is.na(annot[[col]]) & trimws(as.character(annot[[col]])) != "")
        },
        logical(1)
      )
      if (any(non_empty)) {
        groupingVAR <- groupingVAR[non_empty]
      }
      # QC datasets may carry no grouping column at all. Synthesize one so a
      # valid factor is still produced; the NA-coercion below turns the empty
      # column into a single "NA" group (mirrors the all-empty-column case).
      if (length(groupingVAR) < 1) {
        annot[["group"]] <- NA_character_
        groupingVAR <- "group"
      }
      if (any(grepl("^bait", groupingVAR, ignore.case = TRUE))) {
        groupingVAR <- grep(
          "^bait",
          groupingVAR,
          value = TRUE,
          ignore.case = TRUE
        )[1]
      } else {
        groupingVAR <- groupingVAR[1]
      }

      # Coerce missing / blank entries to the literal "NA" so an empty (or
      # partially empty) grouping column still yields a valid factor instead of
      # an all-NA grouping that crashes downstream (e.g. the missingness
      # heatmap, pheatmap: "'gpar' element 'fill' must not be length 0").
      # All-empty -> a single "NA" group; partial -> real groups plus an "NA"
      # group. Mirrors the display-name handling in set_sample_name().
      grouping_vals <- as.character(annot[[groupingVAR]])
      grouping_vals[is.na(grouping_vals) | trimws(grouping_vals) == ""] <- "NA"
      annot[[groupingVAR]] <- grouping_vals

      annot[[groupingVAR]] <- gsub("[[:space:]]", "", annot[[groupingVAR]])
      annot[[groupingVAR]] <- gsub(
        "[-\\+\\/\\*\\(\\)]",
        "_",
        annot[[groupingVAR]]
      )

      if (self$SAINT) {
        atable$factors[["Bait_"]] <- groupingVAR
      } else {
        atable$factors[[self$prefix]] <- groupingVAR
      }

      atable$factor_depth <- 1
      return(annot)
    },

    process_subject_var = function(annot, atable) {
      if (
        sum(grepl(self$subject_pattern, colnames(annot), ignore.case = TRUE)) ==
          1 &
          self$repeated
      ) {
        subvar <- grep(
          self$subject_pattern,
          colnames(annot),
          value = TRUE,
          ignore.case = TRUE
        )
        atable$factors[["Subject_"]] <- subvar
        factor_key <- private$primary_factor_key()

        fct <- dplyr::distinct(annot[, c(
          atable$file_name,
          atable$factors[[factor_key]],
          subvar
        )])
        tmp <- data.frame(table(fct[, c(
          atable$factors[[factor_key]],
          subvar
        )]))
        if (all(tmp$Freq >= 1)) {
          atable$factor_depth <- 2
        }
      }
    },

    set_control_var = function(annot, atable) {
      ctrl <- grep(
        self$control_col_pattern,
        colnames(annot),
        value = TRUE,
        ignore.case = TRUE
      )
      if (length(ctrl) == 1) {
        atable$factors[["CONTROL"]] <- ctrl
        factor_key <- private$primary_factor_key()

        stopifnot(length(setdiff(unique(annot[[ctrl]]), c("C", "T"))) == 0)
        # TODO add check that
        tt <- table(annot[[ctrl]], annot[[atable$factors[[factor_key]]]])
      }
    },

    set_norm_value = function(annot, atable) {
      norm_col <- grep(
        self$norm_value_pattern,
        colnames(annot),
        value = TRUE,
        ignore.case = TRUE
      )
      if (length(norm_col) >= 1) {
        atable$norm_value <- norm_col[1]
        if (length(norm_col) > 1) {
          warning(
            "Multiple normalization value columns found: ",
            paste(norm_col, collapse = ", "),
            ". Using: ",
            norm_col[1]
          )
        }
      }
    },

    get_levels = function(annot, group) {
      factor_key <- private$primary_factor_key()
      levels <- annot |>
        dplyr::select(
          !!factor_key := starts_with(group, ignore.case = TRUE),
          control = starts_with("control", ignore.case = TRUE)
        ) |>
        dplyr::distinct()
      return(levels)
    },

    get_defined_contrasts = function(annot) {
      factor_key <- private$primary_factor_key()
      contr <- annot |>
        dplyr::select(all_of(c("ContrastName", "Contrast"))) |>
        dplyr::filter(nchar(!!rlang::sym("Contrast")) > 0)

      Contrasts <- contr$Contrast
      names(Contrasts) <- contr$ContrastName
      nrpr <- sum(grepl(paste0("\\b", factor_key), Contrasts))
      if (nrpr < 1) {
        stop(
          "Group prefix should be: ",
          factor_key,
          "; but contrasts look like this: ",
          paste(Contrasts, collapse = "\n")
        )
      }
      return(Contrasts)
    },

    generate_contrasts = function(annot, levels, group) {
      factor_key <- private$primary_factor_key()
      if (ncol(levels) != 2) {
        stop(
          "either column ",
          group,
          " or column control are missing. We found only column: ",
          paste(colnames(levels), collapse = " ")
        )
      }
      Contrasts <- character()
      Names <- character()
      ## Generate contrasts from dataset
      if (!is.null(levels$control)) {
        for (i in seq_len(nrow(levels))) {
          for (j in seq_len(nrow(levels))) {
            if (i != j && levels$control[j] == "C") {
              logger::log_info(
                "contrast: {levels[[factor_key]][i]} vs {levels[[factor_key]][j]}"
              )
              Contrasts <- c(
                Contrasts,
                paste0(
                  factor_key,
                  levels[[factor_key]][i],
                  " - ",
                  factor_key,
                  levels[[factor_key]][j]
                )
              )
              Names <- c(
                Names,
                paste0(
                  levels[[factor_key]][i],
                  "_vs_",
                  levels[[factor_key]][j]
                )
              )
            }
          }
        }
        names(Contrasts) <- Names
      }
      return(Contrasts)
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
  res <- AnnotationProcessor$new(
    repeated = repeated,
    SAINT = SAINT,
    prefix = prefix,
    QC = QC,
    shorten_sample_names = shorten_sample_names,
    sample_name_suffix_length = sample_name_suffix_length,
    sample_name_display_column = sample_name_display_column
  )$read_annotation(dsf)
  return(res)
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
