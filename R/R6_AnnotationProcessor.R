# Create a named list of functions
#' read dataset file in csv, tsv or xlsx format
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
#' # af$group <- NULL
#' af$CONTROL <- NULL
#' af$Subject <- NULL
#' ap$check_annotation(af)
#' aa <- ap$read_annotation(af)
#'
#' stopifnot(aa$atable$factor_keys() == "G_")
#' stopifnot(aa$atable$factors == "group")
#' aa <- ap$read_annotation(annot)
#' aa$atable$fileName
#' aa$atable$sampleName
#' as <- annot
#' as$sample <- c("s1","s2","s3","s4")
#' aa <- ap$read_annotation(annot)
#' aa$atable$sampleName
#' stopifnot(is.null(aa$annotation))
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
    file_pattern = "^channel|^Relative|^raw|^file",
    #' @field grouping_pattern colnames grouping variable
    grouping_pattern = "^group|^bait|^Experiment",
    #' @field subject_pattern colnames for pairing variable
    subject_pattern = "^subject|^BioReplicate",
    #' @field control_pattern contrast specification columns
    control_pattern = "ContrastName|Contrast|control",
    #' @field control_col_pattern columns which contains C or T.
    control_col_pattern = "^control",
    #' @field sample_name_pattern sample name column
    sample_name_pattern = "^name|^sample",
    #' @field strict should name check be strict
    strict = FALSE,

    #' @description initialize
    #' @param QC default FALSE
    #' @param prefix default "G_"
    #' @param repeated default TRUE
    #' @param SAINT default FALSE
    initialize = function(QC = FALSE,
                          prefix = "G_",
                          repeated = TRUE,
                          SAINT = FALSE) {
      self$QC <- QC
      self$prefix <- prefix
      self$repeated <- repeated
      self$SAINT <- SAINT
    },
    #' @description
    #' check annotation
    #' @param annot annotation
    check_annotation = function(annot) {
      filename <- grep(self$file_pattern, colnames(annot), ignore.case = TRUE, value = TRUE)
      if (length(filename) < 1) { stop("column starting with :", self$file_pattern , " is missing.") }
      if (length(filename) > 1) { warning("there are more than one column for sample: ", paste0(filename)) }

      samples <- grep(self$sample_name_pattern, colnames(annot), ignore.case = TRUE, value = TRUE)
      if (length(samples) < 1) { warning("column starting with :", self$sample_name_pattern , " is missing.") }
      if (length(samples) > 1) { warning("there are more than one column for sample: ", paste0(samples)) }

      grouping <- grep(self$grouping_pattern, colnames(annot), ignore.case = TRUE, value = TRUE)
      if (length(grouping) < 1) {  stop("column starting with :", self$sample_name_pattern , " is missing.")  }
      if (length(grouping) > 1) { warning("there are more than one column for sample: ", paste0(grouping)) }

      if (!self$QC) {
        contrast <- grep(self$control_pattern, colnames(annot), ignore.case = TRUE, value = TRUE)
        if (length(contrast) < 1) { stop(paste0("you must specify a CONTROL column.")) }
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
      annot <- data.frame(lapply(annot, as.character))
      self$check_annotation(annot)
      res <- private$dataset_set_factors(annot)
      if (!self$QC) {
        contrasts <- self$extract_contrasts(res$annot, group = res$atable$factors[[self$prefix]])
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
      if (!length(levels[[self$prefix]]) > 1) {
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
        annot$ContrastName <- c(names(Contrasts), rep(NA, nrow(annot) - length(Contrasts)))
        annot$Contrast <- c(Contrasts, rep(NA, nrow(annot) - length(Contrasts)))
      } else {
        warning("There are more Contrasts than samples.")
      }
      return(annot)
    }
  ),

  private = list(

    dataset_set_factors = function(annot) {
      atable <- prolfqua::AnalysisTableAnnotation$new()
      annot <- private$set_sample_name(annot, atable)
      private$set_file_name(annot, atable)
      annot <- private$set_grouping_var(annot, atable)
      private$process_subject_var(annot, atable)
      private$set_control_var(annot, atable)
      return(list(atable = atable, annot = annot))
    },

    set_sample_name = function(annot, atable) {
      if (sum(grepl(self$sample_name_pattern, colnames(annot), ignore.case = TRUE)) > 0) {
        atable$sampleName <- grep(self$sample_name_pattern, colnames(annot), value = TRUE, ignore.case = TRUE)[1]
      }
      if (self$strict && any(duplicated(annot[[atable$sampleName]]))) {
        stop("sample Names must be unique.")
      }else if (any(duplicated(annot[[atable$sampleName]]))) {
        annot[[atable$sampleName]] <- data.frame(xx = annot[[atable$sampleName]]) |>
          dplyr::group_by(xx) |>
          dplyr::mutate(count = dplyr::row_number()) |>
          tidyr::unite("name",c("xx","count")) |>
          dplyr::pull("name")
      } else {}
      return(annot)
    },

    set_file_name = function(annot, atable) {
      fileName <- grep(self$file_pattern, colnames(annot), value = TRUE, ignore.case = TRUE)[1]
      atable$fileName <- fileName
      if (any(duplicated(annot[[atable$fileName]]))) {
        stop("file Names must be unique.")
      }
    },

    set_grouping_var = function(annot, atable) {
      groupingVAR <- grep(self$grouping_pattern, colnames(annot), value = TRUE, ignore.case = TRUE)
      if (any(grepl("^bait", groupingVAR, ignore.case = TRUE))) {
        groupingVAR <- grep("^bait", groupingVAR, value = TRUE, ignore.case = TRUE)[1]
      } else {
        groupingVAR <- groupingVAR[1]
      }

      annot[[groupingVAR]] <- gsub("[[:space:]]", "", annot[[groupingVAR]])
      annot[[groupingVAR]] <- gsub("[-\\+\\/\\*\\(\\)]", "_", annot[[groupingVAR]])

      if (self$SAINT) {
        atable$factors[["Bait_"]] <- groupingVAR
      } else {
        atable$factors[[self$prefix]] <- groupingVAR
      }

      atable$factorDepth <- 1
      return(annot)
    },

    process_subject_var = function(annot, atable) {
      if (sum(grepl(self$subject_pattern, colnames(annot), ignore.case = TRUE)) == 1 & self$repeated) {
        subvar <- grep(self$subject_pattern, colnames(annot), value = TRUE, ignore.case = TRUE)
        atable$factors[["Subject_"]] <- subvar

        fct <- dplyr::distinct(annot[, c(atable$fileName, atable$factors[[self$prefix]], subvar)])
        tmp <- data.frame(table(fct[, c(atable$factors[[self$prefix]], subvar)]))
        if (all(tmp$Freq >= 1)) {
          atable$factorDepth <- 2
        }
      }
    },

    set_control_var = function(annot, atable) {
      ctrl <- grep(self$control_col_pattern, colnames(annot), value = TRUE, ignore.case = TRUE)
      if (length(ctrl) == 1) {
        atable$factors[["CONTROL"]] <- ctrl

        stopifnot(length(setdiff(unique(annot[[ctrl]]), c("C", "T"))) == 0)
        # TODO add check that
        tt <- table(annot[[ctrl]], annot[[atable$factors[[self$prefix]]]])
      }
    },


    get_levels = function(annot, group) {
      levels <- annot |>
        dplyr::select(
          !!self$prefix := starts_with(group, ignore.case = TRUE),
          control = starts_with("control", ignore.case = TRUE)) |>
        dplyr::distinct()
      return(levels)
    },

    get_defined_contrasts = function(annot) {
      contr <- annot |>
        dplyr::select(all_of(c("ContrastName", "Contrast"))) |>
        dplyr::filter(nchar(!!rlang::sym("Contrast")) > 0)

      Contrasts <- contr$Contrast
      names(Contrasts) <- contr$ContrastName
      nrpr <- sum(grepl(paste0("\\b", self$prefix), Contrasts))
      if (nrpr < 1) {
        stop("Group prefix should be: ", self$prefix, "; but contrasts look like this: ", paste(Contrasts, collapse = "\n"))
      }
      return(Contrasts)
    },

    generate_contrasts = function(annot, levels, group) {
      if (ncol(levels) != 2) {
        stop("either column ", group, " or column control are missing. We found only column: ", paste(colnames(levels), collapse = " "))
      }
      Contrasts <- character()
      Names <- character()
      ## Generate contrasts from dataset
      if (!is.null(levels$control)) {
        for (i in 1:nrow(levels)) {
          for (j in 1:nrow(levels)) {
            if (i != j && levels$control[j] == "C") {
              cat(levels[[self$prefix]][i], levels[[self$prefix]][j], "\n")
              Contrasts <- c(Contrasts, paste0(self$prefix, levels[[self$prefix]][i], " - ", self$prefix, levels[[self$prefix]][j]))
              Names <- c(Names, paste0(levels[[self$prefix]][i], "_vs_", levels[[self$prefix]][j]))
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
#' @export
#' @examples
#' annot <- data.frame(
#' file = c("a1.raw","a2.raw","a3.raw","a4.raw"),
#' name = c("aa","ba","aa","ba"),
#' group = c("a","a","b","b"))
#' read_annotation(annot, QC = TRUE)
#'
read_annotation <- function(dsf, repeated = TRUE, SAINT = FALSE, prefix = "G_", QC = FALSE){
  res <- AnnotationProcessor$new(repeated = repeated, SAINT = SAINT, prefix = prefix,QC = QC)$read_annotation(dsf)
  return(res)
}

#' extract contrast from annotation file
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
  AnnotationProcessor$new(prefix = prefix)$extract_contrasts(annot, group = group)
}

#' add vector of contrasts to annotation data frame
#' @export
#' @examples
#' annot <- data.frame(Group = rep(c("A","B","C"), each = 3))
#' annot$Name
add_contrasts_vec <- function(xx, Contrasts){
  AnnotationProcessor$new()$add_contrasts_vec(xx, Contrasts)
}


