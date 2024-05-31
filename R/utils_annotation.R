#' check if all required columns in annotation file are there.
check_annotation <- function(annot, QC = FALSE) {
  samples <- grep("^channel|^Relative|^raw|^file", colnames(annot), ignore.case = TRUE, value = TRUE)
  if (length(samples) < 1) { stop("column starting with channel, Relativ, raw, or file is missing.") }
  if (length(samples) > 1) { warning("there are more then one column for sample : ", paste0(samples)) }

  grouping <- grep("^group|^bait|^Experiment", colnames(annot), ignore.case = TRUE, value = TRUE)
  if (length(grouping) < 1) { stop("column with grouping variable (starting with group, bait, Experiment) is missing.") }
  if (length(grouping) > 1) { warning("there are more then one column for sample : ", paste0(grouping)) }
  if (!QC) {
  contrast <- grep("ContrastName|Contrast|CONTROL", colnames(annot), ignore.case = TRUE, value = TRUE)
  if (length(contrast) < 1) { stop(paste0("you must specify a CONTROL column.")) }
  }
  if ("CONTROL" %in%  colnames(annot)) {
    stopifnot(all(c("C","T") %in% annot[["CONTROL"]]))
  }
}


#' read annotation files
#' @return list with annot (annotation table), atable (analtysis table configuration), contrasts list with contrasts.
#' @param dsf annotation table
#' @param REPEATED is this a repeated measurement
#' @param SAINT is this a SAINTexpress analysis
#' @export
read_annotation <- function(dsf, REPEATED = TRUE, SAINT = FALSE, prefix = "G_", QC = FALSE){
  if ("data.frame" %in% class(dsf) ) {
    annot <- dsf
  } else {
    annot <- read.csv(dsf)
  }
  annot <- data.frame(lapply(annot, as.character))
  check_annotation(dsf, QC = QC)
  res <- dataset_set_factors(annot, REPEATED = REPEATED, SAINT = SAINT, prefix = prefix)
  if (!QC) {
    contrasts <- extract_contrasts(res$annot, prefix = prefix, group = res$atable$factors[[prefix]])
    res[["contrasts"]] <- contrasts
  }
  return(res)
}


dataset_set_factors <- function(annot, REPEATED = TRUE, SAINT = FALSE, prefix = "G_") {
  atable <- prolfqua::AnalysisTableAnnotation$new()

  if (sum(grepl("^name|^sample", colnames(annot), ignore.case = TRUE)) > 0) {
    atable$sampleName <- grep("^name|^sample", colnames(annot), value = TRUE, ignore.case = TRUE)[1]
  }
  if (any(duplicated(annot[[atable$sampleName]]))) {
    stop("sample Names must be unique.")
  }

  fileName <- grep("^channel|^Relative|^raw", colnames(annot), value = TRUE, ignore.case = TRUE)[1]
  atable$fileName <- fileName
  if (any(duplicated(annot[[atable$fileName]]))) {
    stop("file Names must be unique.")
  }


  groupingVAR <- grep("^group|^bait|^Experiment", colnames(annot), value = TRUE, ignore.case = TRUE)
  if (any(grepl("^bait", groupingVAR, ignore.case = TRUE))) {
    groupingVAR <- grep("^bait", groupingVAR, value = TRUE, ignore.case = TRUE)[1]
  } else {
    groupingVAR <- groupingVAR[1]
  }

  annot[[groupingVAR]] <- gsub("[[:space:]]", "", annot[[groupingVAR]])
  annot[[groupingVAR]] <- gsub("[-\\+\\/\\*\\(\\)]", "_", annot[[groupingVAR]])

  if (SAINT) {
    atable$factors[["Bait_"]] = groupingVAR
  } else {
    atable$factors[[prefix]] = groupingVAR
  }

  atable$factorDepth <- 1

  if (sum(grepl("^subject|^BioReplicate", colnames(annot), ignore.case = TRUE)) == 1 & REPEATED) {
    subvar <- grep("^subject|^BioReplicate", colnames(annot), value = TRUE, ignore.case = TRUE)
    atable$factors[["Subject_"]] = subvar

    fct <- dplyr::distinct(annot[,c(atable$fileName, groupingVAR, subvar)])
    tmp <- data.frame(table(fct[,c(groupingVAR,subvar)]))
    if (all(tmp$Freq > 1)) {
      atable$factorDepth <- 2
    }
  }
  ctrl <- grep("^control", colnames(annot), value = TRUE, ignore.case = TRUE)
  if (length(ctrl) == 1) {
    atable$factors[["CONTROL"]] = ctrl

    stopifnot(length(setdiff(unique(annot[[ctrl]]),c("C","T"))) == 0)
    # TODO add check that
    tt <- table(annot[[ctrl]], annot[[groupingVAR]])

  }
  return(list(atable = atable , annot = annot))
}


#' extract contrast from annotation file
#' @export
#' @example
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
  levels  <- annot |>
    dplyr::select(
      !!prefix := dplyr::starts_with(group, ignore.case = TRUE),
      control = dplyr::starts_with("control", ignore.case = TRUE)) |>
    dplyr::distinct()
  logger::log_info("levels : ", paste(levels, collapse = " "))
  if (!length(levels[[prefix]]) > 1) {
    logger::log_error("not enough group levels_ to make comparisons.")
  }
  if ( all(c("ContrastName", "Contrast") %in% colnames(annot)) ) {
    contr <- annot |>
      dplyr::select(all_of(c("ContrastName", "Contrast"))) |>
      dplyr::filter(nchar(!!rlang::sym("Contrast")) > 0)

    Contrasts <- contr$Contrast
    names(Contrasts) <- contr$ContrastName
    nrpr <- sum(grepl(paste0("\\b",prefix), Contrasts))
    if ( nrpr < 1 ) {
      stop("Group prefix should be :", prefix, "; but contrasts look like this : ", paste(Contrasts, collapse = "\n"))
    }
    return(Contrasts)
  } else {
    if (ncol(levels) != 2) {
      stop("either column", group, " or column control are missing. We found only column :", paste(colnames(levels), collapse = " "))
    }
    Contrasts <- character()
    ## Generate contrasts from dataset
    if (!is.null(levels$control)) {
      Names <- character()
      for (i in 1:nrow(levels)) {
        for (j in 1:nrow(levels)) {
          if (i != j) {
            if (levels$control[j] == "C") {
              cat(levels[[prefix]][i], levels[[prefix]][j], "\n")
              Contrasts <- c(Contrasts, paste0(prefix,levels[[prefix]][i], " - ", prefix,levels[[prefix]][j]))
              Names <- c(Names, paste0(levels[[prefix]][i], "_vs_", levels[[prefix]][j]))
            }
          }
        }
      }
      names(Contrasts) <- Names
    }
    return(Contrasts)
  }
}

#' add vector of contrasts to annotation data frame
#' @export
#' @examples
#' annot <- data.frame(Group = rep(c("A","B","C"), each = 3))
#' annot$Name
add_contrasts_vec <- function(xx, Contrasts){
  if (length(Contrasts) <= nrow(xx)) {
    xx$CONTROL <- NULL
    xx$ContrastName <- c(names(Contrasts), rep(NA, nrow(xx) - length(Contrasts) ))
    xx$Contrast <- c(Contrasts, rep(NA, nrow(xx) - length(Contrasts) ))
  } else {
    warning("There are more Contrasts then samples.")
  }
  return(xx)
}

