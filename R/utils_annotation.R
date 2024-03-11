#' check if all required columns in annotation file are there.
check_annotation <- function(annot) {
  samples <- grep("^channel|^Relative|^raw", colnames(annot), ignore.case = TRUE, value = TRUE)
  stopifnot(length(samples) >= 1)
  if (length(samples) > 1) {warning("there are more then one column for sample : ", paste0(samples)) }
  grouping <- grep("^group|^bait|^Experiment", colnames(annot), ignore.case = TRUE, value = TRUE)
  stopifnot(length(grouping) >= 1 )
  if (length(grouping) > 1) {warning("there are more then one column for sample : ", paste0(grouping)) }
}


#' read annotation files
#' @return list with annot (annotation table), atable (analtysis table configuration), contrasts list with contrasts.
#' @export
read_annotation <- function(dsf, REPEATED = TRUE, SAINT = FALSE){
  if ("data.frame" %in% class(dsf) ) {
    annot <- dsf
  } else {
    annot <- read.csv(dsf)
  }
  annot <- data.frame(lapply(annot, as.character))
  check_annotation(dsf)
  res <- dataset_set_factors(annot, REPEATED = REPEATED, SAINT = SAINT)
  contrasts <- extract_contrasts(res$annot)
  res[["contrasts"]] <- contrasts
  return(res)
}


dataset_set_factors <- function(annot, REPEATED = TRUE, SAINT = FALSE) {

  atable <- prolfqua::AnalysisTableAnnotation$new()

  if (sum(grepl("^name|^sample", colnames(annot), ignore.case = TRUE)) > 0) {
    atable$sampleName <- grep("^name|^sample", colnames(annot), value = TRUE, ignore.case = TRUE)[1]
  }

  fileName <- grep("^channel|^Relative|^raw", colnames(annot), value = TRUE, ignore.case = TRUE)[1]
  atable$fileName <- fileName

  groupingVAR <- grep("^group|^bait|^Experiment", colnames(annot), value = TRUE, ignore.case = TRUE)
  if (any(grepl("^bait", groupingVAR, ignore.case = TRUE))) {
    groupingVAR <- grep("^bait", groupingVAR, value = TRUE, ignore.case = TRUE)[1]
  } else {
    groupingVAR <- groupingVAR[1]
  }

  annot[[groupingVAR]] <- gsub("[[:space:]]", "", annot[[groupingVAR]])
  annot[[groupingVAR]] <- gsub("[-\\+\\/\\*]", "_", annot[[groupingVAR]])

  if (SAINT) {
    atable$factors[["Bait_"]] = groupingVAR
  } else {
    atable$factors[["Group_"]] = groupingVAR
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
  if (sum(grepl("^control", colnames(annot), ignore.case = TRUE)) == 1) {
    atable$factors[["CONTROL"]] = grep("^control", colnames(annot), value = TRUE, ignore.case = TRUE)
  }
  return(list(atable = atable , annot = annot))
}


# private interface
extract_contrasts <- function(annot) {

  levels  <- annot |>
    dplyr::select(
      Group_ = dplyr::starts_with("group", ignore.case = TRUE),
      control = dplyr::starts_with("control", ignore.case = TRUE)) |>
    dplyr::distinct()
  logger::log_info("levels : ", paste(levels, collapse = " "))
  if (!length(levels$Group_) > 1) {
    logger::log_error("not enough group levels_ to make comparisons.")
  }

  if ( all(c("ContrastName", "Contrast") %in% colnames(annot)) ) {
    contr <- annot |>
      dplyr::select(all_of(c("ContrastName", "Contrast"))) |>
      dplyr::filter(nchar(!!rlang::sym("Contrast")) > 0)
    Contrasts <- contr$Contrast
    names(Contrasts) <- contr$ContrastName

  } else {

    ## Generate contrasts from dataset
    if (!is.null(levels$control)) {


      Contrasts <- character()
      Names <- character()
      for (i in 1:nrow(levels)) {
        for (j in 1:nrow(levels)) {
          if (i != j) {
            if (levels$control[j] == "C") {
              cat(levels$Group_[i], levels$Group_[j], "\n")
              Contrasts <- c(Contrasts, paste0("Group_",levels$Group_[i], " - ", "Group_",levels$Group_[j]))
              Names <- c(Names, paste0(levels$Group_[i], "_vs_", levels$Group_[j]))
            }
          }
        }
      }
    }
    names(Contrasts) <- Names
  }
  return(Contrasts)

}
