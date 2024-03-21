#' set factors and sample names columns
#' @export
#' @examples
#'
#' dataset_csv_examples <- system.file("application/dataset_csv",package = "prolfquapp")
#' files <- dir(dataset_csv_examples, pattern = "*.csv")
#'
#' res <- list()
#' for (i in seq_along(files)) {
#'   print(i)
#'   res[[i]] <- readr::read_csv(file.path(dataset_csv_examples, files[i]))
#' }
#'
#' for (i in seq_along(res)) {
#'   atable <- prolfqua::AnalysisTableAnnotation$new()
#'   #atable$fileName = "channel"
#'   atable$hierarchy[["protein_Id"]] <- c("Protein")
#'   atable$hierarchy[["peptide_Id"]] <- c("Peptide")
#'   tmp <- prolfquapp::dataset_set_factors(atable, res[[i]] )
#'   cat(i, " : " , length(tmp$atable$factors), "factors : ", paste(tmp$atable$factors, collapse = "; "), "\n")
#' }
#'
dataset_set_factors_deprecated <- function(atable, msdata, REPEATED = TRUE, SAINT = FALSE) {
  if (sum(grepl("^name", colnames(msdata), ignore.case = TRUE)) > 0) {
    atable$sampleName <- grep("^name", colnames(msdata), value = TRUE, ignore.case = TRUE)
  }

  stopifnot(sum(grepl("^channel|^Relative|^raw", colnames(msdata), ignore.case = TRUE)) >= 1)
  fileName <- grep("^channel|^Relative|^raw", colnames(msdata), value = TRUE, ignore.case = TRUE)[1]
  atable$fileName <- fileName

  stopifnot(sum(grepl("^group|^bait|^Experiment", colnames(msdata), ignore.case = TRUE)) >= 1)

  groupingVAR <- grep("^group|^bait|^Experiment", colnames(msdata), value = TRUE, ignore.case = TRUE)
  if (any(grepl("^bait", groupingVAR, ignore.case = TRUE))) {
    groupingVAR <- grep("^bait", groupingVAR, value = TRUE, ignore.case = TRUE)[1]
  } else {
    groupingVAR <- groupingVAR[1]
  }

  msdata[[groupingVAR]] <- gsub("[[:space:]]", "", msdata[[groupingVAR]])
  msdata[[groupingVAR]] <- gsub("[-\\+\\/\\*]", "_", msdata[[groupingVAR]])

  if (SAINT) {
    atable$factors[["Bait_"]] = groupingVAR
  } else {
    atable$factors[["Group_"]] = groupingVAR
  }


  atable$factorDepth <- 1

  if (sum(grepl("^subject|^BioReplicate", colnames(msdata), ignore.case = TRUE)) == 1 & REPEATED) {
    subvar <- grep("^subject|^BioReplicate", colnames(msdata), value = TRUE, ignore.case = TRUE)
    atable$factors[["Subject_"]] = subvar

    fct <- dplyr::distinct(msdata[,c(atable$fileName, groupingVAR, subvar)])
    tmp <- data.frame(table(fct[,c(groupingVAR,subvar)]))
    if (all(tmp$Freq > 1)) {
      atable$factorDepth <- 2
    }
  }
  if (sum(grepl("^control", colnames(msdata), ignore.case = TRUE)) == 1) {
    atable$factors[["CONTROL"]] = grep("^control", colnames(msdata), value = TRUE, ignore.case = TRUE)
  }
  return(list(atable = atable , msdata = msdata))
}
