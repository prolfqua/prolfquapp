#' copy Markdown and runscript for FragPipe combined_protein.tsv
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_DEA_FragPipe_DDA <- function(workdir = getwd(), run_script = FALSE) {
  runscripts <- c(
    "application/_Grp2Analysis.Rmd",
    "application/bibliography.bib",
    "application/_DiffExpQC.Rmd",
    if (run_script) {"application/FragPipeDDA/FP_DDA.R"}
  )
  prolfqua:::.scriptCopyHelperVec(runscripts, workdir = workdir, packagename = "prolfquapp")
}

#' copy Markdown and runscript for MaxQuant peptide.txt proteinGroups.txt
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_DEA_MaxQuant <- function(workdir = getwd(), run_script = FALSE) {
  runscripts <- c(
    "application/_Grp2Analysis.Rmd",
    "application/bibliography.bib",
    "application/_DiffExpQC.Rmd",
    if (run_script) {"application/MaxQuantDDA/MQ_DDA.R"}
  )
  prolfqua:::.scriptCopyHelperVec(runscripts, workdir = workdir, packagename = "prolfquapp")
}


#' copy Markdown and runscript for Fragpipe TMT psm.txt
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_DEA_FragPipe_TMT <- function(workdir = getwd(), run_script = FALSE) {
  runscripts <- c(
    "application/_Grp2Analysis.Rmd",
    "application/bibliography.bib",
    "application/_DiffExpQC.Rmd",
    if (run_script) { "application/FragPipeTMT/FP_TMT_V2.R" }
  )
  prolfqua:::.scriptCopyHelperVec(runscripts, workdir = workdir, packagename = "prolfquapp")
}

#' copy Markdown and runscript for DIANN diann-output.tsv
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_DEA_DIANN <- function(workdir = getwd(), run_script = FALSE) {
  runscripts <- c(
    "application/_Grp2Analysis.Rmd",
    "application/bibliography.bib",
    "application/_DiffExpQC.Rmd",
    if (run_script) {"application/DIANN/DIANN.R"}
  )
  prolfqua:::.scriptCopyHelperVec(runscripts, workdir = workdir, packagename = "prolfquapp")
}

