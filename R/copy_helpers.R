#' copy Markdown and runscript for 2 Grp Fragpipe data.
#' @param workdir directory where to copy file - default is current working directory.
#' @keywords internal
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

#' copy Markdown and runscript for 2 Grp Fragpipe data.
#' @param workdir directory where to copy file - default is current working directory.
#' @keywords internal
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


#' copy Markdown and runscript for DEA Fragpipe TMT
#' @param workdir directory where to copy file - default is current working directory.
#' @keywords internal
#' @export
#'
copy_DEA_FragPipe_TMT <- function(workdir = getwd(), run_script = FALSE) {
  runscripts <- c(
    "application/_Grp2Analysis.Rmd",
    "application/bibliography.bib",
    "application/_DiffExpQC.Rmd",
    if (run_script) { "application/FragPipeTMT/FP_TMT.R" }
  )
  prolfqua:::.scriptCopyHelperVec(runscripts, workdir = workdir, packagename = "prolfquapp")
}

#' copy Markdown and runscripts for DEA Fragpipe DIA
#' @param workdir directory where to copy file - default is current working directory.
#' @keywords internal
#' @export
#'
copy_DEA_FragPipe_DIA <- function(workdir = getwd(), run_script = FALSE) {
  runscripts <- c(
    "application/_Grp2Analysis.Rmd",
    "application/bibliography.bib",
    "application/_DiffExpQC.Rmd",
    if (run_script) {"application/FragPipeDIA/FP_DIA.R"}
  )
  prolfqua:::.scriptCopyHelperVec(runscripts, workdir = workdir, packagename = "prolfquapp")
}

