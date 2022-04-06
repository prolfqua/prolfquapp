#' copy Markdown and runscript for 2 Grp Fragpipe data.
#' @param workdir directory where to copy file - default is current working directory.
#' @keywords internal
#' @export
#'
copy_2grp_FragPipe <- function(workdir = getwd()){
  runscripts <- c("application/_Grp2Analysis.Rmd", "application/bibliography.bib", "application/_DiffExpQC.Rmd", "application/FragPipe2Grp/FRG2GRP.R")
  prolfqua:::.scriptCopyHelperVec(runscripts, workdir = workdir, packagename = "prolfquapp")
}


#' copy Markdown and runscript for 2 Grp Fragpipe data.
#' @param workdir directory where to copy file - default is current working directory.
#' @keywords internal
#' @export
#'
copy_2grp_MaxQuant <- function(workdir = getwd()){
  runscripts <- c("application/_Grp2Analysis.Rmd", "application/bibliography.bib", "application/_DiffExpQC.Rmd", "application/FragPipe2Grp/FRG2GRP.R")
  prolfqua:::.scriptCopyHelperVec(runscripts, workdir = workdir, packagename = "prolfquapp")
}


