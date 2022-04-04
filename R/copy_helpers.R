#' copy Markdown for 2 Grp.
#' @param workdir directory where to copy file - default is current working directory.
#' @keywords internal
#' @export
#'
copy_2grp_markdown <- function(workdir = getwd()){
  runscripts <- c("_Grp2Analysis.Rmd", "bibliography.bib", "_DiffExpQC.Rmd", "MQ2GRPALLVSALL.R")
  runscripts <- file.path("application/FragPipe2Grp/",runscripts )
  prolfqua:::.scriptCopyHelperVec(runscripts, workdir = workdir, packagename = "prolfquapp")
}


