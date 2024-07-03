#' copy shellscript to run the DEA app
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_shell_script <- function(workdir = getwd()) {
  runscripts <- c(
    "application/bin/CMD_DEA.sh",
    "application/bin/CMD_DEA.bat"
  )
  prolfqua::scriptCopyHelperVec(runscripts, workdir = workdir, packagename = "prolfquapp")
}


#' copy Markdown and runscripts for DEA
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_DEA_Files <- function(workdir = getwd()) {
  runscripts <- c(
    "application/_Grp2Analysis_V2.Rmd",
    "application/bibliography.bib",
    "application/_DiffExpQC.Rmd"
  )
  prolfqua::scriptCopyHelperVec(runscripts, workdir = workdir, packagename = "prolfquapp")
}

