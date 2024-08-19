#' copy shellscript to run the DEA app
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_shell_script <- function(workdir = getwd()) {
  runscripts <- c(
    "application/bin/CMD_DEA.sh",
    "application/bin/CMD_MAKE_YAML.sh",
    "application/bin/CMD_QUANT_QC.sh",
    "application/bin/CMD_MAKE_DATASET.sh"
  )
  prolfqua::scriptCopyHelperVec(runscripts, workdir = workdir, packagename = "prolfquapp")
}

#' copy bat files to run the DEA app on windows
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_bat_files <- function(workdir = getwd()) {
  runscripts <- c(
    "application/bin/CMD_DEA.bat",
    "application/bin/CMD_MAKE_YAML.bat",
    "application/bin/CMD_QUANT_QC.bat",
    "application/bin/CMD_MAKE_DATASET.bat"
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

