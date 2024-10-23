#' copy dockerfile to run the DEA app
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_docker_script <- function(workdir = getwd() ) {
  runscripts <- c(
    "application/bin/prolfquapp_docker.sh"
  )
  # Check the operating system and add the appropriate extension
  prolfqua::scriptCopyHelperVec(runscripts, workdir = workdir, packagename = "prolfquapp")
}


#' copy shellscript to run the DEA app
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_shell_script <- function(workdir = getwd() ) {
  runscripts <- c(
    "application/bin/prolfqua_dea",
    "application/bin/prolfqua_yaml",
    "application/bin/prolfqua_qc",
    "application/bin/prolfqua_dataset"
  )
  # Check the operating system and add the appropriate extension
  if (.Platform$OS.type == "windows") {
    runscripts <- paste0(runscripts, ".bat")
  } else {
    runscripts <- paste0(runscripts, ".sh")
  }
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


#' copy Markdown and runscripts for DEA
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_DEA_Metabolomics_Files <- function(workdir = getwd()) {
  runscripts <- c(
    "application/_Grp2Analysis_Metabolomics.Rmd",
    "application/bibliography.bib",
    "application/_DiffExpQC.Rmd"
  )
  prolfqua::scriptCopyHelperVec(runscripts, workdir = workdir, packagename = "prolfquapp")
}

