#' copy dockerfile to run the DEA app
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_docker_script <- function(workdir = getwd()) {
  runscripts <- c(
    "application/bin/prolfquapp_docker.sh"
  )
  # Check the operating system and add the appropriate extension
  prolfqua::script_copy_helper_vec(
    runscripts,
    workdir = workdir,
    packagename = "prolfquapp"
  )
}


#' copy shellscript to run the DEA app
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_shell_script <- function(workdir = getwd()) {
  runscripts <- c(
    "application/bin/prolfqua_dea",
    "application/bin/prolfqua_yaml",
    "application/bin/prolfqua_qc",
    "application/bin/prolfqua_dataset",
    "application/bin/prolfqua_contrasts"
  )
  # Check the operating system and add the appropriate extension
  if (.Platform$OS.type == "windows") {
    runscripts <- paste0(runscripts, ".bat")
  } else {
    runscripts <- paste0(runscripts, ".sh")
  }
  prolfqua::script_copy_helper_vec(
    runscripts,
    workdir = workdir,
    packagename = "prolfquapp"
  )
}


#' copy Markdown templates for DEA (R6-based CMD_DEA_V2.R)
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_DEA_R6_Files <- function(workdir = getwd()) {
  runscripts <- c(
    "doc/Grp2Analysis_V2_R6.Rmd",
    "doc/DiffExpQC_R6.Rmd"
  )
  prolfqua::script_copy_helper_vec(
    runscripts,
    workdir = workdir,
    packagename = "prolfquapp"
  )
  # bibliography.bib lives in application/, not doc/
  bib_src <- system.file("application", "bibliography.bib", package = "prolfquapp")
  if (nzchar(bib_src) && file.exists(bib_src)) {
    file.copy(bib_src, file.path(workdir, "bibliography.bib"), overwrite = TRUE)
  }
}


