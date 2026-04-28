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
  files <- c("Grp2Analysis_V2_R6.Rmd", "DiffExpQC_R6.Rmd")
  copied <- character()
  for (file in files) {
    src <- system.file("doc", file, package = "prolfquapp", mustWork = FALSE)
    if (!nzchar(src) || !file.exists(src)) {
      stop(
        "DEA report template not found in installed prolfquapp package: doc/",
        file,
        "\nReinstall prolfquapp with vignettes so installed package doc/ files are available.",
        call. = FALSE
      )
    }
    dest <- file.path(workdir, file)
    message("copy ", src, " to ", dest)
    if (!file.copy(src, dest, overwrite = TRUE)) {
      stop("Could not copy DEA report template from ", src, " to ", dest, call. = FALSE)
    }
    copied <- c(copied, dest)
  }
  message("your working directory now should contain: ", length(copied), " new files :\n")
  # bibliography.bib lives in application/, not doc/
  bib_src <- system.file("application", "bibliography.bib", package = "prolfquapp")
  if (nzchar(bib_src) && file.exists(bib_src)) {
    file.copy(bib_src, file.path(workdir, "bibliography.bib"), overwrite = TRUE)
  }
  invisible(copied)
}

