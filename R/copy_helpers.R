#' copy shellscript to run the DEA app
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_shell_script <- function(workdir = getwd()) {
  scripts <- c(
    "prolfqua_dea",
    "prolfqua_dea_cd",
    "prolfqua_yaml",
    "prolfqua_qc",
    "prolfqua_dataset",
    "prolfqua_contrasts"
  )
  runscripts <- paste0("application/bin/", scripts, if (.Platform$OS.type == "windows") ".bat" else ".sh")
  script_copy_helper_vec(runscripts, workdir = workdir, packagename = "prolfquapp")
}


#' Copy files shipped in a package into a working directory
#'
#' Used by the `copy_*()` helpers of prolfquapp and of the packages built on it (prophosqua, prolfquasaint) to
#' copy run scripts and report templates out of the installed package.
#' @param runscripts paths of the files, relative to the package's installation directory
#' @param workdir directory to copy the files to, default the working directory
#' @param packagename package to copy the files from
#' @param overwrite overwrite existing files in `workdir`, default TRUE
#' @return paths of the copied files; a file that is missing or cannot be copied gives a warning and is left out
#' @export
#' @examples
#' copied <- script_copy_helper_vec("DESCRIPTION", workdir = tempdir(), packagename = "prolfquapp")
#' basename(copied)
script_copy_helper_vec <- function(runscripts, workdir = getwd(), packagename = "prolfquapp", overwrite = TRUE) {
  res <- NULL
  for (scripts in runscripts) {
    src_script <- prolfqua::find_package_file(packagename, scripts)
    dest_script <- file.path(workdir, basename(scripts))
    logger::log_info("copy ", src_script, " to ", dest_script)
    if (is.null(src_script)) {
      warning(sprintf("Could not copy script file %s.", dest_script), call. = FALSE)
    } else if (file.copy(src_script, dest_script, overwrite = overwrite)) {
      res <- c(res, dest_script)
    } else {
      warning(sprintf("Could not copy script file %s to %s.", src_script, dest_script), call. = FALSE)
    }
  }
  logger::log_info(sprintf("your working directory now should contain: %d new files:", length(res)))
  res
}
