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
  prolfqua::script_copy_helper_vec(runscripts, workdir = workdir, packagename = "prolfquapp")
}
