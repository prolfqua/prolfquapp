# Internal Quarto report helpers.

render_quarto_se_report <- function(
  se_file,
  output_dir,
  output_file = "Grp2Analysis_V2_SE.html",
  template = "Grp2Analysis_V2_SE.qmd",
  template_dir = system.file("templates/quarto", package = "prolfquapp"),
  env = character()
) {
  if (!file.exists(se_file)) {
    stop("SummarizedExperiment file not found: ", se_file, call. = FALSE)
  }
  if (!nzchar(template_dir) || !dir.exists(template_dir)) {
    stop("Quarto template directory not found.", call. = FALSE)
  }
  template_file <- file.path(template_dir, template)
  if (!file.exists(template_file)) {
    stop("Quarto template not found: ", template_file, call. = FALSE)
  }
  quarto <- Sys.which("quarto")
  if (!nzchar(quarto)) {
    logger::log_warn("Quarto CLI not found; skipping SE Quarto DEA report.")
    return(NULL)
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  render_dir <- tempfile("prolfquapp_quarto_report_")
  dir.create(render_dir, recursive = TRUE)
  on.exit(unlink(render_dir, recursive = TRUE), add = TRUE)

  support_files <- c("_fgcz-report.yml", "fgcz_header_quarto.html", template)
  missing_files <- support_files[!file.exists(file.path(template_dir, support_files))]
  if (length(missing_files) > 0) {
    stop(
      "Quarto template support file(s) not found: ",
      paste(missing_files, collapse = ", "),
      call. = FALSE
    )
  }
  file.copy(file.path(template_dir, support_files), render_dir, overwrite = TRUE)

  oldwd <- setwd(render_dir)
  on.exit(setwd(oldwd), add = TRUE)

  status <- system2(
    quarto,
    c("render", template, "-P", paste0("se_file:", normalizePath(se_file, mustWork = TRUE))),
    stdout = TRUE,
    stderr = TRUE,
    env = env
  )
  exit_code <- attr(status, "status")
  if (!is.null(exit_code) && exit_code != 0) {
    logger::log_error("Quarto render output:\n", paste(status, collapse = "\n"))
    stop("Quarto SE DEA report rendering failed.", call. = FALSE)
  }

  rendered_file <- file.path(render_dir, sub("[.]qmd$", ".html", template))
  if (!file.exists(rendered_file)) {
    stop("Quarto render did not create expected HTML file: ", rendered_file, call. = FALSE)
  }
  out_file <- file.path(output_dir, output_file)
  file.copy(rendered_file, out_file, overwrite = TRUE)
  out_file
}
