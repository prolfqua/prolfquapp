# Internal Quarto report helpers.

render_quarto_se_report <- function(
  se_file,
  output_dir,
  output_file = "Grp2Analysis_V2_SE_tabset.html",
  template = "Grp2Analysis_V2_SE_tabset.qmd",
  template_dir = system.file("templates/quarto", package = "prolfquapp"),
  prolfquapp_source_path = NULL,
  buttons = TRUE
) {
  if (!file.exists(se_file)) {
    stop("SummarizedExperiment file not found: ", se_file, call. = FALSE)
  }
  se_file <- normalizePath(se_file, mustWork = TRUE)
  if (!is.null(prolfquapp_source_path)) {
    prolfquapp_source_path <- normalizePath(
      prolfquapp_source_path,
      mustWork = TRUE
    )
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
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  render_dir <- tempfile("prolfquapp_quarto_report_")
  dir.create(render_dir, recursive = TRUE)
  on.exit(unlink(render_dir, recursive = TRUE), add = TRUE)

  file.copy(
    template_file,
    render_dir,
    overwrite = TRUE
  )

  oldwd <- setwd(render_dir)
  on.exit(setwd(oldwd), add = TRUE)

  execute_params <- list(se_file = se_file)
  if (!is.null(prolfquapp_source_path)) {
    execute_params$prolfquapp_source_path <- prolfquapp_source_path
  }

  render_result <- tryCatch(
    fgczquartotemplate::fgcz_render(
      template,
      buttons = buttons,
      execute_params = execute_params
    ),
    error = function(e) e
  )
  if (inherits(render_result, "error")) {
    setwd(oldwd)
    logger::log_error("Quarto render failed: ", conditionMessage(render_result))
    stop("Quarto SE DEA report rendering failed.", call. = FALSE)
  }

  rendered_file <- file.path(render_dir, sub("[.]qmd$", ".html", template))
  if (!file.exists(rendered_file)) {
    stop(
      "Quarto render did not create expected HTML file: ",
      rendered_file,
      call. = FALSE
    )
  }
  out_file <- file.path(output_dir, output_file)
  copied <- file.copy(rendered_file, out_file, overwrite = TRUE)
  if (!isTRUE(copied) || !file.exists(out_file)) {
    stop(
      "Could not copy Quarto report to output file: ",
      out_file,
      call. = FALSE
    )
  }
  out_file
}
