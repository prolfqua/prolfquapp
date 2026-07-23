# Internal Quarto report helpers.
#
# The Quarto report sources live in `vignettes/` as plain `format: html` reports
# whose FGCZ styling (SCSS theme, header) comes from a directory-level
# `_metadata.yml`, and whose Find/Download toolbar is wired via
# `include-after-body: fgcz-plot-finder.html`. `data-raw/sync_quarto_assets.R`
# keeps these build-time assets synchronized with `fgczquartotemplate`. The file
# `vignettes/.install_extras` ships the qmd sources into the installed package's
# `doc/` directory, so they are reachable at runtime via `system.file("doc", ...)`.
#
# Rendering copies the qmd out of `doc/` into a temporary directory and calls
# `fgczquartotemplate::fgcz_render()`, which stages the FGCZ assets
# (`_metadata.yml`, `fgcz.scss`, `fgcz_header_quarto.html`, `fgcz-plot-finder.html`)
# next to the qmd from the installed `fgczquartotemplate` package before rendering.
# Runtime rendering therefore receives the same template assets without depending
# on an `_extensions/` tree being shipped into `doc/`. `buttons = FALSE`: the
# toolbar is already wired by the report's own `include-after-body`, so
# `fgcz_render()` must not inject it a second time.
#
# Because the sources are shipped via the vignette machinery, this requires the
# package to have been installed with vignettes built (the default for
# `R CMD build` + `R CMD INSTALL`, which is how `make install` installs). If the
# sources are absent from `doc/`, rendering is skipped with a warning.

.quarto_visual_abstract_names <- c(
  "differential-expression.png",
  "differential-expression-tabset.png",
  "differential-expression-qc.png",
  "protein-abundances.png",
  "quality-control-sample-size.png"
)

.render_quarto_doc_report <- function(
  qmd_name,
  execute_params,
  output_dir,
  output_file
) {
  quarto <- Sys.which("quarto")
  if (!nzchar(quarto)) {
    logger::log_warn("Quarto CLI not found; skipping ", qmd_name, " report.")
    return(NULL)
  }
  qmd_src <- system.file("doc", qmd_name, package = "prolfquapp")
  if (!nzchar(qmd_src) || !file.exists(qmd_src)) {
    logger::log_warn(
      "Quarto report source not found in installed doc/: ",
      qmd_name,
      " (install prolfquapp with vignettes built). Skipping."
    )
    return(NULL)
  }
  bib_src <- system.file("doc/bibliography.bib", package = "prolfquapp")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  render_dir <- tempfile("prolfquapp_quarto_report_")
  dir.create(render_dir, recursive = TRUE)
  on.exit(unlink(render_dir, recursive = TRUE), add = TRUE)

  file.copy(qmd_src, render_dir, overwrite = TRUE)
  if (nzchar(bib_src) && file.exists(bib_src)) {
    file.copy(bib_src, render_dir, overwrite = TRUE)
  }
  oldwd <- setwd(render_dir)
  on.exit(setwd(oldwd), add = TRUE)

  # fgcz_render() stages the FGCZ assets from the installed fgczquartotemplate
  # package next to the qmd, then renders. buttons = FALSE because the report's
  # own include-after-body already wires the Find/Download toolbar.
  render_result <- tryCatch(
    fgczquartotemplate::fgcz_render(
      input = qmd_name,
      buttons = FALSE,
      execute_params = execute_params
    ),
    error = function(e) e
  )
  if (inherits(render_result, "error")) {
    setwd(oldwd)
    logger::log_error("Quarto render failed: ", conditionMessage(render_result))
    stop("Quarto report rendering failed for ", qmd_name, ".", call. = FALSE)
  }

  rendered_file <- file.path(render_dir, sub("[.]qmd$", ".html", qmd_name))
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

# Render the SummarizedExperiment-backed tabbed Quarto report from a serialized
# SummarizedExperiment `.rds`.
render_quarto_se_report <- function(
  se_file,
  output_dir,
  output_file = "Grp2Analysis_V2_SE_tabset.html",
  fdr_threshold = 0.05,
  diff_threshold = 1
) {
  if (!file.exists(se_file)) {
    stop("SummarizedExperiment file not found: ", se_file, call. = FALSE)
  }
  se_file <- normalizePath(se_file, mustWork = TRUE)
  .render_quarto_doc_report(
    qmd_name = "Grp2Analysis_V2_SE_tabset.qmd",
    execute_params = list(
      se_file = se_file,
      fdr_threshold = fdr_threshold,
      diff_threshold = diff_threshold
    ),
    output_dir = output_dir,
    output_file = output_file
  )
}

# Render the DEAnalyse-backed Grp2 Quarto report from a serialized DEAnalyse
# `.rds`.
render_quarto_dea_report <- function(
  deanalyse_file,
  output_dir,
  output_file = "Grp2Analysis_V2_R6.html"
) {
  if (!file.exists(deanalyse_file)) {
    stop("DEAnalyse file not found: ", deanalyse_file, call. = FALSE)
  }
  deanalyse_file <- normalizePath(deanalyse_file, mustWork = TRUE)
  .render_quarto_doc_report(
    qmd_name = "Grp2Analysis_V2_R6.qmd",
    execute_params = list(deanalyse_file = deanalyse_file),
    output_dir = output_dir,
    output_file = output_file
  )
}

# Render the QC & sample-size-estimation Quarto report from a serialized
# `list(data, configuration)` `.rds` (the structure the QCandSSE Quarto report
# expects via its `qc_data_file` parameter). `project_conf` must be a plain list
# of scalars (e.g. project_Id / order_Id / workunit_Id): Quarto serializes
# execute params to YAML, so an R6 project-spec object cannot be passed through.
render_quarto_qc_sse_report <- function(
  qc_data_file,
  output_dir,
  output_file = "QCandSSE_tabset.html",
  project_conf = NULL,
  target_type = "protein",
  plot_density = TRUE,
  plot_sd_vs_mean = FALSE
) {
  if (!file.exists(qc_data_file)) {
    stop("QC data file not found: ", qc_data_file, call. = FALSE)
  }
  qc_data_file <- normalizePath(qc_data_file, mustWork = TRUE)
  .render_quarto_doc_report(
    qmd_name = "QCandSSE_tabset.qmd",
    execute_params = list(
      qc_data_file = qc_data_file,
      project_conf = project_conf,
      target_type = target_type,
      plot_density = plot_density,
      plot_sd_vs_mean = plot_sd_vs_mean
    ),
    output_dir = output_dir,
    output_file = output_file
  )
}

# Render the differential-expression QC (tabbed) Quarto report from a serialized
# DEAnalyse `.rds`.
render_quarto_diffexpqc_report <- function(
  deanalyse_file,
  output_dir,
  output_file = "DiffExpQC_R6_tabset.html"
) {
  if (!file.exists(deanalyse_file)) {
    stop("DEAnalyse file not found: ", deanalyse_file, call. = FALSE)
  }
  deanalyse_file <- normalizePath(deanalyse_file, mustWork = TRUE)
  .render_quarto_doc_report(
    qmd_name = "DiffExpQC_R6_tabset.qmd",
    execute_params = list(deanalyse_file = deanalyse_file),
    output_dir = output_dir,
    output_file = output_file
  )
}

# Render the QC protein-abundances (tabbed) Quarto report from a serialized
# QC_generator `.rds`. `project_info` must be a plain list of scalars (order_Id,
# workunit_Id): Quarto serializes execute params to YAML.
render_quarto_protein_abundances_report <- function(
  pap_file,
  output_dir,
  output_file = "QC_ProteinAbundances_tabset.html",
  project_info = NULL,
  factors = TRUE
) {
  if (!file.exists(pap_file)) {
    stop("QC protein-abundances data file not found: ", pap_file, call. = FALSE)
  }
  pap_file <- normalizePath(pap_file, mustWork = TRUE)
  .render_quarto_doc_report(
    qmd_name = "QC_ProteinAbundances_tabset.qmd",
    execute_params = list(
      pap_file = pap_file,
      project_info = project_info,
      factors = factors
    ),
    output_dir = output_dir,
    output_file = output_file
  )
}

# Coerce a project-spec field to a length-1 character, or NULL when unset, so it
# survives Quarto's YAML execute-param serialization (and the report falls back
# to "n/a").
.as_project_id <- function(x) {
  if (length(x) >= 1 && nzchar(as.character(x)[[1]])) {
    as.character(x)[[1]]
  } else {
    NULL
  }
}

.try_report_step <- function(expr, label) {
  tryCatch(
    expr,
    error = function(error) {
      logger::log_warn(
        "prolfquapp: skipping ",
        label,
        ": ",
        conditionMessage(error)
      )
      NULL
    }
  )
}

.write_dea_report_inputs <- function(reporter) {
  deanalyse_file <- file.path(
    reporter$resultdir,
    "DEAnalyse.rds"
  )
  se_file <- file.path(
    reporter$resultdir,
    "SummarizedExperiment.rds"
  )
  list(
    deanalyse_file = .try_report_step(
      {
        saveRDS(reporter$deanalyse, file = deanalyse_file)
        deanalyse_file
      },
      "DEAnalyse.rds"
    ),
    se_file = .try_report_step(
      {
        saveRDS(
          reporter$make_SummarizedExperiment(),
          file = se_file
        )
        se_file
      },
      "SummarizedExperiment.rds"
    )
  )
}

.deanalyse_supports_qc <- function(deanalyse) {
  tryCatch(
    isTRUE(
      deanalyse$contrast_results[[
        deanalyse$default_model
      ]]$get_config()$supports_dea_qc
    ),
    error = function(error) FALSE
  )
}

.dea_report_project_conf <- function(deanalyse) {
  project_spec <- deanalyse$prolfq_app_config$project_spec
  list(
    project_Id = .as_project_id(project_spec$project_Id),
    project_name = .as_project_id(project_spec$project_name),
    order_Id = .as_project_id(project_spec$order_Id),
    workunit_Id = .as_project_id(project_spec$workunit_Id),
    input_URL = .as_project_id(project_spec$input_URL),
    software = .as_project_id(
      deanalyse$prolfq_app_config$software
    ),
    model = .as_project_id(deanalyse$default_model)
  )
}

.render_dea_sse_report <- function(reporter) {
  qc_data_file <- file.path(
    reporter$resultdir,
    "QC_sampleSizeEstimation.rds"
  )
  saveRDS(
    list(
      data = reporter$deanalyse$lfq_data_raw$data_long(),
      configuration = reporter$deanalyse$lfq_data_raw$get_config()
    ),
    file = qc_data_file
  )
  render_quarto_qc_sse_report(
    qc_data_file = qc_data_file,
    output_dir = reporter$resultdir,
    output_file = "QCandSSE_tabset.html",
    project_conf = .dea_report_project_conf(
      reporter$deanalyse
    ),
    target_type = "protein"
  )
}

# Render the full set of DEA Quarto reports from a DEAReportGenerator and write
# the SummarizedExperiment + DEAnalyse `.rds`. Each report is rendered
# independently, so one failure logs a warning without dropping the others.
# Returns a named list of output paths (NULL for any that could not be produced).
render_dea_reports <- function(reporter) {
  out <- .write_dea_report_inputs(reporter)

  # Primary DEA report (R6 Quarto).
  out$dea_file <- if (!is.null(out$deanalyse_file)) {
    .try_report_step(
      render_quarto_dea_report(
        deanalyse_file = out$deanalyse_file,
        output_dir = reporter$resultdir,
        output_file = "Grp2Analysis_V2_R6.html"
      ),
      "DEA Quarto report"
    )
  }

  # Secondary overview report (SE tabset).
  out$tabset_file <- if (!is.null(out$se_file)) {
    .try_report_step(
      render_quarto_se_report(
        se_file = out$se_file,
        output_dir = reporter$resultdir,
        output_file = "Grp2Analysis_V2_SE_tabset.html",
        fdr_threshold = reporter$deanalyse$FDR_threshold,
        diff_threshold = reporter$deanalyse$diff_threshold
      ),
      "SE tabset report"
    )
  }

  # Differential-expression QC report, when the model supports it.
  supports_qc <- .deanalyse_supports_qc(reporter$deanalyse)
  out$qc_file <- if (!is.null(out$deanalyse_file) && supports_qc) {
    .try_report_step(
      render_quarto_diffexpqc_report(
        deanalyse_file = out$deanalyse_file,
        output_dir = reporter$resultdir,
        output_file = "DiffExpQC_R6_tabset.html"
      ),
      "DEA QC report"
    )
  }

  # Sample-size (SSE) report, built from the raw feature-level data.
  out$sse_file <- .try_report_step(
    .render_dea_sse_report(reporter),
    "sample-size report"
  )

  out
}
