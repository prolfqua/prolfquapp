# Internal Quarto report helpers.
#
# The Quarto report sources live in `vignettes/` as `format: html` reports; FGCZ
# styling comes from a directory-level `_metadata.yml` and the Find/Download
# toolbar from `include-after-body: fgcz-plot-finder.html` (both kept in sync
# with fgczQuartoTemplate by `data-raw/sync_quarto_assets.R`).
# `vignettes/.install_extras` ships the qmd sources into the installed `doc/`,
# so rendering requires prolfquapp installed with vignettes built; without the
# sources, rendering is skipped with a warning. Rendering copies the qmd into a
# temporary directory and calls `fgczQuartoTemplate::fgcz_render()`, which
# stages the FGCZ assets next to it; `buttons = FALSE` because the report wires
# the toolbar itself.

.quarto_visual_abstract_names <- c(
  "differential-expression.png",
  "differential-expression-tabset.png",
  "differential-expression-qc.png",
  "protein-abundances.png",
  "quality-control-sample-size.png"
)

.render_quarto_doc_report <- function(qmd_name, execute_params, output_dir, output_file) {
  if (!nzchar(Sys.which("quarto"))) {
    logger::log_warn("Quarto CLI not found; skipping ", qmd_name, " report.")
    return(NULL)
  }
  qmd_src <- system.file("doc", qmd_name, package = "prolfquapp")
  if (!file.exists(qmd_src)) {
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
  if (file.exists(bib_src)) {
    file.copy(bib_src, render_dir, overwrite = TRUE)
  }
  oldwd <- setwd(render_dir)
  on.exit(setwd(oldwd), add = TRUE)
  tryCatch(
    fgczQuartoTemplate::fgcz_render(input = qmd_name, buttons = FALSE, execute_params = execute_params),
    error = function(e) {
      logger::log_error("Quarto render failed: ", conditionMessage(e))
      stop("Quarto report rendering failed for ", qmd_name, ".", call. = FALSE)
    }
  )

  rendered_file <- file.path(render_dir, sub("[.]qmd$", ".html", qmd_name))
  if (!file.exists(rendered_file)) {
    stop("Quarto render did not create expected HTML file: ", rendered_file, call. = FALSE)
  }
  out_file <- file.path(output_dir, output_file)
  if (!file.copy(rendered_file, out_file, overwrite = TRUE)) {
    stop("Could not copy Quarto report to output file: ", out_file, call. = FALSE)
  }
  out_file
}

# The normalized path of an existing report input file.
.report_input <- function(path, what) {
  if (!file.exists(path)) {
    stop(what, " not found: ", path, call. = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

# Render the tabbed Quarto report from a DEA result artifact -- an `AnnData.h5ad`
# or a serialized `SummarizedExperiment`; the report reads either through
# `DEAResultReader`.
render_quarto_se_report <- function(
  artifact_file,
  output_dir,
  output_file = "Grp2Analysis_V2_SE_tabset.html",
  fdr_threshold = 0.05,
  diff_threshold = 1
) {
  params <- list(
    artifact_file = .report_input(artifact_file, "DEA result artifact"),
    fdr_threshold = fdr_threshold,
    diff_threshold = diff_threshold
  )
  .render_quarto_doc_report("Grp2Analysis_V2_SE_tabset.qmd", params, output_dir, output_file)
}

# Render the DEAnalyse-backed Grp2 Quarto report from a serialized DEAnalyse `.rds`.
render_quarto_dea_report <- function(deanalyse_file, output_dir, output_file = "Grp2Analysis_V2_R6.html") {
  params <- list(deanalyse_file = .report_input(deanalyse_file, "DEAnalyse file"))
  .render_quarto_doc_report("Grp2Analysis_V2_R6.qmd", params, output_dir, output_file)
}

# Render the QC & sample-size-estimation Quarto report from a serialized
# `list(data, configuration)` `.rds`. Quarto serializes execute params to YAML,
# so `project_conf` must be a plain list of scalars, not an R6 project spec.
render_quarto_qc_sse_report <- function(
  qc_data_file,
  output_dir,
  output_file = "QCandSSE_tabset.html",
  project_conf = NULL,
  target_type = "protein",
  plot_density = TRUE,
  plot_sd_vs_mean = FALSE
) {
  params <- list(
    qc_data_file = .report_input(qc_data_file, "QC data file"),
    project_conf = project_conf,
    target_type = target_type,
    plot_density = plot_density,
    plot_sd_vs_mean = plot_sd_vs_mean
  )
  .render_quarto_doc_report("QCandSSE_tabset.qmd", params, output_dir, output_file)
}

# Render the differential-expression QC (tabbed) Quarto report from a serialized DEAnalyse `.rds`.
render_quarto_diffexpqc_report <- function(deanalyse_file, output_dir, output_file = "DiffExpQC_R6_tabset.html") {
  params <- list(deanalyse_file = .report_input(deanalyse_file, "DEAnalyse file"))
  .render_quarto_doc_report("DiffExpQC_R6_tabset.qmd", params, output_dir, output_file)
}

# Render the QC protein-abundances (tabbed) Quarto report from a serialized
# QC_generator `.rds`. `project_info` must be a plain list of scalars.
render_quarto_protein_abundances_report <- function(
  pap_file,
  output_dir,
  output_file = "QC_ProteinAbundances_tabset.html",
  project_info = NULL,
  factors = TRUE
) {
  params <- list(
    pap_file = .report_input(pap_file, "QC protein-abundances data file"),
    project_info = project_info,
    factors = factors
  )
  .render_quarto_doc_report("QC_ProteinAbundances_tabset.qmd", params, output_dir, output_file)
}

.try_report_step <- function(expr, label) {
  tryCatch(expr, error = function(error) {
    logger::log_warn("prolfquapp: skipping ", label, ": ", conditionMessage(error))
    NULL
  })
}

# Project fields as length-1 characters, or NULL when unset, so they survive
# Quarto's YAML execute-param serialization (the report falls back to "n/a").
.dea_report_project_conf <- function(deanalyse) {
  ps <- deanalyse$prolfq_app_config$project_spec
  fields <- list(
    project_Id = ps$project_Id,
    project_name = ps$project_name,
    order_Id = ps$order_Id,
    workunit_Id = ps$workunit_Id,
    input_URL = ps$input_URL,
    software = deanalyse$prolfq_app_config$software,
    model = deanalyse$default_model
  )
  lapply(fields, function(x) if (length(x) >= 1 && nzchar(as.character(x)[[1]])) as.character(x)[[1]])
}

# Write DEAnalyse.rds, SummarizedExperiment.rds and AnnData.h5ad and render the
# DEA Quarto reports. Each step runs independently, so one failure logs a
# warning without dropping the others. Returns a named list of output paths
# (NULL for any that could not be produced).
render_dea_reports <- function(reporter, summarized_experiment = reporter$make_SummarizedExperiment()) {
  dea <- reporter$deanalyse
  resultdir <- reporter$resultdir
  save_rds <- function(object, name) {
    path <- file.path(resultdir, name)
    .try_report_step(
      {
        saveRDS(object, file = path)
        path
      },
      name
    )
  }
  out <- list(
    deanalyse_file = save_rds(dea, "DEAnalyse.rds"),
    se_file = save_rds(summarized_experiment, "SummarizedExperiment.rds"),
    # The tabset report renders from this file, so a DEA run that produces a
    # report has demonstrated that its AnnData can be read back.
    anndata_file = .try_report_step(
      write_summarized_experiment_h5ad(summarized_experiment, file.path(resultdir, "AnnData.h5ad")),
      "AnnData.h5ad"
    )
  )

  out$dea_file <- if (!is.null(out$deanalyse_file)) {
    .try_report_step(
      render_quarto_dea_report(
        deanalyse_file = out$deanalyse_file,
        output_dir = resultdir,
        output_file = "Grp2Analysis_V2_R6.html"
      ),
      "DEA Quarto report"
    )
  }
  out$tabset_file <- if (!is.null(out$anndata_file)) {
    .try_report_step(
      render_quarto_se_report(
        artifact_file = out$anndata_file,
        output_dir = resultdir,
        output_file = "Grp2Analysis_V2_SE_tabset.html",
        fdr_threshold = dea$FDR_threshold,
        diff_threshold = dea$diff_threshold
      ),
      "SE tabset report"
    )
  }
  supports_qc <- tryCatch(
    isTRUE(dea$contrast_results[[dea$default_model]]$get_config()$supports_dea_qc),
    error = function(error) FALSE
  )
  out$qc_file <- if (!is.null(out$deanalyse_file) && supports_qc) {
    .try_report_step(
      render_quarto_diffexpqc_report(
        deanalyse_file = out$deanalyse_file,
        output_dir = resultdir,
        output_file = "DiffExpQC_R6_tabset.html"
      ),
      "DEA QC report"
    )
  }
  # Sample-size (SSE) report, built from the raw feature-level data.
  out$sse_file <- .try_report_step(
    {
      qc_data_file <- file.path(resultdir, "QC_sampleSizeEstimation.rds")
      saveRDS(list(data = dea$lfq_data_raw$data_long(), configuration = dea$lfq_data_raw$get_config()), qc_data_file)
      render_quarto_qc_sse_report(
        qc_data_file = qc_data_file,
        output_dir = resultdir,
        output_file = "QCandSSE_tabset.html",
        project_conf = .dea_report_project_conf(dea),
        target_type = "protein"
      )
    },
    "sample-size report"
  )
  out
}
