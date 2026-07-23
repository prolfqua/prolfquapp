fake_dea_reporter <- function(resultdir, supports_qc = TRUE) {
  lfq_data_raw <- list(
    data_long = function() {
      data.frame(
        protein_Id = "P1",
        sampleName = "S1",
        abundance = 10
      )
    },
    get_config = function() list(hierarchy_depth = 1)
  )
  model_result <- list(
    get_config = function() list(supports_dea_qc = supports_qc)
  )
  deanalyse <- list(
    FDR_threshold = 0.05,
    diff_threshold = 1,
    default_model = "lm",
    contrast_results = list(lm = model_result),
    lfq_data_raw = lfq_data_raw,
    prolfq_app_config = list(
      software = "DIA-NN",
      project_spec = list(
        project_Id = "123",
        project_name = "Example",
        order_Id = "456",
        workunit_Id = "789",
        input_URL = "https://example.org/input"
      )
    )
  )
  list(
    resultdir = resultdir,
    deanalyse = deanalyse,
    make_SummarizedExperiment = function() list(assay = "test")
  )
}

test_that("render_dea_reports writes inputs and dispatches every supported report", {
  resultdir <- tempfile("dea-reports-")
  dir.create(resultdir)
  on.exit(unlink(resultdir, recursive = TRUE), add = TRUE)
  reporter <- fake_dea_reporter(resultdir)
  calls <- new.env(parent = emptyenv())

  record_render <- function(kind) {
    force(kind)
    function(..., output_dir, output_file) {
      calls[[kind]] <- list(...)
      file.path(output_dir, output_file)
    }
  }
  local_mocked_bindings(
    render_quarto_dea_report = record_render("dea"),
    render_quarto_se_report = record_render("tabset"),
    render_quarto_diffexpqc_report = record_render("qc"),
    render_quarto_qc_sse_report = record_render("sse"),
    .package = "prolfquapp"
  )

  result <- prolfquapp:::render_dea_reports(reporter)

  expect_named(
    result,
    c(
      "deanalyse_file",
      "se_file",
      "dea_file",
      "tabset_file",
      "qc_file",
      "sse_file"
    )
  )
  expect_true(file.exists(result$deanalyse_file))
  expect_true(file.exists(result$se_file))
  expect_equal(
    basename(result$dea_file),
    "Grp2Analysis_V2_R6.html"
  )
  expect_equal(
    basename(result$tabset_file),
    "Grp2Analysis_V2_SE_tabset.html"
  )
  expect_equal(
    basename(result$qc_file),
    "DiffExpQC_R6_tabset.html"
  )
  expect_equal(
    basename(result$sse_file),
    "QCandSSE_tabset.html"
  )
  expect_equal(calls$tabset$fdr_threshold, 0.05)
  expect_equal(calls$tabset$diff_threshold, 1)
  expect_equal(calls$sse$project_conf$workunit_Id, "789")
  expect_equal(calls$sse$project_conf$software, "DIA-NN")
  expect_equal(calls$sse$project_conf$model, "lm")
})

test_that("render_dea_reports isolates render failures and skips unsupported QC", {
  resultdir <- tempfile("dea-reports-failure-")
  dir.create(resultdir)
  on.exit(unlink(resultdir, recursive = TRUE), add = TRUE)
  reporter <- fake_dea_reporter(
    resultdir,
    supports_qc = FALSE
  )

  local_mocked_bindings(
    render_quarto_dea_report = function(...) {
      stop("render failed")
    },
    render_quarto_se_report = function(
      ...,
      output_dir,
      output_file
    ) {
      file.path(output_dir, output_file)
    },
    render_quarto_diffexpqc_report = function(...) {
      stop("QC renderer must not be called")
    },
    render_quarto_qc_sse_report = function(
      ...,
      output_dir,
      output_file
    ) {
      file.path(output_dir, output_file)
    },
    .package = "prolfquapp"
  )

  expect_no_error(
    result <- prolfquapp:::render_dea_reports(reporter)
  )

  expect_null(result$dea_file)
  expect_null(result$qc_file)
  expect_true(file.exists(result$deanalyse_file))
  expect_true(file.exists(result$se_file))
  expect_equal(
    basename(result$tabset_file),
    "Grp2Analysis_V2_SE_tabset.html"
  )
  expect_equal(
    basename(result$sse_file),
    "QCandSSE_tabset.html"
  )
})

test_that("Quarto report wrappers validate inputs and forward parameters", {
  input_file <- tempfile(fileext = ".rds")
  saveRDS(list(test = TRUE), input_file)
  on.exit(unlink(input_file), add = TRUE)
  calls <- list()

  local_mocked_bindings(
    .render_quarto_doc_report = function(
      qmd_name,
      execute_params,
      output_dir,
      output_file
    ) {
      calls[[length(calls) + 1]] <<- list(
        qmd_name = qmd_name,
        execute_params = execute_params,
        output_dir = output_dir,
        output_file = output_file
      )
      file.path(output_dir, output_file)
    },
    .package = "prolfquapp"
  )

  prolfquapp:::render_quarto_se_report(
    input_file,
    tempdir(),
    fdr_threshold = 0.1,
    diff_threshold = 2
  )
  prolfquapp:::render_quarto_dea_report(input_file, tempdir())
  prolfquapp:::render_quarto_qc_sse_report(
    input_file,
    tempdir(),
    project_conf = list(workunit_Id = "123"),
    target_type = "peptide",
    plot_density = FALSE,
    plot_sd_vs_mean = TRUE
  )
  prolfquapp:::render_quarto_diffexpqc_report(
    input_file,
    tempdir()
  )
  prolfquapp:::render_quarto_protein_abundances_report(
    input_file,
    tempdir(),
    project_info = list(order_Id = "456"),
    factors = FALSE
  )

  expect_equal(
    vapply(calls, `[[`, character(1), "qmd_name"),
    c(
      "Grp2Analysis_V2_SE_tabset.qmd",
      "Grp2Analysis_V2_R6.qmd",
      "QCandSSE_tabset.qmd",
      "DiffExpQC_R6_tabset.qmd",
      "QC_ProteinAbundances_tabset.qmd"
    )
  )
  expect_equal(calls[[1]]$execute_params$fdr_threshold, 0.1)
  expect_equal(calls[[1]]$execute_params$diff_threshold, 2)
  expect_equal(calls[[3]]$execute_params$target_type, "peptide")
  expect_false(calls[[3]]$execute_params$plot_density)
  expect_true(calls[[3]]$execute_params$plot_sd_vs_mean)
  expect_false(calls[[5]]$execute_params$factors)

  expect_error(
    prolfquapp:::render_quarto_se_report(
      "/does/not/exist.rds",
      tempdir()
    ),
    "SummarizedExperiment file not found"
  )
  expect_error(
    prolfquapp:::render_quarto_protein_abundances_report(
      "/does/not/exist.rds",
      tempdir()
    ),
    "QC protein-abundances data file not found"
  )
})
