test_package_path <- function() {
  source_path <- Sys.getenv("PROLFQUAPP_SOURCE_PATH", unset = "")
  if (nzchar(source_path) && file.exists(file.path(source_path, "DESCRIPTION"))) {
    return(normalizePath(source_path, mustWork = TRUE))
  }
  wd <- normalizePath(getwd(), mustWork = TRUE)
  if (file.exists(file.path(wd, "DESCRIPTION"))) {
    return(wd)
  }
  package_path <- normalizePath(file.path(system.file(package = "prolfquapp"), ".."), mustWork = TRUE)
  if (!file.exists(file.path(package_path, "DESCRIPTION"))) {
    package_path <- normalizePath(system.file(package = "prolfquapp"), mustWork = TRUE)
  }
  package_path
}

test_that("se_report_lfqdata reconstructs LFQData objects from SummarizedExperiment", {
  skip_on_cran()

  dea <- prolfquapp::example_deanalyse(Nprot = 20)
  reporter <- prolfquapp::DEAReportGenerator$new(dea, dea$prolfq_app_config)
  se <- reporter$make_SummarizedExperiment()

  report <- prolfquapp:::se_report_lfqdata(se)

  expect_true("LFQData" %in% class(report$lfq_raw))
  expect_true("LFQData" %in% class(report$lfq_transformed))
  expect_gt(nrow(report$lfq_raw$data_long()), 0)
  expect_gt(nrow(report$lfq_transformed$data_long()), 0)
  expect_gt(nrow(report$lfq_raw$data_wide(as.matrix = TRUE)$data), 0)
  expect_gt(nrow(report$lfq_transformed$data_wide(as.matrix = TRUE)$data), 0)
  expect_true("LFQDataPlotter" %in% class(report$lfq_raw$get_Plotter()))
  expect_true("LFQDataSummariser" %in% class(report$lfq_raw$get_Summariser()))
  expect_true("LFQDataStats" %in% class(report$lfq_transformed$get_Stats()))
  expect_true("protein_Id" %in% colnames(report$contrast_table))
  expect_true("contrast" %in% colnames(report$contrast_table))
})

test_that("SE Quarto templates render with reconstructed LFQData", {
  skip_on_cran()
  skip_if(Sys.which("quarto") == "", "Quarto CLI not installed")
  skip_if_not_installed("devtools")
  skip_if_not_installed("DT")
  skip_if_not_installed("gridExtra")

  dea <- prolfquapp::example_deanalyse(Nprot = 12)
  reporter <- prolfquapp::DEAReportGenerator$new(dea, dea$prolfq_app_config)
  se <- reporter$make_SummarizedExperiment()

  template_dir <- system.file("templates/quarto", package = "prolfquapp")
  skip_if(!nzchar(template_dir), "Quarto templates not installed")

  workdir <- file.path(tempdir(), "quarto_se_report_test")
  unlink(workdir, recursive = TRUE)
  dir.create(workdir, recursive = TRUE)
  on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

  se_file <- file.path(workdir, "SummarizedExperiment.rds")
  saveRDS(se, se_file)

  support_files <- c("_fgcz-report.yml", "fgcz_header_quarto.html")
  report_files <- "Grp2Analysis_V2_SE.qmd"
  file.copy(file.path(template_dir, c(support_files, report_files)), workdir)

  package_path <- test_package_path()

  r_profile <- file.path(workdir, ".Rprofile")
  writeLines(
    c(
      "source_path <- Sys.getenv('PROLFQUAPP_SOURCE_PATH')",
      "if (nzchar(source_path) && file.exists(file.path(source_path, 'DESCRIPTION'))) {",
      "  suppressPackageStartupMessages(devtools::load_all(source_path, quiet = TRUE))",
      "}"
    ),
    r_profile
  )

  oldwd <- setwd(workdir)
  on.exit(setwd(oldwd), add = TRUE)

  for (report_file in report_files) {
    status <- system2(
      "quarto",
      c("render", report_file, "-P", paste0("se_file:", normalizePath(se_file))),
      env = c(
        paste0("PROLFQUAPP_SOURCE_PATH=", package_path),
        paste0("R_PROFILE_USER=", r_profile)
      ),
      stdout = TRUE,
      stderr = TRUE
    )
    exit_code <- attr(status, "status")
    if (!is.null(exit_code) && exit_code != 0) {
      message("quarto render output:\n", paste(status, collapse = "\n"))
    }
    expect_true(is.null(exit_code) || exit_code == 0)
    html_file <- sub("[.]qmd$", ".html", report_file)
    expect_true(file.exists(html_file))
    html <- paste(readLines(html_file, warn = FALSE), collapse = "\n")
    expect_match(html, "fgcz-banner", fixed = TRUE)
    expect_match(html, "panel-tabset", fixed = TRUE)
  }
})

test_that("internal SE Quarto report helper renders HTML", {
  skip_on_cran()
  skip_if(Sys.which("quarto") == "", "Quarto CLI not installed")
  skip_if_not_installed("devtools")
  skip_if_not_installed("DT")
  skip_if_not_installed("gridExtra")

  dea <- prolfquapp::example_deanalyse(Nprot = 12)
  reporter <- prolfquapp::DEAReportGenerator$new(dea, dea$prolfq_app_config)
  se <- reporter$make_SummarizedExperiment()

  workdir <- file.path(tempdir(), "quarto_se_report_helper_test")
  unlink(workdir, recursive = TRUE)
  dir.create(workdir, recursive = TRUE)
  on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

  se_file <- file.path(workdir, "SummarizedExperiment.rds")
  saveRDS(se, se_file)

  package_path <- test_package_path()

  r_profile <- file.path(workdir, ".Rprofile")
  writeLines(
    c(
      "source_path <- Sys.getenv('PROLFQUAPP_SOURCE_PATH')",
      "if (nzchar(source_path) && file.exists(file.path(source_path, 'DESCRIPTION'))) {",
      "  suppressPackageStartupMessages(devtools::load_all(source_path, quiet = TRUE))",
      "}"
    ),
    r_profile
  )

  html_file <- prolfquapp:::render_quarto_se_report(
    se_file = se_file,
    output_dir = workdir,
    output_file = "helper-report.html",
    env = c(
      paste0("PROLFQUAPP_SOURCE_PATH=", package_path),
      paste0("R_PROFILE_USER=", r_profile)
    )
  )

  expect_true(file.exists(html_file))
  html <- paste(readLines(html_file, warn = FALSE), collapse = "\n")
  expect_match(html, "fgcz-banner", fixed = TRUE)
  expect_match(html, "Differential Expression Analysis", fixed = TRUE)
})
