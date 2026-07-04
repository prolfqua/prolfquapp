test_package_path <- function() {
  wd <- normalizePath(getwd(), mustWork = TRUE)
  current <- wd
  repeat {
    description <- file.path(current, "DESCRIPTION")
    if (file.exists(description)) {
      desc <- read.dcf(description)
      if (identical(unname(desc[1, "Package"]), "prolfquapp")) {
        return(current)
      }
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      break
    }
    current <- parent
  }
  package_path <- normalizePath(
    file.path(system.file(package = "prolfquapp"), ".."),
    mustWork = TRUE
  )
  if (!file.exists(file.path(package_path, "DESCRIPTION"))) {
    package_path <- normalizePath(
      system.file(package = "prolfquapp"),
      mustWork = TRUE
    )
  }
  package_path
}

test_template_dir <- function() {
  package_path <- test_package_path()
  source_template_dir <- file.path(package_path, "inst", "templates", "quarto")
  if (dir.exists(source_template_dir)) {
    return(normalizePath(source_template_dir, mustWork = TRUE))
  }
  system.file("templates/quarto", package = "prolfquapp")
}

test_that("se_report_lfqdata reconstructs LFQData objects from SummarizedExperiment", {
  skip_on_cran()

  dea <- prolfquapp::example_deanalyse(Nprot = 20)
  reporter <- prolfquapp::DEAReportGenerator$new(dea, dea$prolfq_app_config)
  se <- reporter$make_SummarizedExperiment()
  expect_contains(SummarizedExperiment::assayNames(se), "nr_children")

  report <- prolfquapp:::se_report_lfqdata(se)

  expect_s3_class(report$lfq_raw, "LFQData")
  expect_s3_class(report$lfq_transformed, "LFQData")
  expect_gt(nrow(report$lfq_raw$data_long()), 0)
  expect_gt(nrow(report$lfq_transformed$data_long()), 0)
  expect_gt(nrow(report$lfq_raw$data_wide(as.matrix = TRUE)$data), 0)
  expect_gt(nrow(report$lfq_transformed$data_wide(as.matrix = TRUE)$data), 0)
  expect_s3_class(report$lfq_raw$get_Plotter(), "LFQDataPlotter")
  expect_s3_class(report$lfq_raw$get_Summariser(), "LFQDataSummariser")
  expect_s3_class(report$lfq_transformed$get_Stats(), "LFQDataStats")
  expect_contains(colnames(report$contrast_table), "protein_Id")
  expect_contains(colnames(report$contrast_table), "contrast")
  expect_contains(colnames(report$feature_annotation), "nrPeptides")
  expect_gt(sum(report$feature_annotation$nrPeptides >= 2, na.rm = TRUE), 0)
  child_counts <- report$lfq_raw$get_Summariser()$hierarchy_counts_sample(
    value = "long",
    nr_children = 2
  )
  expect_gt(nrow(child_counts), 0)
})

test_that("se_report_lfqdata drops padded empty contrast rows", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      rawData = matrix(
        1:6,
        nrow = 3,
        dimnames = list(c("P1", "P2", "P3"), c("S1", "S2"))
      ),
      transformedData = matrix(
        1:6,
        nrow = 3,
        dimnames = list(c("P1", "P2", "P3"), c("S1", "S2"))
      )
    ),
    colData = S4Vectors::DataFrame(
      sampleName = c("S1", "S2"),
      group = c("A", "B")
    )
  )
  SummarizedExperiment::rowData(se)[["constrast_A"]] <- data.frame(
    protein_Id = c("P1", "P2", NA_character_),
    contrast = c("A", "A", NA_character_),
    diff = c(1, -1, NA_real_),
    FDR = c(0.01, 0.2, NA_real_),
    modelName = c("TableTest", "TableTest", NA_character_)
  )

  report <- prolfquapp:::se_report_lfqdata(se)

  expect_equal(nrow(report$contrast_table), 2)
  expect_equal(sum(is.na(report$contrast_table$contrast)), 0)
  expect_equal(sum(is.na(report$contrast_table$protein_Id)), 0)
})

test_that("se_report_lfqdata drops SAINT rows without bait estimates", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      rawData = matrix(
        1:6,
        nrow = 3,
        dimnames = list(c("P1", "P2", "P3"), c("S1", "S2"))
      ),
      transformedData = matrix(
        1:6,
        nrow = 3,
        dimnames = list(c("P1", "P2", "P3"), c("S1", "S2"))
      )
    ),
    colData = S4Vectors::DataFrame(
      sampleName = c("S1", "S2"),
      group = c("A", "B")
    ),
    metadata = list(default_model = "saint")
  )
  SummarizedExperiment::rowData(se)[["constrast_PPE4"]] <- data.frame(
    protein_Id = c("P1", "P2", "P3"),
    modelName = c("ContrastSaint", NA_character_, "ContrastSaint"),
    Bait = c("PPE4", NA_character_, "PPE4"),
    log2_EFCs = c(1, NA_real_, -1),
    SaintScore = c(0.9, NA_real_, 0.2),
    BFDR = c(0.01, NA_real_, 0.8)
  )

  report <- prolfquapp:::se_report_lfqdata(se)

  expect_equal(report$contrast_table$protein_Id, c("P1", "P3"))
  expect_equal(sum(is.na(report$contrast_table$Bait)), 0)
})

test_that("SE Quarto templates render with reconstructed LFQData", {
  skip_on_cran()
  skip_if(Sys.which("quarto") == "", "Quarto CLI not installed")
  skip_if_not_installed("devtools")
  skip_if_not_installed("DT")
  skip_if_not_installed("fgczquartotemplate")
  skip_if_not_installed("gridExtra")
  skip_if_not_installed("plotly")

  dea <- prolfquapp::example_deanalyse(Nprot = 12)
  reporter <- prolfquapp::DEAReportGenerator$new(dea, dea$prolfq_app_config)
  se <- reporter$make_SummarizedExperiment()

  template_dir <- test_template_dir()
  skip_if(!nzchar(template_dir), "Quarto templates not installed")

  workdir <- file.path(tempdir(), "quarto_se_report_test")
  unlink(workdir, recursive = TRUE)
  dir.create(workdir, recursive = TRUE)
  on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

  se_file <- file.path(workdir, "SummarizedExperiment.rds")
  saveRDS(se, se_file)

  report_files <- "Grp2Analysis_V2_SE_tabset.qmd"
  file.copy(file.path(template_dir, report_files), workdir)
  fgczquartotemplate::fgcz_copy_assets(workdir)

  package_path <- test_package_path()

  oldwd <- setwd(workdir)
  on.exit(setwd(oldwd), add = TRUE)

  for (report_file in report_files) {
    status <- system2(
      "quarto",
      c(
        "render",
        report_file,
        "-P",
        paste0("se_file:", normalizePath(se_file)),
        "-P",
        paste0("prolfquapp_source_path:", package_path)
      ),
      stdout = TRUE,
      stderr = TRUE
    )
    exit_code <- attr(status, "status")
    if (!is.null(exit_code) && exit_code != 0) {
      message("quarto render output:\n", paste(status, collapse = "\n"))
    }
    expect_null(exit_code)
    html_file <- sub("[.]qmd$", ".html", report_file)
    expect_equal(file.exists(html_file), TRUE)
    html <- paste(readLines(html_file, warn = FALSE), collapse = "\n")
    expect_match(html, "fgcz-banner", fixed = TRUE)
    expect_match(html, "panel-tabset", fixed = TRUE)
    expect_match(html, "Feature Detection", fixed = TRUE)
    expect_match(html, "Differential Abundance", fixed = TRUE)
    expect_match(html, "Result Table", fixed = TRUE)
    expect_match(html, "at least two peptides in the experiment", fixed = TRUE)
    expect_equal(grepl("Protein Identification", html, fixed = TRUE), FALSE)
  }
})

test_that("internal SE Quarto report helper renders HTML", {
  skip_on_cran()
  skip_if(Sys.which("quarto") == "", "Quarto CLI not installed")
  skip_if_not_installed("devtools")
  skip_if_not_installed("DT")
  skip_if_not_installed("fgczquartotemplate")
  skip_if_not_installed("gridExtra")
  skip_if_not_installed("plotly")

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

  html_file <- prolfquapp:::render_quarto_se_report(
    se_file = se_file,
    output_dir = workdir,
    output_file = "helper-report.html",
    template_dir = test_template_dir(),
    prolfquapp_source_path = package_path,
    buttons = TRUE
  )

  expect_equal(file.exists(html_file), TRUE)
  html <- paste(readLines(html_file, warn = FALSE), collapse = "\n")
  expect_match(html, "fgcz-banner", fixed = TRUE)
  expect_match(html, "Differential Abundance Analysis", fixed = TRUE)
  expect_match(html, "Feature Detection", fixed = TRUE)
  expect_match(html, "Result Table", fixed = TRUE)
  expect_match(html, "fgcz-pf-toolbar", fixed = TRUE)
})
