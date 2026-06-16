test_that("run_dea supports the prolfquasaint-backed saint model", {
  skip_on_cran()
  skip_if_not_installed("prolfquasaint")
  dataset <- system.file(
    "application/sim_test/dataset_sim.csv",
    package = "prolfquapp"
  )
  skip_if(nchar(dataset) == 0, "sim_test fixture not installed")

  config <- prolfquapp::make_DEA_config_R6(
    WORKUNITID = "TEST_SAINT",
    Normalization = "none",
    model = "saint"
  )

  result <- suppressWarnings(prolfquapp::run_dea(
    indir = tempdir(),
    dataset = dataset,
    software = "prolfquapp.SIM",
    config = config
  ))

  expect_equal(result$deanalyse$default_model, "saint")
  # SAINT is reached through the prolfqua facade registry; the stored
  # object is a ContrastsSAINTFacade wrapping a ContrastsSAINTexpress.
  expect_match(
    paste(class(result$deanalyse$contrast_results$saint), collapse = " "),
    "ContrastsSAINTFacade"
  )
  expect_s3_class(
    result$deanalyse$contrast_results$saint,
    "ContrastsInterface"
  )
  saint_extras <- result$deanalyse$contrast_results$saint$extra_artifacts()
  expect_true(all(
    c("saint_inter", "saint_prey", "saint_bait", "saint_list") %in%
      names(saint_extras)
  ))
  expect_gt(nrow(saint_extras$saint_list), 0)
  expect_true(all(
    c("BFDR", "SaintScore", "log2_EFCs") %in%
      colnames(result$deanalyse$annotated_contrasts)
  ))
  expect_setequal(
    intersect(
      c("FDR", "statistic", "diff", "contrast"),
      colnames(result$deanalyse$annotated_contrasts)
    ),
    character()
  )
})

test_that("SAINT report writer adds SAINT sheets and enrichment inputs", {
  skip_on_cran()
  skip_if_not_installed("prolfquasaint")
  dataset <- system.file(
    "application/sim_test/dataset_sim.csv",
    package = "prolfquapp"
  )
  skip_if(nchar(dataset) == 0, "sim_test fixture not installed")

  workdir <- file.path(tempdir(), "saint_report_writer")
  dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

  config <- prolfquapp::make_DEA_config_R6(
    PATH = workdir,
    WORKUNITID = "TEST_SAINT_REPORT",
    Normalization = "none",
    FDR_threshold = 1,
    diff_threshold = 0,
    model = "saint"
  )

  result <- suppressWarnings(prolfquapp::run_dea(
    indir = tempdir(),
    dataset = dataset,
    software = "prolfquapp.SIM",
    config = config
  ))

  reporter <- prolfquapp::DEAReportGenerator$new(result$deanalyse, config)
  files <- reporter$write_DEA(ORA = TRUE, GSEA = TRUE)

  expect_gt(length(files$ora_files), 0)
  expect_gt(length(files$gsea_files), 0)
  expect_match(names(files$ora_files)[1], "^ORA_Bait_")
  expect_match(names(files$gsea_files)[1], "^Bait_")
  expect_equal(all(file.exists(unlist(files$ora_files))), TRUE)
  expect_equal(all(file.exists(unlist(files$gsea_files))), TRUE)
  sheets <- readxl::excel_sheets(files$xlsx_file)
  diff_exp <- readxl::read_xlsx(files$xlsx_file, sheet = "diff_exp_analysis")
  expect_setequal(
    intersect(c("Bait", "BFDR", "SaintScore", "log2_EFCs"), colnames(diff_exp)),
    c("Bait", "BFDR", "SaintScore", "log2_EFCs")
  )
  expect_setequal(
    intersect(c("contrast", "FDR", "statistic", "diff"), colnames(diff_exp)),
    character()
  )
  expect_setequal(
    intersect(
      c("saint_inter", "saint_prey", "saint_bait", "saint_list"),
      sheets
    ),
    c("saint_inter", "saint_prey", "saint_bait", "saint_list")
  )
})

test_that("regular DEA report writer keeps ORA and GSEA exports", {
  skip_on_cran()

  workdir <- file.path(tempdir(), "regular_report_writer")
  dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

  dea <- prolfquapp::example_deanalyse(Nprot = 20)
  dea$prolfq_app_config$path <- workdir
  dea$prolfq_app_config$project_spec$workunit_Id <- "TEST_REGULAR_REPORT"
  dea$prolfq_app_config$set_zipdir_name()
  dea$FDR_threshold <- 1
  dea$diff_threshold <- 0

  reporter <- prolfquapp::DEAReportGenerator$new(
    dea,
    dea$prolfq_app_config
  )
  files <- reporter$write_DEA(ORA = TRUE, GSEA = TRUE)

  expect_gt(length(files$ora_files), 0)
  expect_gt(length(files$gsea_files), 0)
  expect_match(names(files$ora_files)[1], "^ORA_")
  expect_match(names(files$gsea_files)[1], "^GSEA_")
  expect_equal(all(file.exists(unlist(files$ora_files))), TRUE)
  expect_equal(all(file.exists(unlist(files$gsea_files))), TRUE)
})
