dataset <- system.file(
  "application/sim_test/dataset_sim.csv",
  package = "prolfquapp"
)

test_that("run_qc_preprocess returns xd and config with SIM", {
  skip_if(nchar(dataset) == 0, "sim_test fixture not installed")

  result <- prolfquapp::run_qc_preprocess(
    indir = tempdir(),
    dataset = dataset,
    software = "SIM",
    outdir = tempdir(),
    workunit = "TEST_QC"
  )
  expect_type(result, "list")
  expect_true("xd" %in% names(result))
  expect_true("config" %in% names(result))
  expect_true("LFQData" %in% class(result$xd$lfqdata))
  expect_true(
    "ProteinAnnotation" %in% class(result$xd$protein_annotation)
  )
  expect_true(
    "ProlfquAppConfig" %in% class(result$config)
  )
  expect_equal(result$xd$lfqdata$get_config()$hierarchy_depth, 1)
})

test_that("run_qc_preprocess errors on missing annotation", {
  expect_error(
    prolfquapp::run_qc_preprocess(
      indir = tempdir(),
      dataset = "/no/such/file.csv",
      software = "SIM"
    ),
    "No annotation file found"
  )
})

test_that("run_qc_preprocess uses yaml config when provided", {
  skip_if(nchar(dataset) == 0, "sim_test fixture not installed")

  cfg <- prolfquapp::run_make_yaml(workunit = "WU_YAML")
  tmp_yaml <- tempfile(fileext = ".yaml")
  on.exit(unlink(tmp_yaml), add = TRUE)
  yaml::write_yaml(cfg, tmp_yaml)

  result <- prolfquapp::run_qc_preprocess(
    indir = tempdir(),
    dataset = dataset,
    software = "SIM",
    yaml_file = tmp_yaml
  )
  expect_equal(
    result$config$project_spec$workunit_Id, "WU_YAML"
  )
})
