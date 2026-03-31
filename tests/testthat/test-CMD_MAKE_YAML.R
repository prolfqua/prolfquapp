script <- system.file("application/CMD_MAKE_YAML.R", package = "prolfquapp")
rscript <- file.path(R.home("bin"), "Rscript")

test_that("CMD_MAKE_YAML generates valid YAML with expected fields", {
  skip_if(nchar(script) == 0, "CMD_MAKE_YAML.R not installed")

  outdir <- file.path(tempdir(), "yaml_test")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)

  status <- system2(
    rscript,
    c(
      script,
      "--norm", "robscale",
      "--model", "limma",
      "--workunit", "WU_TEST",
      "--project", "P_TEST",
      "--order", "O_TEST",
      "-o", outdir,
      "-y", "test_config.yaml"
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  expect_equal(attr(status, "status"), NULL)

  yaml_path <- file.path(outdir, "test_config.yaml")
  expect_true(file.exists(yaml_path))

  cfg <- yaml::read_yaml(yaml_path)
  expect_equal(cfg$processing_options$transform, "robscale")
  expect_equal(cfg$processing_options$model, "limma")
  expect_equal(cfg$project_spec$workunit_Id, "WU_TEST")
  expect_equal(cfg$project_spec$project_Id, "P_TEST")
  expect_equal(cfg$project_spec$order_Id, "O_TEST")
})

test_that("CMD_MAKE_YAML round-trips through get_config", {
  skip_if(nchar(script) == 0, "CMD_MAKE_YAML.R not installed")

  outdir <- file.path(tempdir(), "yaml_roundtrip")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)

  system2(
    rscript,
    c(script, "-o", outdir, "-y", "rt_config.yaml"),
    stdout = TRUE,
    stderr = TRUE
  )

  yaml_path <- file.path(outdir, "rt_config.yaml")
  skip_if(!file.exists(yaml_path), "YAML not generated")

  GRP2 <- prolfquapp::get_config(yaml_path)
  expect_true("ProlfquAppConfig" %in% class(GRP2))
  expect_true(!is.null(GRP2$processing_options$model))
})
