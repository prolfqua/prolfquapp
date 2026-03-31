script <- system.file("application/CMD_DEA_V2.R", package = "prolfquapp")
dataset <- system.file(
  "application/sim_test/dataset_sim.csv",
  package = "prolfquapp"
)
rscript <- file.path(R.home("bin"), "Rscript")

test_that("CMD_DEA_V2 runs full pipeline with SIM preprocessor", {
  skip_on_cran()
  skip_if(nchar(script) == 0, "CMD_DEA_V2.R not installed")
  skip_if(nchar(dataset) == 0, "sim_test fixture not installed")

  workdir <- file.path(tempdir(), "dea_v2_test")
  dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

  # Generate config.yaml
  yaml_script <- system.file("application/CMD_MAKE_YAML.R", package = "prolfquapp")
  skip_if(nchar(yaml_script) == 0, "CMD_MAKE_YAML.R not installed")

  system2(
    rscript,
    c(yaml_script, "--norm", "robscale", "-o", workdir, "-y", "config.yaml"),
    stdout = TRUE,
    stderr = TRUE
  )
  yaml_path <- file.path(workdir, "config.yaml")
  skip_if(!file.exists(yaml_path), "config.yaml not generated")

  # Run DEA pipeline
  status <- system2(
    rscript,
    c(
      script,
      "-i", workdir,
      "-d", dataset,
      "-y", yaml_path,
      "-s", "prolfquapp.SIM",
      "-w", "TEST_WU",
      "-o", workdir
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  exit_code <- attr(status, "status")
  if (!is.null(exit_code) && exit_code != 0) {
    message("CMD_DEA_V2 stderr:\n", paste(status, collapse = "\n"))
  }
  expect_true(is.null(exit_code) || exit_code == 0)

  # Check output directory was created
  result_dirs <- list.dirs(workdir, recursive = FALSE)
  result_dir <- grep("^DEA_", basename(result_dirs), value = TRUE)
  expect_true(length(result_dir) >= 1)
})
