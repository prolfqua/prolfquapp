script <- system.file("application/CMD_QUANT_QC.R", package = "prolfquapp")
dataset <- system.file(
  "application/sim_test/dataset_sim.csv",
  package = "prolfquapp"
)
rscript <- file.path(R.home("bin"), "Rscript")

test_that("CMD_QUANT_QC runs with SIM preprocessor", {
  skip_on_cran()
  skip_if(nchar(script) == 0, "CMD_QUANT_QC.R not installed")
  skip_if(nchar(dataset) == 0, "sim_test fixture not installed")

  workdir <- file.path(tempdir(), "qc_test")
  dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

  # Generate config.yaml
  yaml_script <- system.file("application/CMD_MAKE_YAML.R", package = "prolfquapp")
  skip_if(nchar(yaml_script) == 0, "CMD_MAKE_YAML.R not installed")

  system2(
    rscript,
    c(yaml_script, "-o", workdir, "-y", "config.yaml"),
    stdout = TRUE,
    stderr = TRUE
  )
  yaml_path <- file.path(workdir, "config.yaml")
  skip_if(!file.exists(yaml_path), "config.yaml not generated")

  # Run QC pipeline
  status <- system2(
    rscript,
    c(
      script,
      "-i", workdir,
      "-d", dataset,
      "-s", "SIM",
      "-w", "TEST_QC",
      "-o", workdir,
      "-y", yaml_path
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  exit_code <- attr(status, "status")
  if (!is.null(exit_code) && exit_code != 0) {
    message("CMD_QUANT_QC stderr:\n", paste(status, collapse = "\n"))
  }
  expect_true(is.null(exit_code) || exit_code == 0)

  # Check QC output directory was created with files
  qc_dirs <- list.dirs(workdir, recursive = FALSE)
  expect_true(length(qc_dirs) >= 1)
})
