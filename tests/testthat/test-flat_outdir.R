dataset <- system.file(
  "application/sim_test/dataset_sim.csv",
  package = "prolfquapp"
)

test_that("get_zipdir respects flat_outdir flag", {
  cfg <- prolfquapp::make_DEA_config_R6(
    PROJECTID = "123",
    ORDERID = "456",
    WORKUNITID = "789",
    prefix = "QC"
  )
  cfg$path <- "myout"

  # default (dynamic): dated subdir under path
  expect_false(cfg$flat_outdir)
  expect_equal(cfg$get_zipdir(), file.path("myout", cfg$zipdir_name))

  # flat (static): outputs go directly into path
  cfg$flat_outdir <- TRUE
  expect_equal(cfg$get_zipdir(), "myout")
  expect_equal(cfg$get_result_dir(), file.path("myout", "Results_WU_789"))
  expect_equal(cfg$get_input_dir(), file.path("myout", "Inputs_WU_789"))
})

test_that("flat_outdir round-trips through the config list", {
  cfg <- prolfquapp::make_DEA_config_R6(WORKUNITID = "WU1")
  cfg$flat_outdir <- TRUE

  lst <- cfg$as_list()
  expect_true(lst$flat_outdir)

  restored <- prolfquapp::list_to_R6_app_config(lst)
  expect_true(restored$flat_outdir)
  expect_equal(restored$get_zipdir(), restored$path)
})

test_that("flat_outdir defaults to FALSE through the config list", {
  cfg <- prolfquapp::make_DEA_config_R6(WORKUNITID = "WU1")
  restored <- prolfquapp::list_to_R6_app_config(cfg$as_list())
  expect_false(restored$flat_outdir)
  expect_equal(restored$get_zipdir(), file.path(restored$path, restored$zipdir_name))
})

test_that("sync_opt_config applies flat_outdir override", {
  cfg <- prolfquapp::make_DEA_config_R6(WORKUNITID = "WU1")
  expect_false(cfg$flat_outdir)

  res <- prolfquapp::sync_opt_config(list(flat_outdir = TRUE), cfg)
  expect_true(res$config$flat_outdir)
  expect_equal(res$config$get_zipdir(), res$config$path)
})

test_that("run_qc_preprocess honors flat_outdir", {
  skip_if(nchar(dataset) == 0, "sim_test fixture not installed")

  out <- tempfile("qc_flat_")
  dir.create(out)
  on.exit(unlink(out, recursive = TRUE), add = TRUE)

  result <- prolfquapp::run_qc_preprocess(
    indir = tempdir(),
    dataset = dataset,
    software = "SIM",
    outdir = out,
    workunit = "TEST_QC",
    flat_outdir = TRUE
  )
  expect_true(result$config$flat_outdir)
  expect_equal(result$config$get_zipdir(), out)
})
