test_that("run_make_yaml returns config list with expected fields", {
  cfg <- prolfquapp::run_make_yaml(
    project = "P_TEST",
    order = "O_TEST",
    workunit = "WU_TEST",
    norm = "robscale",
    model = "limma"
  )
  expect_type(cfg, "list")
  expect_equal(cfg$processing_options$transform, "robscale")
  expect_equal(cfg$processing_options$model, "limma")
  expect_equal(cfg$project_spec$workunit_Id, "WU_TEST")
  expect_equal(cfg$project_spec$project_Id, "P_TEST")
  expect_equal(cfg$project_spec$order_Id, "O_TEST")
})

test_that("run_make_yaml defaults produce valid config", {
  cfg <- prolfquapp::run_make_yaml()
  expect_type(cfg, "list")
  expect_equal(cfg$processing_options$transform, "vsn")
  expect_equal(cfg$processing_options$model, "lm_missing")
})

test_that("run_make_yaml reorders fields (internal sections at bottom)", {
  cfg <- prolfquapp::run_make_yaml()
  nms <- names(cfg)
  bottom_fields <- c("ext_reader", "group")
  present <- intersect(bottom_fields, nms)
  if (length(present) > 0) {
    positions <- match(present, nms)
    # These fields should be at the end
    expect_true(all(positions > length(nms) - length(present)))
  }
})

test_that("run_make_yaml round-trips through get_config", {
  cfg <- prolfquapp::run_make_yaml(workunit = "WU_RT")
  tmp <- tempfile(fileext = ".yaml")
  on.exit(unlink(tmp), add = TRUE)
  yaml::write_yaml(cfg, tmp)

  GRP2 <- prolfquapp::get_config(tmp)
  expect_true("ProlfquAppConfig" %in% class(GRP2))
  expect_equal(GRP2$project_spec$workunit_Id, "WU_RT")
})
