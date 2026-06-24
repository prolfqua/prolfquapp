test_that("set_list_to_R6 sets scalar and nested fields on an R6 config", {
  cfg <- prolfquapp::ProlfquAppConfig$new(
    prolfquapp::ProcessingOptions$new(),
    prolfquapp::ProjectSpec$new(),
    prolfquapp::ExternalReader$new()
  )

  config_list <- list(
    software = "DIANN",
    prefix = "DEA",
    processing_options = list(
      transform = "robscale",
      FDR_threshold = 0.05
    ),
    project_spec = list(
      project_Id = "1234",
      workunit_Id = "WU99"
    )
  )

  prolfquapp::set_list_to_R6(config_list, cfg)

  expect_equal(cfg$software, "DIANN")
  expect_equal(cfg$prefix, "DEA")
  expect_equal(cfg$processing_options$transform, "robscale")
  expect_equal(cfg$processing_options$FDR_threshold, 0.05)
  expect_equal(cfg$project_spec$project_Id, "1234")
  expect_equal(cfg$project_spec$workunit_Id, "WU99")
})

test_that("set_list_to_R6 returns the (mutated) R6 object invisibly", {
  cfg <- prolfquapp::ProlfquAppConfig$new(
    prolfquapp::ProcessingOptions$new(),
    prolfquapp::ProjectSpec$new(),
    prolfquapp::ExternalReader$new()
  )

  out <- prolfquapp::set_list_to_R6(list(software = "FP_TMT"), cfg)

  expect_identical(out, cfg)
  expect_equal(cfg$software, "FP_TMT")
})
