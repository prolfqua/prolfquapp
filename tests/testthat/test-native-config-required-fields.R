complete_native_config <- function() {
  list(
    project_spec = list(project_Id = 1, order_Id = 1, workunit_Id = 1),
    processing_options = list(transform = "vsn", aggregate = "medpolish"),
    prefix = "DEA",
    group = "G_"
  )
}

test_that("list_to_R6_app_config accepts a complete native config", {
  cfg <- prolfquapp::list_to_R6_app_config(complete_native_config())
  expect_equal(cfg$group, "G_")
  expect_equal(cfg$prefix, "DEA")
})

test_that("list_to_R6_app_config errors clearly on a missing required field", {
  no_group <- complete_native_config()
  no_group$group <- NULL
  expect_error(prolfquapp::list_to_R6_app_config(no_group), "group")

  no_prefix <- complete_native_config()
  no_prefix$prefix <- NULL
  expect_error(prolfquapp::list_to_R6_app_config(no_prefix), "prefix")

  no_workunit <- complete_native_config()
  no_workunit$project_spec$workunit_Id <- NULL
  expect_error(prolfquapp::list_to_R6_app_config(no_workunit), "workunit_Id")
})
