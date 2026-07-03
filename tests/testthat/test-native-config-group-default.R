# Regression: native configs (e.g. from the A414 app-runner) omit `group`.
# list_to_R6_app_config must NOT clobber the "G_" default with NULL, and
# read_annotation must tolerate a NULL/empty prefix -- otherwise set_grouping_var
# crashes with "attempt to select less than one element in OneIndex".

test_that("native config without `group` keeps the G_ default (not NULL)", {
  dd <- list(
    project_spec = list(project_Id = 1, order_Id = 1, workunit_Id = 1),
    processing_options = list(transform = "vsn", aggregate = "medpolish"),
    prefix = "DEA"
  )
  cfg <- prolfquapp::list_to_R6_app_config(dd)
  expect_false(is.null(cfg$group))
  expect_equal(cfg$group, "G_")
})

test_that("read_annotation tolerates a NULL / empty group prefix (falls back to G_)", {
  annot <- data.frame(
    "Relative Path" = c("a.raw", "b.raw", "c.raw", "d.raw"),
    "Name" = c("s1", "s2", "s3", "s4"),
    "Grouping Var" = c("a", "a", "b", "b"),
    "CONTROL" = c("T", "T", "C", "C"),
    check.names = FALSE
  )

  res_null <- prolfquapp::read_annotation(annot, prefix = NULL)
  expect_true("G_" %in% names(res_null$atable$factors))
  expect_equal(length(res_null$contrasts), 1L)

  res_empty <- prolfquapp::read_annotation(annot, prefix = "")
  expect_true("G_" %in% names(res_empty$atable$factors))
})
