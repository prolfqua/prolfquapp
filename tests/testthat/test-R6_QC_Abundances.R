test_that("render_sample_size_QC reads project_spec from self, not a free global GRP2", {
  # Static check: the method must not reference a free variable named `GRP2`.
  # Before the fix, `project_conf = GRP2$project_spec` referenced a global
  # GRP2, so the report rendered with a stray global's metadata (or errored).
  fn <- prolfquapp::QC_generator$public_methods$render_sample_size_QC
  free_vars <- codetools::findGlobals(fn, merge = FALSE)$variables

  expect_false("GRP2" %in% free_vars)
})
