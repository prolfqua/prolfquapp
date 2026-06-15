test_that("copy_DEA_R6_Files copies R6 report templates and bibliography", {
  workdir <- tempfile("dea-r6-files-")
  dir.create(workdir)
  on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

  files <- c("Grp2Analysis_V2_R6.Rmd", "DiffExpQC_R6.Rmd")

  # Assert on what the function produces, not on whether the vignettes happen
  # to be built into inst/doc. copy_DEA_R6_Files resolves templates from the
  # installed doc/ and falls back to the source tree (doc/ or vignettes/), so
  # this test must not depend on `make build-vignettes` having run first.
  copied <- prolfquapp::copy_DEA_R6_Files(workdir = workdir)

  expect_setequal(basename(copied), files)
  expect_true(all(file.exists(file.path(workdir, files))))
  expect_true(file.exists(file.path(workdir, "bibliography.bib")))
})
