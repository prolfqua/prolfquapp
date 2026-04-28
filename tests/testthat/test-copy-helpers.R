test_that("copy_DEA_R6_Files copies installed R6 report templates", {
  workdir <- tempfile("dea-r6-files-")
  dir.create(workdir)
  on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

  files <- c("Grp2Analysis_V2_R6.Rmd", "DiffExpQC_R6.Rmd")
  src <- system.file("doc", files, package = "prolfquapp")

  expect_true(all(nzchar(src)))
  expect_true(all(file.exists(src)))

  copied <- prolfquapp::copy_DEA_R6_Files(workdir = workdir)

  expect_setequal(basename(copied), files)
  expect_true(all(file.exists(file.path(workdir, files))))
  expect_true(file.exists(file.path(workdir, "bibliography.bib")))
})
