test_that("copy_DEA_FragPipe_DDA", {
  tmpdir = tempdir()
  dir(tmpdir, recursive = TRUE)
  unlink(paste0(tmpdir,"/*"), recursive = TRUE)
  tmp <- copy_DEA_FragPipe_DDA(workdir = tempdir())
  testthat::expect_equal(length(dir(tmpdir)), 4)
})

test_that("copy_DEA_MaxQuant", {
  tmpdir = tempdir()
  dir(tmpdir, recursive = TRUE)
  unlink(paste0(tmpdir,"/*"), recursive = TRUE)
  tmp <- length(dir(tmpdir, recursive = TRUE))
  tmp <- copy_DEA_MaxQuant(workdir = tempdir())
  testthat::expect_equal(length(dir(tmpdir)), 4)
})

test_that("copy_DEA_FragPipe_TMT", {
  tmpdir = tempdir()
  dir(tmpdir, recursive = TRUE)
  unlink(paste0(tmpdir,"/*"), recursive = TRUE)
  tmp <- length(dir(tmpdir, recursive = TRUE))
  tmp <- copy_DEA_MaxQuant(workdir = tempdir())
  testthat::expect_equal(length(dir(tmpdir)), 4)
})

test_that("copy_DEA_DIANN", {
  tmpdir = tempdir()
  dir(tmpdir, recursive = TRUE)
  unlink(paste0(tmpdir,"/*"), recursive = TRUE)
  tmp <- length(dir(tmpdir, recursive = TRUE))
  tmp <- copy_DEA_MaxQuant(workdir = tempdir())
  testthat::expect_equal(length(dir(tmpdir)), 4)
})
