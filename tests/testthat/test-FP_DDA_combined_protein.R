test_that("FP_DDA_combined_protein", {
  tmp <- system.file("application/FragPipeDDA/TESTDATA.zip",package = "prolfquapp")
  tmpdir <- tempdir()
  unlink(paste0(tmpdir,"/*"), recursive = TRUE)
  file.copy(tmp, tmpdir)
  file <- dir(tmpdir,full.names = TRUE)
  unzip(file, exdir = tmpdir)

  xd <- prolfquapp::copy_DEA_FragPipe_DDA(workdir = file.path(tmpdir,"TESTDATA"))
  curdir <- getwd()
  setwd(file.path(tmpdir,"TESTDATA"))
  source(xd[4])
  tmp <- dir(file.path(tmpdir,"TESTDATA","C28350WU277930","DEA_"))
  testthat::expect_true(any(grepl(".rnk",tmp)))
  testthat::expect_true(any(grepl(".html",tmp)))
  testthat::expect_true(any(grepl(".xlsx",tmp)))

  testthat::expect_equal(length(tmp),7)
  setwd(curdir)

  })
