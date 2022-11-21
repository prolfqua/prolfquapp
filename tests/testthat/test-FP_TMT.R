test_that("FP_DDA_peptide.tsc", {

  tmp <- system.file("application/FragPipeTMT/TESTDATA.zip", package = "prolfquapp")
  tmpdir <- tempdir()
  testthat::expect_true(file.exists(tmp))

  unlink(paste0(tmpdir,"/*"), recursive = TRUE)
  testthat::expect_equal(0,length(dir(tmpdir)))

  file.copy(tmp, tmpdir)
  file <- dir(tmpdir,full.names = TRUE)
  unzip(file, exdir = tmpdir)
  xd <- prolfquapp::copy_DEA_FragPipe_TMT(workdir = file.path(tmpdir,"TESTDATA"))

  curdir <- getwd()
  setwd(file.path(tmpdir,"TESTDATA"))
  source(xd[4])

  #tmp <- dir(file.path(tmpdir,"TESTDATA","C29460WU282164","DE_Groups_vs_Controls"))
  testthat::expect_true(any(grepl(".rnk",tmp)))
  testthat::expect_true(any(grepl(".html",tmp)))
  testthat::expect_true(any(grepl(".xlsx",tmp)))
  testthat::expect_equal(length(tmp),9)
  unlink(paste0(tmpdir,"/*"), recursive = TRUE)
  setwd(curdir)
})
