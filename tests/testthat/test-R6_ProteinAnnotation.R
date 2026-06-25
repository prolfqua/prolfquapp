test_that("ProteinAnnotation$clean guards and filters map contaminants->CON, decoys->REV", {
  res <- suppressWarnings(suppressMessages(
    prolfquapp::sim_data_protAnnot(Nprot = 100, PROTEIN = TRUE)
  ))
  pa <- res$pannot
  n <- nrow(pa$row_annot)
  testthat::skip_if(n < 20)

  # known decoy (REV) / contaminant (CON) layout, non-overlapping
  pa$row_annot$REV <- c(rep(TRUE, 10), rep(FALSE, n - 10))
  pa$row_annot$CON <- c(rep(FALSE, 10), rep(TRUE, 5), rep(FALSE, n - 15))

  # clean() must agree with nr_clean() for every flag combination
  expect_equal(nrow(pa$clean()), pa$nr_clean())
  expect_equal(nrow(pa$clean(contaminants = TRUE, decoys = FALSE)),
    pa$nr_clean(contaminants = TRUE, decoys = FALSE))
  expect_equal(nrow(pa$clean(contaminants = FALSE, decoys = TRUE)),
    pa$nr_clean(contaminants = FALSE, decoys = TRUE))

  # contaminants filtering removes CON only; decoys filtering removes REV only
  expect_equal(nrow(pa$clean()), n - 15)
  expect_equal(nrow(pa$clean(contaminants = TRUE, decoys = FALSE)), n - 5)
  expect_equal(nrow(pa$clean(contaminants = FALSE, decoys = TRUE)), n - 10)

  # the guard must check the column the requested filter needs, not its partner:
  # contaminants filtering needs CON present (REV may be absent)
  pa_no_rev <- res$pannot
  pa_no_rev$row_annot$CON <- c(rep(TRUE, 5), rep(FALSE, n - 5))
  pa_no_rev$row_annot$REV <- NULL
  expect_error(pa_no_rev$clean(contaminants = TRUE, decoys = FALSE), NA)

  # decoys filtering needs REV present (CON may be absent)
  res2 <- suppressWarnings(suppressMessages(
    prolfquapp::sim_data_protAnnot(Nprot = 100, PROTEIN = TRUE)
  ))
  pa_no_con <- res2$pannot
  m <- nrow(pa_no_con$row_annot)
  pa_no_con$row_annot$REV <- c(rep(TRUE, 5), rep(FALSE, m - 5))
  pa_no_con$row_annot$CON <- NULL
  expect_error(pa_no_con$clean(contaminants = FALSE, decoys = TRUE), NA)
})
