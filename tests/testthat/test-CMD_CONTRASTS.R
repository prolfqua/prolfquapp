script <- system.file("application/CMD_CONTRASTS.R", package = "prolfquapp")
csv_s1 <- system.file(
  "application/contrasts/scenario1_single_factor.csv",
  package = "prolfquapp"
)
csv_s2 <- system.file(
  "application/contrasts/scenario2_two_factor.csv",
  package = "prolfquapp"
)
rscript <- file.path(R.home("bin"), "Rscript")

test_that("scenario 1 - single factor: writes output with ContrastName/Contrast columns", {
  skip_if(nchar(script) == 0, "CMD_CONTRASTS.R not installed")
  skip_if(nchar(csv_s1) == 0, "scenario1 CSV not installed")

  out <- file.path(tempdir(), "s1_out.csv")
  status <- system2(
    rscript,
    c(script, csv_s1, "--control", "WT", "-o", out),
    stdout = TRUE,
    stderr = TRUE
  )
  expect_equal(attr(status, "status"), NULL) # exit 0 → no status attribute
  expect_true(file.exists(out))

  result <- readr::read_csv(out, show_col_types = FALSE)
  expect_true("CONTROL" %in% colnames(result))
  expect_true(all(result$CONTROL %in% c("C", "T")))
  expect_true(any(result$CONTROL == "C"))
  expect_true(any(result$CONTROL == "T"))
  # WT rows should be C, others T
  expect_true(all(result$CONTROL[result$group == "WT"] == "C"))
  expect_true(all(result$CONTROL[result$group != "WT"] == "T"))
})

test_that("scenario 2 - two factor: writes output with ContrastName/Contrast columns", {
  skip_if(nchar(script) == 0, "CMD_CONTRASTS.R not installed")
  skip_if(nchar(csv_s2) == 0, "scenario2 CSV not installed")

  out <- file.path(tempdir(), "s2_out.csv")
  status <- system2(
    rscript,
    c(script, csv_s2, "--f1", "treatment", "--f2", "time", "-o", out),
    stdout = TRUE,
    stderr = TRUE
  )
  expect_equal(attr(status, "status"), NULL)

  expect_true(file.exists(out))

  result <- readr::read_csv(out, show_col_types = FALSE)
  expect_true("ContrastName" %in% colnames(result))
  expect_true("Contrast" %in% colnames(result))

  ct <- dplyr::distinct(result[
    !is.na(result$ContrastName),
    c("ContrastName", "Contrast")
  ])
  expect_gt(nrow(ct), 0)
  # two-factor: expect main effect, level-specific, and interaction contrasts
  expect_true(any(grepl("interaction", ct$ContrastName)))
})
