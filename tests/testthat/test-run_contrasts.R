csv_s1 <- system.file(
  "application/contrasts/scenario1_single_factor.csv",
  package = "prolfquapp"
)
csv_s2 <- system.file(
  "application/contrasts/scenario2_two_factor.csv",
  package = "prolfquapp"
)

# --- run_contrasts_single -----------------------------------------------------

test_that("run_contrasts_single adds CONTROL column", {
  skip_if(nchar(csv_s1) == 0, "scenario1 CSV not installed")

  result <- prolfquapp::run_contrasts_single(csv_s1, control = "WT")
  expect_true("CONTROL" %in% colnames(result))
  expect_true(all(result$CONTROL[result$group == "WT"] == "C"))
  expect_true(all(result$CONTROL[result$group != "WT"] == "T"))
})

test_that("run_contrasts_single with explicit group column", {
  skip_if(nchar(csv_s1) == 0, "scenario1 CSV not installed")

  result <- prolfquapp::run_contrasts_single(
    csv_s1, control = "WT", group = "group"
  )
  expect_true("CONTROL" %in% colnames(result))
  expect_equal(sum(result$CONTROL == "C"), 3)
})

test_that("run_contrasts_single errors on invalid control level", {
  skip_if(nchar(csv_s1) == 0, "scenario1 CSV not installed")

  expect_error(
    prolfquapp::run_contrasts_single(csv_s1, control = "NONEXISTENT"),
    "not found"
  )
})

test_that("run_contrasts_single errors on missing file", {
  expect_error(
    prolfquapp::run_contrasts_single("/no/such/file.csv", control = "WT")
  )
})

# --- run_contrasts_twofactor --------------------------------------------------

test_that("run_contrasts_twofactor adds contrast columns", {
  skip_if(nchar(csv_s2) == 0, "scenario2 CSV not installed")

  result <- prolfquapp::run_contrasts_twofactor(
    csv_s2, f1 = "treatment", f2 = "time"
  )
  expect_true("ContrastName" %in% colnames(result))
  expect_true("Contrast" %in% colnames(result))

  ct <- dplyr::distinct(result[
    !is.na(result$ContrastName),
    c("ContrastName", "Contrast")
  ])
  expect_gt(nrow(ct), 0)
  expect_true(any(grepl("interaction", ct$ContrastName, ignore.case = TRUE)))
})

test_that("run_contrasts_twofactor without interactions", {
  skip_if(nchar(csv_s2) == 0, "scenario2 CSV not installed")

  result <- prolfquapp::run_contrasts_twofactor(
    csv_s2, f1 = "treatment", f2 = "time", interactions = FALSE
  )
  ct <- dplyr::distinct(result[
    !is.na(result$ContrastName),
    c("ContrastName", "Contrast")
  ])
  expect_false(any(grepl("interaction", ct$ContrastName, ignore.case = TRUE)))
})

test_that("run_contrasts_twofactor errors on missing columns", {
  skip_if(nchar(csv_s1) == 0, "scenario1 CSV not installed")

  expect_error(
    prolfquapp::run_contrasts_twofactor(csv_s1, f1 = "nope", f2 = "nope2"),
    "not found"
  )
})
