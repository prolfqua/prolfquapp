dataset <- system.file(
  "application/sim_test/dataset_sim.csv",
  package = "prolfquapp"
)

test_that("run_dea returns deanalyse and supporting data with SIM", {
  skip_on_cran()
  skip_if(nchar(dataset) == 0, "sim_test fixture not installed")

  GRP2 <- prolfquapp::make_DEA_config_R6(
    WORKUNITID = "TEST_DEA",
    Normalization = "robscale"
  )

  result <- prolfquapp::run_dea(
    indir = tempdir(),
    dataset = dataset,
    software = "prolfquapp.SIM",
    config = GRP2
  )

  expect_type(result, "list")
  expect_true("DEAnalyse" %in% class(result$deanalyse))
  expect_true("LFQData" %in% class(result$xd$lfqdata))
  expect_true(
    "ProteinAnnotation" %in% class(result$xd$protein_annotation)
  )
  expect_true(!is.null(result$annotation$contrasts))
  expect_true(length(result$annotation$contrasts) > 0)
  expect_equal(result$software, "prolfquapp.SIM")
  expect_equal(result$requested_software, "prolfquapp.SIM")
  expect_equal(GRP2$software, "prolfquapp.SIM")
})

test_that("run_dea maps legacy prolfqua model to current facade names", {
  skip_on_cran()
  skip_if(nchar(dataset) == 0, "sim_test fixture not installed")

  GRP2 <- prolfquapp::make_DEA_config_R6(
    WORKUNITID = "TEST_DEA_LEGACY_MODEL",
    Normalization = "robscale"
  )
  GRP2$processing_options$model <- "prolfqua"
  GRP2$processing_options$model_missing <- FALSE

  result <- prolfquapp::run_dea(
    indir = tempdir(),
    dataset = dataset,
    software = "prolfquapp.SIM",
    config = GRP2
  )

  expect_equal(result$deanalyse$default_model, "lm")

  GRP2$processing_options$model_missing <- TRUE
  result <- prolfquapp::run_dea(
    indir = tempdir(),
    dataset = dataset,
    software = "prolfquapp.SIM",
    config = GRP2
  )

  expect_equal(result$deanalyse$default_model, "lm_impute")
})

test_that(".resolve_nested_reader switches a nested facade to the peptide reader", {
  available <- c("prolfquapp.DIANN", "prolfquapp.DIANN_PEPTIDE")
  out <- prolfquapp:::.resolve_nested_reader(
    "prolfquapp.DIANN",
    is_nested = TRUE,
    available = available,
    facade = "firth_nested"
  )
  expect_equal(out, "prolfquapp.DIANN_PEPTIDE")
})

test_that(".resolve_nested_reader leaves non-nested and peptide readers unchanged", {
  available <- c("prolfquapp.DIANN", "prolfquapp.DIANN_PEPTIDE")
  # non-nested facade: unchanged
  expect_equal(
    prolfquapp:::.resolve_nested_reader(
      "prolfquapp.DIANN",
      is_nested = FALSE,
      available = available
    ),
    "prolfquapp.DIANN"
  )
  # already a peptide reader: unchanged
  expect_equal(
    prolfquapp:::.resolve_nested_reader(
      "prolfquapp.DIANN_PEPTIDE",
      is_nested = TRUE,
      available = available
    ),
    "prolfquapp.DIANN_PEPTIDE"
  )
})

test_that(".resolve_nested_reader errors when no peptide counterpart is registered", {
  expect_error(
    prolfquapp:::.resolve_nested_reader(
      "prolfquapp.SIM",
      is_nested = TRUE,
      available = c("prolfquapp.SIM"),
      facade = "firth_nested"
    ),
    "no peptide-level counterpart"
  )
})

test_that("run_dea errors when no peptide-level reader exists for a nested facade", {
  skip_if(nchar(dataset) == 0, "sim_test fixture not installed")

  GRP2 <- prolfquapp::make_DEA_config_R6(WORKUNITID = "TEST_DEA_NESTED_NOPEP")
  GRP2$processing_options$model <- "firth_nested"
  GRP2$processing_options$model_missing <- FALSE

  expect_error(
    prolfquapp::run_dea(
      indir = tempdir(),
      dataset = dataset,
      software = "prolfquapp.SIM",
      config = GRP2
    ),
    "no peptide-level counterpart"
  )
})

test_that("run_dea errors on missing annotation file", {
  GRP2 <- prolfquapp::make_DEA_config_R6()
  expect_error(
    prolfquapp::run_dea(
      indir = tempdir(),
      dataset = "/no/such/file.csv",
      software = "prolfquapp.SIM",
      config = GRP2
    ),
    "Annotation file not found"
  )
})

test_that("run_dea errors on unknown software", {
  skip_if(nchar(dataset) == 0, "sim_test fixture not installed")

  GRP2 <- prolfquapp::make_DEA_config_R6()
  expect_error(
    prolfquapp::run_dea(
      indir = tempdir(),
      dataset = dataset,
      software = "NONEXISTENT",
      config = GRP2
    ),
    "not found"
  )
})
