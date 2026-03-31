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
