test_that("ProteinDataPrep aggregate returns unchanged site-level data", {
  sim <- prolfqua::sim_lfq_data_peptide_config(n_proteins = 10)
  lfq <- prolfqua::LFQData$new(sim$data, sim$config)
  lfq$set_config_value("hierarchy_depth", length(lfq$hierarchy_keys()))

  row_annot <- data.frame(
    protein_Id = unique(lfq$data_long()$protein_Id),
    description = unique(lfq$data_long()$protein_Id),
    nr_peptides = 1
  )
  pannot <- prolfquapp::ProteinAnnotation$new(
    lfq,
    row_annot = row_annot,
    description = "description",
    exp_nr_children = "nr_peptides"
  )
  config <- prolfquapp::make_DEA_config_R6(aggregation = "medpolish")

  data_prep <- prolfquapp::ProteinDataPrep$new(lfq, pannot, config)

  expect_warning(
    result <- data_prep$aggregate(),
    "nothing to aggregate from, returning unchanged data."
  )

  expect_true("LFQData" %in% class(result))
  expect_identical(data_prep$lfq_data, lfq)
  expect_equal(data_prep$lfq_data$data_long(), lfq$data_long())
  expect_equal(
    data_prep$lfq_data$get_config()$hierarchy_depth,
    length(data_prep$lfq_data$hierarchy_keys())
  )
  expect_null(data_prep$aggregator)
})

test_that("ProteinDataPrep aggregates with every aggregation make_DEA_config_R6 offers", {
  sim <- prolfqua::sim_lfq_data_peptide_config(n_proteins = 10)
  lfq <- prolfqua::LFQData$new(sim$data, sim$config)
  row_annot <- data.frame(
    protein_Id = unique(lfq$data_long()$protein_Id),
    description = unique(lfq$data_long()$protein_Id),
    nr_peptides = 1
  )
  pannot <- prolfquapp::ProteinAnnotation$new(
    lfq,
    row_annot = row_annot,
    description = "description",
    exp_nr_children = "nr_peptides"
  )

  for (agg in eval(formals(prolfquapp::make_DEA_config_R6)$aggregation)) {
    config <- prolfquapp::make_DEA_config_R6(aggregation = agg)
    data_prep <- prolfquapp::ProteinDataPrep$new(lfq, pannot, config)
    if (agg == "rlm") {
      # rlm does not converge for two simulated proteins; they are reported, not dropped
      expect_warning(data_prep$aggregate(), "failed to converge", class = "prolfqua_fit_failed")
    } else {
      data_prep$aggregate()
    }
    expect_true("LFQData" %in% class(data_prep$lfq_data), info = agg)
    expect_equal(data_prep$lfq_data$hierarchy_keys(), "protein_Id", info = agg)
  }
})

make_calibration_prep <- function() {
  samples <- c("run_B", "run_A", "run_C")
  data <- data.frame(
    protein_Id = rep(c("std1", "std2", "target"), each = 3),
    sampleName = rep(samples, 3),
    raw.file = paste0(rep(samples, 3), ".raw"),
    intensity = 2^c(2, 3, 4, 4, 5, 6, 9, 9, 9),
    nr_children = 1
  )
  cfg <- prolfqua::AnalysisConfiguration$new()
  cfg$file_name <- "raw.file"
  cfg$hierarchy <- list(protein_Id = "protein")
  cfg$set_response("intensity")
  lfq <- prolfqua::LFQData$new(data, cfg)
  row_annot <- data.frame(protein_Id = c("std1", "std2", "target"), description = "protein", nr_peptides = 1)
  annotation <- prolfquapp::ProteinAnnotation$new(
    lfq,
    row_annot = row_annot,
    description = "description",
    exp_nr_children = "nr_peptides"
  )
  config <- prolfquapp::make_DEA_config_R6(Normalization = "none", model = "lm")
  config$processing_options$internal <- c("std2", "std1")
  prep <- prolfquapp::ProteinDataPrep$new(lfq, annotation, config)
  prep$lfq_data <- lfq
  prep
}

test_that("protein log2 calibration subtracts the per-sample standard medians", {
  for (method in c("none", "log2")) {
    prep <- make_calibration_prep()
    prep$prolfq_app_config$processing_options$transform <- method
    calibrated <- prep$transform_data()
    expect_identical(calibrated$response(), "normalized_abundance")
    expect_equal(calibrated$data_long()$normalized_abundance, c(-1, -1, -1, 1, 1, 1, 6, 5, 4))
    expect_false("centered_abundance_by_mean" %in% names(calibrated$data_long()))
    expect_identical(prep$lfq_data$response(), "abundance")
  }
})

test_that("protein log2 calibration names samples without a usable standard", {
  prep <- make_calibration_prep()
  data <- prep$lfq_data$data_long()
  standards <- data$protein_Id %in% c("std1", "std2")
  data$intensity[standards & data$sampleName == "run_B"] <- NA_real_
  data <- data[!(standards & data$sampleName == "run_C"), ]
  prep$lfq_data$set_data(data)

  expect_error(prep$transform_data(), "run_B, run_C")
  expect_identical(prep$lfq_data$response(), "intensity")
  expect_null(prep$lfq_data_transformed)

  prep$prolfq_app_config$processing_options$internal <- "unknown"
  expect_error(prep$transform_data(), "run_B, run_A, run_C")
})

test_that("nested protein report preparation retains its existing calibration", {
  prep <- make_calibration_prep()
  prep$prolfq_app_config$processing_options$model <- "ropeca_nested"
  data <- prep$lfq_data$data_long()
  data$intensity[data$sampleName == "run_B" & data$protein_Id != "target"] <- NA_real_
  prep$lfq_data$set_data(data)

  calibrated <- prep$transform_data()$data_long()
  expect_true(all(is.na(calibrated$normalized_abundance[calibrated$sampleName == "run_B"])))
  expect_true("centered_abundance_by_mean" %in% names(calibrated))
})
