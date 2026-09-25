test_that("lm_impute DEA carries a complete imputedData assay and its imputation block", {
  dea <- example_deanalyse(Nprot = 30)
  se <- DEAReportGenerator$new(dea, dea$prolfq_app_config, name = "")$make_SummarizedExperiment()

  expect_true("imputedData" %in% SummarizedExperiment::assayNames(se))
  imputed <- SummarizedExperiment::assay(se, "imputedData")
  transformed <- SummarizedExperiment::assay(se, "transformedData")
  expect_equal(dim(imputed), dim(transformed))
  expect_false(anyNA(imputed))
  observed <- !is.na(transformed)
  expect_equal(imputed[observed], transformed[observed])

  block <- SummarizedExperiment::rowData(se)[["imputation"]]
  expect_equal(
    colnames(block),
    c(S4Vectors::metadata(se)$feature_keys, "n_observed", "n_imputed", "route")
  )
  expect_equal(rownames(block), rownames(imputed))
  expect_equal(block$n_imputed, unname(rowSums(!observed)))
  expect_equal(block$n_observed, unname(rowSums(observed)))
  expect_true(all(block$route %in% c("complete", "fitted", "lod_refit")))
  expect_equal(S4Vectors::metadata(se)$schema_version, "2.1.0")
})

test_that("imputedData fills a missing cell with the protein model's prediction", {
  dea <- example_deanalyse(Nprot = 30)
  se <- DEAReportGenerator$new(dea, dea$prolfq_app_config, name = "")$make_SummarizedExperiment()
  transformed <- SummarizedExperiment::assay(se, "transformedData")
  imputed <- SummarizedExperiment::assay(se, "imputedData")
  block <- SummarizedExperiment::rowData(se)[["imputation"]]
  protein <- rownames(block)[block$route == "fitted"][1]
  skip_if(is.na(protein), "no protein took the fitted route")

  facade <- dea$contrast_results[[dea$default_model]]
  model_df <- facade$model$model_df
  fit <- model_df$linear_model[[which(model_df$protein_Id == protein)]]
  sample <- colnames(transformed)[is.na(transformed[protein, ])][1]
  annotation <- dplyr::distinct(dplyr::select(
    facade$.lfqdata$data_long(),
    dplyr::all_of(facade$.lfqdata$get_config()$annotation_vars())
  ))
  newdata <- annotation[annotation[[facade$.lfqdata$sample_name()]] == sample, ]
  expect_equal(imputed[protein, sample], unname(stats::predict(fit, newdata = newdata)))
})

test_that("imputedData is written only for the lm_impute default model", {
  dea <- example_deanalyse(Nprot = 30)
  dea$default_model <- "lm"
  dea$build_default()
  se <- DEAReportGenerator$new(dea, dea$prolfq_app_config, name = "")$make_SummarizedExperiment()
  expect_false("imputedData" %in% SummarizedExperiment::assayNames(se))
  expect_null(SummarizedExperiment::rowData(se)[["imputation"]])
})

test_that("imputedData and the imputation block survive the h5ad round trip", {
  dea <- example_deanalyse(Nprot = 30)
  se <- DEAReportGenerator$new(dea, dea$prolfq_app_config, name = "")$make_SummarizedExperiment()
  path <- tempfile(fileext = ".h5ad")
  write_summarized_experiment_h5ad(se, path)
  adata <- anndataR::read_h5ad(path)

  expect_true("imputedData" %in% adata$layers_keys())
  expect_true("imputation" %in% adata$varm_keys())
  expect_equal(
    unname(as.matrix(adata$layers[["imputedData"]])),
    unname(t(SummarizedExperiment::assay(se, "imputedData")))
  )
  expect_false(any(c("n_observed", "n_imputed", "route") %in% colnames(adata$var)))
})

test_that("every rowData frame carries the feature keys", {
  dea <- example_deanalyse(Nprot = 30)
  se <- DEAReportGenerator$new(dea, dea$prolfq_app_config, name = "")$make_SummarizedExperiment()
  keys <- S4Vectors::metadata(se)$feature_keys
  row_data <- SummarizedExperiment::rowData(se)
  for (name in names(row_data)) {
    frame <- as.data.frame(row_data[[name]])
    expect_true(all(keys %in% colnames(frame)), info = name)
    present <- !is.na(frame[[keys[[1]]]])
    expect_equal(frame[[keys[[1]]]][present], rownames(se)[present], info = name)
  }
})
