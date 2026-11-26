test_that("set_list_to_R6 sets scalar and nested fields on an R6 config", {
  cfg <- prolfquapp::ProlfquAppConfig$new(
    prolfquapp::ProcessingOptions$new(),
    prolfquapp::ProjectSpec$new(),
    prolfquapp::ExternalReader$new()
  )

  config_list <- list(
    software = "DIANN",
    prefix = "DEA",
    processing_options = list(
      transform = "robscale",
      FDR_threshold = 0.05
    ),
    project_spec = list(
      project_Id = "1234",
      workunit_Id = "WU99"
    )
  )

  prolfquapp::set_list_to_R6(config_list, cfg)

  expect_equal(cfg$software, "DIANN")
  expect_equal(cfg$prefix, "DEA")
  expect_equal(cfg$processing_options$transform, "robscale")
  expect_equal(cfg$processing_options$FDR_threshold, 0.05)
  expect_equal(cfg$project_spec$project_Id, "1234")
  expect_equal(cfg$project_spec$workunit_Id, "WU99")
})

test_that("set_list_to_R6 returns the (mutated) R6 object invisibly", {
  cfg <- prolfquapp::ProlfquAppConfig$new(
    prolfquapp::ProcessingOptions$new(),
    prolfquapp::ProjectSpec$new(),
    prolfquapp::ExternalReader$new()
  )

  out <- prolfquapp::set_list_to_R6(list(software = "FP_TMT"), cfg)

  expect_identical(out, cfg)
  expect_equal(cfg$software, "FP_TMT")
})

test_that("make_DEA_config_R6 offers only aggregations and normalizations the pipeline runs", {
  expect_true(all(
    eval(formals(prolfquapp::make_DEA_config_R6)$Normalization) %in%
      eval(formals(prolfquapp::transform_lfqdata)$method)
  ))
  cfg <- prolfquapp::make_DEA_config_R6(aggregation = "topN")
  expect_equal(cfg$processing_options$aggregate, "topN")
  expect_error(prolfquapp::make_DEA_config_R6(aggregation = "top3"), "should be one of")
  expect_error(prolfquapp::make_DEA_config_R6(Normalization = "quantile"), "should be one of")
})

test_that("use_identification_confidence reaches the YAML list, the report table and the modelling config", {
  cfg <- prolfquapp::ProlfquAppConfig$new(
    prolfquapp::ProcessingOptions$new(),
    prolfquapp::ProjectSpec$new(),
    prolfquapp::ExternalReader$new()
  )
  expect_true(cfg$processing_options$use_identification_confidence)
  prolfquapp::set_list_to_R6(list(processing_options = list(use_identification_confidence = FALSE)), cfg)
  expect_false(prolfqua::r6_extract_values(cfg)$processing_options$use_identification_confidence)
  table <- cfg$processing_options$parameters_table(software = "DIANN")
  expect_equal(table$Value[table$Parameter == "Identification confidence as model weight"], "no")

  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 5)
  lfq <- prolfqua::LFQData$new(istar$data, istar$config)
  dea <- prolfquapp::DEAnalyse$new(lfq, NULL, cfg, c(A_vs_Ctrl = "group_A - group_Ctrl"), default_model = "lm")
  expect_false(dea$lfq_data$get_config()$use_identification_confidence)
})
