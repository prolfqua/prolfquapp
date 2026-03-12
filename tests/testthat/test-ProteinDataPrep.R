test_that("ProteinDataPrep aggregate returns unchanged site-level data", {
  sim <- prolfqua::sim_lfq_data_peptide_config(Nprot = 10)
  lfq <- prolfqua::LFQData$new(sim$data, sim$config)
  lfq$config$hierarchyDepth <- length(lfq$config$hierarchy_keys())

  row_annot <- data.frame(
    protein_Id = unique(lfq$data$protein_Id),
    description = unique(lfq$data$protein_Id),
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
  expect_equal(data_prep$lfq_data$data, lfq$data)
  expect_equal(
    data_prep$lfq_data$config$hierarchyDepth,
    length(data_prep$lfq_data$config$hierarchy_keys())
  )
  expect_null(data_prep$aggregator)
})
