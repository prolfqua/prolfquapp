test_that("ProteinDataPrep aggregate returns unchanged site-level data", {
  sim <- prolfqua::sim_lfq_data_peptide_config(Nprot = 10)
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
  sim <- prolfqua::sim_lfq_data_peptide_config(Nprot = 10)
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
    data_prep$aggregate()
    expect_true("LFQData" %in% class(data_prep$lfq_data), info = agg)
    expect_equal(data_prep$lfq_data$hierarchy_keys(), "protein_Id", info = agg)
  }
})
