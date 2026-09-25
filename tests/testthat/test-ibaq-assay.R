test_that("the ibaq assay lines up with the features and reaches the AnnData", {
  dea <- example_deanalyse(Nprot = 30)
  annotation <- dea$rowAnnot$clone()
  annotation$row_annot$protein_length <- 300
  annotation$row_annot$nr_tryptic_peptides <- 10
  # A protein without annotation gets no IBAQ value.
  dropped <- annotation$row_annot[[annotation$pID]][[1]]
  annotation$row_annot <- annotation$row_annot[-1, ]
  # The same simulated peptides example_deanalyse() starts from.
  pep <- prolfqua::sim_lfq_data_peptide_config(Nprot = 30)
  ibaq <- compute_IBAQ_values(prolfqua::LFQData$new(pep$data, pep$config), annotation)

  se <- DEAReportGenerator$new(dea, dea$prolfq_app_config, name = "")$make_SummarizedExperiment(ibaq = ibaq)
  values <- SummarizedExperiment::assay(se, "ibaq")
  expect_equal(dimnames(values), dimnames(SummarizedExperiment::assay(se, "rawData")))
  expect_true(all(is.na(values[dropped, ])))
  expected <- strip_rownames(ibaq$data_wide(as.matrix = TRUE)$data, "~lfq~light")
  shared <- setdiff(intersect(rownames(values), rownames(expected)), dropped)
  expect_gt(length(shared), 0)
  expect_equal(values[shared, ], expected[shared, colnames(values)])

  path <- file.path(withr::local_tempdir(), "AnnData.h5ad")
  restored <- anndataR::read_h5ad(prolfquapp:::write_summarized_experiment_h5ad(se, path))
  expect_equal(as.matrix(restored$layers[["ibaq"]]), t(values), ignore_attr = TRUE)
})

test_that("without ibaq the SummarizedExperiment has no ibaq assay", {
  dea <- example_deanalyse(Nprot = 30)
  se <- DEAReportGenerator$new(dea, dea$prolfq_app_config, name = "")$make_SummarizedExperiment()
  expect_false("ibaq" %in% SummarizedExperiment::assayNames(se))
})
