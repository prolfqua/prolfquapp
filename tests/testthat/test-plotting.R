paired_peptide_lfq <- function() {
  sim <- prolfqua::sim_lfq_data_peptide_config(Nprot = 3, with_missing = FALSE)
  sim$data$subject_ <- sub("^[^_]*_", "", sim$data$sample)
  sim$config$factors$subject_ <- "subject"
  prolfqua::LFQData$new(sim$data, sim$config)
}

test_that("writeLinesPaired draws one plot per protein, peptides as facets", {
  lfq <- paired_peptide_lfq()
  expect_identical(lfq$hierarchy_keys(), c("protein_Id", "peptide_Id"))

  plots <- prolfquapp::writeLinesPaired(lfq)

  proteins <- unique(lfq$data_long()$protein_Id)
  expect_length(plots, length(proteins))
  for (p in plots) {
    built <- ggplot2::ggplot_build(p)
    expect_s3_class(p$facet, "FacetWrap")
    expect_true(p$labels$title %in% proteins)
    expect_gt(nrow(built$data[[1]]), 0)
  }
})

test_that("writeLinesPaired draws one plot per protein for protein-level data", {
  sim <- prolfqua::sim_lfq_data_protein_config(Nprot = 3, with_missing = FALSE, paired = TRUE)
  lfq <- prolfqua::LFQData$new(sim$data, sim$config)

  plots <- prolfquapp::writeLinesPaired(lfq)

  expect_length(plots, 3)
  for (p in plots) {
    expect_s3_class(p$facet, "FacetNull")
    expect_gt(nrow(ggplot2::ggplot_build(p)$data[[1]]), 0)
  }
})

test_that("writeLinesPaired writes all plots to boxplots_paired.pdf", {
  outdir <- withr::local_tempdir()
  prolfquapp::writeLinesPaired(paired_peptide_lfq(), outdir)
  expect_true(file.exists(file.path(outdir, "boxplots_paired.pdf")))
})
