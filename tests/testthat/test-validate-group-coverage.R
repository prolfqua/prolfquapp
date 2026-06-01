# Regression test for the informative error raised when a contrast references a
# group that has no samples after the annotation is matched to the quant data.
# Previously this failed deep inside prolfqua's linfct construction with a
# cryptic "subscript out of bounds".

make_prep <- function(lfq) {
  pA <- data.frame(protein_Id = unique(lfq$data_long()$protein_Id))
  pA$fasta.annot <- paste0(pA$protein_Id, "_desc")
  pA <- prolfquapp::ProteinAnnotation$new(lfq, row_annot = pA, description = "fasta.annot")
  GRP2 <- prolfquapp::make_DEA_config_R6(Normalization = "robscale")
  prep <- prolfquapp::ProteinDataPrep$new(lfq, pA, GRP2)
  prep$aggregate()
  prep$transform_data()
  prep
}

test_that("missing contrast group aborts with an informative error", {
  skip_on_cran()
  pep <- prolfqua::sim_lfq_data_peptide_config(Nprot = 30)
  pep <- prolfqua::LFQData$new(pep$data, pep$config)
  fk <- pep$relevant_factor_keys()[1]

  # Drop one whole group, as would happen if its raw files are absent from the
  # quantification report.
  kept <- dplyr::filter(pep$data_long(), .data[[fk]] != "Ctrl")
  pep_missing <- prolfqua::LFQData$new(kept, pep$get_config())

  dea <- make_prep(pep_missing)$build_deanalyse(
    c("AVsC" = "group_A - group_Ctrl", "BVsA" = "group_B - group_A")
  )

  expect_error(
    dea$build_default(),
    "group 'group_Ctrl' has 0 samples",
    fixed = TRUE
  )
})

test_that("full group coverage does not trigger a false positive", {
  skip_on_cran()
  pep <- prolfqua::sim_lfq_data_peptide_config(Nprot = 30)
  pep <- prolfqua::LFQData$new(pep$data, pep$config)

  dea <- make_prep(pep)$build_deanalyse(
    c("AVsC" = "group_A - group_Ctrl", "BVsC" = "group_B - group_Ctrl")
  )

  expect_no_error(dea$build_default())
})
