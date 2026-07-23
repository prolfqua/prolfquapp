# Single explicit quant filtering path (see prolfqua/TODO/TODO_revpattern_handling.md):
#   - contaminants: KEPT everywhere + labelled (never removed);
#   - decoys: KEPT through aggregation/normalization, dropped ONLY at the fit,
#     preserved in the raw/abundance data for export (NA contrast stats).
# No get_subset(clean()) inner-join filter; one join (the export right_join).

make_fixture <- function(Nprot = 40) {
  istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = Nprot)
  lfq <- prolfqua::LFQData$new(istar$data, istar$config)
  d <- lfq$data_long()
  d$protein_Id <- prolfquapp::add_RevCon(d$protein_Id) # inject REV_ (~10%) + zz (~5%)
  lfq$set_data(d)
  ids <- unique(lfq$data_long()$protein_Id)
  addannot <- data.frame(protein_Id = ids, description = paste0(ids, "_desc"))
  pA <- prolfquapp::ProteinAnnotation$new(
    lfq,
    addannot,
    description = "description",
    pattern_contaminants = "^zz",
    pattern_decoys = "^REV"
  )
  GRP2 <- prolfquapp::make_DEA_config_R6()
  GRP2$processing_options$transform <- "robscale"
  GRP2$processing_options$aggregate <- "medpolish"
  GRP2$processing_options$model <- "lm"
  GRP2$processing_options$pattern_decoys <- "^REV"
  GRP2$processing_options$pattern_contaminants <- "^CON|^zz"
  list(lfq = lfq, pA = pA, GRP2 = GRP2, ids = ids)
}

test_that("cont_decoy_summary does not filter; reports decoy + contaminant QC", {
  f <- make_fixture()
  dp <- prolfquapp::ProteinDataPrep$new(f$lfq, f$pA, f$GRP2)
  n_before <- length(unique(dp$lfq_data_peptide$data_long()$protein_Id))
  s <- dp$cont_decoy_summary()
  # QC surfaced (decoys kept, proportion > 0), no removal happened
  expect_true("percentOfDecoys" %in% names(s))
  expect_gt(s$percentOfDecoys, 0)
  expect_equal(
    length(unique(dp$lfq_data_peptide$data_long()$protein_Id)),
    n_before
  )
  # decoys AND contaminants still present in the quant data
  ids <- unique(dp$lfq_data_peptide$data_long()$protein_Id)
  expect_true(any(prolfqua::is_decoy(ids, "^REV")))
  expect_true(any(prolfqua::is_contaminant(ids, "^zz")))
})

test_that("decoys + contaminants survive aggregation and normalization", {
  f <- make_fixture()
  dp <- prolfquapp::ProteinDataPrep$new(f$lfq, f$pA, f$GRP2)
  dp$cont_decoy_summary()
  dp$aggregate()
  dp$transform_data()
  ids <- unique(dp$lfq_data_transformed$data_long()$protein_Id)
  expect_true(any(prolfqua::is_decoy(ids, "^REV")))
  expect_true(any(prolfqua::is_contaminant(ids, "^zz")))
})

test_that("decoys dropped ONLY at the fit; contaminants kept; raw keeps decoys", {
  f <- make_fixture()
  dp <- prolfquapp::ProteinDataPrep$new(f$lfq, f$pA, f$GRP2)
  dp$cont_decoy_summary()
  dp$aggregate()
  dp$transform_data()
  contrasts <- c("AVsC" = "group_A - group_Ctrl")
  dea <- dp$build_deanalyse(contrasts)
  dea$build_default()
  contr <- dea$contrast_results[[dea$default_model]]$get_contrasts()
  cids <- unique(contr$protein_Id)
  # targets-only fit: no decoy in the contrasts
  expect_false(any(prolfqua::is_decoy(cids, "^REV")))
  # contaminants are kept in the fit / contrasts (real proteins)
  expect_true(any(prolfqua::is_contaminant(cids, "^zz")))
  # decoys are preserved in the export-source (raw) LFQData
  raw_ids <- unique(dea$lfq_data_raw$data_long()$protein_Id)
  expect_true(any(prolfqua::is_decoy(raw_ids, "^REV")))
})

test_that("pattern_decoys NULL -> decoy machinery off, decoys kept in the fit", {
  # Opt out is NULL (matches prolfqua::build_contrast_analysis; the app maps an
  # empty REVpattern to NULL upstream). A non-NULL "" would instead opt IN
  # (defaults only), consistent with the core detector.
  f <- make_fixture()
  f$GRP2$processing_options$pattern_decoys <- NULL
  dp <- prolfquapp::ProteinDataPrep$new(f$lfq, f$pA, f$GRP2)
  dp$cont_decoy_summary()
  dp$aggregate()
  dp$transform_data()
  dea <- dp$build_deanalyse(c("AVsC" = "group_A - group_Ctrl"))
  dea$build_default()
  contr <- dea$contrast_results[[dea$default_model]]$get_contrasts()
  expect_true(any(prolfqua::is_decoy(unique(contr$protein_Id), "^REV")))
})

test_that("pattern_decoys '' -> opt in (defaults), decoys dropped at the fit", {
  f <- make_fixture()
  f$GRP2$processing_options$pattern_decoys <- ""
  dp <- prolfquapp::ProteinDataPrep$new(f$lfq, f$pA, f$GRP2)
  dp$cont_decoy_summary()
  dp$aggregate()
  dp$transform_data()
  dea <- dp$build_deanalyse(c("AVsC" = "group_A - group_Ctrl"))
  dea$build_default()
  contr <- dea$contrast_results[[dea$default_model]]$get_contrasts()
  expect_false(any(prolfqua::is_decoy(unique(contr$protein_Id), "^REV")))
})
