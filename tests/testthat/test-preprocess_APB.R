# testdata/apb_diann.h5ad is apb-export's prolfqua export of its synthetic
# DIA-NN report (consumers/export_examples.py), with --fasta testdata/apb_diann.fasta:
# 4 runs, 6 precursors, one per peptide and protein, with the qValue, pg_qValue,
# pg_qValue_experiment and pep layers. The FASTA holds P1, P2, P3 and P7.

apb_annotation <- function() {
  prolfquapp::read_annotation(
    data.frame(
      file = c("run_A1", "run_A2", "run_B1", "run_B2"),
      name = c("a1", "a2", "b1", "b2"),
      group = c("A", "A", "B", "B")
    ),
    QC = TRUE
  )
}

apb_fixture <- function() {
  testthat::test_path("testdata", "apb_diann.h5ad")
}

# the fixture, rewritten with `var` changed by `edit_var`
apb_edited_fixture <- function(edit_var) {
  adata <- anndataR::read_h5ad(apb_fixture())
  edited <- anndataR::AnnData(
    X = adata$X,
    obs = as.data.frame(adata$obs),
    var = edit_var(as.data.frame(adata$var)),
    layers = sapply(adata$layers_keys(), function(key) adata$layers[[key]], simplify = FALSE),
    uns = adata$uns
  )
  path <- file.path(withr::local_tempdir(.local_envir = parent.frame()), "apb.h5ad")
  prolfquapp::write_h5ad_atomic(edited, path)
  path
}

test_that("APB reader builds protein- and peptide-level LFQData from an apb-export file", {
  skip_if_not_installed("anndataR")
  protein <- suppressWarnings(prolfquapp::preprocess_APB(apb_fixture(), character(), apb_annotation()))
  peptide <- suppressWarnings(
    prolfquapp::preprocess_APB(apb_fixture(), character(), apb_annotation(), hierarchy_depth = 2)
  )

  expect_equal(protein$lfqdata$hierarchy_keys(), c("protein_Id", "peptide_Id"))
  expect_equal(protein$lfqdata$get_config()$hierarchy_depth, 1)
  expect_equal(peptide$lfqdata$get_config()$hierarchy_depth, 2)
  expect_equal(protein$lfqdata$hierarchy_counts()$protein_Id, 6)
  expect_equal(protein$lfqdata$response(), "intensity")
  annot <- protein$protein_annotation$row_annot
  expect_equal(annot$nrPeptides, rep(1L, 6))
  expect_true(all(c("IDcolumn", "fasta.id", "description", "gene_name") %in% colnames(annot)))
  annot <- annot[order(annot$protein_Id), ]
  expect_equal(annot$protein_length, c(NA, 22, 22, NA, 29, 22))
  expect_equal(annot$nr_tryptic_peptides, c(NA, 1, 2, NA, 2, 2))
})

test_that("APB reader drops observations above the threshold of each q-value layer present", {
  skip_if_not_installed("anndataR")
  # qValue is 0.001 to 0.006 by precursor; pg_qValue_experiment 0.001 to 0.006 by protein
  strict <- suppressWarnings(
    prolfquapp::preprocess_APB(apb_fixture(), character(), apb_annotation(), q_values = list(qValue = 0.0035))
  )
  expect_equal(sort(unique(strict$lfqdata$data_long()$peptide_Id)), c("AACLLK", "PEPTIDEK", "YEASTPEPK"))

  absent <- suppressWarnings(
    prolfquapp::preprocess_APB(apb_fixture(), character(), apb_annotation(), q_values = list(no_such_layer = 0))
  )
  expect_equal(absent$lfqdata$hierarchy_counts()$protein_Id, 6)

  expect_error(
    prolfquapp::preprocess_APB(apb_fixture(), character(), apb_annotation(), q_values = list(pg_qValue = 0)),
    "no observations below the q-value thresholds"
  )
})

test_that("APB reader sums the precursors of a peptide and counts them in nr_children", {
  skip_if_not_installed("anndataR")
  # the AACLLK precursor becomes a second precursor of P1's PEPTIDEK
  path <- apb_edited_fixture(function(var) {
    var$protein_Id[var$peptide_Id == "AACLLK"] <- "P1"
    var$peptide_Id[var$peptide_Id == "AACLLK"] <- "PEPTIDEK"
    var
  })
  xd <- suppressWarnings(prolfquapp::preprocess_APB(path, character(), apb_annotation()))
  p1 <- dplyr::filter(xd$lfqdata$data_long(), .data$peptide_Id == "PEPTIDEK", .data$raw.file == "run_A1")
  expect_equal(p1$intensity, 900 + 1800)
  expect_equal(p1$nr_children, 2)
  expect_equal(p1$qValue, 0.001)
})

test_that("APB files and dataset template come from the folder's one h5ad file", {
  skip_if_not_installed("anndataR")
  files <- prolfquapp::get_APB_files(dirname(apb_fixture()))
  expect_equal(basename(files$data), "apb_diann.h5ad")
  expect_equal(prolfquapp::dataset_template_APB(files)$raw.file, c("run_A1", "run_A2", "run_B1", "run_B2"))
  expect_true(all(c("APB", "APB_PEPTIDE") %in% names(prolfquapp::prolfqua_preprocess_functions)))
})

test_that("LFQData_from_anndata restores the hierarchy order from hierarchy_keys", {
  skip_if_not_installed("anndataR")
  restored <- prolfquapp::LFQData_from_anndata(anndataR::read_h5ad(apb_fixture()))
  expect_equal(restored$lfqdata$hierarchy_keys(), c("protein_Id", "peptide_Id", "precursor_Id"))

  res <- prolfquapp::sim_data_protAnnot(Nprot = 10, PROTEIN = FALSE)
  path <- file.path(withr::local_tempdir(), "lfqdata.h5ad")
  prolfquapp::write_h5ad_atomic(prolfquapp::preprocess_anndata_from_lfq(res$lfqdata, res$pannot), path)
  adata <- anndataR::read_h5ad(path)
  back <- prolfquapp::LFQData_from_anndata(adata)
  expect_equal(back$lfqdata$hierarchy_keys(), res$lfqdata$hierarchy_keys())

  adata$uns$prolfquapp$analysis_configuration$hierarchy_keys <- NULL
  expect_error(prolfquapp::LFQData_from_anndata(adata), "no 'hierarchy_keys'")
})

test_that("a DEA of an APB file writes the input and its results to MuData.h5mu", {
  skip_if_not_installed("anndataR")
  skip_on_cran()
  workdir <- withr::local_tempdir()
  dir.create(file.path(workdir, "in"))
  file.copy(apb_fixture(), file.path(workdir, "in"))
  dataset <- file.path(workdir, "dataset.csv")
  utils::write.csv(
    data.frame(
      file = c("run_A1", "run_A2", "run_B1", "run_B2"),
      name = c("a1", "a2", "b1", "b2"),
      group = c("A", "A", "B", "B"),
      CONTROL = c("C", "C", "T", "T")
    ),
    dataset,
    row.names = FALSE
  )
  config <- prolfquapp::make_DEA_config_R6(
    PATH = file.path(workdir, "out"),
    WORKUNITID = "WU1",
    application = "prolfquapp.APB",
    Normalization = "robscale"
  )
  dir.create(config$get_zipdir(), recursive = TRUE)
  ymlfile <- file.path(workdir, "config.yaml")
  yaml::write_yaml(prolfqua::R6_extract_values(config), ymlfile)
  opt <- list(workunit = "WU1", software = "prolfquapp.APB", dataset = dataset)
  outdir <- withr::with_dir(
    workdir,
    suppressWarnings({
      result <- prolfquapp::run_dea(file.path(workdir, "in"), dataset, "prolfquapp.APB", config)
      prolfquapp:::write_dea_run_outputs(result, config, opt, ymlfile)
    })
  )

  mudata <- prolfquapp::read_h5mu(outdir$data_files$mudata_file)
  expect_equal(names(mudata$modalities), c("lfqdata", "dea"))
  expect_equal(rownames(mudata$obs), c("a1", "a2", "b1", "b2"))
  input <- mudata$modalities$lfqdata
  expect_equal(input$obs_names, c("a1", "a2", "b1", "b2"))
  expect_equal(as.data.frame(input$obs)$raw.file, c("run_A1", "run_A2", "run_B1", "run_B2"))
  expect_true(all(c("qValue", "pg_qValue") %in% input$layers_keys()))
  expect_equal(input$uns$apb$source_level, anndataR::read_h5ad(apb_fixture())$uns$apb$source_level)
  reader <- prolfquapp::DEAResultReader$new(outdir$data_files$mudata_file)
  expect_setequal(unique(reader$contrast_table$protein_Id), mudata$modalities$dea$var_names)
})
