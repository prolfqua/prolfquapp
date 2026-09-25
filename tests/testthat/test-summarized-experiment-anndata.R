test_that("DEA SummarizedExperiment maps to typed AnnData slots", {
  se <- make_dea_summarized_experiment()

  adata <- prolfquapp:::summarized_experiment_to_anndata(se)

  expect_equal(adata$n_obs(), ncol(se))
  expect_equal(adata$n_vars(), nrow(se))
  expect_equal(rownames(as.data.frame(adata$obs)), colnames(se))
  expect_equal(rownames(as.data.frame(adata$var)), rownames(se))
  expect_equal(
    names(as.data.frame(adata$var)),
    c("protein_Id", "site", "SequenceWindow")
  )
  expect_equal(
    as.matrix(adata$X),
    t(SummarizedExperiment::assay(se, "transformedData"))
  )
  expect_equal(
    as.matrix(adata$layers[["rawData"]]),
    t(SummarizedExperiment::assay(se, "rawData"))
  )
  expect_setequal(
    adata$layers_keys(),
    c("rawData", "transformedData", "nr_children")
  )
  expect_setequal(
    adata$varm_keys(),
    c(
      "constrast_A%2FB",
      "constrast_A%252FB",
      "stats_normalized_wide",
      "stats_raw_wide"
    )
  )
  dea <- as.data.frame(adata$varm[["constrast_A%2FB"]])
  # A contrast frame carries the feature keys, not the rest of the annotation.
  expect_true(all(c("protein_Id", "site") %in% names(dea)))
  expect_false("SequenceWindow" %in% names(dea))
  expect_equal(dea$diff, c(1, -1))
  expect_equal(adata$uns$prolfquapp$artifact_type, "dea_results")
  expect_equal(adata$uns$prolfquapp$schema_version, "2.1.0")
  expect_equal(adata$uns$prolfquapp$source_software, "DIANN")
  expect_s3_class(adata$varm[["constrast_A%2FB"]], "data.frame")
  expect_equal(dea$modelName, c("lm", "lm"))
  expect_no_error(prolfquapp::validate_prolfquapp_anndata(adata))
  expect_error(
    prolfquapp::LFQData_from_anndata(adata),
    "requires a prolfquapp LFQData artifact"
  )
})

test_that("DEA AnnData writes and reads without changing axes or values", {
  se <- make_dea_summarized_experiment()
  output_dir <- tempfile("dea-anndata-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)
  path <- file.path(output_dir, "AnnData.h5ad")

  result <- prolfquapp:::write_summarized_experiment_h5ad(se, path)
  restored <- anndataR::read_h5ad(result)

  expect_equal(normalizePath(path), result)
  expect_true(file.exists(path))
  expect_equal(rownames(as.data.frame(restored$obs)), colnames(se))
  expect_equal(rownames(as.data.frame(restored$var)), rownames(se))
  expect_equal(
    as.matrix(restored$layers[["rawData"]]),
    t(SummarizedExperiment::assay(se, "rawData"))
  )
  expect_equal(restored$varm[["constrast_A%2FB"]]$statistic, c(2, -2))
  # The varm data frame index is the feature axis, as Python anndata requires.
  expect_equal(
    as.vector(rhdf5::h5read(path, "varm/constrast_A%2FB/_index")),
    rownames(se)
  )
  expect_equal(restored$uns$prolfquapp$provenance$workunit_Id, 42)
  expect_length(
    list.files(output_dir, pattern = "^[.]AnnData-", all.files = TRUE),
    0
  )
})

test_that("DEA AnnData rejects missing assays and misaligned result tables", {
  se <- make_dea_summarized_experiment()
  SummarizedExperiment::assays(se)[["rawData"]] <- NULL
  expect_error(
    prolfquapp:::summarized_experiment_to_anndata(se),
    "missing required assay.*rawData"
  )

  se <- make_dea_summarized_experiment()
  result <- SummarizedExperiment::rowData(se)[["constrast_A/B"]]
  rownames(result) <- rev(rownames(result))
  SummarizedExperiment::rowData(se)[["constrast_A/B"]] <- result
  expect_error(
    prolfquapp:::summarized_experiment_to_anndata(se),
    "not aligned to the feature axis"
  )
})

test_that("AnnData converts back to an equivalent SummarizedExperiment", {
  se <- make_dea_summarized_experiment()

  adata <- prolfquapp:::summarized_experiment_to_anndata(se)
  back <- prolfquapp:::anndata_to_summarized_experiment(adata)

  expect_equal(dim(back), dim(se))
  expect_equal(rownames(back), rownames(se))
  expect_equal(colnames(back), colnames(se))
  expect_setequal(
    SummarizedExperiment::assayNames(back),
    SummarizedExperiment::assayNames(se)
  )
  for (assay_name in SummarizedExperiment::assayNames(se)) {
    expect_equal(
      SummarizedExperiment::assay(back, assay_name),
      SummarizedExperiment::assay(se, assay_name),
      info = assay_name
    )
  }
  expect_equal(
    as.data.frame(SummarizedExperiment::colData(back)),
    as.data.frame(SummarizedExperiment::colData(se))
  )
  # Percent-encoded varm keys are decoded back to the original frame names.
  expect_setequal(
    colnames(SummarizedExperiment::rowData(back)),
    colnames(SummarizedExperiment::rowData(se))
  )
  for (frame_name in colnames(SummarizedExperiment::rowData(se))) {
    expected <- as.data.frame(SummarizedExperiment::rowData(se)[[frame_name]])
    restored <- as.data.frame(SummarizedExperiment::rowData(back)[[frame_name]])
    expect_setequal(names(restored), names(expected))
    expect_equal(restored[names(expected)], expected, info = frame_name)
  }
  # The reshape keys describe the AnnData layout, not the artifact.
  expect_equal(
    S4Vectors::metadata(back)$contrast_configuration,
    S4Vectors::metadata(se)$contrast_configuration
  )
  expect_setequal(
    names(S4Vectors::metadata(back)),
    names(S4Vectors::metadata(se))
  )
})

test_that("DEAResultReader reads a written h5ad file", {
  skip_on_cran()

  dea <- prolfquapp::example_deanalyse(Nprot = 12)
  reporter <- prolfquapp::DEAReportGenerator$new(dea, dea$prolfq_app_config)
  se <- reporter$make_SummarizedExperiment()
  output_dir <- tempfile("dea-anndata-reader-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)
  path <- prolfquapp:::write_summarized_experiment_h5ad(
    se,
    file.path(output_dir, "AnnData.h5ad")
  )

  from_h5ad <- prolfquapp::DEAResultReader$new(path)
  from_se <- prolfquapp::DEAResultReader$new(se)

  expect_equal(from_h5ad$subject_id, from_se$subject_id)
  expect_equal(
    from_h5ad$contrast_config$effect_col,
    from_se$contrast_config$effect_col
  )
  expect_equal(nrow(from_h5ad$contrast_table), nrow(from_se$contrast_table))
  expect_setequal(
    colnames(from_h5ad$contrast_table),
    colnames(from_se$contrast_table)
  )
  expect_equal(
    from_h5ad$lfq_transformed$data_long(),
    from_se$lfq_transformed$data_long()
  )
  expect_equal(
    nrow(from_h5ad$significant(FDR_threshold = 0.25, diff_threshold = 0.5)),
    nrow(from_se$significant(FDR_threshold = 0.25, diff_threshold = 0.5))
  )
})

test_that("DEA AnnData requires an existing output directory", {
  se <- make_dea_summarized_experiment()
  missing_dir <- tempfile("missing-anndata-dir-")

  expect_error(
    prolfquapp:::write_summarized_experiment_h5ad(
      se,
      file.path(missing_dir, "AnnData.h5ad")
    ),
    "output directory does not exist"
  )
})

test_that("the AnnData and the SummarizedExperiment read back the same", {
  se <- make_dea_summarized_experiment()
  output_dir <- tempfile("dea-artifacts-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)
  rds <- file.path(output_dir, "SummarizedExperiment.rds")
  saveRDS(se, rds)
  h5ad <- prolfquapp:::write_summarized_experiment_h5ad(
    se,
    file.path(output_dir, "AnnData.h5ad")
  )

  from_rds <- prolfquapp::DEAResultReader$new(rds)
  from_h5ad <- prolfquapp::DEAResultReader$new(h5ad)

  expect_equal(from_h5ad$contrast_table, from_rds$contrast_table)
  expect_equal(from_h5ad$lfq_raw$data_long(), from_rds$lfq_raw$data_long())
  expect_equal(
    from_h5ad$lfq_transformed$data_long(),
    from_rds$lfq_transformed$data_long()
  )
  expect_equal(
    from_h5ad$significant(0.5, 0.5),
    from_rds$significant(0.5, 0.5)
  )
  expect_equal(
    prolfqua::R6_extract_values(from_h5ad$contrast_config),
    prolfqua::R6_extract_values(from_rds$contrast_config)
  )
  # Metadata the report reads by name, including the frames that would be lost
  # if uns flattened them to lists of columns.
  expect_s3_class(from_h5ad$metadata$contrasts, "data.frame")
  expect_equal(
    tibble::as_tibble(from_h5ad$metadata$contrasts),
    tibble::as_tibble(from_rds$metadata$contrasts)
  )
  expect_equal(
    from_h5ad$metadata$provenance$workunit_Id,
    from_rds$metadata$provenance$workunit_Id
  )
  expect_equal(from_h5ad$metadata$feature_keys, c("protein_Id", "site"))
  expect_equal(from_h5ad$metadata$sample_key, "sampleName")
  # HDF5 hands list names back in its own order; the recorded order restores
  # them, at the top level and in nested lists.
  list_names <- function(x) {
    if (!is.list(x) || is.data.frame(x)) {
      return(NULL)
    }
    c(list(names(x)), lapply(unname(x), list_names))
  }
  expect_identical(list_names(from_h5ad$metadata), list_names(from_rds$metadata))
})

test_that("DEAResultReader keys every table by the artifact's feature keys", {
  se <- make_dea_summarized_experiment()
  imputed <- SummarizedExperiment::assay(se, "transformedData")
  imputed[is.na(imputed)] <- 3
  SummarizedExperiment::assays(se)[["imputedData"]] <- imputed
  SummarizedExperiment::rowData(se)[["imputation"]] <- data.frame(
    protein_Id = c("P1", "P2"),
    site = c("S10", "S20"),
    n_observed = c(3, 2),
    n_imputed = c(0, 1),
    route = c("complete", "fitted"),
    row.names = rownames(se)
  )
  output_dir <- tempfile("dea-reader-keys-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)
  h5ad <- prolfquapp::write_summarized_experiment_h5ad(
    se,
    file.path(output_dir, "AnnData.h5ad")
  )

  for (reader in list(
    prolfquapp::DEAResultReader$new(se),
    prolfquapp::DEAResultReader$new(h5ad)
  )) {
    expect_equal(reader$subject_id, c("protein_Id", "site"))
    for (lfq in list(reader$lfq_raw, reader$lfq_transformed, reader$lfq_imputed)) {
      expect_equal(lfq$hierarchy_keys(), c("protein_Id", "site"))
      expect_setequal(lfq$data_long()$site, c("S10", "S20"))
    }
    abundances <- dplyr::inner_join(
      reader$contrast_table,
      reader$lfq_transformed$data_long(),
      by = reader$subject_id,
      relationship = "many-to-many"
    )
    expect_equal(nrow(abundances), nrow(reader$contrast_table) * ncol(se))
    imputed_long <- reader$lfq_imputed$data_long()
    expect_false(anyNA(imputed_long[[reader$lfq_imputed$response()]]))
    expect_equal(
      reader$imputation$route[reader$imputation$site == "S20"],
      "fitted"
    )
    expect_equal(
      reader$annotation$SequenceWindow,
      c("AAAAASAAAAA", "BBBBBSBBBBB")
    )
    expect_equal(reader$samples$sampleName, colnames(se))
  }
})

test_that("uns rejects a list AnnData cannot store", {
  se <- make_dea_summarized_experiment()
  S4Vectors::metadata(se)$feature_keys <- list("protein_Id", "site")

  expect_error(
    prolfquapp:::summarized_experiment_to_anndata(se),
    "unnamed list at metadata\\$feature_keys"
  )
})

test_that("a single-row metadata table survives the h5ad round-trip", {
  se <- make_dea_summarized_experiment()
  # A single-contrast analysis writes one-row tables. anndataR encodes a
  # one-row data frame's columns as HDF5 scalars, which anndata in Python
  # refuses to read, so the tables are stored column by column instead.
  S4Vectors::metadata(se)$contrasts <- data.frame(
    contrast_name = "A/B",
    contrast = "A - B"
  )
  output_dir <- tempfile("single-contrast-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)

  path <- prolfquapp:::write_summarized_experiment_h5ad(
    se,
    file.path(output_dir, "AnnData.h5ad")
  )
  restored <- prolfquapp::DEAResultReader$new(path)

  expect_s3_class(restored$metadata$contrasts, "data.frame")
  expect_equal(nrow(restored$metadata$contrasts), 1L)
  expect_equal(
    names(restored$metadata$contrasts),
    c("contrast_name", "contrast")
  )
  expect_equal(restored$metadata$contrasts$contrast, "A - B")
  expect_equal(restored$metadata$formula$formula, "abundance ~ group")
})

test_that("uns rejects a table nested inside a list", {
  se <- make_dea_summarized_experiment()
  S4Vectors::metadata(se)$provenance$tables <- data.frame(a = 1:2)

  expect_error(
    prolfquapp:::summarized_experiment_to_anndata(se),
    "nested table at metadata\\$provenance\\$tables"
  )
})
