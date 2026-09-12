make_dea_summarized_experiment <- function() {
  feature_names <- c("P1~S10", "P2~S20")
  sample_names <- c("S1", "S2", "S3")
  raw <- matrix(
    c(10, NA, 30, 20, 40, 60),
    nrow = 2,
    dimnames = list(feature_names, sample_names)
  )
  transformed <- log2(raw)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      rawData = raw,
      transformedData = transformed,
      nr_children = matrix(
        c(2, 3, 2, 3, 2, 3),
        nrow = 2,
        dimnames = list(feature_names, sample_names)
      )
    ),
    colData = S4Vectors::DataFrame(
      sampleName = sample_names,
      group = c("A", "A", "B"),
      row.names = sample_names
    ),
    metadata = list(
      artifact_type = "dea_results",
      schema_version = "2.0.0",
      source_software = "DIANN",
      feature_keys = list("protein_Id", "site"),
      sample_key = "sampleName",
      bfabric_urls = list(projectURL = "https://example.org/project/1"),
      provenance = list(software = "DIANN", workunit_Id = 42),
      contrasts = data.frame(
        contrast_name = c("A/B", "A%2FB"),
        contrast = c("A - B", "A + B")
      ),
      formula = data.frame(formula = "abundance ~ group"),
      default_model = "lm",
      analysis_configuration_raw = list(
        sample_name = "sampleName",
        hierarchy = list(protein_Id = "protein", site = "site")
      ),
      analysis_configuration = list(
        sample_name = "sampleName",
        hierarchy = list(protein_Id = "protein", site = "site")
      ),
      contrast_configuration = list(
        subject_id = "protein_Id",
        model_name_col = "modelName",
        contrast_col = "contrast",
        effect_col = "diff",
        score_col = "statistic",
        pvalue_col = "p.value",
        fdr_col = "FDR",
        avg_abundance_col = "avgAbd"
      )
    )
  )
  annotation <- data.frame(
    protein_Id = c("P1", "P2"),
    site = c("S10", "S20"),
    SequenceWindow = c("AAAAASAAAAA", "BBBBBSBBBBB"),
    row.names = feature_names
  )
  contrast_result <- function(contrast, diff) {
    data.frame(
      modelName = "lm",
      estimate_type = "observed",
      contrast = contrast,
      diff = diff,
      statistic = diff * 2,
      p.value = c(0.01, 0.2),
      FDR = c(0.02, 0.2),
      row.names = feature_names
    )
  }
  SummarizedExperiment::rowData(se)[["annotation"]] <- annotation
  SummarizedExperiment::rowData(se)[["constrast_A/B"]] <-
    contrast_result("A/B", c(1, -1))
  SummarizedExperiment::rowData(se)[["constrast_A%2FB"]] <-
    contrast_result("A%2FB", c(2, -2))
  SummarizedExperiment::rowData(se)[["stats_normalized_wide"]] <-
    data.frame(
      protein_Id = c("P1", "P2"),
      site = c("S10", "S20"),
      meanAbundance_All = c(5, 6),
      row.names = feature_names
    )
  SummarizedExperiment::rowData(se)[["stats_raw_wide"]] <-
    data.frame(
      protein_Id = c("P1", "P2"),
      site = c("S10", "S20"),
      meanAbundance_All = c(25, 40),
      row.names = feature_names
    )
  se
}

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
  expect_false(any(c("protein_Id", "site", "SequenceWindow") %in% names(dea)))
  expect_equal(dea$diff, c(1, -1))
  expect_equal(adata$uns$prolfquapp$artifact_type, "dea_results")
  expect_equal(adata$uns$prolfquapp$schema_version, "2.0.0")
  expect_equal(adata$uns$prolfquapp$source_software, "DIANN")
  expect_equal(
    adata$uns$prolfquapp$varm_columns[["constrast_A%2FB"]],
    c("diff", "statistic", "p.value", "FDR")
  )
  expect_equal(
    adata$uns$prolfquapp$varm_annotations[["constrast_A%2FB"]]$modelName,
    c("lm", "lm")
  )
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
  columns <- restored$uns$prolfquapp$varm_columns[["constrast_A%2FB"]]
  statistic_column <- match("statistic", columns)
  expect_equal(
    unname(as.matrix(restored$varm[["constrast_A%2FB"]])[, statistic_column]),
    c(2, -2)
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
