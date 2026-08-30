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
      bfabric_urls = list(projectURL = "https://example.org/project/1"),
      report_provenance = list(software = "DIANN", workunit_Id = 42),
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
      analysis_configuration_transformed = list(
        sample_name = "sampleName",
        hierarchy = list(protein_Id = "protein", site = "site")
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
      annotation,
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
    as.matrix(adata$layers[["raw"]]),
    t(SummarizedExperiment::assay(se, "rawData"))
  )
  expect_setequal(
    adata$layers_keys(),
    c("raw", "transformed", "nr_children")
  )
  expect_setequal(
    adata$varm_keys(),
    c(
      "dea__A%2FB",
      "dea__A%252FB",
      "stats_normalized_wide",
      "stats_raw_wide"
    )
  )
  dea <- as.data.frame(adata$varm[["dea__A%2FB"]])
  expect_false(any(c("protein_Id", "site", "SequenceWindow") %in% names(dea)))
  expect_equal(dea$diff, c(1, -1))
  expect_equal(adata$uns$prolfquapp$artifact_type, "dea_results")
  expect_equal(adata$uns$prolfquapp$schema_version, "1.0.0")
  expect_equal(adata$uns$prolfquapp$source_software, "DIANN")
  expect_equal(
    adata$uns$prolfquapp$varm_columns[["dea__A%2FB"]],
    c("diff", "statistic", "p.value", "FDR")
  )
  expect_equal(
    adata$uns$prolfquapp$varm_annotations[["dea__A%2FB"]]$modelName,
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
    as.matrix(restored$layers[["raw"]]),
    t(SummarizedExperiment::assay(se, "rawData"))
  )
  columns <- restored$uns$prolfquapp$varm_columns[["dea__A%2FB"]]
  statistic_column <- match("statistic", columns)
  expect_equal(
    unname(as.matrix(restored$varm[["dea__A%2FB"]])[, statistic_column]),
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
