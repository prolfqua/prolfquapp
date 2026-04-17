test_that("round-trip peptide-level (depth=2) preserves data", {
  skip_if_not_installed("anndataR")

  res <- prolfquapp::sim_data_protAnnot(Nprot = 20, PROTEIN = FALSE)
  lfqdata <- res$lfqdata
  pannot <- res$pannot

  # Forward: LFQData → AnnData
  adata <- prolfquapp::preprocess_anndata_from_lfq(
    lfqdata,
    pannot,
    source_software = "simulated"
  )

  expect_true(inherits(adata, "AbstractAnnData"))
  expect_true(!is.null(adata$uns[["prolfquapp"]]))
  expect_equal(adata$uns[["prolfquapp"]]$schema_version, "1.0.0")
  expect_equal(adata$uns[["prolfquapp"]]$source_software, "simulated")

  # AnnData dimensions: obs = samples, var = peptides
  n_samples <- nrow(lfqdata$factors())
  n_features <- nrow(lfqdata$data_wide(as.matrix = TRUE)$data)
  expect_equal(nrow(adata$obs), n_samples)
  expect_equal(nrow(adata$var), n_features)

  # Backward: AnnData → LFQData
  back <- prolfquapp::LFQData_from_anndata(adata)

  expect_true(inherits(back$lfqdata, "LFQData"))
  expect_true(inherits(back$protein_annotation, "ProteinAnnotation"))

  # Config fields match
  expect_equal(back$lfqdata$get_config()$sep, lfqdata$get_config()$sep)
  expect_equal(back$lfqdata$file_name(), lfqdata$file_name())
  expect_equal(back$lfqdata$sample_name(), lfqdata$sample_name())
  expect_equal(back$lfqdata$isotope_label(), lfqdata$isotope_label())
  expect_equal(back$lfqdata$config$hierarchy, lfqdata$config$hierarchy)
  expect_equal(
    back$lfqdata$get_config()$hierarchy_depth,
    lfqdata$get_config()$hierarchy_depth
  )
  expect_equal(back$lfqdata$get_config()$factors, lfqdata$get_config()$factors)
  expect_equal(
    back$lfqdata$response(),
    lfqdata$response()
  )

  # Hierarchy counts match
  orig_counts <- lfqdata$hierarchy_counts()
  back_counts <- back$lfqdata$hierarchy_counts()
  expect_equal(nrow(orig_counts), nrow(back_counts))

  # Intensity values match (round-trip fidelity)
  orig_wide <- lfqdata$data_wide(as.matrix = TRUE)$data
  back_wide <- back$lfqdata$data_wide(as.matrix = TRUE)$data
  expect_equal(orig_wide, back_wide, tolerance = 1e-10)

  # ProteinAnnotation fields match
  expect_equal(back$protein_annotation$pID, pannot$pID)
  expect_equal(back$protein_annotation$description, pannot$description)
  expect_equal(
    back$protein_annotation$pattern_contaminants,
    pannot$pattern_contaminants
  )
  expect_equal(back$protein_annotation$pattern_decoys, pannot$pattern_decoys)
})


test_that("round-trip protein-level (depth=1) preserves data", {
  skip_if_not_installed("anndataR")

  res <- prolfquapp::sim_data_protAnnot(Nprot = 20, PROTEIN = TRUE)
  lfqdata <- res$lfqdata
  pannot <- res$pannot

  adata <- prolfquapp::preprocess_anndata_from_lfq(
    lfqdata,
    pannot,
    source_software = "simulated"
  )

  # AnnData dimensions: obs = samples, var = proteins
  n_samples <- nrow(lfqdata$factors())
  n_features <- nrow(lfqdata$data_wide(as.matrix = TRUE)$data)
  expect_equal(nrow(adata$obs), n_samples)
  expect_equal(nrow(adata$var), n_features)

  back <- prolfquapp::LFQData_from_anndata(adata)

  orig_counts <- lfqdata$hierarchy_counts()
  back_counts <- back$lfqdata$hierarchy_counts()
  expect_equal(nrow(orig_counts), nrow(back_counts))

  # Intensity values match (round-trip fidelity)
  orig_wide <- lfqdata$data_wide(as.matrix = TRUE)$data
  back_wide <- back$lfqdata$data_wide(as.matrix = TRUE)$data
  expect_equal(orig_wide, back_wide, tolerance = 1e-10)
})


test_that("NA handling: missing values survive round-trip", {
  skip_if_not_installed("anndataR")

  res <- prolfquapp::sim_data_protAnnot(Nprot = 20, PROTEIN = TRUE)
  lfqdata <- res$lfqdata
  pannot <- res$pannot

  # Inject some NAs
  resp <- lfqdata$response()
  tmp <- lfqdata$data_long()
  na_idx <- sample(seq_len(nrow(tmp)), size = 10)
  tmp[[resp]][na_idx] <- NA
  lfqdata$set_data(tmp)

  adata <- prolfquapp::preprocess_anndata_from_lfq(lfqdata, pannot)
  back <- prolfquapp::LFQData_from_anndata(adata)

  # Count NAs in original and round-tripped
  orig_na <- sum(is.na(lfqdata$data_long()[[resp]]))
  back_na <- sum(is.na(back$lfqdata$data_long()[[resp]]))
  expect_equal(orig_na, back_na)
})


test_that("validation rejects bare AnnData without prolfquapp metadata", {
  skip_if_not_installed("anndataR")

  bare <- anndataR::AnnData(
    X = matrix(1:6, nrow = 2),
    obs = data.frame(s = c("s1", "s2"), row.names = c("s1", "s2")),
    var = data.frame(p = c("p1", "p2", "p3"), row.names = c("p1", "p2", "p3"))
  )

  expect_error(
    prolfquapp::validate_prolfquapp_anndata(bare),
    "uns"
  )

  expect_error(
    prolfquapp::LFQData_from_anndata(bare),
    "uns"
  )
})


test_that("uns schema has all required namespaces", {
  skip_if_not_installed("anndataR")

  res <- prolfquapp::sim_data_protAnnot(Nprot = 10, PROTEIN = TRUE)

  adata <- prolfquapp::preprocess_anndata_from_lfq(
    res$lfqdata,
    res$pannot,
    source_software = "test"
  )

  expect_true("X_layer_name" %in% names(adata$uns))
  expect_true("exploreDE" %in% names(adata$uns))
  expect_true("prolfquapp" %in% names(adata$uns))

  # exploreDE structure
  cr <- adata$uns$exploreDE$column_roles
  expect_true("var" %in% names(cr))
  expect_true("obs" %in% names(cr))
  expect_true("description" %in% names(cr$var))
  expect_true("factor" %in% names(cr$obs))

  # prolfquapp structure
  pm <- adata$uns$prolfquapp
  expect_equal(pm$schema_version, "1.0.0")
  expect_true("analysis_configuration" %in% names(pm))
  expect_true("protein_annotation" %in% names(pm))
  expect_true("layer_names" %in% names(pm))
})
