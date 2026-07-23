dataset <- system.file(
  "application/sim_test/dataset_sim.csv",
  package = "prolfquapp"
)

test_that("run_dea returns deanalyse and supporting data with SIM", {
  skip_on_cran()
  skip_if(nchar(dataset) == 0, "sim_test fixture not installed")

  GRP2 <- prolfquapp::make_DEA_config_R6(
    WORKUNITID = "TEST_DEA",
    Normalization = "robscale"
  )

  result <- prolfquapp::run_dea(
    indir = tempdir(),
    dataset = dataset,
    software = "prolfquapp.SIM",
    config = GRP2
  )

  expect_type(result, "list")
  expect_true("DEAnalyse" %in% class(result$deanalyse))
  expect_true("LFQData" %in% class(result$xd$lfqdata))
  expect_true(
    "ProteinAnnotation" %in% class(result$xd$protein_annotation)
  )
  expect_true(!is.null(result$annotation$contrasts))
  expect_true(length(result$annotation$contrasts) > 0)
  expect_equal(result$software, "prolfquapp.SIM")
  expect_equal(result$requested_software, "prolfquapp.SIM")
  expect_equal(GRP2$software, "prolfquapp.SIM")
})

test_that("run_dea maps legacy prolfqua model to current facade names", {
  skip_on_cran()
  skip_if(nchar(dataset) == 0, "sim_test fixture not installed")

  GRP2 <- prolfquapp::make_DEA_config_R6(
    WORKUNITID = "TEST_DEA_LEGACY_MODEL",
    Normalization = "robscale"
  )
  GRP2$processing_options$model <- "prolfqua"
  GRP2$processing_options$model_missing <- FALSE

  result <- prolfquapp::run_dea(
    indir = tempdir(),
    dataset = dataset,
    software = "prolfquapp.SIM",
    config = GRP2
  )

  expect_equal(result$deanalyse$default_model, "lm")

  GRP2$processing_options$model_missing <- TRUE
  result <- prolfquapp::run_dea(
    indir = tempdir(),
    dataset = dataset,
    software = "prolfquapp.SIM",
    config = GRP2
  )

  expect_equal(result$deanalyse$default_model, "lm_impute")
})

test_that(".resolve_nested_reader switches a nested facade to the peptide reader", {
  available <- c("prolfquapp.DIANN", "prolfquapp.DIANN_PEPTIDE")
  out <- prolfquapp:::.resolve_nested_reader(
    "prolfquapp.DIANN",
    is_nested = TRUE,
    available = available,
    facade = "firth_nested"
  )
  expect_equal(out, "prolfquapp.DIANN_PEPTIDE")
})

test_that(".resolve_nested_reader leaves non-nested and peptide readers unchanged", {
  available <- c("prolfquapp.DIANN", "prolfquapp.DIANN_PEPTIDE")
  # non-nested facade: unchanged
  expect_equal(
    prolfquapp:::.resolve_nested_reader(
      "prolfquapp.DIANN",
      is_nested = FALSE,
      available = available
    ),
    "prolfquapp.DIANN"
  )
  # already a peptide reader: unchanged
  expect_equal(
    prolfquapp:::.resolve_nested_reader(
      "prolfquapp.DIANN_PEPTIDE",
      is_nested = TRUE,
      available = available
    ),
    "prolfquapp.DIANN_PEPTIDE"
  )
})

test_that(".resolve_nested_reader errors when no peptide counterpart is registered", {
  expect_error(
    prolfquapp:::.resolve_nested_reader(
      "prolfquapp.SIM",
      is_nested = TRUE,
      available = c("prolfquapp.SIM"),
      facade = "firth_nested"
    ),
    "no peptide-level counterpart"
  )
})

test_that("run_dea errors when no peptide-level reader exists for a nested facade", {
  skip_if(nchar(dataset) == 0, "sim_test fixture not installed")

  GRP2 <- prolfquapp::make_DEA_config_R6(WORKUNITID = "TEST_DEA_NESTED_NOPEP")
  GRP2$processing_options$model <- "firth_nested"
  GRP2$processing_options$model_missing <- FALSE

  expect_error(
    prolfquapp::run_dea(
      indir = tempdir(),
      dataset = dataset,
      software = "prolfquapp.SIM",
      config = GRP2
    ),
    "no peptide-level counterpart"
  )
})

test_that("run_dea errors on missing annotation file", {
  GRP2 <- prolfquapp::make_DEA_config_R6()
  expect_error(
    prolfquapp::run_dea(
      indir = tempdir(),
      dataset = "/no/such/file.csv",
      software = "prolfquapp.SIM",
      config = GRP2
    ),
    "Annotation file not found"
  )
})

test_that("run_dea errors on unknown software", {
  skip_if(nchar(dataset) == 0, "sim_test fixture not installed")

  GRP2 <- prolfquapp::make_DEA_config_R6()
  expect_error(
    prolfquapp::run_dea(
      indir = tempdir(),
      dataset = dataset,
      software = "NONEXISTENT",
      config = GRP2
    ),
    "not found"
  )
})

test_that("run_dea orchestrates the non-nested workflow and forwards options", {
  dataset_file <- tempfile(fileext = ".csv")
  writeLines("file,group\nsample.raw,A", dataset_file)
  on.exit(unlink(dataset_file), add = TRUE)
  config <- prolfquapp::make_DEA_config_R6(
    model = "lm",
    nr_peptides = 3
  )
  config$processing_options$pattern_contaminants <- "^CONT_"
  config$processing_options$pattern_decoys <- "^DECOY_"
  events <- new.env(parent = emptyenv())
  events$pipeline <- character()
  events$preprocess_args <- NULL

  deanalyse <- new.env(parent = emptyenv())
  deanalyse$build_default <- function() {
    events$pipeline <- c(events$pipeline, "build_default")
  }
  deanalyse$get_annotated_contrasts <- function() {
    events$pipeline <- c(
      events$pipeline,
      "get_annotated_contrasts"
    )
  }
  data_prep <- new.env(parent = emptyenv())
  data_prep$cont_decoy_summary <- function() {
    events$pipeline <- c(
      events$pipeline,
      "cont_decoy_summary"
    )
  }
  data_prep$aggregate <- function() {
    events$pipeline <- c(events$pipeline, "aggregate")
  }
  data_prep$transform_data <- function() {
    events$pipeline <- c(events$pipeline, "transform_data")
  }
  data_prep$build_deanalyse <- function(contrasts) {
    events$contrasts <- contrasts
    deanalyse
  }
  protein_data_prep <- list(
    new = function(...) data_prep
  )

  withr::local_options(prolfqua.progress = NULL)
  local_mocked_bindings(
    get_procfuncs = function() {
      list("prolfquapp.SIM" = list(reader = "sim"))
    },
    read_table_data = function(path) {
      data.frame(file = "sample.raw", group = "A")
    },
    read_annotation = function(data, prefix, SAINT) {
      events$annotation_args <- list(
        prefix = prefix,
        SAINT = SAINT
      )
      list(contrasts = c(A_vs_B = "group_A - group_B"))
    },
    preprocess_software = function(
      indir,
      annotation,
      preprocess_functions,
      pattern_contaminants,
      pattern_decoys,
      nr_peptides
    ) {
      events$preprocess_args <- list(
        indir = indir,
        annotation = annotation,
        preprocess_functions = preprocess_functions,
        pattern_contaminants = pattern_contaminants,
        pattern_decoys = pattern_decoys,
        nr_peptides = nr_peptides
      )
      list(
        xd = list(
          lfqdata = "lfqdata",
          protein_annotation = "protein annotation"
        ),
        files = list(data = "quant.csv", fasta = "database.fasta")
      )
    },
    ProteinDataPrep = protein_data_prep,
    .package = "prolfquapp"
  )

  result <- prolfquapp::run_dea(
    indir = "/analysis",
    dataset = dataset_file,
    software = "prolfquapp.SIM",
    config = config
  )

  expect_identical(result$deanalyse, deanalyse)
  expect_equal(
    events$pipeline,
    c(
      "cont_decoy_summary",
      "aggregate",
      "transform_data",
      "build_default",
      "get_annotated_contrasts"
    )
  )
  expect_equal(events$preprocess_args$nr_peptides, 3)
  expect_equal(
    events$preprocess_args$pattern_contaminants,
    "^CONT_"
  )
  expect_equal(
    events$preprocess_args$pattern_decoys,
    "^DECOY_"
  )
  expect_false(events$annotation_args$SAINT)
  expect_equal(result$software, "prolfquapp.SIM")
  expect_equal(result$requested_software, "prolfquapp.SIM")
  expect_type(getOption("prolfqua.progress"), "closure")
})

test_that("run_dea orchestrates nested peptide models and report aggregation", {
  dataset_file <- tempfile(fileext = ".csv")
  writeLines("file,group\nsample.raw,A", dataset_file)
  on.exit(unlink(dataset_file), add = TRUE)
  config <- prolfquapp::make_DEA_config_R6(model = "firth_nested")
  events <- new.env(parent = emptyenv())
  events$pipeline <- character()
  events$constructor_args <- NULL

  make_lfq <- function(name) {
    object <- new.env(parent = emptyenv())
    object$get_copy <- function() make_lfq(paste0(name, "_copy"))
    object$set_config_value <- function(key, value) {
      events$pipeline <- c(
        events$pipeline,
        paste(name, "set", key, value)
      )
    }
    object$rename_response <- function(response) {
      events$pipeline <- c(
        events$pipeline,
        paste(name, "rename", response)
      )
    }
    object
  }
  peptide_data <- make_lfq("peptide")
  transformed_data <- make_lfq("transformed")
  report_raw <- make_lfq("report_raw")
  report_transformed <- make_lfq("report_transformed")

  data_prep <- new.env(parent = emptyenv())
  data_prep$lfq_data_peptide <- peptide_data
  data_prep$rowAnnot <- "row annotation"
  data_prep$summary <- "summary"
  data_prep$cont_decoy_summary <- function() {
    events$pipeline <- c(
      events$pipeline,
      "cont_decoy_summary"
    )
  }
  report_prep <- new.env(parent = emptyenv())
  report_prep$lfq_data <- report_raw
  report_prep$lfq_data_transformed <- report_transformed
  report_prep$aggregate <- function() {
    events$pipeline <- c(
      events$pipeline,
      "report aggregate"
    )
  }
  report_prep$transform_data <- function() {
    events$pipeline <- c(
      events$pipeline,
      "report transform"
    )
  }
  constructor_count <- 0L
  protein_data_prep <- list(
    new = function(...) {
      constructor_count <<- constructor_count + 1L
      if (constructor_count == 1L) data_prep else report_prep
    }
  )

  deanalyse <- new.env(parent = emptyenv())
  deanalyse$build_default <- function() {
    events$pipeline <- c(events$pipeline, "build_default")
  }
  deanalyse$get_annotated_contrasts <- function() {
    events$pipeline <- c(
      events$pipeline,
      "get_annotated_contrasts"
    )
  }
  deanalyse_peptide_to_protein <- list(
    new = function(...) {
      events$constructor_args <- list(...)
      deanalyse
    }
  )

  local_mocked_bindings(
    get_procfuncs = function() {
      list(
        "prolfquapp.DIANN" = list(reader = "protein"),
        "prolfquapp.DIANN_PEPTIDE" = list(reader = "peptide")
      )
    },
    read_table_data = function(path) {
      data.frame(file = "sample.raw", group = "A")
    },
    read_annotation = function(data, prefix, SAINT) {
      list(contrasts = c(A_vs_B = "group_A - group_B"))
    },
    preprocess_software = function(
      indir,
      annotation,
      preprocess_functions,
      ...
    ) {
      events$selected_reader <- preprocess_functions$reader
      list(
        xd = list(
          lfqdata = "lfqdata",
          protein_annotation = "protein annotation"
        ),
        files = list(data = "quant.csv", fasta = "database.fasta")
      )
    },
    transform_lfqdata = function(lfqdata, method) {
      events$transform_method <- method
      transformed_data
    },
    ProteinDataPrep = protein_data_prep,
    DEAnalysePeptideToProtein = deanalyse_peptide_to_protein,
    .package = "prolfquapp"
  )

  result <- prolfquapp::run_dea(
    indir = "/analysis",
    dataset = dataset_file,
    software = "prolfquapp.DIANN",
    config = config
  )

  expect_equal(
    result$software,
    "prolfquapp.DIANN_PEPTIDE"
  )
  expect_equal(result$requested_software, "prolfquapp.DIANN")
  expect_equal(events$selected_reader, "peptide")
  expect_equal(
    events$constructor_args$default_model,
    "firth_nested"
  )
  expect_identical(
    events$constructor_args$lfq_data,
    transformed_data
  )
  expect_identical(deanalyse$lfq_data, report_transformed)
  expect_identical(deanalyse$lfq_data_raw, report_raw)
  expect_true("report aggregate" %in% events$pipeline)
  expect_true("report transform" %in% events$pipeline)
  expect_true("build_default" %in% events$pipeline)
  expect_true(
    "get_annotated_contrasts" %in% events$pipeline
  )
})

test_that("run_dea rejects an unknown facade before preprocessing", {
  dataset_file <- tempfile(fileext = ".csv")
  writeLines("file,group\nsample.raw,A", dataset_file)
  on.exit(unlink(dataset_file), add = TRUE)
  config <- prolfquapp::make_DEA_config_R6(model = "not-a-facade")

  expect_error(
    prolfquapp::run_dea(
      indir = tempdir(),
      dataset = dataset_file,
      software = "prolfquapp.SIM",
      config = config
    ),
    "Unknown facade: not-a-facade",
    fixed = TRUE
  )
})
