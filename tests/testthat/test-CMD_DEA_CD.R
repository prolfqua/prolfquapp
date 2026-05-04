make_test_cd_zip <- function(path, n_features = 40) {
  samples <- data.frame(
    Sample = paste0("sample_", seq_len(6)),
    Name = paste0("sample_", seq_len(6)),
    Group = rep(c("control", "treated"), each = 3),
    Control = rep(c("C", "T"), each = 3)
  )

  features <- sprintf("Feature%03d_mz%.4f_rt%.2f", seq_len(n_features), 100 + seq_len(n_features), seq_len(n_features) / 10)
  long <- expand.grid(
    Feature_ID = features,
    Sample = samples$Sample,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  long <- merge(long, samples[, c("Sample", "Group")], by = "Sample", sort = FALSE)
  long$Intensity <- 100000 + seq_len(nrow(long)) * 10
  long$Intensity[long$Group == "treated"] <- long$Intensity[long$Group == "treated"] * 1.5

  dup <- long[long$Feature_ID == features[[5]], ]
  dup$Intensity <- dup$Intensity * 0.5
  long <- rbind(long, dup)
  long$Curated <- ifelse(long$Feature_ID %in% features[1:10], "yes", NA)
  long$Filtered <- ifelse(long$Feature_ID %in% features[6:20], "yes", NA)
  long <- long[, c("Feature_ID", "Sample", "Intensity", "Group", "Curated", "Filtered")]

  workdir <- tempfile("cd_zip_")
  dir.create(workdir)
  long_file <- file.path(workdir, "test_long.csv")
  sample_file <- file.path(workdir, "test_prolfqua_samples.csv")
  readr::write_csv(long, long_file)
  readr::write_csv(samples, sample_file)

  oldwd <- setwd(workdir)
  on.exit({
    setwd(oldwd)
    unlink(workdir, recursive = TRUE)
  }, add = TRUE)
  utils::zip(path, files = c(basename(long_file), basename(sample_file)))
  path
}

test_that("make_CD_duplicate_features_unique suffixes duplicated feature/sample rows", {
  input <- data.frame(
    Feature_ID = c("A_mz1_rt1", "A_mz1_rt1", "A_mz1_rt1", "B_mz2_rt2"),
    Sample = c("s1", "s1", "s2", "s1"),
    Intensity = c(1, 2, 3, 4)
  )

  result <- prolfquapp::make_CD_duplicate_features_unique(input)

  expect_equal(
    result$Feature_ID,
    c("A_mz1_rt1 [duplicate 1]", "A_mz1_rt1 [duplicate 2]", "A_mz1_rt1 [duplicate 1]", "B_mz2_rt2")
  )
  expect_false(any(duplicated(result[, c("Feature_ID", "Sample")])))
  expect_equal(result$Feature_ID_original, c("A_mz1_rt1", "A_mz1_rt1", "A_mz1_rt1", "B_mz2_rt2"))
})

test_that("run_dea_cd builds a DEAnalyse object from a CD ZIP export", {
  skip_on_cran()

  zipfile <- file.path(tempdir(), "cd_test_export.zip")
  make_test_cd_zip(zipfile)
  on.exit(unlink(zipfile), add = TRUE)

  cfg <- prolfquapp::make_DEA_config_R6(
    PATH = tempdir(),
    WORKUNITID = "CD_TEST",
    Normalization = "none",
    application = "CompoundDiscoverer"
  )

  result <- prolfquapp::run_dea_cd(zipfile, cfg)
  on.exit(unlink(result$files$tempdir, recursive = TRUE), add = TRUE)

  expect_s3_class(result$deanalyse, "DEAnalyse")
  expect_s3_class(result$xd$lfqdata, "LFQData")
  expect_s3_class(result$xd$protein_annotation, "ProteinAnnotation")
  expect_false(any(duplicated(result$xd$lfqdata$data_long()[, c("metabolite_feature_Id", "Name")])))
  expect_true(any(grepl("\\[duplicate 2\\]", result$xd$lfqdata$data_long()$metabolite_feature_Id)))
  expect_true(length(result$annotation$contrasts) >= 1)
  expect_equal(
    prolfquapp::get_CD_subset_columns(result$files$data),
    c("Curated", "Filtered")
  )

  curated <- prolfquapp::run_dea_cd(
    config = cfg,
    files = result$files,
    subset_column = "Curated"
  )
  curated_ids <- unique(curated$xd$lfqdata$data_long()$metabolite_feature_Id)
  expect_true(length(curated_ids) < length(unique(result$xd$lfqdata$data_long()$metabolite_feature_Id)))
  expect_true(all(grepl("^Feature0(0[1-9]|10)_", curated_ids) | grepl("\\[duplicate ", curated_ids)))

  datax <- result$deanalyse$annotated_contrasts |>
    tidyr::unite("ID", metabolite_feature_Id, sep = "~")
  default_palette <- c(
    Linear_Model_moderated = "black",
    Imputed_Mean_moderated = "orange",
    WaldTest_moderated = "black"
  )
  model_levels <- sort(unique(as.character(stats::na.omit(datax$modelName))))
  palette <- default_palette[model_levels]
  missing_palette <- is.na(palette)
  if (any(missing_palette)) {
    palette[missing_palette] <- grDevices::hcl.colors(
      sum(missing_palette),
      palette = "Dark 3"
    )
  }
  volcano <- prolfqua::volcano_plotly(
    datax,
    proteinID = "ID",
    effect = "diff",
    significance = "FDR",
    contrast = "contrast",
    color = "modelName",
    palette = palette
  )
  built <- plotly::plotly_build(volcano[[1]])
  expect_false(any(vapply(
    built$x$data,
    function(trace) identical(trace$marker$color, "transparent"),
    logical(1)
  )))
})

test_that("CMD_DEA_CD runs the full output pipeline", {
  skip_on_cran()

  script <- system.file("application/CMD_DEA_CD.R", package = "prolfquapp")
  skip_if(nchar(script) == 0, "CMD_DEA_CD.R not installed")
  source_root <- normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)
  skip_if(
    startsWith(normalizePath(script), source_root),
    "CMD_DEA_CD.R subprocess test requires an installed package with current exports"
  )

  zipfile <- file.path(tempdir(), "cd_test_export.zip")
  make_test_cd_zip(zipfile)
  on.exit(unlink(zipfile), add = TRUE)

  workdir <- file.path(tempdir(), "dea_cd_test")
  dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

  rscript <- file.path(R.home("bin"), "Rscript")
  status <- system2(
    rscript,
    c(
      script,
      "-i", zipfile,
      "-n", "none",
      "-w", "CD_TEST",
      "-o", workdir,
      "--subset-columns", "none"
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  exit_code <- attr(status, "status")
  if (!is.null(exit_code) && exit_code != 0) {
    message("CMD_DEA_CD stderr:\n", paste(status, collapse = "\n"))
  }
  expect_true(is.null(exit_code) || exit_code == 0)

  result_dirs <- list.dirs(workdir, recursive = FALSE)
  zipdir <- grep("^DEA_", basename(result_dirs), value = TRUE)
  expect_true(length(zipdir) >= 1)
  expect_true(any(grepl("_none$", zipdir)))

  result_dir <- list.dirs(file.path(workdir, zipdir[[1]]), recursive = FALSE)
  result_dir <- grep("Results_WU_CD_TEST$", result_dir, value = TRUE)
  expect_true(length(result_dir) == 1)
  expect_true(file.exists(file.path(workdir, zipdir[[1]], "index.html")))
  expect_true(file.exists(file.path(result_dir, "lfqdata_normalized.parquet")))
  expect_true(file.exists(file.path(result_dir, "lfqdata.yaml")))
  expect_true(file.exists(file.path(result_dir, "SummarizedExperiment.rds")))
  if (nzchar(Sys.which("quarto"))) {
    quarto_file <- file.path(result_dir, "DE_WUCD_TEST_quarto.html")
    expect_true(file.exists(quarto_file))
    index_html <- paste(readLines(file.path(workdir, zipdir[[1]], "index.html"), warn = FALSE), collapse = "\n")
    expect_match(index_html, "Quarto DEA report", fixed = TRUE)
  }
  expect_true(length(list.files(result_dir, pattern = "[.]xlsx$", full.names = TRUE)) >= 1)
})
