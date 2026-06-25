# Tests for reading native DIA-NN 2.x output (parquet with bare `Run`, no `File.Name`).
# See TODO/TODO_diann2x_native_output.md.

# minimal DIA-NN report rows that pass the q-value filter in diann_read_output()
make_report_rows <- function(run_values, run_col = "Run") {
  n <- length(run_values)
  df <- data.frame(
    Protein.Group = rep("P11111", n),
    Protein.Names = rep("PROT1_TEST", n),
    Stripped.Sequence = rep("SAMPLEPEPTIDEK", n),
    Precursor.Quantity = rep(1000, n),
    Precursor.Normalised = rep(1000, n),
    PEP = rep(0.001, n),
    PG.Quantity = rep(5000, n),
    PG.Q.Value = rep(0.001, n),
    Lib.PG.Q.Value = rep(0.001, n),
    stringsAsFactors = FALSE
  )
  df[[run_col]] <- run_values
  df
}

test_that(".normalize_raw_file strips expected raw-file wrappers", {
  raw_files <- c(
    "xsample.raw",
    "C:\\runs\\sample.d",
    "/mnt/data/control.mzML",
    "fraction.d.zip",
    "plain_sample"
  )

  expect_identical(
    prolfquapp:::.normalize_raw_file(raw_files),
    c("sample", "sample", "control", "fraction", "plain_sample")
  )
})

test_that(".diagnose_sample_join reports annotation-only files but allows quant-only files", {
  log_file <- withr::local_tempfile()
  logger::log_appender(logger::appender_file(log_file))
  withr::defer(logger::log_appender(logger::appender_console))

  result <- prolfquapp:::.diagnose_sample_join(
    annotation_keys = c("sample", "missing"),
    quant_keys = c("sample", "extra"),
    matched_keys = "sample",
    context = "unit-test"
  )

  expect_equal(result$annotation_missing, "missing")
  expect_equal(result$quant_only, "extra")
  expect_match(
    paste(readLines(log_file, warn = FALSE), collapse = "\n"),
    "unit-test: annotated files not found in quantification data: missing"
  )
})

test_that(".diagnose_sample_join does not warn for quant-only files", {
  log_file <- withr::local_tempfile()
  logger::log_appender(logger::appender_file(log_file))
  withr::defer(logger::log_appender(logger::appender_console))

  result <- prolfquapp:::.diagnose_sample_join(
    annotation_keys = "sample",
    quant_keys = c("sample", "extra"),
    matched_keys = "sample",
    context = "unit-test"
  )

  expect_equal(result$annotation_missing, character())
  expect_equal(result$quant_only, "extra")
  expect_false(file.exists(log_file))
})

test_that("diann_read_output derives identical raw.file from native `Run` and renamed `File.Name`", {
  runs <- c(
    "20260623_010_C42222_S1172811_Plate_7001_H1",
    "20260623_011_C42222_S1172812_Plate_7001_H2"
  )
  # native DIA-NN 2.x: bare `Run`, no `File.Name`
  native <- prolfquapp::diann_read_output(make_report_rows(
    runs,
    run_col = "Run"
  ))
  # legacy shim: `Run` renamed to `File.Name`, value is still the bare basename
  renamed <- prolfquapp::diann_read_output(make_report_rows(
    runs,
    run_col = "File.Name"
  ))

  expect_identical(native$raw.file, runs)
  expect_identical(native$raw.file, renamed$raw.file)
})

test_that("diann_read_output falls back to DIA-NN 1.x `File.Name` full paths", {
  file_names <- c(
    "D:\\data\\runs\\20260623_010_C42222_S1172811_Plate_7001_H1.raw",
    "/mnt/data/20260623_011_C42222_S1172812_Plate_7001_H2.d"
  )
  out <- prolfquapp::diann_read_output(make_report_rows(
    file_names,
    run_col = "File.Name"
  ))
  expect_identical(
    out$raw.file,
    c(
      "20260623_010_C42222_S1172811_Plate_7001_H1",
      "20260623_011_C42222_S1172812_Plate_7001_H2"
    )
  )
})

test_that("diann_read_output errors when neither `Run` nor `File.Name` is present", {
  rows <- make_report_rows("x", run_col = "Run")
  rows$Run <- NULL
  expect_error(
    prolfquapp::diann_read_output(rows),
    "neither 'Run' nor 'File.Name'"
  )
})

test_that("get_DIANN_files discovers and prefers the native parquet", {
  dir <- withr::local_tempdir()
  parquet <- file.path(dir, "report.parquet")
  tsv <- file.path(dir, "report.tsv")
  arrow::write_parquet(make_report_rows("sample", run_col = "Run"), parquet)
  readr::write_tsv(make_report_rows("sample", run_col = "Run"), tsv)
  writeLines(
    c(">sp|P11111|PROT1_TEST one", "PEPTIDEK"),
    file.path(dir, "database.fasta")
  )

  files <- prolfquapp::get_DIANN_files(dir)
  # only the parquet is returned when both report.parquet and report.tsv exist
  expect_equal(basename(files$data), "report.parquet")
})

test_that("preprocess_DIANN builds an LFQData from a native DIA-NN 2.x parquet", {
  dir <- withr::local_tempdir()
  runs <- c("sample", "control")
  proteins <- c("P11111", "P22222")
  peptides <- list(
    P11111 = c("SAMPLEPEPTIDEK", "ANOTHERPEPTIDER"),
    P22222 = c("THIRDPEPTIDEK", "FOURTHPEPTIDER")
  )

  rows <- do.call(
    rbind,
    lapply(runs, function(run) {
      do.call(
        rbind,
        lapply(proteins, function(pg) {
          data.frame(
            Run = run,
            Run.Index = match(run, runs) - 1L,
            Protein.Group = pg,
            Protein.Names = paste0(pg, "_TEST"),
            Stripped.Sequence = peptides[[pg]],
            Precursor.Quantity = c(2000, 3000),
            Precursor.Normalised = c(2000, 3000),
            PEP = 0.001,
            PG.Quantity = 5000,
            PG.Q.Value = 0.001,
            Lib.PG.Q.Value = 0.001,
            stringsAsFactors = FALSE
          )
        })
      )
    })
  )

  parquet <- file.path(dir, "report.parquet")
  arrow::write_parquet(rows, parquet)
  expect_false("File.Name" %in% names(arrow::read_parquet(parquet)))

  fasta <- file.path(dir, "database.fasta")
  writeLines(
    c(
      ">sp|P11111|PROT1_TEST Protein one OS=Test OX=1 GN=GENE1 PE=1 SV=1",
      "MSAMPLEPEPTIDEKANOTHERPEPTIDERAAAAAAAAAA",
      ">sp|P22222|PROT2_TEST Protein two OS=Test OX=1 GN=GENE2 PE=1 SV=1",
      "MTHIRDPEPTIDEKFOURTHPEPTIDERBBBBBBBBBB"
    ),
    fasta
  )

  annotation <- prolfquapp::read_annotation(
    data.frame(
      file = c("C:\\runs\\sample.d", "/mnt/data/control.mzML"),
      name = c("sample", "control"),
      group = c("bait", "control"),
      CONTROL = c("T", "C")
    )
  )

  result <- suppressWarnings(suppressMessages(
    prolfquapp::preprocess_DIANN(parquet, fasta, annotation)
  ))

  expect_s3_class(result$lfqdata, "LFQData")
  # raw.file (the file_name column) is derived from the parquet's bare `Run`
  expect_equal(result$lfqdata$file_name(), "raw.file")
  expect_setequal(unique(result$lfqdata$data_long()[["raw.file"]]), runs)
})
