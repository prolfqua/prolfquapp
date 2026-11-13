psm_annotation_columns <- function() {
  tibble::tibble(
    Spectrum = c("run.100.100.2", "run.200.200.2"),
    `Spectrum File` = "run.mzML",
    Peptide = c("PEPTIDEA", "PEPTIDEB"),
    `Modified Peptide` = NA_character_,
    Charge = 2,
    Retention = c(100, 200),
    `Calibrated Observed Mass` = c(1000, 1100),
    Probability = 0.99,
    Intensity = 1e7,
    `Assigned Modifications` = NA_character_,
    Purity = 0.9,
    Protein = "sp|P1|PROT1_HUMAN",
    `Protein Description` = "Protein 1",
    `Mapped Genes` = NA_character_,
    `Mapped Proteins` = NA_character_
  )
}

write_psm <- function(psm) {
  path <- withr::local_tempfile(fileext = ".tsv", .local_envir = parent.frame())
  readr::write_tsv(psm, path, na = "")
  path
}

reporter_values <- list(s1 = c(10, 30), s2 = c(20, 40))

expect_reporter_channels <- function(result) {
  data <- result$data
  expect_setequal(unique(data$channel), c("s1", "s2"))
  expect_equal(sum(data$abundance), 100)
  expect_true(all(data$abundance < 1e7))
}

test_that("FragPipe psm.tsv with Quan Usage before the channels", {
  psm <- dplyr::bind_cols(
    psm_annotation_columns(),
    tibble::tibble(`Quan Usage` = "Use"),
    tibble::as_tibble(reporter_values)
  )

  expect_reporter_channels(prolfquapp::tidy_FragPipe_psm(write_psm(psm)))
})

test_that("FragPipe 22 psm.tsv with the channels after Mapped Proteins", {
  psm <- dplyr::bind_cols(psm_annotation_columns(), tibble::as_tibble(reporter_values))

  expect_reporter_channels(prolfquapp::tidy_FragPipe_psm(write_psm(psm)))
})

test_that("FragPipe 24 psm.tsv reads the Intensity <channel> columns", {
  trace <- tibble::tibble(
    `Parent Scan Number` = 1,
    `Apex Retention Time` = 100,
    `Apex Scan Number` = 1,
    `Retention Time Start` = 99,
    `Retention Time End` = 101,
    `Retention Time FWHM` = 1,
    `Traced Scans` = 5
  )
  per_channel <- function(prefix, values) {
    stats::setNames(tibble::as_tibble(values), paste(prefix, names(values)))
  }
  psm <- dplyr::bind_cols(
    psm_annotation_columns(),
    trace,
    per_channel("Intensity", reporter_values),
    per_channel("Resolution", list(s1 = c(6e4, 6e4), s2 = c(6e4, 6e4))),
    per_channel("SNR", list(s1 = c(500, 500), s2 = c(500, 500)))
  )

  expect_reporter_channels(prolfquapp::tidy_FragPipe_psm(write_psm(psm)))
})

test_that("get_FP_PSM_files finds a FASTA file ending in .fas", {
  dir <- withr::local_tempdir()
  file.create(file.path(dir, c("psm.tsv", "database.fas")))

  files <- prolfquapp::get_FP_PSM_files(dir)

  expect_equal(basename(files$fasta), "database.fas")
  expect_equal(basename(files$data), "psm.tsv")
})
