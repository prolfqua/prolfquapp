reader_test_annotation <- function(
  files = c("sample1.raw", "sample2.raw"),
  file_column = "file"
) {
  data <- data.frame(
    file = files,
    name = c("sample1", "sample2"),
    group = c("control", "treated")
  )
  names(data)[names(data) == "file"] <- file_column
  prolfquapp::read_annotation(data, QC = TRUE)
}

reader_test_fasta_annotation <- function() {
  data.frame(
    fasta.id = "P1",
    fasta.header = "Protein 1",
    proteinname = "P1",
    protein_length = 100,
    nr_tryptic_peptides = 5
  )
}

expect_reader_result <- function(result, hierarchy_depth) {
  expect_s3_class(result$lfqdata, "LFQData")
  expect_s3_class(
    result$protein_annotation,
    "ProteinAnnotation"
  )
  expect_equal(
    result$lfqdata$get_config()$hierarchy_depth,
    hierarchy_depth
  )
  expect_gt(nrow(result$lfqdata$data_long()), 0)
}

test_that("MaxQuant reader builds protein- and peptide-level LFQData", {
  peptide <- tidyr::expand_grid(
    leading.razor.protein = "P1",
    sequence = c("PEPTIDEA", "PEPTIDEB"),
    raw.file = c("sample1", "sample2")
  ) |>
    dplyr::mutate(
      pep = 0.99,
      peptide.intensity = seq_len(dplyr::n())
    )
  local_mocked_bindings(
    tidyMQ_Peptides = function(...) peptide,
    get_annot_from_fasta = function(...) {
      reader_test_fasta_annotation()
    },
    .package = "prolfquapp"
  )
  annotation <- reader_test_annotation()

  protein <- suppressWarnings(
    prolfquapp::preprocess_MQ_peptide(
      "peptides.txt",
      "database.fasta",
      annotation,
      hierarchy_depth = 1
    )
  )
  peptide_level <- suppressWarnings(
    prolfquapp::preprocess_MQ_peptide(
      "peptides.txt",
      "database.fasta",
      annotation,
      hierarchy_depth = 2
    )
  )

  expect_reader_result(protein, 1)
  expect_reader_result(peptide_level, 2)
})

test_that("FragPipe TMT reader preserves the configured hierarchy depth", {
  psm <- tidyr::expand_grid(
    Protein = "P1",
    Peptide = c("PEPTIDEA", "PEPTIDEB"),
    channel = c("sample1", "sample2")
  ) |>
    dplyr::mutate(
      Probability = 0.99,
      Modified.Peptide = .data$Peptide,
      Assigned.Modifications = "",
      abundance = seq_len(dplyr::n()),
      nr_psm = 1
    )
  parse_fun <- function(...) {
    list(
      data = psm,
      nrPeptides = data.frame(
        Protein = "P1",
        nrPeptides = 2
      )
    )
  }
  local_mocked_bindings(
    get_annot_from_fasta = function(...) {
      reader_test_fasta_annotation()
    },
    .package = "prolfquapp"
  )

  result <- suppressWarnings(
    prolfquapp::preprocess_FP_PSM(
      "psm.tsv",
      "database.fasta",
      reader_test_annotation(file_column = "raw.file"),
      hierarchy_depth = 2,
      parse_fun = parse_fun
    )
  )

  expect_reader_result(result, 2)
  expect_equal(
    result$lfqdata$get_config()$nr_children,
    "nr_psm"
  )
})

test_that("MSstats readers build their documented LFQData variants", {
  peptide <- tidyr::expand_grid(
    ProteinName = "P1",
    PeptideSequence = c("PEPTIDEA", "PEPTIDEB"),
    Run = c("sample1", "sample2")
  ) |>
    dplyr::mutate(
      IsotopeLabelType = "L",
      nr_children = 1,
      Intensity = seq_len(dplyr::n())
    )
  local_mocked_bindings(
    read_msstats = function(...) peptide,
    get_annot_from_fasta = function(...) {
      reader_test_fasta_annotation()
    },
    .package = "prolfquapp"
  )
  annotation <- reader_test_annotation()

  generic <- suppressWarnings(
    prolfquapp::preprocess_MSstats(
      "msstats.csv",
      "database.fasta",
      annotation,
      hierarchy_depth = 1
    )
  )
  fragpipe_dia <- suppressWarnings(
    prolfquapp::preprocess_MSstats_FPDIA(
      "msstats.csv",
      "database.fasta",
      annotation,
      hierarchy_depth = 2
    )
  )

  expect_reader_result(generic, 1)
  expect_reader_result(fragpipe_dia, 2)
})

test_that("BGS reader maps Spectronaut hierarchy columns", {
  report <- tidyr::expand_grid(
    PG.ProteinGroups = "P1",
    PEP.GroupingKey = c("PEPTIDEA", "PEPTIDEB"),
    R.FileName = c("sample1", "sample2")
  ) |>
    dplyr::mutate(
      FG.Qvalue = 0.001,
      EG.ModifiedSequence = .data$PEP.GroupingKey,
      FG.Charge = 2,
      FG.Quantity = seq_len(dplyr::n())
    )
  local_mocked_bindings(
    read_BGS = function(...) report,
    get_annot_from_fasta = function(...) {
      reader_test_fasta_annotation()
    },
    .package = "prolfquapp"
  )

  result <- suppressWarnings(
    prolfquapp::preprocess_BGS(
      "report.tsv",
      "database.fasta",
      reader_test_annotation(),
      hierarchy_depth = 2
    )
  )

  expect_reader_result(result, 2)
  expect_contains(
    result$lfqdata$get_config()$hierarchy_keys(),
    "elution_group"
  )
})

test_that("MzMine reader distinguishes annotated and unannotated features", {
  quant_file <- tempfile(fileext = ".csv")
  feature_file <- tempfile(fileext = ".csv")
  writeLines("placeholder\n1", quant_file)
  writeLines("placeholder\n1", feature_file)
  on.exit(unlink(c(quant_file, feature_file)), add = TRUE)
  features <- data.frame(
    feature_id = c("F1", "F2"),
    datafile = c("sample1.csv", "sample2.csv"),
    metabolite_feature_Id = c("F1", "F2"),
    area = c(100, 200)
  )
  feature_annotation <- data.frame(
    id = "F1",
    feature_rt = 1,
    feature_mz = 100,
    feature_charge = 1,
    description = "Known feature"
  )
  local_mocked_bindings(
    tidy_mzMineFeatures = function(...) features,
    make_feature_annotation = function(...) feature_annotation,
    .package = "prolfquapp"
  )
  annotation <- reader_test_annotation(
    files = c("sample1.csv", "sample2.csv"),
    file_column = "relative_path"
  )

  unannotated <- suppressWarnings(
    prolfquapp::preprocess_mzMine(
      quant_file,
      feature_file,
      annotation,
      annotated = FALSE
    )
  )
  annotated <- suppressWarnings(
    prolfquapp::preprocess_mzMine(
      quant_file,
      feature_file,
      annotation,
      annotated = TRUE
    )
  )

  expect_reader_result(unannotated, 1)
  expect_reader_result(annotated, 1)
  expect_gt(
    nrow(unannotated$lfqdata$data_long()),
    nrow(annotated$lfqdata$data_long())
  )
})
