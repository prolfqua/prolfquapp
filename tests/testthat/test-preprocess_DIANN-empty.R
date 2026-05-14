test_that("preprocess_DIANN reports empty filtered DIA-NN output clearly", {
  report <- tempfile(fileext = ".tsv")
  fasta <- tempfile(fileext = ".fasta")
  writeLines(
    paste(
      c(
        "File.Name",
        "Run",
        "Protein.Group",
        "Protein.Ids",
        "Protein.Names",
        "Genes",
        "PG.Quantity",
        "PG.Q.Value",
        "Lib.PG.Q.Value",
        "Stripped.Sequence",
        "Precursor.Quantity",
        "Precursor.Normalised",
        "PEP"
      ),
      collapse = "\t"
    ),
    report
  )
  writeLines(c(">sp|P1|P1_TEST Test protein", "PEPTIDE"), fasta)
  annotation <- prolfquapp::read_annotation(
    data.frame(
      file = c("sample.raw", "control.raw"),
      name = c("sample", "control"),
      group = c("bait", "control"),
      CONTROL = c("T", "C")
    )
  )

  expect_error(
    prolfquapp::preprocess_DIANN(report, fasta, annotation),
    "DIA-NN report contains no rows after filtering"
  )
})
