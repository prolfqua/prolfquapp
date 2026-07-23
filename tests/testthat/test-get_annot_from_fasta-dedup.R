# get_annot_from_fasta() is a pure parser: it returns every FASTA record as read.
# Decoy removal and protein-ID uniqueness are now resolved downstream by
# ProteinAnnotation, so this reader no longer dedups or filters decoys.

test_that("get_annot_from_fasta returns every record (no dedup of accessions)", {
  # Two overlapping FASTAs both reducing to accession "ACC1" (the WU347806 case):
  # a pure parser must return BOTH rows, not collapse them.
  seq1 <- paste(rep("A", 28), collapse = "")
  file_a <- tempfile(fileext = ".fasta")
  file_b <- tempfile(fileext = ".fasta")
  on.exit(unlink(c(file_a, file_b)), add = TRUE)
  writeLines(
    c(
      ">tr|ACC1|ACC1_TEST Protein one OS=Test OX=1 GN=GENEA PE=4 SV=1",
      seq1,
      ">sp|ACC2|ACC2_TEST Protein two OS=Test OX=1 GN=GENE2 PE=1 SV=1",
      seq1
    ),
    file_a
  )
  writeLines(
    c(
      ">sp|ACC1|ACC1_TEST Protein one reviewed OS=Test OX=1 GN=GENEB PE=1 SV=2",
      seq1
    ),
    file_b
  )

  fa <- suppressWarnings(suppressMessages(
    prolfquapp:::get_annot_from_fasta(c(file_a, file_b))
  ))
  # 3 records in, 3 records out; ACC1 appears twice (duplicate retained)
  expect_equal(nrow(fa), 3)
  expect_equal(sum(fa$proteinname == "ACC1"), 2)
})

test_that("get_annot_from_fasta does not filter decoys (deferred to ProteinAnnotation)", {
  # Even with a decoy pattern supplied, the parser keeps REV_ records.
  seq1 <- paste(rep("A", 28), collapse = "")
  fasta <- tempfile(fileext = ".fasta")
  on.exit(unlink(fasta), add = TRUE)
  writeLines(
    c(
      ">sp|ACC1|ACC1_TEST forward OS=Test OX=1 GN=GENE1 PE=1 SV=1",
      seq1,
      ">REV_sp|ACC1|ACC1_TEST decoy OS=Test OX=1 GN=GENE1 PE=1 SV=1",
      seq1
    ),
    fasta
  )

  fa <- suppressWarnings(suppressMessages(
    prolfquapp:::get_annot_from_fasta(fasta, pattern_decoys = "^REV_|^rev_")
  ))
  expect_equal(nrow(fa), 2)
  expect_true(any(grepl("^REV_", fa$fasta.id)))
})
