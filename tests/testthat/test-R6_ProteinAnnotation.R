# Tests for ProteinAnnotation (single uniqueness + decoy authority).
# The REV-column API (annotate_decoys, clean(decoys=), nr_clean(decoys=)) is gone;
# uniqueness + decoy resolution happen at construction, decoy removal from the
# annotation is pattern-gated in clean().

# ---- .detect_decoy_ids: defaults vs configured pattern --------------------
test_that(".detect_decoy_ids uses built-in defaults for NULL/empty/'a^'", {
  ids <- c("REV_sp|P1|X", "sp|P2|X", "decoy_P3", "DECOY_P4", "normalProtein")
  for (pat in list(NULL, "", "a^")) {
    res <- prolfquapp:::.detect_decoy_ids(ids, pattern = pat)
    expect_equal(res, c(TRUE, FALSE, TRUE, TRUE, FALSE))
  }
})

test_that(".detect_decoy_ids unions a configured pattern with the defaults", {
  ids <- c("CUSTOM_P1", "REV_sp|P2|X", "sp|P3|X")
  res <- prolfquapp:::.detect_decoy_ids(ids, pattern = "^CUSTOM_")
  expect_equal(res, c(TRUE, TRUE, FALSE))
})

# ---- .resolve_unique_protein_ids: the resolution chain --------------------
test_that(".resolve_unique_protein_ids leaves unique input untouched", {
  df <- data.frame(protein_Id = c("P1", "P2"),
    fid = c("sp|P1|G", "tr|P2|G"), stringsAsFactors = FALSE)
  res <- prolfquapp:::.resolve_unique_protein_ids(df, "protein_Id", "fid")
  expect_equal(nrow(res), 2)
})

test_that(".resolve_unique_protein_ids drops the decoy, keeps the forward", {
  df <- data.frame(
    protein_Id = c("P1", "P1", "P2"),
    fid = c("sp|P1|G", "REV_sp|P1|G", "sp|P2|G"),
    stringsAsFactors = FALSE
  )
  res <- suppressWarnings(
    prolfquapp:::.resolve_unique_protein_ids(df, "protein_Id", "fid", "^REV")
  )
  expect_equal(nrow(res), 2)
  expect_equal(res$fid[res$protein_Id == "P1"], "sp|P1|G")
})

test_that(".resolve_unique_protein_ids prefers sp| over tr| (tiebreak)", {
  df <- data.frame(
    protein_Id = c("P1", "P1"),
    fid = c("tr|P1|G", "sp|P1|G"),
    stringsAsFactors = FALSE
  )
  res <- suppressWarnings(
    prolfquapp:::.resolve_unique_protein_ids(df, "protein_Id", "fid")
  )
  expect_equal(nrow(res), 1)
  expect_equal(res$fid, "sp|P1|G")
})

test_that(".resolve_unique_protein_ids falls back to keep-first", {
  df <- data.frame(
    protein_Id = c("P1", "P1"),
    fid = c("tr|P1|A", "tr|P1|B"),
    stringsAsFactors = FALSE
  )
  res <- suppressWarnings(
    prolfquapp:::.resolve_unique_protein_ids(df, "protein_Id", "fid")
  )
  expect_equal(nrow(res), 1)
  expect_equal(res$fid, "tr|P1|A")
})

test_that(".resolve_unique_protein_ids keeps one row for an all-decoy group", {
  df <- data.frame(
    protein_Id = c("P1", "P1"),
    fid = c("REV_sp|P1|G", "REV_tr|P1|G"),
    stringsAsFactors = FALSE
  )
  res <- suppressWarnings(
    prolfquapp:::.resolve_unique_protein_ids(df, "protein_Id", "fid", "^REV")
  )
  expect_equal(nrow(res), 1)
})

# ---- Integration via the real constructor (sim data) ----------------------
test_that("ProteinAnnotation summary/clean: contaminants + pattern-gated decoys", {
  res <- suppressWarnings(suppressMessages(
    sim_data_protAnnot(Nprot = 100, PROTEIN = TRUE)
  ))
  pa <- res$pannot
  testthat::skip_if(nrow(pa$row_annot) != 100)

  # one row per protein id; 5 contaminants (zz), 10 decoys (REV) in sim data
  expect_equal(nrow(pa$row_annot), 100)
  expect_equal(pa$annotate_contaminants(), 5)

  # configured decoy pattern is exposed
  expect_equal(pa$get_rev_pattern(), "^REV")

  # clean() removes the 5 contaminants AND the 10 decoys (pattern-gated)
  cleaned <- pa$clean()
  expect_equal(nrow(cleaned), 85)
  expect_equal(pa$nr_clean(), 85)

  # summary has no decoy fields anymore
  s <- pa$get_summary()
  expect_true(all(c("totalNrOfProteins", "percentOfContaminants") %in% names(s)))
  expect_false("percentOfFalsePositives" %in% names(s))
  expect_false("NrOfProteinsNoDecoys" %in% names(s))
  expect_equal(s$totalNrOfProteins, 100)

  # the old REV API is gone
  expect_false("annotate_decoys" %in% names(pa))
})

test_that("clean() keeps decoys when no decoy pattern is configured", {
  res <- suppressWarnings(suppressMessages(
    sim_data_protAnnot(Nprot = 100, PROTEIN = TRUE)
  ))
  lfqdata <- res$lfqdata
  pids <- unique(lfqdata$data_long()$protein_Id)
  addannot <- data.frame(
    protein_Id = pids, cleanID = pids, description = "d",
    nr_peptides = 2, stringsAsFactors = FALSE
  )
  pa <- suppressWarnings(suppressMessages(ProteinAnnotation$new(
    lfqdata, addannot, description = "description", cleaned_ids = "cleanID",
    exp_nr_children = "nr_peptides",
    pattern_contaminants = "^zz", pattern_decoys = NULL
  )))
  expect_null(pa$get_rev_pattern())
  pa$annotate_contaminants()
  cleaned <- pa$clean()
  # only contaminants removed; REV decoys remain (no pattern configured)
  expect_true(any(grepl("^REV", cleaned$protein_Id)))
})

test_that("ProteinAnnotation$new resolves a duplicate id (keeps forward)", {
  res <- suppressWarnings(suppressMessages(
    sim_data_protAnnot(Nprot = 30, PROTEIN = TRUE)
  ))
  lfqdata <- res$lfqdata
  fwd <- grep("^zz|^REV", unique(lfqdata$data_long()$protein_Id),
    value = TRUE, invert = TRUE)
  testthat::skip_if(length(fwd) < 2)
  pick <- fwd[1]
  ra <- data.frame(
    protein_Id = c(fwd, pick),
    cleanID = c(fwd, pick),
    description = "d",
    fasta.id = c(paste0("sp|", fwd, "|G"), paste0("REV_sp|", pick, "|G")),
    nr_peptides = 2,
    stringsAsFactors = FALSE
  )
  pa <- suppressWarnings(suppressMessages(ProteinAnnotation$new(
    lfqdata, ra, description = "description", cleaned_ids = "cleanID",
    full_id = "fasta.id", exp_nr_children = "nr_peptides",
    pattern_contaminants = "^zz", pattern_decoys = "^REV"
  )))
  expect_equal(sum(pa$row_annot$protein_Id == pick), 1L)
  expect_true(grepl("^sp\\|",
    pa$row_annot$fasta.id[pa$row_annot$protein_Id == pick]))
})

test_that("annotate_contaminants with empty pattern does NOT flag all (no match-all)", {
  res <- suppressWarnings(suppressMessages(
    sim_data_protAnnot(Nprot = 40, PROTEIN = TRUE)
  ))
  lfqdata <- res$lfqdata
  ids <- unique(lfqdata$data_long()$protein_Id)
  ra <- data.frame(
    protein_Id = ids, cleanID = ids, description = "d",
    fasta.id = ids, nr_peptides = 2, stringsAsFactors = FALSE
  )
  # pattern_contaminants = "" must fall back to defaults, never grepl("", x)
  pa <- suppressWarnings(suppressMessages(ProteinAnnotation$new(
    lfqdata, ra, description = "description", cleaned_ids = "cleanID",
    full_id = "fasta.id", exp_nr_children = "nr_peptides",
    pattern_contaminants = "", pattern_decoys = "^REV"
  )))
  n_con <- pa$annotate_contaminants()
  expect_lt(n_con, nrow(pa$row_annot)) # NOT everything flagged
  expect_gt(n_con, 0) # zz contaminants still flagged (defaults)
  expect_true(all(prolfqua::is_contaminant(pa$row_annot$fasta.id[pa$row_annot$CON])))
})
