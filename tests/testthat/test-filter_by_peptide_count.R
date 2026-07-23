# Tests for the reader-local minimum-peptides-per-protein filter and its wiring.

test_that("filter_by_peptide_count drops proteins below the threshold", {
  d <- data.frame(
    prot = c("A", "A", "A", "B", "C", "C", "C"),
    pep = c("p1", "p2", "p2", "p1", "p1", "p2", "p3"),
    x = 1:7,
    stringsAsFactors = FALSE
  )
  # distinct peptides per protein: A = 2, B = 1, C = 3
  expect_equal(nrow(filter_by_peptide_count(d, "prot", "pep", 1)), 7) # no-op
  expect_equal(nrow(filter_by_peptide_count(d, "prot", "pep", 0)), 7) # no-op

  r2 <- filter_by_peptide_count(d, "prot", "pep", 2)
  expect_setequal(unique(r2$prot), c("A", "C")) # B dropped

  r3 <- filter_by_peptide_count(d, "prot", "pep", 3)
  expect_setequal(unique(r3$prot), "C") # only C survives
})

test_that("filter_by_peptide_count counts DISTINCT peptides (repeats don't inflate)", {
  d <- data.frame(
    prot = c("A", "A", "A", "A"),
    pep = c("p1", "p1", "p1", "p1"), # 1 distinct peptide across 4 rows
    stringsAsFactors = FALSE
  )
  expect_equal(nrow(filter_by_peptide_count(d, "prot", "pep", 2)), 0)
  expect_equal(nrow(filter_by_peptide_count(d, "prot", "pep", 1)), 4)
})

test_that("filter_by_peptide_count is a no-op for NULL nr_peptides", {
  d <- data.frame(prot = c("A", "B"), pep = c("p1", "p1"))
  expect_equal(nrow(filter_by_peptide_count(d, "prot", "pep", NULL)), 2)
})

test_that(".validate_nr_peptides accepts whole numbers and coerces cleanly", {
  expect_identical(prolfquapp:::.validate_nr_peptides(NULL), 1L)
  expect_identical(prolfquapp:::.validate_nr_peptides(1), 1L)
  expect_identical(prolfquapp:::.validate_nr_peptides(3), 3L)
  expect_identical(prolfquapp:::.validate_nr_peptides(2.0), 2L)
  expect_identical(prolfquapp:::.validate_nr_peptides("2"), 2L)
})

test_that(".validate_nr_peptides rejects invalid thresholds", {
  expect_error(prolfquapp:::.validate_nr_peptides(0))
  expect_error(prolfquapp:::.validate_nr_peptides(-1))
  expect_error(prolfquapp:::.validate_nr_peptides(1.5))
  expect_error(prolfquapp:::.validate_nr_peptides(NA))
  expect_error(prolfquapp:::.validate_nr_peptides("abc"))
  expect_error(prolfquapp:::.validate_nr_peptides(c(1, 2)))
  # non-finite / out-of-range must be rejected (else coercion -> NA -> crash)
  expect_error(prolfquapp:::.validate_nr_peptides(Inf))
  expect_error(prolfquapp:::.validate_nr_peptides(-Inf))
  expect_error(prolfquapp:::.validate_nr_peptides(NaN))
  expect_error(prolfquapp:::.validate_nr_peptides(1e30))
})

test_that(".forward_nr_peptides adds the arg only for readers that declare it", {
  supports <- function(quant_data, nr_peptides = 1) NULL
  no_support <- function(quant_data, pattern_decoys = NULL) NULL

  # supported -> arg added
  a <- prolfquapp:::.forward_nr_peptides(list(quant_data = 1), supports, 3)
  expect_equal(a$nr_peptides, 3)

  # unsupported + threshold 1 -> not added, no warning
  b <- prolfquapp:::.forward_nr_peptides(list(quant_data = 1), no_support, 1)
  expect_false("nr_peptides" %in% names(b))

  # unsupported + threshold > 1 -> not added, warns (explicit skip)
  expect_warning(
    d <- prolfquapp:::.forward_nr_peptides(list(quant_data = 1), no_support, 3),
    "does not support nr_peptides"
  )
  expect_false("nr_peptides" %in% names(d))
})

test_that("every reader declares nr_peptides (MZMine/dummy accept + ignore)", {
  readers <- c(
    "preprocess_DIANN",
    "preprocess_MQ_peptide",
    "preprocess_FP_PSM",
    "preprocess_MSstats",
    "preprocess_MSstats_FPDIA",
    "preprocess_BGS",
    "preprocess_SIM",
    "preprocess_mzMine",
    "preprocess_dummy"
  )
  for (r in readers) {
    fn <- utils::getFromNamespace(r, "prolfquapp")
    expect_true("nr_peptides" %in% names(formals(fn)), info = r)
  }
})

test_that("ProteinAnnotation follows the filtered quant (annotation is not a filter)", {
  istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 30)
  long <- istar$data

  filtered <- filter_by_peptide_count(long, "protein_Id", "peptide_Id", 3)
  lfq <- prolfqua::LFQData$new(filtered, istar$config)

  addannot <- long |>
    dplyr::distinct(dplyr::across(dplyr::all_of("protein_Id"))) |>
    dplyr::mutate(description = "d", nr_peptides = 1)

  pannot <- prolfquapp::ProteinAnnotation$new(
    lfq,
    addannot,
    description = "description"
  )

  # annotation holds exactly the proteins that survived the quant filter
  expect_setequal(
    pannot$row_annot[[pannot$pID]],
    unique(filtered$protein_Id)
  )
  # and every surviving protein has >= 3 distinct peptides
  counts <- filtered |>
    dplyr::distinct(protein_Id, peptide_Id) |>
    dplyr::count(protein_Id)
  expect_true(all(counts$n >= 3))
})

test_that("preprocess_software forwards nr_peptides to a supporting reader (DUMMY)", {
  annot <- data.frame(
    file = c("a1.raw", "a2.raw", "a3.raw", "a4.raw"),
    name = c("aa", "ba", "aa", "ba"),
    group = c("a", "a", "b", "b")
  )
  annot <- prolfquapp::read_annotation(annot, QC = TRUE)
  funcs <- prolfquapp::prolfqua_preprocess_functions[["DUMMY"]]

  # DUMMY accepts + ignores nr_peptides: no "unused argument" error, no warning
  res <- prolfquapp::preprocess_software(".", annot, funcs, nr_peptides = 3)
  expect_false(is.null(res$xd))
})

test_that("run_dea applies nr_peptides via the SIM reader and keeps xd consistent", {
  skip_on_cran()
  dataset <- system.file(
    "application/sim_test/dataset_sim.csv",
    package = "prolfquapp"
  )
  skip_if(nchar(dataset) == 0, "sim_test fixture not installed")

  GRP2 <- prolfquapp::make_DEA_config_R6(
    WORKUNITID = "TEST_NRPEP",
    Normalization = "robscale"
  )
  GRP2$processing_options$nr_peptides <- 2L

  result <- prolfquapp::run_dea(
    indir = tempdir(),
    dataset = dataset,
    software = "prolfquapp.SIM",
    config = GRP2
  )

  lfq <- result$xd$lfqdata
  hk <- lfq$hierarchy_keys()
  counts <- lfq$data_long() |>
    dplyr::distinct(dplyr::across(dplyr::all_of(hk[1:2]))) |>
    dplyr::count(dplyr::across(dplyr::all_of(hk[1])), name = "n")

  # every surviving protein meets the >= 2 distinct-peptides threshold
  expect_true(all(counts$n >= 2))

  # finding #1: protein_annotation follows the filtered quant, no drift
  pid <- result$xd$protein_annotation$pID
  expect_setequal(
    unique(result$xd$protein_annotation$row_annot[[pid]]),
    unique(lfq$data_long()[[hk[1]]])
  )
})

test_that("run_make_yaml writes the chosen nr_peptides (validated)", {
  cfg <- prolfquapp::run_make_yaml(workunit = "WU1", nr_peptides = 3)
  expect_equal(cfg$processing_options$nr_peptides, 3L)
})

test_that("sync_opt_config applies and validates a CLI nr_peptides override", {
  res <- prolfquapp::sync_opt_config(
    list(nr_peptides = 2),
    prolfquapp::make_DEA_config_R6()
  )
  expect_equal(res$config$processing_options$nr_peptides, 2L)

  expect_error(
    prolfquapp::sync_opt_config(
      list(nr_peptides = 0),
      prolfquapp::make_DEA_config_R6()
    )
  )
})
