annotation_fixture <- function(names) {
  data.frame(
    raw.file = paste0("file_", seq_along(names), ".raw"),
    Name = names,
    group = c(rep("control", length(names) - 1L), "treat"),
    CONTROL = c(rep("C", length(names) - 1L), "T"),
    Subject = paste0("subject_", seq_along(names)),
    stringsAsFactors = FALSE
  )
}

test_that("read_annotation keeps short unique sample names unchanged", {
  annotation <- annotation_fixture(c("control_1", "control_2", "treat_1"))

  result <- prolfquapp::read_annotation(annotation, repeated = FALSE)

  expect_equal(result$atable$sample_name, "Name")
  expect_disjoint(colnames(result$annot), "sampleName")
  expect_equal(result$annot$Name, annotation$Name)
})

test_that("read_annotation derives suffix display names for long sample names", {
  annotation <- annotation_fixture(c(
    "very_long_control_sample_001",
    "very_long_control_sample_002",
    "very_long_treatment_sample_001"
  ))
  expected_names <- substring(
    annotation$Name,
    nchar(annotation$Name) - 14L + 1L,
    nchar(annotation$Name)
  )

  result <- prolfquapp::read_annotation(annotation, repeated = FALSE)

  expect_equal(result$atable$sample_name, "sampleName")
  expect_contains(colnames(result$annot), "sampleName")
  expect_equal(result$annot$Name, annotation$Name)
  expect_equal(result$annot$sampleName, expected_names)
  expect_lte(max(nchar(result$annot$sampleName)), 14L)
})

test_that("read_annotation makes duplicate derived suffix display names unique", {
  annotation <- annotation_fixture(c(
    "alpha_prefix_SHARED_SUFFIX",
    "beta_prefix_SHARED_SUFFIX",
    "gamma_prefix_OTHER_SUFFIX"
  ))

  result <- prolfquapp::read_annotation(annotation, repeated = FALSE)

  expect_equal(result$atable$sample_name, "sampleName")
  expect_equal(result$annot$Name, annotation$Name)
  expect_length(unique(result$annot$sampleName), nrow(annotation))
  expect_match(result$annot$sampleName[2], "_1$")
})

test_that("read_annotation does not overwrite duplicate short sample names", {
  annotation <- annotation_fixture(c("control", "control", "treat"))

  result <- prolfquapp::read_annotation(annotation, repeated = FALSE)

  expect_equal(result$atable$sample_name, "sampleName")
  expect_equal(result$annot$Name, annotation$Name)
  expect_equal(result$annot$sampleName, c("control", "control_1", "treat"))
})

test_that("empty bait column does not hijack a populated grouping variable", {
  # Mirrors A386 QC datasets: a populated "Grouping Var" alongside an empty
  # "Bait ID". The grouping pattern matches both; the empty bait column must
  # not win, otherwise the grouping factor is all-NA and the QC missingness
  # heatmap crashes in pheatmap (gpar fill length 0).
  annotation <- data.frame(
    raw.file = paste0("file_", 1:3, ".raw"),
    Name = c("s1", "s2", "s3"),
    "Grouping Var" = c("A", "A", "B"),
    "Bait ID" = c("", "", ""),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  result <- prolfquapp::read_annotation(annotation, repeated = FALSE, QC = TRUE)

  grouping_col <- result$atable$factors[[result$atable$factor_keys()[1]]]
  expect_equal(grouping_col, "Grouping Var")
  expect_true(all(!is.na(result$annot[[grouping_col]])))
  expect_setequal(unique(result$annot[[grouping_col]]), c("A", "B"))
})

test_that("a fully empty grouping/bait column collapses to a single NA group", {
  # Both candidate columns are empty: nothing to prefer, so instead of an
  # all-NA factor (which crashes the QC heatmap) every sample lands in one
  # group named "NA".
  annotation <- data.frame(
    raw.file = paste0("file_", 1:3, ".raw"),
    Name = c("s1", "s2", "s3"),
    "Grouping Var" = c("", "", ""),
    "Bait ID" = c(NA, NA, NA),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  result <- prolfquapp::read_annotation(annotation, repeated = FALSE, QC = TRUE)

  grouping_col <- result$atable$factors[[result$atable$factor_keys()[1]]]
  expect_true(all(!is.na(result$annot[[grouping_col]])))
  expect_equal(unique(result$annot[[grouping_col]]), "NA")
})

test_that("partially empty grouping column keeps real groups plus an NA group", {
  annotation <- data.frame(
    raw.file = paste0("file_", 1:4, ".raw"),
    Name = c("s1", "s2", "s3", "s4"),
    "Grouping Var" = c("A", "A", "", NA),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  result <- prolfquapp::read_annotation(annotation, repeated = FALSE, QC = TRUE)

  grouping_col <- result$atable$factors[[result$atable$factor_keys()[1]]]
  expect_equal(grouping_col, "Grouping Var")
  expect_true(all(!is.na(result$annot[[grouping_col]])))
  expect_setequal(unique(result$annot[[grouping_col]]), c("A", "NA"))
})

test_that("QC injects a single dummy group when no grouping column exists", {
  # SUSHI/apprunner QC datasets may ship only file + name columns. QC must not
  # error (unlike DEA); set_grouping_var synthesizes a single "NA" group.
  annotation <- data.frame(
    raw.file = paste0("file_", 1:3, ".raw"),
    Name = c("s1", "s2", "s3"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  expect_warning(
    result <- prolfquapp::read_annotation(
      annotation,
      repeated = FALSE,
      QC = TRUE
    ),
    "grouping column"
  )

  grouping_col <- result$atable$factors[[result$atable$factor_keys()[1]]]
  expect_equal(grouping_col, "group")
  expect_true(all(!is.na(result$annot[[grouping_col]])))
  expect_equal(unique(result$annot[[grouping_col]]), "NA")
})

test_that("non-QC errors naming the grouping pattern when grouping is missing", {
  # DEA still requires a grouping column; the message must reference the
  # grouping pattern (^group...), not the sample-name pattern (^name).
  annotation <- data.frame(
    raw.file = paste0("file_", 1:3, ".raw"),
    Name = c("s1", "s2", "s3"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  expect_error(
    prolfquapp::read_annotation(annotation, repeated = FALSE, QC = FALSE),
    "\\^group"
  )
})
