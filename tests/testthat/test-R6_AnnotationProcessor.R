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
