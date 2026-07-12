test_that("report provenance supplies stable compact fields", {
  provenance <- prolfquapp:::.report_provenance(
    project_spec = list(
      workunit_Id = "348267",
      order_Id = "23078",
      project_Id = "12345",
      project_name = "Example project",
      input_URL = "https://example.org/input"
    ),
    software = "DIA-NN",
    model = "lm_impute"
  )

  expect_equal(provenance$workunit_id, "348267")
  expect_equal(provenance$order_id, "23078")
  expect_equal(provenance$input_data, "https://example.org/input")
  expect_equal(provenance$software, "DIA-NN")
  expect_equal(provenance$model, "lm_impute")

  table <- prolfquapp:::.report_provenance_table(provenance)
  expect_named(table, c("Field", "Value"))
  expect_equal(table$Field[[1]], "Workunit ID")
  expect_equal(table$Value[[7]], "https://example.org/input")
})
