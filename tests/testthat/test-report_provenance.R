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

test_that("report provenance uses project defaults only for missing overrides", {
  project_spec <- list(
    input_URL = "project input",
    software = "project software",
    model = "project model"
  )

  defaults <- prolfquapp:::.report_provenance(
    project_spec = project_spec,
    input_data = character(),
    software = NA_character_
  )
  expect_equal(defaults$input_data, "project input")
  expect_equal(defaults$software, "project software")
  expect_equal(defaults$model, "project model")

  overrides <- prolfquapp:::.report_provenance(
    project_spec = project_spec,
    input_data = "explicit input",
    software = "explicit software",
    model = "explicit model"
  )
  expect_equal(overrides$input_data, "explicit input")
  expect_equal(overrides$software, "explicit software")
  expect_equal(overrides$model, "explicit model")
})

test_that("overview cards summarize samples, groups, and quantified proteins", {
  dea <- prolfquapp::example_deanalyse(Nprot = 12)
  lfq_data <- dea$lfq_data_raw

  overview <- prolfquapp:::.report_overview_summary(lfq_data)
  factor_columns <- lfq_data$relevant_factor_keys()
  expected_group_count <- nrow(unique(
    lfq_data$factors()[, factor_columns, drop = FALSE]
  ))

  expect_equal(overview$label, c("Samples", "Groups", "Proteins"))
  expect_equal(overview$count[[1]], nrow(lfq_data$factors()))
  expect_equal(overview$count[[2]], expected_group_count)
  expect_equal(overview$count[[3]], nrow(lfq_data$hierarchy()))

  cards <- prolfquapp:::.report_overview_cards(lfq_data)
  expect_match(cards, "prolfquapp-overview-cards", fixed = TRUE)
  expect_match(cards, "Samples", fixed = TRUE)
  expect_match(cards, "Groups", fixed = TRUE)
  expect_match(cards, "Proteins", fixed = TRUE)
})
