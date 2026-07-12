# Internal helpers for compact report provenance in the Quarto reports.

.report_provenance_scalar <- function(x, fallback = "n/a") {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) {
    return(fallback)
  }

  value <- as.character(x[[1]])
  if (is.na(value) || !nzchar(value)) {
    fallback
  } else {
    value
  }
}

.report_provenance_field <- function(x, name) {
  if (is.null(x)) {
    return(NULL)
  }

  if (is.environment(x)) {
    return(tryCatch(x[[name]], error = function(e) NULL))
  }

  if (is.list(x)) {
    return(x[[name]])
  }

  NULL
}

.report_provenance_creator <- function() {
  candidates <- c(
    Sys.getenv("BFABRIC_USER", unset = ""),
    Sys.getenv("SUSHI_USER", unset = ""),
    Sys.getenv("USER", unset = ""),
    Sys.getenv("USERNAME", unset = ""),
    Sys.info()[["user"]]
  )
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  .report_provenance_scalar(candidates)
}

.report_provenance <- function(
  project_spec = NULL,
  input_data = NULL,
  software = NULL,
  model = NULL
) {
  if (
    is.null(input_data) || length(input_data) == 0 || all(is.na(input_data))
  ) {
    input_data <- .report_provenance_field(project_spec, "input_URL")
  }
  if (is.null(software) || length(software) == 0 || all(is.na(software))) {
    software <- .report_provenance_field(project_spec, "software")
  }
  if (is.null(model) || length(model) == 0 || all(is.na(model))) {
    model <- .report_provenance_field(project_spec, "model")
  }

  list(
    workunit_id = .report_provenance_scalar(.report_provenance_field(
      project_spec,
      "workunit_Id"
    )),
    order_id = .report_provenance_scalar(.report_provenance_field(
      project_spec,
      "order_Id"
    )),
    project_id = .report_provenance_scalar(.report_provenance_field(
      project_spec,
      "project_Id"
    )),
    project_name = .report_provenance_scalar(.report_provenance_field(
      project_spec,
      "project_name"
    )),
    creator = .report_provenance_creator(),
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    input_data = .report_provenance_scalar(input_data),
    software = .report_provenance_scalar(software),
    model = .report_provenance_scalar(model),
    prolfquapp_version = as.character(utils::packageVersion("prolfquapp"))
  )
}

.report_provenance_table <- function(provenance) {
  fields <- c(
    "Workunit ID",
    "Order ID",
    "Project ID",
    "Project name",
    "Creator",
    "Created at",
    "Input data",
    "Quantification software",
    "Model",
    "prolfquapp version"
  )
  values <- c(
    provenance$workunit_id,
    provenance$order_id,
    provenance$project_id,
    provenance$project_name,
    provenance$creator,
    provenance$created_at,
    provenance$input_data,
    provenance$software,
    provenance$model,
    provenance$prolfquapp_version
  )

  data.frame(Field = fields, Value = values, stringsAsFactors = FALSE)
}
