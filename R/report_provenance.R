# Internal helpers for compact report provenance in the Quarto reports.

.report_provenance_scalar <- function(x, fallback = "n/a") {
  if (length(x) == 0 || all(is.na(x))) {
    return(fallback)
  }
  value <- as.character(x[[1]])
  if (is.na(value) || !nzchar(value)) fallback else value
}

.report_provenance_field <- function(x, name) {
  if (is.environment(x)) {
    tryCatch(x[[name]], error = function(e) NULL)
  } else if (is.list(x)) {
    x[[name]]
  }
}

# Explicit value, or the project_spec field when the value is missing.
.report_provenance_value <- function(value, project_spec, name) {
  if (length(value) == 0 || all(is.na(value))) .report_provenance_field(project_spec, name) else value
}

.report_overview_summary <- function(lfq_data, feature_label = "Proteins") {
  factors <- lfq_data$factors()
  factor_columns <- intersect(lfq_data$relevant_factor_keys(), colnames(factors))
  factor_columns <- grep("^control", factor_columns, ignore.case = TRUE, invert = TRUE, value = TRUE)
  group_count <- if (nrow(factors) == 0) {
    0L
  } else if (length(factor_columns) == 0) {
    1L
  } else {
    nrow(unique(factors[, factor_columns, drop = FALSE]))
  }
  data.frame(
    label = c("Samples", "Groups", feature_label),
    count = c(nrow(factors), group_count, nrow(lfq_data$hierarchy())),
    stringsAsFactors = FALSE
  )
}

.report_overview_html_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub('"', "&quot;", x, fixed = TRUE)
}

.report_overview_cards <- function(lfq_data, feature_label = "Proteins") {
  overview <- .report_overview_summary(lfq_data, feature_label)
  cards <- sprintf(
    paste0(
      '<div class="prolfquapp-overview-card">',
      '<span class="prolfquapp-overview-card-title">%s</span>',
      '<strong class="prolfquapp-overview-card-value">%s</strong>',
      "</div>"
    ),
    .report_overview_html_escape(overview$label),
    format(overview$count, big.mark = ",", scientific = FALSE, trim = TRUE)
  )

  paste0(
    "<style>\n",
    ".prolfquapp-overview-cards{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:.75rem;margin:0 0 1rem}",
    ".prolfquapp-overview-card{",
    "display:flex;align-items:baseline;justify-content:space-between;",
    "gap:.5rem;padding:.75rem 1rem;border:1px solid #d7e0e8;",
    "border-radius:.35rem;background:#f4f8fb}",
    ".prolfquapp-overview-card-title{font-size:.95rem;font-weight:600;color:#425466}",
    ".prolfquapp-overview-card-value{font-size:1.5rem;line-height:1;color:#356da3}",
    "@media (max-width:600px){.prolfquapp-overview-cards{grid-template-columns:1fr}}\n",
    "</style>\n",
    '<div class="prolfquapp-overview-cards">',
    paste(cards, collapse = ""),
    "</div>\n"
  )
}

.report_provenance <- function(project_spec = NULL, input_data = NULL, software = NULL, model = NULL) {
  users <- c(Sys.getenv(c("BFABRIC_USER", "SUSHI_USER", "USER", "USERNAME"), unset = ""), Sys.info()[["user"]])
  list(
    workunit_id = .report_provenance_scalar(.report_provenance_field(project_spec, "workunit_Id")),
    order_id = .report_provenance_scalar(.report_provenance_field(project_spec, "order_Id")),
    project_id = .report_provenance_scalar(.report_provenance_field(project_spec, "project_Id")),
    project_name = .report_provenance_scalar(.report_provenance_field(project_spec, "project_name")),
    creator = .report_provenance_scalar(users[!is.na(users) & nzchar(users)]),
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    input_data = .report_provenance_scalar(.report_provenance_value(input_data, project_spec, "input_URL")),
    software = .report_provenance_scalar(.report_provenance_value(software, project_spec, "software")),
    model = .report_provenance_scalar(.report_provenance_value(model, project_spec, "model")),
    prolfquapp_version = as.character(utils::packageVersion("prolfquapp"))
  )
}

.report_provenance_table <- function(provenance) {
  fields <- c(
    workunit_id = "Workunit ID",
    order_id = "Order ID",
    project_id = "Project ID",
    project_name = "Project name",
    creator = "Creator",
    created_at = "Created at",
    input_data = "Input data",
    software = "Quantification software",
    model = "Model",
    prolfquapp_version = "prolfquapp version"
  )
  data.frame(
    Field = unname(fields),
    Value = unlist(provenance[names(fields)], use.names = FALSE),
    stringsAsFactors = FALSE
  )
}
