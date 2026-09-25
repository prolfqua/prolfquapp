#' convert tibble to data.frame with rownames
#' @param .data a tibble or data.frame
#' @param var name of the column with new row.names
#' @param sep separator for uniting columns
#' @return a data.frame with rownames
#' @export
#' @examples
#' ind <- tibble::tibble(a = 1:3, rowname = letters[1:3])
#' column_to_rownames(ind)
column_to_rownames <- function(.data, var = "rowname", sep = "~lfq~") {
  res <- as.data.frame(.data)
  rownames(res) <- .data[, var, drop = FALSE] |>
    tidyr::unite("id", tidyselect::everything(), sep = sep) |>
    dplyr::pull("id")
  res
}

#' Strip pattern from row names of a matrix or data.frame
#'
#' @param .data a matrix or data.frame with row names
#' @param strip regex pattern to remove from row names
#' @return the input with cleaned row names
#' @export
strip_rownames <- function(.data, strip = "~lfq~light$") {
  rownames(.data) <- gsub(strip, "", rownames(.data))
  .data
}

# Enrich a quant/result table `x` with the protein annotation: a right join on
# the annotation's key `protein_id` (one column for protein-level analyses,
# protein_Id and site for sites), so every row of `x` is kept, never multiplied.
.join_annotation <- function(annotation, x, protein_id) {
  if (
    length(protein_id) < 1L ||
      !all(protein_id %in% colnames(annotation)) ||
      !all(protein_id %in% colnames(x))
  ) {
    stop("internal: protein annotation and result must share the key column(s): ", paste(protein_id, collapse = ", "))
  }
  if (anyDuplicated(annotation[, protein_id, drop = FALSE]) > 0) {
    stop(
      "internal: protein annotation is not unique on '",
      paste(protein_id, collapse = " + "),
      "' before the annotation join."
    )
  }
  dplyr::right_join(annotation, x, by = protein_id, multiple = "all")
}


#' build bfabric urls
#' @param project_spec ProjectSpec R6 object with project, order, workunit IDs
#' @export
#' @examples
#' ps <- ProjectSpec$new()
#' ps$project_Id <- 32258
#' ps$order_Id <- 34628
#' ps$workunit_Id <- 302212
#' bfabric_url_builder(ps)
#'
#' ps <- ProjectSpec$new()
#' ps$order_Id <- 34628
#' ps$workunit_Id <- 302212
#' bfabric_url_builder(ps)
#'
bfabric_url_builder <- function(project_spec) {
  url <- function(entity, id) {
    id <- suppressWarnings(as.numeric(id))
    if ((length(id) > 0) && !is.na(id)) {
      paste0("https://fgcz-bfabric.uzh.ch/bfabric/", entity, "/show.html?id=", id, "&tab=details")
    }
  }
  list(
    orderURL = url("order", project_spec$order_Id),
    projectURL = url("project", project_spec$project_Id),
    workunitURL = url("workunit", project_spec$workunit_Id)
  )
}

.test_links <- local({
  in_dir <- function(x) {
    file.path("./DEA_20250704_PI35298_O38953_WUtotal_proteome_none/Results_WU_total_proteome", x)
  }
  ora <- c("ORA_Treated_vs_Control_down_WUtotal_proteome.txt", "ORA_Treated_vs_Control_up_WUtotal_proteome.txt")
  gsea <- "GSEA_Treated_vs_Control_WUtotal_proteome.rnk"
  list(
    dea_file = in_dir("DE_WUtotal_proteome.html"),
    qc_file = in_dir("QC_WUtotal_proteome.html"),
    data_files = list(
      xlsx_file = in_dir("DE_WUtotal_proteome.xlsx"),
      ora_files = as.list(stats::setNames(in_dir(ora), ora)),
      gsea_files = as.list(stats::setNames(in_dir(gsea), gsea)),
      ibaq_file = in_dir("IBAQ_total_proteome.xlsx")
    )
  )
})

.path_to_url_path <- function(path) {
  gsub("\\\\", "/", path)
}

.index_relative_href <- function(path, result_dir) {
  if (file.exists(path) && dir.exists(result_dir)) {
    path <- normalizePath(path, mustWork = TRUE)
    result_dir <- normalizePath(result_dir, mustWork = TRUE)
  }
  path_url <- .path_to_url_path(path)
  result_prefix <- paste0(sub("/+$", "", .path_to_url_path(result_dir)), "/")
  rel <- if (startsWith(tolower(path_url), tolower(result_prefix))) {
    substring(path_url, nchar(result_prefix) + 1)
  } else {
    basename(path_url)
  }
  paste0("./", paste(utils::URLencode(strsplit(rel, "/", fixed = TRUE)[[1]], reserved = TRUE), collapse = "/"))
}

#' write index.html file with links to all relevant files:
#' @param file_path_list named list of output file paths
#' @param result_dir directory for the index.html output
#' @export
#' @examples
#' .resdir <- "."
#' \dontrun{
#' write_index_html(prolfquapp:::.test_links,tempdir())
#' }
write_index_html <- function(file_path_list, result_dir) {
  dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
  result_dir <- normalizePath(result_dir, mustWork = TRUE)
  index_file <- file.path(result_dir, "index.html")
  render_dir <- tempfile("prolfquapp_index_")
  dir.create(render_dir, recursive = TRUE)
  on.exit(unlink(render_dir, recursive = TRUE), add = TRUE)

  index_data_file <- file.path(render_dir, "index_data.rds")
  saveRDS(.index_deliverables(file_path_list, result_dir), index_data_file)
  index_template <- system.file("templates", "dea_index.qmd", package = "prolfquapp", mustWork = TRUE)
  file.copy(index_template, file.path(render_dir, "index.qmd"), overwrite = TRUE)

  oldwd <- setwd(render_dir)
  on.exit(setwd(oldwd), add = TRUE)
  fgczQuartoTemplate::fgcz_render(
    input = "index.qmd",
    buttons = FALSE,
    execute_params = list(index_data_file = normalizePath(index_data_file))
  )

  rendered_file <- file.path(render_dir, "index.html")
  if (!file.exists(rendered_file)) {
    stop("Quarto render did not create expected HTML file: ", rendered_file, call. = FALSE)
  }
  if (!file.copy(rendered_file, index_file, overwrite = TRUE)) {
    stop("Could not copy Quarto index to output file: ", index_file, call. = FALSE)
  }
  message("Wrote HTML index to: ", index_file)
  invisible(index_file)
}

.index_deliverables <- function(file_path_list, result_dir) {
  reports <- c(file_path_list$dea_file, file_path_list$qc_file, file_path_list$quarto_file, file_path_list$sse_file)
  reports <- reports[.index_has_value(reports)]
  topdir_name <- basename(if (length(reports) > 0) dirname(reports[1]) else result_dir)

  list(
    workunit = sub("^Results_WU_", "", topdir_name),
    result_name = topdir_name,
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    prolfquapp_version = as.character(packageVersion("prolfquapp")),
    reports = .index_entry_table(
      result_dir,
      list(
        list(
          file = file_path_list$dea_file,
          label = "DEA Report (read first)",
          description = paste(
            "Main differential-expression report with analysis settings,",
            "quality-control plots, volcano plots and per-contrast result tables."
          )
        ),
        list(
          file = file_path_list$qc_file,
          label = "Differential-expression QC report",
          description = paste(
            "Model diagnostics covering missing values, protein variance,",
            "fold-change and p-value distributions, and MA plots."
          )
        ),
        list(
          file = file_path_list$sse_file,
          label = "Sample-size estimation report",
          description = paste(
            "Feature-level variability and two-sample t-test sample-size /",
            "power estimates for follow-up experiment planning."
          )
        ),
        list(
          file = file_path_list$quarto_file,
          label = "Overview report (SummarizedExperiment tabset)",
          description = paste(
            "Tabbed overview of settings, feature detection, quality control,",
            "differential abundance and result tables."
          )
        )
      )
    ),
    spreadsheets = .index_entry_table(
      result_dir,
      list(
        list(
          file = file_path_list$data_files$xlsx_file,
          label = "Differential-expression results (XLSX)",
          description = "Workbook with the differential-expression result tables.",
          contents = paste(
            "Feature identifiers and annotation, model estimates, log2 fold",
            "changes, p-values, FDR values, significance flags and model metadata",
            "where available."
          )
        ),
        list(
          file = file_path_list$data_files$ibaq_file,
          label = "Protein abundances / iBAQ (XLSX)",
          description = "Workbook with protein abundance and iBAQ summaries.",
          contents = paste(
            "Protein-level abundance summaries per sample or group for downstream",
            "review, filtering and quality-control interpretation."
          )
        )
      ),
      include_contents = TRUE
    ),
    ora = .index_file_table(file_path_list$data_files$ora_files, result_dir),
    gsea = .index_file_table(file_path_list$data_files$gsea_files, result_dir)
  )
}

.index_has_value <- function(x) {
  if (is.null(x)) {
    return(logical(0))
  }
  !is.na(x) & nzchar(as.character(x))
}

.index_entry_table <- function(result_dir, entries, include_contents = FALSE) {
  entries <- Filter(function(entry) isTRUE(.index_has_value(entry$file)[1]), entries)
  columns <- c("File", "Description", if (include_contents) "Contents", "Size")
  if (length(entries) == 0) {
    return(.index_empty_table(columns))
  }

  rows <- lapply(entries, function(entry) {
    row <- data.frame(File = .index_file_link(entry$file, result_dir, entry$label), Description = entry$description)
    if (include_contents) {
      row$Contents <- entry$contents
    }
    row$Size <- .index_file_size(entry$file)
    row
  })
  data.frame(do.call(rbind, rows), row.names = NULL, check.names = FALSE)
}

.index_file_table <- function(paths, result_dir) {
  paths <- unlist(paths, use.names = TRUE)
  paths <- paths[.index_has_value(paths)]
  if (length(paths) == 0) {
    return(.index_empty_table(c("File", "Size")))
  }
  data.frame(
    File = vapply(
      seq_along(paths),
      function(i) .index_file_link(paths[[i]], result_dir, names(paths)[i]),
      character(1)
    ),
    Size = vapply(paths, .index_file_size, character(1), USE.NAMES = FALSE),
    row.names = NULL,
    check.names = FALSE
  )
}

.index_empty_table <- function(columns) {
  stats::setNames(data.frame(matrix(ncol = length(columns), nrow = 0)), columns)
}

.index_file_link <- function(path, result_dir, label = NULL) {
  label <- if (isTRUE(.index_has_value(label)[1])) label else basename(path)
  sprintf(
    "<a href='%s'>%s</a>",
    .index_relative_href(path, result_dir),
    .index_html_escape(label)
  )
}

.index_file_size <- function(path) {
  size <- file.info(path)$size
  if (is.na(size)) {
    ""
  } else if (size >= 1e6) {
    sprintf("%.1f MB", size / 1e6)
  } else if (size >= 1e3) {
    sprintf("%.0f KB", size / 1e3)
  } else {
    sprintf("%d B", size)
  }
}

.index_html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}
