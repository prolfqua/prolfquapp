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
  rownames(res) <- .data[, var] |>
    tidyr::unite("id", tidyselect::everything(), sep = sep) |>
    dplyr::pull("id")
  return(res)
}

#' Strip pattern from row names of a matrix or data.frame
#'
#' @param .data a matrix or data.frame with row names
#' @param strip regex pattern to remove from row names
#' @return the input with cleaned row names
#' @export
strip_rownames <- function(.data, strip = "~lfq~light$") {
  newrnames <- gsub(strip, "", rownames(.data))
  rownames(.data) <- newrnames
  return(.data)
}

#' Enrich a quant/result table with the protein annotation
#'
#' Right-joins so every row of \code{x} (the quant or result table) is preserved
#' and never multiplied.
#'
#' Joins on the **hierarchy keys** shared by both frames — the config's
#' feature-identity columns (e.g. \code{protein_Id}; plus deeper keys like
#' \code{site} for a PTM \code{protein_Id} + \code{site} hierarchy) — never on a
#' coincidentally-shared value column. \code{hierarchy_keys} is intersected with
#' the columns actually present in both frames, so a protein-level annotation
#' joined to protein-level contrasts uses \code{protein_Id}, while a site-level
#' analysis uses \code{protein_Id} + \code{site}. Joining on only \code{protein_Id}
#' when both carry \code{site} would suffix it to \code{site.x}/\code{site.y} and
#' drop the bare \code{site} the report needs (`unite(all_of(hierarchy_keys))`).
#' Keeps the row-preserving \code{right_join} and a uniqueness guard on the
#' resolved key(s).
#' @param annotation protein annotation data frame
#' @param x quant/result table to annotate
#' @param hierarchy_keys config hierarchy keys (e.g. \code{lfqdata$hierarchy_keys()});
#'   the join uses the subset present in both frames
#' @return \code{x} enriched with annotation columns; one row per row of \code{x}
#' @keywords internal
#' @noRd
.join_annotation <- function(annotation, x, hierarchy_keys) {
  join_keys <- intersect(
    hierarchy_keys,
    intersect(colnames(annotation), colnames(x))
  )
  if (length(join_keys) == 0) {
    stop("internal: no shared hierarchy key to join annotation on.")
  }
  if (anyDuplicated(annotation[, join_keys, drop = FALSE]) > 0) {
    stop(
      "internal: protein annotation is not unique on the hierarchy key(s) '",
      paste(join_keys, collapse = "', '"),
      "' before the annotation join."
    )
  }
  dplyr::right_join(annotation, x, by = join_keys, multiple = "all")
}


#' build bfabric urls
#' @param project_spec ProjectSpec R6 object with project, order, workunit IDs
#' @export
#' @examples
#'
#'
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
  as_bfabric_id <- function(x) {
    suppressWarnings(as.numeric(x))
  }
  orderURL <- NULL
  workunitURL <- NULL
  projectURL <- NULL
  orderID <- as_bfabric_id(project_spec$order_Id)
  if ((length(orderID) > 0) && !is.na(orderID)) {
    orderURL <- paste0(
      "https://fgcz-bfabric.uzh.ch/bfabric/order/show.html?id=",
      orderID,
      "&tab=details"
    )
  }
  workunitID <- as_bfabric_id(project_spec$workunit_Id)
  if ((length(workunitID) > 0) && !is.na(workunitID)) {
    workunitURL <- paste0(
      "https://fgcz-bfabric.uzh.ch/bfabric/workunit/show.html?id=",
      workunitID,
      "&tab=details"
    )
  }
  projectID <- as_bfabric_id(project_spec$project_Id)
  if ((length(projectID) > 0) && !is.na(projectID)) {
    projectURL <- paste0(
      "https://fgcz-bfabric.uzh.ch/bfabric/project/show.html?id=",
      projectID,
      "&tab=details"
    )
  }

  return(list(
    orderURL = orderURL,
    projectURL = projectURL,
    workunitURL = workunitURL
  ))
}


.base_dir <- paste0(
  "./DEA_20250704_PI35298_O38953_WUtotal_proteome_none/",
  "Results_WU_total_proteome"
)
.test_links <- list(
  dea_file = file.path(.base_dir, "DE_WUtotal_proteome.html"),
  qc_file = file.path(.base_dir, "QC_WUtotal_proteome.html"),
  data_files = list(
    xlsx_file = file.path(.base_dir, "DE_WUtotal_proteome.xlsx"),
    ora_files = list(
      ORA_Treated_vs_Control_down_WUtotal_proteome.txt = file.path(
        .base_dir,
        "ORA_Treated_vs_Control_down_WUtotal_proteome.txt"
      ),
      ORA_Treated_vs_Control_up_WUtotal_proteome.txt = file.path(
        .base_dir,
        "ORA_Treated_vs_Control_up_WUtotal_proteome.txt"
      )
    ),
    gsea_files = list(
      `GSEA_Treated_vs_Control_WUtotal_proteome.rnk` = file.path(
        .base_dir,
        "GSEA_Treated_vs_Control_WUtotal_proteome.rnk"
      )
    ),
    ibaq_file = file.path(.base_dir, "IBAQ_total_proteome.xlsx")
  )
)

.path_to_url_path <- function(path) {
  gsub("\\\\", "/", path)
}

.encode_url_path <- function(path) {
  parts <- strsplit(path, "/", fixed = TRUE)[[1]]
  paste(utils::URLencode(parts, reserved = TRUE), collapse = "/")
}

.index_relative_href <- function(path, result_dir) {
  path_url <- .path_to_url_path(path)
  result_url <- sub("/+$", "", .path_to_url_path(result_dir))

  if (file.exists(path) && dir.exists(result_dir)) {
    path_url <- .path_to_url_path(normalizePath(path, mustWork = TRUE))
    result_url <- sub(
      "/+$",
      "",
      .path_to_url_path(normalizePath(result_dir, mustWork = TRUE))
    )
  }

  result_prefix <- paste0(result_url, "/")

  if (startsWith(tolower(path_url), tolower(result_prefix))) {
    rel <- substring(path_url, nchar(result_prefix) + 1)
  } else {
    rel <- basename(path_url)
  }

  paste0("./", .encode_url_path(rel))
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
  setwd(oldwd)

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
  candidate_paths <- c(
    file_path_list$dea_file,
    file_path_list$qc_file,
    file_path_list$quarto_file,
    file_path_list$sse_file
  )
  candidate_paths <- candidate_paths[.index_has_value(candidate_paths)]
  topdir_path <- if (length(candidate_paths) > 0) {
    dirname(candidate_paths[1])
  } else {
    result_dir
  }
  topdir_name <- basename(topdir_path)

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
    values <- data.frame(
      File = .index_file_link(entry$file, result_dir, entry$label),
      Description = entry$description,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    if (include_contents) {
      values$Contents <- entry$contents
    }
    values$Size <- .index_file_size(entry$file)
    values
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
      function(index) {
        label <- names(paths)[index]
        .index_file_link(paths[[index]], result_dir, label)
      },
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
  if (!file.exists(path)) {
    return("")
  }
  size <- file.info(path)$size
  if (is.na(size)) {
    return("")
  }
  if (size >= 1e6) {
    return(sprintf("%.1f MB", size / 1e6))
  }
  if (size >= 1e3) {
    return(sprintf("%.0f KB", size / 1e3))
  }
  sprintf("%d B", size)
}

.index_html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}
