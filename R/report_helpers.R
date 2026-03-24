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
  orderURL <- NULL
  workunitURL <- NULL
  projectURL <- NULL
  orderID <- as.numeric(project_spec$order_Id)
  if ((length(orderID) > 0) && !is.na(orderID)) {
    orderURL <- paste0(
      "https://fgcz-bfabric.uzh.ch/bfabric/order/show.html?id=",
      orderID,
      "&tab=details"
    )
  }
  workunitID <- as.numeric(project_spec$workunit_Id)
  if ((length(workunitID) > 0) && !is.na(workunitID)) {
    workunitURL <- paste0(
      "https://fgcz-bfabric.uzh.ch/bfabric/workunit/show.html?id=",
      orderID,
      "&tab=details"
    )
  }
  projectID <- as.numeric(project_spec$project_Id)
  if ((length(projectID) > 0) && !is.na(projectID)) {
    projectURL <- paste0(
      "https://fgcz-bfabric.uzh.ch/bfabric/project/show.html?id=",
      orderID,
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

#' write index.html file with links to all relevant files:
#' @param file_path_list named list of output file paths
#' @param result_dir directory for the index.html output
#' @export
#' @examples
#' .resdir <- "."
#' write_index_html(prolfquapp:::.test_links,tempdir())
write_index_html <- function(file_path_list, result_dir) {
  # Determine top-level directory and name
  dea_path <- file_path_list$dea_file
  topdir_path <- dirname(dea_path)
  topdir_name <- basename(topdir_path)

  # Helper to make a section
  make_section <- function(title, paths) {
    lines <- c(
      sprintf("  <h2>%s</h2>", title),
      "  <ul>"
    )
    for (i in seq_along(paths)) {
      p <- paths[i]
      name <- names(paths)[i]
      text <- if (!is.null(name) && nzchar(name)) name else basename(p)
      rel <- gsub(result_dir, "", p)
      rel <- paste0("./", rel)
      lines <- c(
        lines,
        sprintf('    <li><a href="%s">%s</a></li>', rel, text)
      )
    }
    c(lines, "  </ul>")
  }

  style_block <- R"(
  <style>
    .footer {
      text-align: center;
      margin-top: 40px;
      padding-top: 20px;
      border-top: 1px solid #ecf0f1;
      color: #7f8c8d;
    }
    h2 {
      text-align: left;
      margin-top: 20px;
      padding-bottom: 10px;
      border-bottom: 1px solid #ecf0f1;
      color: #2c3e50;
    }

    h3 {
      text-align: center;
      margin-top: 10px;
      color: #7f8c8d;
    }
    body {
      font-family: Arial, sans-serif;
    }

    h1 {
      text-align: center;
    }
  </style>
  )"

  # Build html
  title <- paste0("Differential Expression Analysis (DEA)")
  html_lines <- c(
    "<!DOCTYPE html>",
    "<html>",
    "<head>",
    "  <meta charset='UTF-8'>",
    sprintf("  <title>%s results</title>", title),

    style_block,
    "</head>",
    "<body>",
    sprintf("  <h1>%s results</h1>", title),
    sprintf(" <h3>for: %s</h3>", topdir_name),
    "<h3>Functional Genomics Center Zurich</h3>"
  )

  # 1. Top-level files: dea + qc
  html_lines <- c(
    html_lines,
    make_section(
      "Overview files",
      c(
        "<b>DEA Report (read first)</b>" = file_path_list$dea_file,
        "QC report" = file_path_list$qc_file
      )
    ),
    make_section(
      "Spreadsheet exports",
      c(
        "DEA results (XLSX)" = file_path_list$data_files$xlsx_file,
        "Intensity Based Absolute Quantitation (XLSX)" = file_path_list$data_files$ibaq_file
      )
    ),
    make_section("ORA inputs:", unlist(file_path_list$data_files$ora_files)),
    make_section(
      "GSEA ranklists:",
      unlist(file_path_list$data_files$gsea_files)
    )
  )

  analysis_date <- format(Sys.time(), "%B %d, %Y")
  version <- as.character(packageVersion("prolfquapp"))

  # Close out
  html_lines <- c(
    html_lines,
    "  <div class='footer'>",
    sprintf("    <p><strong>Analysis Date:</strong> %s</p>", analysis_date),
    "    <p><strong>Project:</strong> Differential Expression analysis</p>",
    sprintf(
      "    <p><strong>Generated by:</strong> prolfquapp package v%s</p>",
      version
    ),
    "    <p><strong>Publications (please cite):</strong></p>",
    paste0(
      "    <p><a href=\"https://doi.org/10.1021/acs.jproteome.2c00441\"",
      " target=\"_blank\">prolfqua - Wolski et al.,",
      " J Proteome Res. 2023;22(4):1092-1104</a></p>"
    ),
    paste0(
      "    <p><a href=\"https://doi.org/10.1021/acs.jproteome.4c00911\"",
      " target=\"_blank\">prolfquapp - Wolski et al.,",
      " J Proteome Res. 2025;24(2):955-965</a></p>"
    ),
    "  </div>",
    "</body>",
    "</html>"
  )
  # Write it
  index_file <- file.path(result_dir, "index.html")
  writeLines(html_lines, con = index_file)
  message("Wrote HTML index to: ", index_file)
  invisible(html_lines)
}
