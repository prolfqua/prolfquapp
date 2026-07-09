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

# Colour palette keyed on the estimate_type column. Observed rows are black;
# rescued rows (LOD imputation or missing-group fallback) get a distinct colour
# so they stand out in tables/plots. Unknown levels fall back to an auto palette.
.estimate_type_palette <- function(estimate_types) {
  default_palette <- c(
    observed = "black",
    lod_imputed = "orange",
    missing_fallback = "gray60"
  )
  levels_ <- sort(unique(as.character(stats::na.omit(estimate_types))))
  palette <- default_palette[levels_]
  missing_palette <- is.na(palette)
  if (any(missing_palette)) {
    palette[missing_palette] <- grDevices::hcl.colors(
      sum(missing_palette),
      palette = "Dark 3"
    )
  }
  names(palette) <- levels_
  palette
}

#' write index.html file with links to all relevant files:
#' @param file_path_list named list of output file paths
#' @param result_dir directory for the index.html output
#' @export
#' @examples
#' .resdir <- "."
#' write_index_html(prolfquapp:::.test_links,tempdir())
write_index_html <- function(file_path_list, result_dir) {
  # Determine the top-level directory/name. Prefer the primary DEA report, then
  # any other available report, and fall back to result_dir, so the index still
  # builds when a report is missing.
  candidate_paths <- c(
    file_path_list$dea_file,
    file_path_list$qc_file,
    file_path_list$quarto_file,
    file_path_list$sse_file
  )
  candidate_paths <- candidate_paths[!is.na(candidate_paths) & nzchar(candidate_paths)]
  topdir_path <- if (length(candidate_paths) > 0) {
    dirname(candidate_paths[1])
  } else {
    result_dir
  }
  topdir_name <- basename(topdir_path)

  # Build link entries from a (possibly named) vector of paths, e.g. the ORA /
  # GSEA file sets. Empty/NA paths are dropped; the shared `desc` (if any) is
  # attached to every entry.
  paths_to_entries <- function(paths, desc = NULL) {
    paths <- unlist(paths, use.names = TRUE)
    paths <- paths[!is.na(paths) & nzchar(paths)]
    lapply(seq_along(paths), function(i) {
      nm <- names(paths)[i]
      list(
        label = if (!is.null(nm) && nzchar(nm)) nm else basename(paths[[i]]),
        path = unname(paths[[i]]),
        desc = desc
      )
    })
  }

  # Render a titled section. `entries` is a list of list(label=, path=, desc=);
  # entries with an empty path are dropped and the section is omitted if none
  # remain. `intro` is optional explanatory text shown above the links.
  make_section <- function(title, entries, intro = NULL) {
    entries <- Filter(
      function(e) length(e$path) == 1 && !is.na(e$path) && nzchar(e$path),
      entries
    )
    if (length(entries) == 0) {
      return(character())
    }
    lines <- sprintf("  <h2>%s</h2>", title)
    if (!is.null(intro) && nzchar(intro)) {
      lines <- c(lines, sprintf("  <p class='intro'>%s</p>", intro))
    }
    lines <- c(lines, "  <ul>")
    for (e in entries) {
      rel <- .index_relative_href(e$path, result_dir)
      label <- if (!is.null(e$label) && nzchar(e$label)) e$label else basename(e$path)
      desc <- if (!is.null(e$desc) && nzchar(e$desc)) {
        sprintf("<br><span class='desc'>%s</span>", e$desc)
      } else {
        ""
      }
      lines <- c(lines, sprintf('    <li><a href="%s">%s</a>%s</li>', rel, label, desc))
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
    .desc {
      color: #7f8c8d;
      font-size: 0.9em;
    }
    .intro {
      color: #34495e;
      margin: 6px 0 10px;
      max-width: 70em;
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

  # Report links, each with a one-line description of the document.
  html_lines <- c(
    html_lines,
    make_section(
      "Reports",
      list(
        list(
          label = "<b>DEA report (read first)</b>",
          path = file_path_list$dea_file,
          desc = paste(
            "Main differential-expression report: analysis settings,",
            "quality-control plots, volcano plots and the per-contrast results",
            "tables."
          )
        ),
        list(
          label = "Differential-expression QC report",
          path = file_path_list$qc_file,
          desc = paste(
            "Model diagnostics: missing values, protein variance, and the",
            "fold-change / p-value distributions and MA plots used to judge the",
            "model fit."
          )
        ),
        list(
          label = "Sample-size estimation report",
          path = file_path_list$sse_file,
          desc = paste(
            "Feature-level variability and a two-sample t-test sample-size /",
            "power estimate to help plan follow-up experiments."
          )
        ),
        list(
          label = "Overview report (SummarizedExperiment tabset)",
          path = file_path_list$quarto_file,
          desc = paste(
            "Tabbed overview built from the SummarizedExperiment: settings,",
            "feature detection, quality control, differential abundance and",
            "result tables."
          )
        )
      )
    ),
    make_section(
      "Spreadsheets",
      list(
        list(
          label = "Differential-expression results (XLSX)",
          path = file_path_list$data_files$xlsx_file,
          desc = paste(
            "All contrasts with estimates, log2 fold-changes, p-values and FDR",
            "(one sheet per contrast)."
          )
        ),
        list(
          label = "Protein abundances / iBAQ (XLSX)",
          path = file_path_list$data_files$ibaq_file,
          desc = paste(
            "Intensity-based absolute quantification (iBAQ) per protein and",
            "sample."
          )
        )
      )
    ),
    make_section(
      "Over-representation analysis (ORA) input",
      paths_to_entries(file_path_list$data_files$ora_files),
      intro = paste(
        "Gene lists of the significantly regulated proteins, one file per",
        "contrast and direction (up / down). Paste a list into an",
        "over-representation tool such as STRING, g:Profiler or WebGestalt to",
        "find which pathways or GO terms are enriched among the regulated",
        "proteins."
      )
    ),
    make_section(
      "GSEA rank files",
      paths_to_entries(file_path_list$data_files$gsea_files),
      intro = paste(
        "Pre-ranked protein lists (.rnk), one per contrast: every quantified",
        "protein ranked by its signed differential-expression statistic. Use",
        "these with pre-ranked gene-set enrichment analysis (e.g. GSEA Preranked",
        "or fgsea), which detects gene sets shifted toward the top or bottom of",
        "the ranking."
      )
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
