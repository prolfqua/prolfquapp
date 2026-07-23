test_package_path <- function() {
  wd <- normalizePath(getwd(), mustWork = TRUE)
  current <- wd
  repeat {
    description <- file.path(current, "DESCRIPTION")
    if (file.exists(description)) {
      desc <- read.dcf(description)
      if (identical(unname(desc[1, "Package"]), "prolfquapp")) {
        return(current)
      }
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      break
    }
    current <- parent
  }
  package_path <- normalizePath(
    file.path(system.file(package = "prolfquapp"), ".."),
    mustWork = TRUE
  )
  if (!file.exists(file.path(package_path, "DESCRIPTION"))) {
    package_path <- normalizePath(
      system.file(package = "prolfquapp"),
      mustWork = TRUE
    )
  }
  package_path
}

# Directory holding the Quarto report sources and static visual-abstract assets.
# Prefers the source `vignettes/` (dev / load_all); falls back to the installed
# `doc/`, where `vignettes/.install_extras` ships the same sources and assets.
test_report_source_dir <- function() {
  package_path <- test_package_path()
  source_vignette_dir <- file.path(package_path, "vignettes")
  marker <- "Grp2Analysis_V2_SE_tabset.qmd"
  if (file.exists(file.path(source_vignette_dir, marker))) {
    return(normalizePath(source_vignette_dir, mustWork = TRUE))
  }
  installed_doc <- system.file("doc", package = "prolfquapp")
  if (nzchar(installed_doc) && file.exists(file.path(installed_doc, marker))) {
    return(normalizePath(installed_doc, mustWork = TRUE))
  }
  ""
}

# The current package as a loadable *source* tree, or "" if only an installed
# copy is reachable (e.g. under R CMD check). Used to render the report against
# the working-tree code via the qmd's `prolfquapp_source_path` parameter.
test_source_tree <- function() {
  package_path <- test_package_path()
  r_dir <- file.path(package_path, "R")
  if (dir.exists(r_dir) && length(list.files(r_dir, pattern = "[.]R$")) > 0) {
    return(package_path)
  }
  ""
}

test_that("se_report_lfqdata reconstructs LFQData objects from SummarizedExperiment", {
  skip_on_cran()

  dea <- prolfquapp::example_deanalyse(Nprot = 20)
  reporter <- prolfquapp::DEAReportGenerator$new(dea, dea$prolfq_app_config)
  se <- reporter$make_SummarizedExperiment()
  expect_contains(SummarizedExperiment::assayNames(se), "nr_children")

  report <- prolfquapp:::se_report_lfqdata(se)

  expect_s3_class(report$lfq_raw, "LFQData")
  expect_s3_class(report$lfq_transformed, "LFQData")
  expect_gt(nrow(report$lfq_raw$data_long()), 0)
  expect_gt(nrow(report$lfq_transformed$data_long()), 0)
  expect_gt(nrow(report$lfq_raw$data_wide(as.matrix = TRUE)$data), 0)
  expect_gt(nrow(report$lfq_transformed$data_wide(as.matrix = TRUE)$data), 0)
  expect_s3_class(report$lfq_raw$get_Plotter(), "LFQDataPlotter")
  expect_s3_class(report$lfq_raw$get_Summariser(), "LFQDataSummariser")
  expect_s3_class(report$lfq_transformed$get_Stats(), "LFQDataStats")
  expect_contains(colnames(report$contrast_table), "protein_Id")
  expect_contains(colnames(report$contrast_table), "contrast")
  expect_contains(colnames(report$feature_annotation), "nrPeptides")
  expect_gt(sum(report$feature_annotation$nrPeptides >= 2, na.rm = TRUE), 0)
  child_counts <- report$lfq_raw$get_Summariser()$hierarchy_counts_sample(
    value = "long",
    nr_children = 2
  )
  expect_gt(nrow(child_counts), 0)
})

test_that("se_report_lfqdata drops padded empty contrast rows", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      rawData = matrix(
        1:6,
        nrow = 3,
        dimnames = list(c("P1", "P2", "P3"), c("S1", "S2"))
      ),
      transformedData = matrix(
        1:6,
        nrow = 3,
        dimnames = list(c("P1", "P2", "P3"), c("S1", "S2"))
      )
    ),
    colData = S4Vectors::DataFrame(
      sampleName = c("S1", "S2"),
      group = c("A", "B")
    )
  )
  SummarizedExperiment::rowData(se)[["constrast_A"]] <- data.frame(
    protein_Id = c("P1", "P2", NA_character_),
    contrast = c("A", "A", NA_character_),
    diff = c(1, -1, NA_real_),
    FDR = c(0.01, 0.2, NA_real_),
    modelName = c("TableTest", "TableTest", NA_character_)
  )

  report <- prolfquapp:::se_report_lfqdata(se)

  expect_equal(nrow(report$contrast_table), 2)
  expect_equal(sum(is.na(report$contrast_table$contrast)), 0)
  expect_equal(sum(is.na(report$contrast_table$protein_Id)), 0)
})

test_that("se_report_lfqdata drops SAINT rows without bait estimates", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      rawData = matrix(
        1:6,
        nrow = 3,
        dimnames = list(c("P1", "P2", "P3"), c("S1", "S2"))
      ),
      transformedData = matrix(
        1:6,
        nrow = 3,
        dimnames = list(c("P1", "P2", "P3"), c("S1", "S2"))
      )
    ),
    colData = S4Vectors::DataFrame(
      sampleName = c("S1", "S2"),
      group = c("A", "B")
    ),
    metadata = list(default_model = "saint")
  )
  SummarizedExperiment::rowData(se)[["constrast_PPE4"]] <- data.frame(
    protein_Id = c("P1", "P2", "P3"),
    modelName = c("ContrastSaint", NA_character_, "ContrastSaint"),
    Bait = c("PPE4", NA_character_, "PPE4"),
    log2_EFCs = c(1, NA_real_, -1),
    SaintScore = c(0.9, NA_real_, 0.2),
    BFDR = c(0.01, NA_real_, 0.8)
  )

  report <- prolfquapp:::se_report_lfqdata(se)

  expect_equal(report$contrast_table$protein_Id, c("P1", "P3"))
  expect_equal(sum(is.na(report$contrast_table$Bait)), 0)
})

test_that("SE Quarto tabset report renders with reconstructed LFQData", {
  skip_on_cran()
  skip_if(Sys.which("quarto") == "", "Quarto CLI not installed")
  skip_if_not_installed("DT")
  skip_if_not_installed("gridExtra")
  skip_if_not_installed("plotly")

  dea <- prolfquapp::example_deanalyse(Nprot = 12)
  reporter <- prolfquapp::DEAReportGenerator$new(dea, dea$prolfq_app_config)
  se <- reporter$make_SummarizedExperiment()

  report_source_dir <- test_report_source_dir()
  skip_if(!nzchar(report_source_dir), "Quarto report sources not available")

  workdir <- file.path(tempdir(), "quarto_se_report_test")
  unlink(workdir, recursive = TRUE)
  dir.create(workdir, recursive = TRUE)
  on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

  se_file <- file.path(workdir, "SummarizedExperiment.rds")
  saveRDS(se, se_file)

  report_file <- "Grp2Analysis_V2_SE_tabset.qmd"
  file.copy(
    file.path(report_source_dir, report_file),
    workdir,
    overwrite = TRUE
  )
  file.copy(
    file.path(report_source_dir, prolfquapp:::.quarto_visual_abstract_names),
    workdir,
    overwrite = TRUE
  )

  # Stage the FGCZ assets and render via fgcz_render (the same path the runtime
  # report helper uses): fgcz_render copies _metadata.yml / fgcz.scss /
  # fgcz_header_quarto.html / fgcz-plot-finder.html from the installed
  # fgczquartotemplate package next to the qmd, so no _extensions/ tree is needed.
  execute_params <- list(se_file = normalizePath(se_file))
  source_tree <- test_source_tree()
  if (nzchar(source_tree)) {
    skip_if_not_installed("devtools")
    execute_params$prolfquapp_source_path <- source_tree
  }

  oldwd <- setwd(workdir)
  on.exit(setwd(oldwd), add = TRUE)

  render_result <- tryCatch(
    fgczquartotemplate::fgcz_render(
      input = report_file,
      buttons = FALSE,
      execute_params = execute_params
    ),
    error = function(e) e
  )
  if (inherits(render_result, "error")) {
    message("fgcz_render failed:\n", conditionMessage(render_result))
  }
  expect_false(inherits(render_result, "error"))
  html_file <- sub("[.]qmd$", ".html", report_file)
  expect_equal(file.exists(html_file), TRUE)
  html <- paste(readLines(html_file, warn = FALSE), collapse = "\n")
  expect_match(html, "fgcz-banner", fixed = TRUE)
  expect_match(html, "panel-tabset", fixed = TRUE)
  expect_match(html, "Feature Detection", fixed = TRUE)
  expect_match(html, "Differential Abundance", fixed = TRUE)
  expect_match(html, "Result Table", fixed = TRUE)
  expect_match(html, "prolfquapp-overview-cards", fixed = TRUE)
  expect_match(html, "Samples", fixed = TRUE)
  expect_match(html, "Groups", fixed = TRUE)
  expect_match(html, "Proteins", fixed = TRUE)
  expect_match(html, "at least two peptides in the experiment", fixed = TRUE)
  expect_equal(grepl("Protein Identification", html, fixed = TRUE), FALSE)
})

# The render helpers render from the installed doc/ (shipped via
# vignettes/.install_extras), which is absent under devtools::load_all / test,
# so these skip unless the package was installed with vignettes built.
test_that("internal SE Quarto report helper renders HTML", {
  skip_on_cran()
  skip_if(Sys.which("quarto") == "", "Quarto CLI not installed")
  skip_if_not_installed("DT")
  skip_if_not_installed("gridExtra")
  skip_if_not_installed("plotly")
  skip_if(
    system.file("doc/Grp2Analysis_V2_SE_tabset.qmd", package = "prolfquapp") == "",
    "Report sources not installed in doc/ (needs a vignette-built install)."
  )

  dea <- prolfquapp::example_deanalyse(Nprot = 12)
  reporter <- prolfquapp::DEAReportGenerator$new(dea, dea$prolfq_app_config)
  se <- reporter$make_SummarizedExperiment()

  workdir <- file.path(tempdir(), "quarto_se_report_helper_test")
  unlink(workdir, recursive = TRUE)
  dir.create(workdir, recursive = TRUE)
  on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

  se_file <- file.path(workdir, "SummarizedExperiment.rds")
  saveRDS(se, se_file)

  html_file <- prolfquapp:::render_quarto_se_report(
    se_file = se_file,
    output_dir = workdir,
    output_file = "helper-report.html"
  )

  expect_equal(file.exists(html_file), TRUE)
  html <- paste(readLines(html_file, warn = FALSE), collapse = "\n")
  expect_match(html, "fgcz-banner", fixed = TRUE)
  expect_match(
    html,
    "Differential Expression Analysis (Tabbed Report)",
    fixed = TRUE
  )
  expect_match(html, "Feature Detection", fixed = TRUE)
  expect_match(html, "Result Table", fixed = TRUE)
  expect_match(html, "fgcz-pf-toolbar", fixed = TRUE)
})

test_that("internal DEA Quarto report helper renders HTML", {
  skip_on_cran()
  skip_if(Sys.which("quarto") == "", "Quarto CLI not installed")
  skip_if_not_installed("DT")
  skip_if_not_installed("plotly")
  skip_if(
    system.file("doc/Grp2Analysis_V2_R6.qmd", package = "prolfquapp") == "",
    "Report sources not installed in doc/ (needs a vignette-built install)."
  )

  dea <- prolfquapp::example_deanalyse(Nprot = 12)

  workdir <- file.path(tempdir(), "quarto_dea_report_helper_test")
  unlink(workdir, recursive = TRUE)
  dir.create(workdir, recursive = TRUE)
  on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

  deanalyse_file <- file.path(workdir, "DEAnalyse.rds")
  saveRDS(dea, deanalyse_file)

  html_file <- prolfquapp:::render_quarto_dea_report(
    deanalyse_file = deanalyse_file,
    output_dir = workdir,
    output_file = "helper-dea-report.html"
  )

  expect_equal(file.exists(html_file), TRUE)
  html <- paste(readLines(html_file, warn = FALSE), collapse = "\n")
  expect_match(html, "fgcz-banner", fixed = TRUE)
  expect_match(html, "Differential Expression Analysis", fixed = TRUE)
  expect_match(html, "fgcz-pf-toolbar", fixed = TRUE)
})

test_that("DiffExpQC Quarto report declares report provenance metadata", {
  report_source_dir <- test_report_source_dir()
  skip_if(!nzchar(report_source_dir), "Quarto report sources not available")

  qmd <- paste(
    readLines(
      file.path(report_source_dir, "DiffExpQC_R6_tabset.qmd"),
      warn = FALSE
    ),
    collapse = "\n"
  )

  expect_match(qmd, "fgcz-report-metadata", fixed = TRUE)
  expect_match(qmd, "data-order-id", fixed = TRUE)
  expect_match(qmd, "data-workunit-id", fixed = TRUE)
  expect_match(qmd, "# Overview", fixed = TRUE)
  expect_match(qmd, "# Session Info", fixed = TRUE)
  expect_match(qmd, "## Report provenance", fixed = TRUE)
  expect_match(qmd, "## R session info", fixed = TRUE)
  expect_match(qmd, ".report_provenance_table", fixed = TRUE)
  expect_match(qmd, "sessionInfo()", fixed = TRUE)
})

test_that("all Quarto reports use the Overview and Session Info layout", {
  report_source_dir <- test_report_source_dir()
  skip_if(!nzchar(report_source_dir), "Quarto report sources not available")

  report_visual_abstracts <- c(
    "Grp2Analysis_V2_R6.qmd" = "differential-expression.png",
    "Grp2Analysis_V2_SE_tabset.qmd" = "differential-expression-tabset.png",
    "DiffExpQC_R6_tabset.qmd" = "differential-expression-qc.png",
    "QC_ProteinAbundances_tabset.qmd" = "protein-abundances.png",
    "QCandSSE_tabset.qmd" = "quality-control-sample-size.png"
  )
  for (report_file in names(report_visual_abstracts)) {
    qmd <- paste(
      readLines(file.path(report_source_dir, report_file), warn = FALSE),
      collapse = "\n"
    )
    expect_match(qmd, "# Overview", fixed = TRUE, info = report_file)
    expect_match(
      qmd,
      ".report_overview_cards",
      fixed = TRUE,
      info = report_file
    )
    expect_match(
      qmd,
      sprintf('<img src="%s"', report_visual_abstracts[[report_file]]),
      fixed = TRUE,
      info = report_file
    )
    expect_match(qmd, "# Session Info", fixed = TRUE, info = report_file)
    expect_match(qmd, "## Report provenance", fixed = TRUE, info = report_file)
    expect_match(qmd, "## R session info", fixed = TRUE, info = report_file)
  }
})

test_that("visual abstracts are individual vignette extras", {
  source_vignette_dir <- file.path(test_package_path(), "vignettes")
  skip_if_not(
    file.exists(file.path(source_vignette_dir, ".install_extras")),
    "Source vignette extras are not available"
  )

  extra_files <- readLines(
    file.path(source_vignette_dir, ".install_extras"),
    warn = FALSE
  )
  expected_extra_files <- c(
    "differential-expression\\.png$",
    "differential-expression-tabset\\.png$",
    "differential-expression-qc\\.png$",
    "protein-abundances\\.png$",
    "quality-control-sample-size\\.png$"
  )

  expect_equal(
    file.exists(file.path(
      source_vignette_dir,
      prolfquapp:::.quarto_visual_abstract_names
    )),
    rep(TRUE, length(prolfquapp:::.quarto_visual_abstract_names))
  )
  expect_equal("visual_abstracts" %in% extra_files, FALSE)
  expect_equal(
    tail(extra_files, length(expected_extra_files)),
    expected_extra_files
  )
})

test_that("vendored FGCZ Quarto assets match the installed template", {
  source_vignette_dir <- file.path(test_package_path(), "vignettes")
  asset_names <- c(
    "_metadata.yml",
    "fgcz.scss",
    "fgcz_header_quarto.html",
    "fgcz-plot-finder.html"
  )
  skip_if_not(
    all(file.exists(file.path(source_vignette_dir, asset_names))),
    "Source Quarto assets are not available"
  )

  source_hashes <- unname(tools::md5sum(
    fgczquartotemplate::fgcz_quarto_dir(asset_names)
  ))
  vendored_hashes <- unname(tools::md5sum(
    file.path(source_vignette_dir, asset_names)
  ))

  expect_equal(vendored_hashes, source_hashes)
})
