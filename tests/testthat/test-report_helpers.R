test_that("index href helper creates relative links for Windows paths", {
  result_dir <- "E:\\CD_output\\prolfqua_dea_results\\run\\DEA_20260507"
  report_file <- paste0(
    result_dir,
    "\\Results WU\\DE WUp31230_o37713_Lipidomics_pos_export.html"
  )

  expect_equal(
    prolfquapp:::.index_relative_href(report_file, result_dir),
    "./Results%20WU/DE%20WUp31230_o37713_Lipidomics_pos_export.html"
  )

  mixed_separator_file <- paste0(
    "E:/CD_output/prolfqua_dea_results/run/DEA_20260507/",
    "Results_WU/QC.html"
  )
  expect_equal(
    prolfquapp:::.index_relative_href(mixed_separator_file, result_dir),
    "./Results_WU/QC.html"
  )
})

test_that("write_index_html does not emit absolute local paths", {
  result_dir <- file.path(tempdir(), "DEA index test")
  result_subdir <- file.path(result_dir, "Results WU")
  dir.create(result_subdir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(result_dir, recursive = TRUE), add = TRUE)

  links <- list(
    dea_file = file.path(result_subdir, "DE WU.html"),
    qc_file = file.path(result_subdir, "QC WU.html"),
    quarto_file = file.path(result_subdir, "DE WU quarto.html"),
    data_files = list(
      xlsx_file = file.path(result_subdir, "DE WU.xlsx"),
      ibaq_file = "",
      ora_files = list(),
      gsea_files = list()
    )
  )

  prolfquapp::write_index_html(links, result_dir)
  html <- paste(
    readLines(file.path(result_dir, "index.html"), warn = FALSE),
    collapse = "\n"
  )

  expect_match(html, 'href="./Results%20WU/DE%20WU.html"', fixed = TRUE)
  expect_equal(
    grepl(prolfquapp:::.path_to_url_path(result_dir), html, fixed = TRUE),
    FALSE
  )
})
