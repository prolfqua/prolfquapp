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

test_that("index href helper handles mixed absolute and relative paths", {
  result_dir <- file.path(tempdir(), "DEA index mixed paths")
  result_subdir <- file.path(result_dir, "Results_WU_1")
  dir.create(result_subdir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(result_dir, recursive = TRUE), add = TRUE)

  report_file <- file.path(result_subdir, "DE_WU1_quarto.html")
  file.create(report_file)

  oldwd <- setwd(dirname(result_dir))
  on.exit(setwd(oldwd), add = TRUE)

  expect_equal(
    prolfquapp:::.index_relative_href(
      normalizePath(report_file, mustWork = TRUE),
      basename(result_dir)
    ),
    "./Results_WU_1/DE_WU1_quarto.html"
  )
})

test_that("estimate_type palette fixes known colours and fills unknown levels", {
  palette <- prolfquapp:::.estimate_type_palette(c(
    "observed",
    "lod_imputed",
    "missing_fallback",
    "some_future_type"
  ))

  # sorted alphabetically by level name
  expect_named(
    palette,
    c("lod_imputed", "missing_fallback", "observed", "some_future_type")
  )
  expect_equal(palette[["observed"]], "black")
  expect_equal(palette[["lod_imputed"]], "orange")
  expect_false(any(is.na(names(palette))))
  expect_false(any(is.na(palette)))
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

test_that(".join_annotation preserves every quant/result row and never multiplies", {
  annotation <- data.frame(
    protein_Id = c("P1", "P2", "P3"),
    description = c("a", "b", "c"),
    stringsAsFactors = FALSE
  )
  # x has several rows per protein (e.g. per-contrast / long format)
  x <- data.frame(
    protein_Id = c("P1", "P1", "P2", "P3", "P3"),
    contrast = c("c1", "c2", "c1", "c1", "c2"),
    diff = 1:5,
    stringsAsFactors = FALSE
  )
  res <- prolfquapp:::.join_annotation(annotation, x, "protein_Id")
  expect_equal(nrow(res), nrow(x)) # no row dropped, no row multiplied
  expect_true(all(c("description", "contrast", "diff") %in% colnames(res)))
  expect_equal(sort(res$diff), 1:5)
})

test_that(".join_annotation keeps result rows lacking annotation (NA enrich)", {
  annotation <- data.frame(protein_Id = "P1", description = "a",
    stringsAsFactors = FALSE)
  x <- data.frame(protein_Id = c("P1", "P2"), v = c(10, 20),
    stringsAsFactors = FALSE)
  res <- prolfquapp:::.join_annotation(annotation, x, "protein_Id")
  expect_equal(nrow(res), 2)
  expect_true(is.na(res$description[res$protein_Id == "P2"]))
})

test_that(".join_annotation errors when annotation is not unique (safety net)", {
  annotation <- data.frame(protein_Id = c("P1", "P1"), description = c("a", "b"),
    stringsAsFactors = FALSE)
  x <- data.frame(protein_Id = "P1", v = 1, stringsAsFactors = FALSE)
  expect_error(
    prolfquapp:::.join_annotation(annotation, x, "protein_Id"),
    "not unique"
  )
})

test_that(".join_annotation joins on the shared hierarchy keys (PTM: protein_Id + site)", {
  # Both frames carry `site`; joining only on protein_Id would yield
  # site.x/site.y and drop the bare `site` the report's unite(all_of(hkey)) needs.
  annotation <- data.frame(
    protein_Id = c("P1", "P1", "P2"), site = c("s1", "s2", "s1"),
    description = c("a", "a", "b"), stringsAsFactors = FALSE
  )
  x <- data.frame(
    protein_Id = c("P1", "P1", "P2"), site = c("s1", "s2", "s1"),
    contrast = "c1", diff = c(1, 2, 3), stringsAsFactors = FALSE
  )
  res <- prolfquapp:::.join_annotation(annotation, x, c("protein_Id", "site"))
  expect_true("site" %in% colnames(res))
  expect_false(any(grepl("^site\\.", colnames(res)))) # no site.x / site.y
  expect_equal(nrow(res), nrow(x))
  expect_true(all(c("description", "diff") %in% colnames(res)))
})

test_that(".join_annotation never joins on a coincidentally-shared value column", {
  # `avgAbd` is shared by name but is NOT a hierarchy key -> must not be a join
  # key (else the mismatch would break enrichment). Joins on protein_Id only.
  annotation <- data.frame(
    protein_Id = c("P1", "P2"), avgAbd = c(99, 99),
    description = c("a", "b"), stringsAsFactors = FALSE
  )
  x <- data.frame(
    protein_Id = c("P1", "P2"), avgAbd = c(1, 2), diff = c(5, 6),
    stringsAsFactors = FALSE
  )
  res <- prolfquapp:::.join_annotation(annotation, x, c("protein_Id", "site"))
  expect_equal(nrow(res), 2)
  # enrichment succeeds despite avgAbd differing (avgAbd was not a join key)
  expect_equal(res$description[res$protein_Id == "P1"], "a")
})
