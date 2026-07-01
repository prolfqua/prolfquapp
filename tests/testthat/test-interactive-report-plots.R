# Interactive (ggiraph) report plots: volcano with contaminant marking + density.

skip_if_no_ggiraph <- function() {
  testthat::skip_if_not_installed("ggiraph")
}

make_contrasts <- function() {
  data.frame(
    protein_Id = c("sp|P1|X", "sp|P2|X", "zz|CON1|X", "REV_P9", "sp|P4|X"),
    contrast = "A_vs_B",
    diff = c(1.2, -0.3, 2.1, NA, 0.1),
    FDR = c(0.001, 0.5, 0.2, NA, 0.9),
    estimate_type = c("observed", "observed", "observed", "observed", "lod_imputed"),
    CON = c(FALSE, FALSE, TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
}

test_that("volcano_ggiraph returns a girafe widget and drops NA-stat (decoy) rows", {
  skip_if_no_ggiraph()
  g <- volcano_ggiraph(make_contrasts(), id_cols = "protein_Id")
  expect_true("girafe" %in% class(g))
})

test_that("volcano_ggiraph works without a CON column (no contaminants marked)", {
  skip_if_no_ggiraph()
  d <- make_contrasts()
  d$CON <- NULL
  expect_true("girafe" %in% class(volcano_ggiraph(d, id_cols = "protein_Id")))
})

test_that("volcano_ggiraph works without estimate_type (single colour)", {
  skip_if_no_ggiraph()
  d <- make_contrasts()
  d$estimate_type <- NULL
  expect_true("girafe" %in% class(volcano_ggiraph(d, id_cols = "protein_Id")))
})

test_that("volcano_ggiraph handles saint-style effect/score columns", {
  skip_if_no_ggiraph()
  d <- data.frame(
    protein_Id = c("sp|P1|X", "zz|CON1|X"),
    Bait = "bait1",
    log2_EFCs = c(2.0, 1.0),
    BFDR = c(0.01, 0.3),
    CON = c(FALSE, TRUE)
  )
  g <- volcano_ggiraph(
    d,
    effect = "log2_EFCs", score = "BFDR", contrast = "Bait",
    id_cols = "protein_Id"
  )
  expect_true("girafe" %in% class(g))
})

test_that("volcano_ggiraph errors when required columns are missing", {
  skip_if_no_ggiraph()
  expect_error(
    volcano_ggiraph(data.frame(x = 1), id_cols = "x"),
    "effect"
  )
})

test_that("intensity_density_ggiraph returns a girafe for a single LFQData", {
  skip_if_no_ggiraph()
  istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 15)
  lfq <- prolfqua::LFQData$new(istar$data, istar$config)
  expect_true("girafe" %in% class(intensity_density_ggiraph(lfq)))
})

test_that("intensity_density_ggiraph facets a named list of LFQData", {
  skip_if_no_ggiraph()
  istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 15)
  lfq <- prolfqua::LFQData$new(istar$data, istar$config)
  g <- intensity_density_ggiraph(
    list("empirical" = lfq, "normalized" = lfq),
    legend = FALSE
  )
  expect_true("girafe" %in% class(g))
})
