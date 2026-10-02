volcano_test_contrasts <- function() {
  data.frame(
    feature = paste0("P", 1:8),
    gene = c("A", "A", "B", "C", "D", "E", "F", "G"),
    contrast = rep(c("B_vs_A", "C_vs_A"), each = 4),
    diff = c(3, 2.5, -2, 0.1, 2, -3, 0.2, NA),
    FDR = c(0.001, 0.01, 0.02, 0.9, 0, 0.03, 0.5, 0.01)
  )
}

test_that("volcano_static colours significant features by direction and counts them per contrast", {
  p <- prolfquapp::volcano_static(volcano_test_contrasts(), label = "gene", fc_threshold = 1, fdr_threshold = 0.05)
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  significant <- p$layers[[2]]$data
  expect_equal(sort(as.character(significant$call)), c("decreased", "decreased", "increased", "increased", "increased"))
  expect_true(all(is.finite(significant$neg_log10_score)))
  expect_setequal(
    levels(significant$contrast),
    c("B_vs_A  (1 down, 2 up)", "C_vs_A  (1 down, 1 up)")
  )
  labels <- p$layers[[length(p$layers)]]$data$gene
  expect_equal(sum(labels == "A"), 1)
  expect_gt(nrow(built$data[[1]]), 0)
})

test_that("volcano_static directional counts only increased effects", {
  p <- prolfquapp::volcano_static(volcano_test_contrasts(), fc_threshold = 1, fdr_threshold = 0.05, directional = TRUE)
  expect_true(all(p$layers[[2]]$data$call == "increased"))
})
