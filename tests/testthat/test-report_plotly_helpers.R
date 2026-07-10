test_that("plotly_ggplot_subplot de-duplicates legend entries across panels", {
  testthat::skip_if_not_installed("plotly")

  p1 <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, colour = factor(cyl))) +
    ggplot2::geom_density()
  p2 <- ggplot2::ggplot(mtcars, ggplot2::aes(disp, colour = factor(cyl))) +
    ggplot2::geom_density()

  widget <- plotly_ggplot_subplot(p1, p2)
  legend_data <- lapply(widget$x$data, function(trace) {
    list(
      name = trace$name,
      showlegend = trace$showlegend,
      legendgroup = trace$legendgroup
    )
  })
  visible_names <- vapply(
    legend_data,
    function(trace) if (isTRUE(trace$showlegend)) trace$name else NA_character_,
    character(1)
  )
  visible_names <- stats::na.omit(visible_names)

  expect_equal(sort(visible_names), sort(unique(visible_names)))
  legend_groups <- vapply(legend_data, `[[`, character(1), "legendgroup")
  expect_equal(sort(visible_names), sort(unique(legend_groups)))
})

test_that("plotly_ggplot_subplot fades keyed traces on hover", {
  testthat::skip_if_not_installed("plotly")

  keyed_mtcars <- plotly::highlight_key(
    mtcars,
    key = as.character(mtcars$cyl),
    group = "test-density"
  )
  p <- ggplot2::ggplot(
    keyed_mtcars,
    ggplot2::aes(mpg, colour = factor(cyl), group = factor(cyl))
  ) +
    ggplot2::geom_density()

  widget <- plotly_ggplot_subplot(p)

  expect_equal(widget$x$highlight$on, "plotly_hover")
  expect_equal(widget$x$highlight$opacityDim, 0.15)
  expect_equal(as.character(widget$x$highlight$ctGroups), "test-density")
})
