#' Combine ggplots into a Plotly subplot
#'
#' Converts ggplot objects with `plotly::ggplotly()`, combines them with
#' `plotly::subplot()`, and keeps only one legend entry per trace name across
#' panels. Traces with the same name are assigned the same legend group so one
#' legend click toggles the corresponding traces in all panels.
#'
#' @param ... ggplot or Plotly objects to combine.
#' @param showlegend Logical. If `FALSE`, hide all subplot legend entries.
#' @param nrows Number of subplot rows passed to `plotly::subplot()`.
#' @param titleX,titleY Logical values passed to `plotly::subplot()`.
#' @param margin Numeric margin passed to `plotly::subplot()`.
#' @param legend_groupclick Plotly legend group click behaviour. Use `NULL` to
#'   leave the Plotly default unchanged.
#' @return A Plotly htmlwidget.
#' @export
#' @examples
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   p1 <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, colour = factor(cyl))) +
#'     ggplot2::geom_density()
#'   p2 <- ggplot2::ggplot(mtcars, ggplot2::aes(disp, colour = factor(cyl))) +
#'     ggplot2::geom_density()
#'   plotly_ggplot_subplot(p1, p2)
#' }
plotly_ggplot_subplot <- function(
    ...,
    showlegend = TRUE,
    nrows = 1,
    titleX = TRUE,
    titleY = TRUE,
    margin = 0.05,
    legend_groupclick = "togglegroup") {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop(
      "Package 'plotly' is required to build this interactive subplot.",
      call. = FALSE
    )
  }

  plots <- list(...)
  if (length(plots) == 0) {
    stop("At least one plot must be supplied.", call. = FALSE)
  }

  plotly_plots <- lapply(plots, function(plot) {
    if (inherits(plot, "plotly")) {
      plot
    } else {
      plotly::ggplotly(plot)
    }
  })

  plotly_obj <- do.call(
    plotly::subplot,
    c(
      plotly_plots,
      list(nrows = nrows, titleX = titleX, titleY = titleY, margin = margin)
    )
  )
  plotly_obj <- .deduplicate_plotly_legend(plotly_obj, showlegend = showlegend)

  if (is.null(legend_groupclick)) {
    return(plotly_obj)
  }
  plotly::layout(plotly_obj, legend = list(groupclick = legend_groupclick))
}

.deduplicate_plotly_legend <- function(plotly_obj, showlegend = TRUE) {
  legend_names <- character()
  for (i in seq_along(plotly_obj$x$data)) {
    trace_name <- plotly_obj$x$data[[i]]$name
    if (is.null(trace_name) || length(trace_name) != 1 || is.na(trace_name)) {
      plotly_obj$x$data[[i]]$showlegend <- FALSE
      next
    }
    trace_name <- as.character(trace_name)
    if (!nzchar(trace_name)) {
      plotly_obj$x$data[[i]]$showlegend <- FALSE
      next
    }
    plotly_obj$x$data[[i]]$legendgroup <- trace_name
    plotly_obj$x$data[[i]]$showlegend <-
      isTRUE(showlegend) && !(trace_name %in% legend_names)
    legend_names <- c(legend_names, trace_name)
  }
  plotly_obj
}
