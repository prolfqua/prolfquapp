#' nice plot for paired analysis
#' @param bb LFQData object with paired data
#' @param outpath output directory for plots
#' @export
#' @examples
#'
#'
#' xd <- prolfqua::sim_lfq_data_protein_config(with_missing = FALSE, paired = TRUE)
#' xd <- prolfqua::LFQData$new(xd$data, xd$config)
#' xa <- prolfquapp::writeLinesPaired(xd)
#' xa[[1]]
#' xa[[2]]
#' xa[[3]]
#' xa[[4]]
#'
writeLinesPaired <- function(bb, outpath = NULL) {
  nested <- bb$data_long() |>
    dplyr::ungroup() |>
    dplyr::group_by(!!!rlang::syms(bb$hierarchy_keys())) |>
    tidyr::nest()
  tr <- nested$data[[1]] # nolint object_usage_linter. used inside plotL
  plotL <- function(tr, pid) {
    ggplot2::ggplot(
      tr,
      ggplot2::aes(
        x = !!sym(bb$factor_keys()[1]),
        y = !!sym(bb$response()),
        group = !!sym(bb$factor_keys()[2]),
        colour = !!sym(bb$factor_keys()[2])
      )
    ) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::labs(title = pid) +
      ggplot2::theme_bw()
  }
  plots <- purrr::map2(nested[[2]], nested[[1]], plotL)
  if (!is.null(outpath)) {
    fn <- file.path(outpath, "boxplots_paired.pdf")
    grDevices::pdf(fn, width = 8, height = 6)
    tryCatch(
      for (p in plots) {
        print(p)
      },
      finally = grDevices::dev.off()
    )
  }
  return(plots)
}
