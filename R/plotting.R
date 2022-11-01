#' nice plot for paired analysis
#' @export
#'
writeLinesPaired <- function(bb, outpath) {
  nested <- bb$data |> ungroup() |> dplyr::group_by(!!!syms(bb$config$table$hierarchy_keys())) |> tidyr::nest()
  tr <- nested$data[[1]]
  plotL <- function(tr, pid){
    ggplot2::ggplot(tr, ggplot2::aes_string(x = bb$config$table$factor_keys()[1],
                          y = bb$config$table$get_response(),
                          group = bb$config$table$factor_keys()[2],
                          colour = bb$config$table$factor_keys()[2] )) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::labs(title = pid)  +
      ggplot2::theme_bw()
  }
  plots <- purrr::map2(nested[[2]], nested[[1]], plotL)
  pdf(file.path(outpath, "lines_paired.pdf"))
  lapply(plots, print)
  dev.off()

}
