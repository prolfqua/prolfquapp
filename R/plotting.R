#' nice plot for paired analysis
#'
#' One plot per element of the modelled hierarchy level (e.g. protein); a finer level
#' (e.g. peptide) is drawn as facets.
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
  rhk <- bb$relevant_hierarchy_keys()
  hk <- bb$hierarchy_keys()
  child <- if (length(rhk) < length(hk)) hk[length(rhk) + 1] else NULL
  nested <- bb$data_long() |>
    tidyr::nest(.by = dplyr::all_of(rhk))
  titles <- do.call(paste, c(nested[rhk], sep = " ~ "))
  treatment <- bb$factor_keys()[1]
  pairing <- bb$factor_keys()[2]
  plotL <- function(tr, pid) {
    p <- ggplot2::ggplot(
      tr,
      ggplot2::aes(
        x = .data[[treatment]],
        y = .data[[bb$response()]],
        group = .data[[pairing]],
        colour = .data[[pairing]]
      )
    ) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::labs(title = pid) +
      prolfqua::theme_prolfqua(grid = TRUE)
    if (!is.null(child)) {
      p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[child]]))
    }
    p
  }
  plots <- purrr::map2(nested$data, titles, plotL)
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
