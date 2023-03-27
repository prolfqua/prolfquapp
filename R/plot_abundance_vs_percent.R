#' Plot relative protein abundance as a function of rank by abundance
#' @export
#' @return ggplot2
#' @examples
#' library(prolfqua)
#' istar <- prolfqua_data('data_ionstar')
#' istar$config <- old2new(istar$config)
#' data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(data, istar$config)
#' sr <- lfqdata$get_Summariser()
#'
#' plot_abundance_vs_percent(sr, top_percent = 20, factors = FALSE, cumulative = TRUE)
#' plot_abundance_vs_percent(sr, top_percent = 3, factors = FALSE, cumulative = FALSE)
#' plot_abundance_vs_percent(sr, top_percent = 20, factors = TRUE, cumulative = TRUE)
#' plot_abundance_vs_percent(sr, top_percent = 20, factors = TRUE, cumulative = FALSE)
#'plot_abundance_vs_percent(sr, top_percent = NULL, factors = TRUE, cumulative = FALSE)
plot_abundance_vs_percent <- function(sr,
                                      top_percent = 5,
                                      top_N = 10,
                                      factors = TRUE ,
                                      log = FALSE,
                                      colors = c("^REV_" =  "red",
                                                 "^CON_" = "orange"),
                                      cumulative = TRUE) {
  columnAb <- if (cumulative) {"abundance_percent_cumulative"} else {"abundance_percent"}
  tmp <- sr$percentage_abundance()
  protID <- sr$lfq$config$table$hierarchy_keys_depth()
  if (!factors) {
    tmp <- tmp |> dplyr::filter(!!rlang::sym(sr$lfq$config$table$factor_keys_depth()[1]) == "ALL")
  }
  colorV <- rep("black", nrow(tmp))

  for (i in seq_along(colors)) {
    colorV[grepl(names(colors)[i], tmp[[protID]])] <- colors[i]
  }
  tmp$color <- colorV

  if (!is.null(top_percent)) {
      topN <- tmp |>
        dplyr::group_by(dplyr::across(sr$lfq$config$table$factor_keys_depth())) |>
        dplyr::filter(!!rlang::sym(columnAb) > top_percent)

  } else {
    topN <- tmp |>
      dplyr::group_by(dplyr::across(sr$lfq$config$table$factor_keys_depth())) |>
      dplyr::slice_max(order_by = !!rlang::sym(columnAb), n = 5)

  }

  ggplot <- ggplot(tmp, aes(x = !!rlang::sym("percent_prot"),
                            y = !!rlang::sym(columnAb))) +
    geom_point(color = tmp$color) +
    facet_wrap(as.formula(paste0("~", paste(sr$lfq$config$table$factor_keys_depth(), sep = " + ")))) +
    ggrepel::geom_label_repel(data = topN,  aes(label = !!rlang::sym(protID)), size = 3, max.overlaps = 100) +
    if (log) {scale_y_continuous(trans = 'log10')} else {NULL}
  return(ggplot)
}
