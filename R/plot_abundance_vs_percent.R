#' Plot relative protein abundance as a function of rank by abundance
#' @export
#' @param sr LFQDataSummarizer
#' @param top_percent deprecate
#' @return ggplot2
#' @examples
#' library(prolfqua)
#' istar <- prolfqua_data('data_ionstar')
#' istar$config <- old2new(istar$config)
#' data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(data, istar$config)
#' sr <- lfqdata$get_Summariser()
#'
#' plot_abundance_vs_percent(sr, top_N = 6, factors = FALSE, cumulative = TRUE)
#' pd <- plot_abundance_vs_percent(sr, top_N = NULL, factors = FALSE, cumulative = FALSE)
#'
#' plot_abundance_vs_percent(sr, top_N = 4, factors = TRUE, cumulative = FALSE)
#' plot_abundance_vs_percent(sr, top_N = NULL, factors = TRUE, cumulative = TRUE)
#'
#'
plot_abundance_vs_percent <- function(
    percInfo,
    cfg_table,
    top_N = 10,
    factors = TRUE ,
    colors = c("^REV_" =  "red",
               "^CON_" = "orange"),
    cumulative = TRUE,
    group = "BB") {

  columnAb <- if (cumulative) {"abundance_percent_cumulative"} else {"abundance_percent"}
  protID <- cfg_table$hierarchy_keys_depth()

  #Select relevant columns
  percInfo <- percInfo |>
    dplyr::select(dplyr::all_of(c(protID,"percent_prot", columnAb, cfg_table$factor_keys_depth())))

  if (!factors) {
    percInfo <- percInfo |> dplyr::filter(!!rlang::sym(cfg_table$factor_keys_depth()[1]) == "ALL")
  }
  colorV <- rep("black", nrow(percInfo))

  for (i in seq_along(colors)) {
    colorV[grepl(names(colors)[i], percInfo[[protID]])] <- colors[i]
  }
  percInfo$color <- colorV


  if (!is.null(top_N)) {
    topN <- percInfo |>
      dplyr::group_by(dplyr::across(cfg_table$factor_keys_depth())) |>
      dplyr::slice_max(order_by = !!rlang::sym(columnAb), n = top_N)
  } else {
    message("creating shared data with key : ", paste(" ~ ", protID ))
    percInfo <- crosstalk::SharedData$new(as.data.frame(percInfo) , key = as.formula(paste(" ~ ", protID )),
                                          group = group)
  }

  myplot <- ggplot(percInfo, aes(x = !!rlang::sym("percent_prot"),
                            y = !!rlang::sym(columnAb),
                            label = !!rlang::sym(protID))) +
    geom_point(color = colorV) +
    facet_wrap(as.formula(paste0("~", paste(cfg_table$factor_keys_depth(), sep = " + "))))

  if (!is.null(top_N) ) {
    myplot <- myplot + ggrepel::geom_label_repel(
      data = topN,
      aes(label = !!rlang::sym(protID)), size = 3, max.overlaps = 100)
  }
  return( myplot)
}
