#' Plot relative protein abundance as a function of rank by abundance
#' @export
#' @param percInfo data frame with percentage abundance info
#' @param cfg_config AnalysisConfiguration object
#' @param top_N number of top proteins to label
#' @param factors if TRUE facet by factors
#' @param colors named vector of colors for special proteins
#' @param columnAb column name for abundance values
#' @param group crosstalk group identifier
#' @param alpha point transparency for regular (black) proteins
#' @param highlight_alpha point transparency for highlighted (coloured) proteins
#' @param logY if TRUE use log10 y-axis
#' @return ggplot2
#' @examples
#'
#' library(prolfqua)
#' istar <- prolfqua::prolfqua_data('data_ionstar')$filtered()
#' data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- prolfqua::LFQData$new(data, istar$config)
#' sr <- lfqdata$get_Summariser()
#' undebug(plot_abundance_vs_percent)
#' plot_abundance_vs_percent(sr$percentage_abundance(),
#'  lfqdata$get_config(),
#'  top_N = 6, factors = FALSE, logY = TRUE)
#'
#' pd <- plot_abundance_vs_percent(
#'   sr$percentage_abundance(),
#'   lfqdata$get_config(), top_N = NULL, factors = FALSE)
#' plot_abundance_vs_percent(
#'   sr$percentage_abundance(),
#'   lfqdata$get_config(), top_N = 4, factors = TRUE)
#' plot_abundance_vs_percent(
#'   sr$percentage_abundance(),
#'   lfqdata$get_config(), top_N = NULL, factors = TRUE)
#'
#'
plot_abundance_vs_percent <- function(
  percInfo,
  cfg_config,
  top_N = 10,
  factors = TRUE,
  colors = c("^REV_" = "red", "^CON_" = "orange"),
  columnAb = "abundance_percent",
  group = "BB",
  alpha = 1,
  highlight_alpha = 1,
  logY = TRUE
) {
  protID <- cfg_config$hierarchy_keys_depth()
  #Select relevant columns
  percInfo <- percInfo |>
    dplyr::select(dplyr::all_of(c(
      protID,
      "percent_prot",
      columnAb,
      cfg_config$factor_keys_depth()
    )))

  if (!factors) {
    percInfo <- percInfo |>
      dplyr::filter(!!rlang::sym(cfg_config$factor_keys_depth()[1]) == "All")
  }
  colorV <- rep("black", nrow(percInfo))

  for (i in seq_along(colors)) {
    colorV[grepl(names(colors)[i], percInfo[[protID]])] <- colors[i]
  }
  percInfo$color <- colorV

  # Draw the highlighted proteins (contaminants / decoys) last so they are
  # plotted on top of the regular (black) points instead of being hidden behind
  # them when points overlap. order() is stable, so the abundance-rank ordering
  # within each colour group is preserved. colorV is reordered in lockstep with
  # percInfo because geom_point() receives it as a positional colour vector.
  draw_order <- order(colorV != "black")
  percInfo <- percInfo[draw_order, , drop = FALSE]
  colorV <- colorV[draw_order]

  # Highlighted proteins get their own (typically higher) opacity so they stand
  # out from the faded background points. Positional, aligned with colorV.
  alphaV <- ifelse(colorV == "black", alpha, highlight_alpha)

  if (!is.null(top_N)) {
    topN <- percInfo |>
      dplyr::group_by(dplyr::across(cfg_config$factor_keys_depth())) |>
      dplyr::slice_max(order_by = !!rlang::sym(columnAb), n = top_N)
  } else {
    message("creating shared data with key : ", paste0(" ~ ", protID))
    percInfo <- crosstalk::SharedData$new(
      as.data.frame(percInfo),
      key = as.formula(paste(" ~ ", protID)),
      group = group
    )
  }

  myplot <- ggplot(
    percInfo,
    aes(
      x = !!rlang::sym("percent_prot"),
      y = !!rlang::sym(columnAb),
      label = !!rlang::sym(protID)
    )
  ) +
    geom_point(color = colorV, alpha = alphaV) +
    facet_wrap(as.formula(paste0(
      " ~ ",
      paste(cfg_config$factor_keys_depth(), collapse = " + ")
    ))) +
    if (logY) {
      ggplot2::scale_y_log10()
    } else {
      NULL
    }

  if (!is.null(top_N)) {
    myplot <- myplot +
      ggrepel::geom_label_repel(
        data = topN,
        aes(label = !!rlang::sym(protID)),
        size = 3,
        max.overlaps = 100
      )
  }
  return(myplot)
}
