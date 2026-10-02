#' Static volcano plot of differential-abundance results
#'
#' One panel per contrast: the fold change against \eqn{-\log_{10}} of the
#' significance score. Features passing both thresholds (score below
#' \code{fdr_threshold}, absolute effect above \code{fc_threshold}, the rule of
#' \code{filter_significant()}) are coloured by direction and the most
#' significant of them labelled; all other features are drawn small and grey.
#' The plot is a plain ggplot, so it renders as a static image however many
#' features the table holds; interactive exploration is left to exploreDE.
#'
#' @param data contrast table, one row per feature and contrast
#' @param effect column with the log2 fold change
#' @param significance column with the score plotted on the y-axis (e.g. FDR)
#' @param contrast column with the contrast name, one panel each
#' @param label column naming the features to label, or \code{NULL} for no labels
#' @param fc_threshold fold-change threshold, drawn as vertical lines
#' @param fdr_threshold significance threshold, drawn as a horizontal line
#' @param directional if \code{TRUE}, only an effect above \code{fc_threshold}
#'   counts as significant (as for SAINT), mirroring the backend's
#'   \code{significance_directional}
#' @param n_labels number of labelled features per contrast and direction
#' @param base_size base font size
#' @return a ggplot object
#' @export
#' @examples
#' set.seed(1)
#' contrasts <- data.frame(
#'   feature = paste0("P", 1:2000),
#'   contrast = rep(c("B_vs_A", "C_vs_A"), each = 1000),
#'   diff = rnorm(2000, sd = 1.2)
#' )
#' contrasts$FDR <- p.adjust(2 * pnorm(-abs(contrasts$diff), sd = 0.5), "BH")
#' volcano_static(contrasts, label = "feature")
volcano_static <- function(
  data,
  effect = "diff",
  significance = "FDR",
  contrast = "contrast",
  label = NULL,
  fc_threshold = 1,
  fdr_threshold = 0.1,
  directional = FALSE,
  n_labels = 5,
  base_size = 12
) {
  data <- as.data.frame(data)
  data <- data[!is.na(data[[effect]]) & !is.na(data[[significance]]), , drop = FALSE]
  score <- data[[significance]]
  # A score of exactly 0 would sit at -log10(0) = Inf; floor it at the smallest
  # positive score so it is drawn at the top of the panel instead of dropped.
  positive <- score[score > 0]
  floor_score <- if (length(positive) > 0) min(positive) else .Machine$double.xmin
  data$neg_log10_score <- -log10(pmax(score, floor_score))

  effect_passes <- if (directional) data[[effect]] > fc_threshold else abs(data[[effect]]) > fc_threshold
  passes <- score < fdr_threshold & effect_passes
  data$call <- factor(
    ifelse(!passes, "not significant", ifelse(data[[effect]] > 0, "increased", "decreased")),
    levels = c("decreased", "not significant", "increased")
  )
  call_colours <- c(decreased = "#2166AC", `not significant` = "grey75", increased = "#B2182B")
  background <- data[data$call == "not significant", , drop = FALSE]
  significant <- data[data$call != "not significant", , drop = FALSE]

  # The panel title carries the significant counts, so they never collide
  # with the point labels.
  counts <- table(factor(significant[[contrast]], levels = unique(data[[contrast]])), significant$call)
  strip_titles <- sprintf(
    "%s  (%d down, %d up)",
    rownames(counts),
    counts[, "decreased"],
    counts[, "increased"]
  )
  background[[contrast]] <- factor(background[[contrast]], levels = rownames(counts), labels = strip_titles)
  significant[[contrast]] <- factor(significant[[contrast]], levels = rownames(counts), labels = strip_titles)

  y_label <- paste0("-log10(", significance, ")")
  p <- ggplot2::ggplot(mapping = ggplot2::aes(x = .data[[effect]], y = .data$neg_log10_score, colour = .data$call)) +
    ggplot2::geom_point(data = background, size = 0.5, alpha = 0.5, stroke = 0) +
    ggplot2::geom_point(data = significant, size = 1.2, alpha = 0.8, stroke = 0) +
    ggplot2::geom_vline(
      xintercept = if (directional) fc_threshold else c(-fc_threshold, fc_threshold),
      linetype = "dashed",
      colour = "grey40",
      linewidth = 0.3
    ) +
    ggplot2::geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    ggplot2::scale_colour_manual(values = call_colours, drop = FALSE, name = NULL) +
    ggplot2::facet_wrap(ggplot2::vars(.data[[contrast]])) +
    ggplot2::labs(x = paste0("log2 fold change (", effect, ")"), y = y_label) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95"),
      legend.position = "bottom"
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1)))

  if (!is.null(label) && n_labels > 0 && nrow(significant) > 0) {
    # One label per name: on peptide or site level a protein's features crowd
    # the top, and repeating its gene name would hide every other hit.
    top <- significant |>
      dplyr::arrange(dplyr::desc(.data$neg_log10_score)) |>
      dplyr::distinct(.data[[contrast]], .data$call, .data[[label]], .keep_all = TRUE) |>
      dplyr::group_by(.data[[contrast]], .data$call) |>
      dplyr::slice_max(.data$neg_log10_score, n = n_labels, with_ties = FALSE) |>
      dplyr::ungroup()
    p <- p +
      ggrepel::geom_text_repel(
        data = top,
        ggplot2::aes(label = .data[[label]]),
        size = (base_size - 3) / ggplot2::.pt,
        max.overlaps = Inf,
        min.segment.length = 0,
        segment.size = 0.2,
        show.legend = FALSE
      )
  }
  p
}
