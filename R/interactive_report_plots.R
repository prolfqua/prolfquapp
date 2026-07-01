# Interactive (ggiraph) report plots.
#
# Kept in prolfquapp (not prolfqua core) because they are report-presentation
# helpers and depend on the report-only `ggiraph` package (Suggests). They also
# consume prolfquapp-specific annotation columns (the contaminant `CON` flag).

.need_ggiraph <- function() {
  if (!requireNamespace("ggiraph", quietly = TRUE)) {
    stop(
      "Package 'ggiraph' is required for interactive report plots. ",
      "Install it with install.packages('ggiraph').",
      call. = FALSE
    )
  }
}

#' Interactive volcano plot with contaminants marked (ggiraph)
#'
#' Builds an interactive volcano from an annotated contrast table (as produced by
#' \code{DEAnalyse$get_annotated_contrasts()}). Points are coloured by
#' \code{colour_col} (e.g. \code{estimate_type}); contaminants (rows where
#' \code{contaminant_col} is \code{TRUE}) are drawn as triangles so they are
#' visually flagged while kept in the plot. Rows with a missing effect or score
#' (e.g. decoys, which carry NA statistics) are dropped, so decoys never appear.
#' Hovering a point shows its identifier and statistics.
#'
#' @param data annotated contrast data.frame
#' @param effect name of the effect-size column (e.g. \code{"diff"} or
#'   \code{"log2_EFCs"})
#' @param score name of the score column plotted as \code{-log10} (e.g.
#'   \code{"FDR"} or \code{"BFDR"})
#' @param contrast name of the contrast (facet) column
#' @param id_cols hierarchy-key column(s) united into the point label / tooltip
#' @param fc_threshold effect-size threshold (green vertical ablines); \code{NULL}
#'   to omit
#' @param score_threshold score threshold (red horizontal abline at
#'   \code{-log10(score_threshold)}); \code{NULL} to omit
#' @param contaminant_col name of the logical contaminant flag column
#'   (default \code{"CON"}); ignored if absent
#' @param colour_col name of the point colour column (default
#'   \code{"estimate_type"}); ignored if absent (single colour)
#' @param width_svg,height_svg girafe canvas size in inches
#' @return a \code{ggiraph::girafe} interactive htmlwidget
#' @export
#' @family reporting
#' @examples
#' contr <- data.frame(
#'   protein_Id = c("sp|P1|X", "sp|P2|X", "zz|CON1|X"),
#'   contrast = "A_vs_B",
#'   diff = c(1.2, -0.3, 2.1),
#'   FDR = c(0.001, 0.5, 0.2),
#'   estimate_type = "observed",
#'   CON = c(FALSE, FALSE, TRUE)
#' )
#' if (requireNamespace("ggiraph", quietly = TRUE)) {
#'   g <- volcano_ggiraph(contr, id_cols = "protein_Id")
#'   stopifnot("girafe" %in% class(g))
#' }
volcano_ggiraph <- function(
  data,
  effect = "diff",
  score = "FDR",
  contrast = "contrast",
  id_cols = "protein_Id",
  fc_threshold = 1,
  score_threshold = 0.05,
  contaminant_col = "CON",
  colour_col = "estimate_type",
  width_svg = 9,
  height_svg = 7
) {
  .need_ggiraph()
  data <- as.data.frame(data)
  stopifnot(all(c(effect, score, contrast) %in% colnames(data)))
  id_cols <- id_cols[id_cols %in% colnames(data)]
  if (length(id_cols) == 0) {
    stop("none of id_cols present in data")
  }

  label <- do.call(paste, c(lapply(id_cols, function(x) as.character(data[[x]])), sep = "~"))
  is_con <- if (contaminant_col %in% colnames(data)) {
    isTRUE_vec <- data[[contaminant_col]]
    !is.na(isTRUE_vec) & as.logical(isTRUE_vec)
  } else {
    rep(FALSE, nrow(data))
  }
  colour <- if (colour_col %in% colnames(data)) {
    as.character(data[[colour_col]])
  } else {
    rep("feature", nrow(data))
  }
  colour[is.na(colour)] <- "NA"

  plt <- data.frame(
    .effect = suppressWarnings(as.numeric(data[[effect]])),
    .score = suppressWarnings(as.numeric(data[[score]])),
    .contrast = as.character(data[[contrast]]),
    .colour = colour,
    .contaminant = ifelse(is_con, "contaminant", "target"),
    .label = label,
    stringsAsFactors = FALSE
  )
  plt$.y <- -log10(plt$.score)
  plt <- plt[is.finite(plt$.effect) & is.finite(plt$.y), , drop = FALSE]
  plt$.tooltip <- paste0(
    plt$.label, "\n", plt$.contrast,
    "\n", effect, " = ", round(plt$.effect, 2),
    "\n", score, " = ", signif(plt$.score, 2),
    ifelse(plt$.contaminant == "contaminant", "\n(contaminant)", "")
  )

  p <- ggplot2::ggplot(
    plt,
    ggplot2::aes(x = .data$.effect, y = .data$.y)
  ) +
    ggiraph::geom_point_interactive(
      ggplot2::aes(
        colour = .data$.colour,
        shape = .data$.contaminant,
        tooltip = .data$.tooltip,
        data_id = .data$.label
      ),
      alpha = 0.7,
      size = 1.9
    ) +
    ggplot2::scale_shape_manual(
      values = c(target = 16, contaminant = 17),
      name = NULL
    ) +
    ggplot2::facet_wrap(ggplot2::vars(.data$.contrast)) +
    ggplot2::labs(
      x = effect,
      y = paste0("-log10(", score, ")"),
      colour = colour_col
    ) +
    ggplot2::theme_bw()

  if (!is.null(score_threshold) && score_threshold > 0) {
    p <- p + ggplot2::geom_hline(
      yintercept = -log10(score_threshold),
      colour = "red", linetype = "dashed"
    )
  }
  if (!is.null(fc_threshold) && fc_threshold > 0) {
    p <- p + ggplot2::geom_vline(
      xintercept = c(-fc_threshold, fc_threshold),
      colour = "darkgreen", linetype = "dashed"
    )
  }
  ggiraph::girafe(ggobj = p, width_svg = width_svg, height_svg = height_svg)
}


#' Interactive per-sample intensity density plot (ggiraph)
#'
#' Interactive kernel-density curves of the working intensity, one line per
#' sample; hovering a curve highlights it and shows the sample name. Accepts one
#' or more \code{LFQData} objects as a named list, rendered as facets (e.g.
#' empirical vs normalized abundance).
#'
#' @param lfqdata_list a single \code{LFQData} or a named list of \code{LFQData}
#'   (names become facet panel titles)
#' @param legend show the (sample) colour legend; forced off above
#'   \code{max_legend_samples}
#' @param max_legend_samples suppress the legend above this many samples
#' @param width_svg,height_svg girafe canvas size in inches
#' @return a \code{ggiraph::girafe} interactive htmlwidget
#' @export
#' @family reporting
#' @examples
#' istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 20)
#' lfq <- prolfqua::LFQData$new(istar$data, istar$config)
#' if (requireNamespace("ggiraph", quietly = TRUE)) {
#'   g <- intensity_density_ggiraph(lfq)
#'   stopifnot("girafe" %in% class(g))
#' }
intensity_density_ggiraph <- function(
  lfqdata_list,
  legend = TRUE,
  max_legend_samples = 16,
  width_svg = 8,
  height_svg = 4
) {
  .need_ggiraph()
  if (inherits(lfqdata_list, "LFQData")) {
    lfqdata_list <- list(intensity = lfqdata_list)
  }
  if (is.null(names(lfqdata_list))) {
    names(lfqdata_list) <- paste0("panel", seq_along(lfqdata_list))
  }

  dens_one <- function(lfq, panel) {
    resp <- lfq$response()
    samp <- lfq$sample_name()
    d <- lfq$data_long()
    parts <- split(d[[resp]], d[[samp]])
    out <- lapply(names(parts), function(s) {
      vals <- parts[[s]][is.finite(parts[[s]])]
      if (length(vals) < 2) {
        return(NULL)
      }
      dd <- stats::density(vals)
      data.frame(x = dd$x, y = dd$y, sample = s, panel = panel, stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
  }

  dens <- do.call(rbind, mapply(
    dens_one, lfqdata_list, names(lfqdata_list),
    SIMPLIFY = FALSE
  ))
  n_samples <- length(unique(dens$sample))

  p <- ggplot2::ggplot(
    dens,
    ggplot2::aes(x = .data$x, y = .data$y, colour = .data$sample, group = .data$sample)
  ) +
    ggiraph::geom_line_interactive(
      ggplot2::aes(tooltip = .data$sample, data_id = .data$sample)
    ) +
    ggplot2::facet_wrap(ggplot2::vars(.data$panel), scales = "free") +
    ggplot2::labs(x = "intensity", y = "density", colour = "sample") +
    ggplot2::theme_bw()

  if (!isTRUE(legend) || n_samples > max_legend_samples) {
    p <- p + ggplot2::guides(colour = "none")
  }
  ggiraph::girafe(ggobj = p, width_svg = width_svg, height_svg = height_svg)
}
