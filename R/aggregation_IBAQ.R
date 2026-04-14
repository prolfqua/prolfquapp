#' dataset transform data
#'
#' @param lfqdata LFQData object
#' @param agg_method aggregation method
#' @param N number of top peptides for topN aggregation
#' @export
#' @examples
#'
#' xx <- prolfqua::sim_lfq_data_peptide_config()
#' lfqdata <- prolfqua::LFQData$new(xx$data, xx$config)
#' aggregated <- aggregate_data(lfqdata, agg_method = "medpolish")
#' aggregated$response()
#' aggregated <- aggregate_data(lfqdata, agg_method = "rlm")
#' aggregated$response()
#' aggregated <- aggregate_data(lfqdata, agg_method = "topN")
#' aggregated$response()
#'
aggregate_data <- function(
  lfqdata,
  agg_method = c("medpolish", "rlm", "topN"),
  N = 3
) {
  agg_method <- match.arg(agg_method)
  if (
    length(lfqdata$hierarchy_keys()) == lfqdata$get_config()$hierarchy_depth
  ) {
    warning("nothing to aggregate from, returning unchanged data.")
    return(lfqdata)
  }

  if (agg_method == "topN") {
    aggregator <- lfqdata$get_Aggregator("topN", N = N)
    aggregator$aggregate()
    lfqdata <- aggregator$lfq_agg
  } else if (agg_method == "rlm" || agg_method == "medpolish") {
    transformed <- lfqdata$get_Transformer()$intensity_array(log)$lfq
    aggregator <- transformed$get_Aggregator(agg_method)
    aggregator$aggregate()
    ag <- aggregator$lfq_agg
    tr <- ag$get_Transformer()
    tr <- tr$intensity_array(exp, force = TRUE)
    lfqdata <- tr$lfq
    lfqdata$is_transformed(FALSE)
  } else {
    logger::log_warn("no such aggregator {agg_method}.")
  }
  return(lfqdata)
}

#' compute IBAQ values
#' @param lfqdata LFQData object
#' @param protein_annotation ProteinAnnotation object
#' @param protein_length column name for protein length
#' @param nr_tryptic_peptides column name for number of tryptic peptides
#' @export
#'
#' @examples
#' pAlf <- sim_data_protAnnot()
#' xd <- compute_IBAQ_values(pAlf$lfqdata, pAlf$pannot)
#' xd$response()
#' pAlf <- sim_data_protAnnot(PROTEIN = TRUE)
#' xd <- compute_IBAQ_values(pAlf$lfqdata, pAlf$pannot)
#' xd$response()
#'
compute_IBAQ_values <- function(
  lfqdata,
  protein_annotation,
  protein_length = "protein_length",
  nr_tryptic_peptides = "nr_tryptic_peptides"
) {
  required <- c(protein_length, nr_tryptic_peptides)
  stopifnot(all(required %in% colnames(protein_annotation$row_annot)))
  rel_annot <- dplyr::select(
    protein_annotation$row_annot,
    c(protein_annotation$pID, required)
  )
  lfqdata$set_config_value("hierarchy_depth", 1) # you want to roll up to protein
  lfqdataProtTotal <- prolfquapp::aggregate_data(
    lfqdata,
    agg_method = "topN",
    N = 10000
  )
  lfqdataProtTotal$set_data(dplyr::inner_join(
    lfqdataProtTotal$data_long(),
    rel_annot,
    by = protein_annotation$pID
  ))
  lfqdataProtTotal$set_data(lfqdataProtTotal$data_long() |>
    dplyr::mutate(
      IBAQValue_proteinLength = !!sym(lfqdataProtTotal$response()) /
        !!sym(protein_length)
    ))
  lfqdataProtTotal$set_data(lfqdataProtTotal$data_long() |>
    dplyr::mutate(
      IBAQValue = !!sym(lfqdataProtTotal$response()) /
        ifelse(!!sym(nr_tryptic_peptides) > 0, !!sym(nr_tryptic_peptides), 1)
    ))
  lfqdataProtTotal$get_config()$set_response("IBAQValue")
  return(lfqdataProtTotal)
}
