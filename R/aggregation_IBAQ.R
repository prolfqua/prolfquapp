#' dataset transform data
#'
#' @export
#' @examples
#'
#' xx <- prolfqua::sim_lfq_data_peptide_config()
#' lfqdata <- prolfqua::LFQData$new(xx$data, xx$config)
#' aggregated <- aggregate_data(lfqdata, agg_method = "medpolish")
#' aggregated$response()
#' aggregated <- aggregate_data(lfqdata, agg_method = "lmrob")
#' aggregated$response()
#' aggregated <- aggregate_data(lfqdata, agg_method = "topN")
#' aggregated$response()
#'
aggregate_data <- function(lfqdata,
                           agg_method = c("medpolish", "lmrob", "topN"), N = 3) {
  agg_method <- match.arg(agg_method)
  if (length(lfqdata$config$table$hierarchy_keys()) == lfqdata$config$table$hierarchyDepth) {
    warning('nothing to aggregate from, returning unchanged data.')
    return(lfqdata)
  }

  if (agg_method == "topN") {
    aggregator <- lfqdata$get_Aggregator()
    aggregator$sum_topN(N = N)
    lfqdata <- aggregator$lfq_agg
  } else if (agg_method == "lmrob" || agg_method == "medpolish") {
    transformed <- lfqdata$get_Transformer()$intensity_array(log)$lfq
    aggregator <- transformed$get_Aggregator()
    if (agg_method == "lmrob" ) {
      aggregator$lmrob()
    } else if (agg_method == "medpolish") {
      aggregator$medpolish()
    }
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
compute_IBAQ_values <- function(lfqdata,
                                protein_annotation,
                                protein_length = "protein_length",
                                nr_tryptic_peptides = "nr_tryptic_peptides") {
  required <- c(protein_length, nr_tryptic_peptides)
  stopifnot(all(required %in% colnames(protein_annotation$row_annot)))
  rel_annot <- dplyr::select(protein_annotation$row_annot, c(protein_annotation$pID, required))
  lfqdata$config$table$hierarchyDepth <- 1 # you want to roll up to portein
  lfqdataProtTotal <- prolfquapp::aggregate_data(lfqdata, agg_method = "topN", N = 10000)
  lfqdataProtTotal$data <- dplyr::inner_join(lfqdataProtTotal$data , rel_annot, by = protein_annotation$pID)
  lfqdataProtTotal$data <- lfqdataProtTotal$data |>
    dplyr::mutate(IBAQValue_proteinLength = !!sym(lfqdataProtTotal$response()) / !!sym(protein_length))
  lfqdataProtTotal$data <- lfqdataProtTotal$data |>
    dplyr::mutate(IBAQValue = !!sym(lfqdataProtTotal$response())  / ifelse(!!sym(nr_tryptic_peptides) > 0,!!sym(nr_tryptic_peptides), 1 ))
  lfqdataProtTotal$config$table$set_response("IBAQValue")
  return(lfqdataProtTotal)
}

