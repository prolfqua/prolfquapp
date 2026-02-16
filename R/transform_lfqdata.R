#' transform lfq data using robscale, vsn or log2, Assumes that data is not transformed (still needs log2 transformation)
#'
#' It will also run internal but then robscale must be used.
#' @param lfqdata \code{\link{LFQData}}
#' @param method normalization method to use
#' @param internal a data.frame with protein ids to be used for internal calibration, column name must be the same as
#' @export
#' @examples
#' istar <- prolfqua::prolfqua_data('data_ionstar')$filtered()
#' config <- prolfqua:::old2new(istar$config)
#' tmp <- prolfqua::LFQData$new(istar$data, config)
#' d <- istar$d
#' internal <- dplyr::filter(d, protein_Id %in% sample(unique(d$protein_Id), 3 )) |>
#'   dplyr::select(all_of(tmp$config$table$hierarchy_keys()[1])) |> dplyr::distinct()
#' tmp2 <- transform_lfqdata(tmp, internal = internal)
#' tmp2 <- transform_lfqdata(tmp)
#'
transform_lfqdata <- function(lfqdata, method = c("robscale", "vsn", "none", "log2"), internal = NULL) {
  method <- match.arg(method)
  lt <- lfqdata$get_Transformer()
  if (method == "robscale") {
    logger::log_info("Transforming using robscale.")
    transformed <- lt$log2()$robscale()$lfq
    if (!is.null(internal)) {
      logger::log_info("Transforming using robscale,")
      logger::log_info("Transforming Attempt of internal calibration.")
      subset <- lfqdata$get_subset(internal)$get_Transformer()$log2()$lfq
      transformed <- lfqdata$get_Transformer()$log2()$robscale_subset(subset)$lfq
    }
  } else if (method == "vsn") {
    n_samples <- length(unique(lfqdata$data[[lfqdata$config$table$sampleName]]))
    if (n_samples < 2) {
      logger::log_warn("vsn requires >= 2 samples, falling back to log2 transformation.")
      transformed <- lt$log2()$lfq
    } else {
      logger::log_info("Transforming using vsn::justvsn")
      transformed <- lt$intensity_matrix( .func = vsn::justvsn)$lfq
    }
  } else if (method == "none" || method == "log2") {
    logger::log_info("Transforming using log2")
    transformed <- lt$log2()$lfq
  } else {
    logger::log_warn("Transforming no such transformaton : {method}")
    return(NULL)
  }
  logger::log_info("Transforming data : {method}.")
  return(transformed)
}

#' transform lfq data with x^2 - apply if non log data is needed
#' @param lfqTrans transformed LFQData
#' @export
exp2 <- function(lfqTrans ) {
  if (!lfqTrans$config$table$is_response_transformed) {
    warning("Data not transformed.")
  }
  tr <- lfqTrans$get_Transformer()
  .exp2 <- function(x){
    2^x
  }
  tr$intensity_array(.exp2, force = TRUE)
  tr$lfq$config$table$is_response_transformed <- FALSE
  lfqdataProt <- tr$lfq
  return(lfqdataProt)
}





