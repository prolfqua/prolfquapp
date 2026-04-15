#' Transform lfq data using robscale, vsn or log2
#'
#' Assumes that data is not transformed (still needs log2 transformation)
#'
#' @param lfqdata \code{\link[prolfqua]{LFQData}}
#' @param method normalization method to use
#' @export
#' @examples
#' istar <- prolfqua::prolfqua_data('data_ionstar')$filtered()
#' tmp <- prolfqua::LFQData$new(istar$data, istar$config)
#' tmp2 <- transform_lfqdata(tmp)
#'
transform_lfqdata <- function(
  lfqdata,
  method = c("robscale", "vsn", "none", "log2")
) {
  method <- match.arg(method)
  lt <- lfqdata$get_Transformer()
  if (method == "robscale") {
    logger::log_info("Transforming using robscale.")
    transformed <- lt$log2()$robscale()$lfq
  } else if (method == "vsn") {
    n_samples <- length(unique(lfqdata$data_long()[[lfqdata$sample_name()]]))
    if (n_samples < 2) {
      logger::log_warn(
        "vsn requires >= 2 samples, falling back to log2 transformation."
      )
      transformed <- lt$log2()$lfq
    } else {
      logger::log_info("Transforming using vsn::justvsn")
      transformed <- lt$intensity_matrix(.func = vsn::justvsn)$lfq
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
exp2 <- function(lfqTrans) {
  if (!lfqTrans$is_transformed()) {
    warning("Data not transformed.")
  }
  tr <- lfqTrans$get_Transformer()
  .exp2 <- function(x) {
    2^x
  }
  tr$intensity_array(.exp2, force = TRUE)
  tr$lfq$is_transformed(FALSE)
  lfqdataProt <- tr$lfq
  return(lfqdataProt)
}
