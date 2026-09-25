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
      logger::log_warn("vsn requires >= 2 samples, falling back to log2 transformation.")
      transformed <- lt$log2()$lfq
    } else {
      logger::log_info("Transforming using vsn::justvsn")
      transformed <- lt$intensity_matrix(.func = vsn::justvsn)$lfq
    }
  } else {
    logger::log_info("Transforming using log2")
    transformed <- lt$log2()$lfq
  }
  logger::log_info("Transforming data : {method}.")
  return(transformed)
}
