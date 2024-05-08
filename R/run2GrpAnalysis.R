#' transform lfq data using robscale, vsn or non,
#'
#' It will also run internal but then robscale must be used.
#' @param lfqdata \code{\link{LFQData}}
#' @param method normalization method to use
#' @param internal a data.frame with protein ids to be used for internal calibration, column name must be the same as
#' @export
#' @examples
#' istar <- prolfqua_data('data_ionstar')$filtered()
#' config <- prolfqua:::old2new(istar$config)
#' tmp <- prolfqua::LFQData$new(istar$data, config)
#' d <- istar$d
#' internal <- dplyr::filter(d, protein_Id %in% sample(unique(d$protein_Id), 3 )) |>
#'   dplyr::select(all_of(tmp$config$table$hierarchy_keys()[1])) |> dplyr::distinct()
#' tmp2 <- transform_lfqdata(tmp, internal = internal)
#' tmp2 <- transform_lfqdata(tmp)
#'
transform_lfqdata <- function(lfqdata, method = c("robscale", "vsn", "none"), internal = NULL) {
  method <- match.arg(method)
  lt <- lfqdata$get_Transformer()
  if (method == "robscale") {
    transformed <- lt$log2()$robscale()$lfq
    if (!is.null(internal)) {
      logger::log_info("Attempt of internal calibration.")
      subset <- lfqdata$get_subset(internal)$get_Transformer()$log2()$lfq
      transformed <- lfqdata$get_Transformer()$log2()$robscale_subset(subset)$lfq
    }
  } else if (method == "vsn") {
    transformed <- lt$intensity_matrix( .func = vsn::justvsn)$lfq
  } else if (method == "none") {
    transformed <- lt$log2()$lfq
  } else {
    logger::log_warn("no such transformaton : {method}")
    return(NULL)
  }
  logger::log_info("data transformed : {method}.")
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



#' convert tibble to data.frame with rownames
#' @param .data a tibble or data.frame
#' @param var name of the column with new row.names
#' @return a data.frame with rownames
#' @export
#' @examples
#' ind <- tibble::tibble(a = 1:3, rowname = letters[1:3])
#' column_to_rownames(ind)
column_to_rownames <- function(.data, var = "rowname"){
  res <- as.data.frame(.data)
  rownames(res) <- .data[[var]]
  return(res)
}

strip_rownames <- function(.data, strip="~lfq~light$"){
  newrnames <- gsub(strip, "", rownames(.data))
  rownames(.data) <- newrnames
  return(.data)
}


#' build bfabric urls
#' @export
#' @examples
#'
#'
#' ps <- ProjectSpec$new()
#' ps$project_Id <- 32258
#' ps$order_Id <- 34628
#' ps$workunit_Id <- 302212
#' bfabric_url_builder(ps)
#'
#' ps <- ProjectSpec$new()
#' ps$order_Id <- 34628
#' ps$workunit_Id <- 302212
#' bfabric_url_builder(ps)
#'
bfabric_url_builder <- function(project_spec){

  orderURL <- NULL
  workunitURL <- NULL
  projectURL <- NULL
  orderID <- as.numeric(project_spec$order_Id)
  if ((length(orderID) > 0) && !is.na(orderID)) {
    orderURL <- paste0("https://fgcz-bfabric.uzh.ch/bfabric/order/show.html?id=", orderID, "&tab=details")
  }
  workunitID <- as.numeric(project_spec$workunit_Id)
  if ((length(workunitID) > 0) && !is.na(workunitID)) {
    workunitURL <- paste0("https://fgcz-bfabric.uzh.ch/bfabric/workunit/show.html?id=", orderID, "&tab=details")
  }
  projectID <- as.numeric(project_spec$project_Id)
  if ((length(projectID) > 0) && !is.na(projectID)) {
    projectURL <- paste0("https://fgcz-bfabric.uzh.ch/bfabric/project/show.html?id=", orderID, "&tab=details")
  }

  return(list(orderURL = orderURL, projectURL = projectURL, workunitURL = workunitURL))
}

#' Convert prolfqua differential expression analysis results to SummarizedExperiment
#'
#' @rdname make_DEA_report
#' @param GRP2 return value of \code{\link{make_DEA_report}}
#' @return SummarizedExperiment
#' @export
#' @family workflow
#'
#'
#'
make_SummarizedExperiment <- function(GRP2, colname = NULL, rowname = NULL, strip="~lfq~light", .url_builder = bfabric_url_builder){
  if (is.null(colname)) {
    colname <- GRP2$RES$lfqData$config$table$sampleName
  }
  if (is.null(rowname)) {
    rowname <- GRP2$RES$lfqData$config$table$hierarchyKeys()
  }
  resTables <- write_DEA(GRP2,".", write = FALSE)
  matTr <- GRP2$RES$transformedlfqData$to_wide(as.matrix = TRUE)
  matRaw <- GRP2$RES$transformedlfqData$to_wide(as.matrix = TRUE)

  mat.raw <- strip_rownames(matRaw$data, strip)
  mat.trans <- strip_rownames(matTr$data, strip)
  col.data <- column_to_rownames(matRaw$annotation, var = colname)
  col.data <- col.data[colnames(mat.raw),]
  x <- SummarizedExperiment::SummarizedExperiment(
    assays = list(rawData = mat.raw, transformedData = mat.trans),
    colData = col.data, metadata = list(bfabric_urls = .url_builder(GRP2$project_spec), contrasts = resTables$contrasts, formula = resTables$formula )
  )

  diffbyContrast <- split(resTables$diff_exp_analysis, resTables$diff_exp_analysis$contrast)
  for (i in names(diffbyContrast)) {
    row.data <- column_to_rownames(diffbyContrast[[i]], var = rowname)
    row.data <- row.data[rownames(mat.raw),]

    SummarizedExperiment::rowData(x)[[paste0("constrast_",i)]] <- row.data
  }

  SummarizedExperiment::rowData(x)[["stats_normalized_wide"]] <- column_to_rownames(resTables$stats_normalized_wide, var = rowname)[rownames(mat.raw),]
  SummarizedExperiment::rowData(x)[["stats_raw_wide"]] <- column_to_rownames(resTables$stats_raw_wide, var = rowname)[rownames(mat.raw),]
  return(x)
}







