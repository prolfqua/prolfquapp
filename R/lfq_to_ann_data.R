#' convert lfqdata to anndata
#' @export
#'
#' @examples
#' # example code
#' library(prolfqua)
#' library(prolfquapp)
#' lfqdata <- prolfqua::sim_lfq_data_2Factor_config()
#' lfqdata <- LFQData$new(lfqdata$data, lfqdata$config)
#' lfqdata$data$protein_Id <- add_RevCon(lfqdata$data$protein_Id)
#' pids <- grep("^zz|^REV", unique(lfqdata$data$protein_Id), value = TRUE, invert = TRUE)
#' addannot <- data.frame(
#'   protein_Id = pids,
#'   description = stringi::stri_rand_strings(length(pids), 13)
#' )
#'
#' addannot <- addannot |> tidyr::separate(protein_Id, c("cleanID", NA), remove = FALSE)
#' pannot <- ProteinAnnotation$new(lfqdata,
#'   addannot,
#'   description = "description",
#'   cleaned_ids = "cleanID",
#'   pattern_contaminants = "^zz",
#'   pattern_decoys = "^REV"
#' )
#'
#' #debug(anndata_from_LFQData)
#' anndata_from_LFQData(lfqdata, pannot)
#' #anndataR::write_h5ad(adata, path = "test.h5ad", mode = "w")
anndata_from_LFQData <- function(lfqdata, pannot) {
  stopifnot(inherits(lfqdata, "LFQData"))
  stopifnot(inherits(pannot, "ProteinAnnotation"))


  layers <- list()
  message("converting to layers: ", paste(lfqdata$config$table$value_vars(), collapse = ", ") )
  for (val in lfqdata$config$table$value_vars()) {
    X <- lfqdata$to_wide(as.matrix = TRUE, value = val)$data
    rownames(X) <- gsub("~lfq~light", "", rownames(X))
    layers[[val]] <- t(X)
  }

  o <- as.data.frame(lfqdata$factors())
  rownames(o) <- o[,lfqdata$config$table$sampleName]

  v <- as.data.frame(pannot$row_annot)
  rownames(v) <- v[,pannot$full_id]
  X <- layers[[1]]
  v <- v[colnames(X), ]
  o <- o[rownames(X), ]
  # vars are proteins or metabolites or genes
  # obs are cells, or samples,

  adata <- anndataR::AnnData(
    X = X,
    var = v,
    obs = o,
    layers = layers
  )
  return(adata)
}


