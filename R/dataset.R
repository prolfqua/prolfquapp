#' extect contrasts from dataset
#' @export
#'
dataset_extract_contrasts <- function(annot, GRP2) {
  if ( all(c("ContrastName", "Contrast") %in% colnames(annot)) ) {
    contr <- annot |>
      dplyr::select(rlang::.data$ContrastName, rlang::.data$Contrast) |>
      dplyr::filter(nchar(rlang::.data$Contrast) > 0)
    Contrasts <- contr$Contrast
    names(Contrasts) <- contr$ContrastName
    GRP2$pop$Contrasts <- Contrasts
  }
  return(GRP2)
}

#' set factors and sample names columns
#' @export
#'
#'
dataset_set_factors <- function(atable, msdata, REPEATED = TRUE) {
  atable$hierarchyDepth <- 1
  if (sum(grepl("^name", colnames(msdata), ignore.case = TRUE)) > 0) {
    atable$sampleName <- grep("^name", colnames(msdata), value = TRUE, ignore.case = TRUE)
  }

  stopifnot(sum(grepl("^group|^bait", colnames(msdata), ignore.case = TRUE)) == 1)
  groupingVAR <- grep("^group|^bait", colnames(msdata), value = TRUE, ignore.case = TRUE)
  msdata[[groupingVAR]] <- gsub("[[:space:]]", "", msdata[[groupingVAR]])
  msdata[[groupingVAR]] <- gsub("[-\\+\\/\\*]", "_", msdata[[groupingVAR]])

  atable$factors[["Group_"]] = groupingVAR
  atable$factorDepth <- 1

  if (sum(grepl("^subject", colnames(msdata), ignore.case = TRUE)) == 1 & REPEATED) {
    subvar <- grep("^subject", colnames(msdata), value = TRUE, ignore.case = TRUE)
    atable$factors[["Subject"]] = subvar
    tmp <- data.frame(table(dplyr::distinct(msdata[,c(groupingVAR,subvar)])) )
    if (all(tmp$Freq > 1)) {
      atable$factorDepth <- 2
    }
  }
  if (sum(grepl("^control", colnames(msdata), ignore.case = TRUE)) == 1) {
    atable$factors[["CONTROL"]] = grep("^control", colnames(msdata), value = TRUE, ignore.case = TRUE)
  }
  return(list(atable = atable , msdata = msdata))
}

#' dataset protein annot
#'
#' @export
#' @param idcolName
#' @param
#'
dataset_protein_annot <- function(
    msdata,
    idcolName,
    proteinID_column = "protein_Id",
    protein_annot = "fasta.header",
    more_columns = c("nrPeptides", "fasta.id")) {
  msdata <- dplyr::rename(msdata, !!proteinID_column := idcolName )
  prot_annot <- dplyr::select(
    msdata ,
    dplyr::all_of(c( proteinID_column, protein_annot, more_columns))) |>
    dplyr::distinct()
  prot_annot <- dplyr::rename(prot_annot, description = !!rlang::sym(protein_annot))
  # figure out if this is an uniprot database.

  UNIPROT <- mean(grepl("^sp\\||^tr\\|", prot_annot[[proteinID_column]])) > 0.8
  message("uniprot database : ", UNIPROT)

  if (UNIPROT) {
    prot_annot <- prolfqua::get_UniprotID_from_fasta_header(prot_annot, idcolumn = proteinID_column)
    prot_annot <- prot_annot |> dplyr::rename(!!"IDcolumn" := !!sym("UniprotID"))
  } else {
    prot_annot$IDcolumn <- prot_annot[[proteinID_column]]
  }
  return(prot_annot)
}
#' dataset transform data
#'
#' @export
#'
aggregate_data <- function(lfqdata,
                           agg_method = c("medpolish", "lmrob", "topN")) {
  agg_method <- match.arg(agg_method)

  if (agg_method == "topN") {
    aggregator <- lfqdata$get_Aggregator()
    aggregator$sum_topN()
    lfqdata <- aggregator$lfq_agg
  } else if (agg_method == "lmrob" || agg_method == "medpolish") {
    transformed <- lfqdata$get_Transformer()$log2()$lfq
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


#' extract fasta header from fasta file.
#' @export
#' @param fasta \code{\link{prozor::readPeptideFasta}}
#'
fasta_to_header <- function(fasta){
  fasta_annot <- as_tibble(
    data.frame(annot = sapply(fasta, seqinr::getAnnot)))

  fasta_annot <- fasta_annot |> tidyr::separate(.data$annot,
                                                c("proteinname","fasta.header"),
                                                sep = "\\s", extra = "merge")
  fasta_annot <- fasta_annot |> dplyr::mutate(proteinname = gsub(">","", .data$proteinname) )
  # remove duplicated id's
  fasta_annot <- fasta_annot[!duplicated(fasta_annot$proteinname),]
  proteinID <- "proteinname"
  UNIPROT <- mean(grepl("^sp\\||^tr\\|", fasta_annot[[proteinID]])) > 0.8
  message("uniprot database : ", UNIPROT)
  if (UNIPROT) {
    fasta_annot <- prolfqua::get_UniprotID_from_fasta_header(fasta_annot, idcolumn = proteinID)
    fasta_annot <- fasta_annot |> dplyr::rename(IDcolumn = UniprotID)
  } else {
    fasta_annot$IDcolumn <- fasta_annot[[proteinID]]
  }

}
