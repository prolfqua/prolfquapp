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
dataset_set_factors <- function(atable, peptide, REPEATED = TRUE) {
  atable$hierarchyDepth <- 1
  if (sum(grepl("^name", colnames(peptide), ignore.case = TRUE)) > 0) {
    atable$sampleName <- grep("^name", colnames(peptide), value = TRUE, ignore.case = TRUE)
  }

  stopifnot(sum(grepl("^group|^bait", colnames(peptide), ignore.case = TRUE)) == 1)
  groupingVAR <- grep("^group|^bait", colnames(peptide), value = TRUE, ignore.case = TRUE)
  peptide[[groupingVAR]] <- gsub("[[:space:]]", "", peptide[[groupingVAR]])
  peptide[[groupingVAR]] <- gsub("[-\\+\\/\\*]", "_", peptide[[groupingVAR]])

  atable$factors[["Group_"]] = groupingVAR
  atable$factorDepth <- 1

  if (sum(grepl("^subject", colnames(peptide), ignore.case = TRUE)) == 1 & REPEATED) {
    subvar <- grep("^subject", colnames(peptide), value = TRUE, ignore.case = TRUE)
    atable$factors[["Subject"]] = subvar
    tmp <- data.frame(table(dplyr::distinct(peptide[,c(groupingVAR,subvar)])) )
    if (all(tmp$Freq > 1)) {
      atable$factorDepth <- 2
    }
  }
  if (sum(grepl("^control", colnames(peptide), ignore.case = TRUE)) == 1) {
    atable$factors[["CONTROL"]] = grep("^control", colnames(peptide), value = TRUE, ignore.case = TRUE)
  }
  return(list(atable = atable , peptide = peptide))
}

#' dataset protein annot
#'
#' @export
#'
dataset_protein_annot <- function(
    peptide,
    atable,
    protein_annot = "Protein.Description") {
  proteinID <- atable$hkeysDepth()
  prot_annot <- dplyr::select(
    peptide ,
    dplyr::all_of(c( atable$hierarchy[[proteinID]], protein_annot))) |>
    dplyr::distinct()
  prot_annot <- dplyr::rename(prot_annot, description = !!rlang::sym(protein_annot))
  prot_annot <- dplyr::rename(prot_annot, !!proteinID := (!!atable$hierarchy[[proteinID]]))

  # figure out if this is an uniprot database.

  UNIPROT <- mean(grepl("^sp\\||^tr\\|", prot_annot[[proteinID]])) > 0.8
  message("uniprot database : ", UNIPROT)

  if (UNIPROT) {
    prot_annot <- prolfqua::get_UniprotID_from_fasta_header(prot_annot, idcolumn = proteinID)
    prot_annot <- prot_annot |> dplyr::rename(IDcolumn = UniprotID)
  } else {
    prot_annot$IDcolumn <- prot_annot[[proteinID]]
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
