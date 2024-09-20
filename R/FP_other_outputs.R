#' get report.tsv and fasta file location in folder
#' @return list with paths to data and fasta
#' @export
get_FP_multi_site_files <- function(path){
  psm_file <- dir(path = path, pattern = "abundance_multi-site_None.tsv", recursive = TRUE, full.names = TRUE)
  fasta.files <- grep("*.fasta$", dir(path = path, recursive = TRUE,full.names = TRUE), value = TRUE)
  if (any(grepl("database[0-9]*.fasta$", fasta.files))) {
    fasta.files <- grep("database[0-9]*.fasta$", fasta.files, value = TRUE)
  }
  if (length(fasta.files) == 0) {
    logger::log_error("No fasta file found!")
    stop()
  }
  return(list(data = psm_file, fasta = fasta.files))
}


#' preprocess FP multisite, filter by purity_threshold and PeptideProphetProb
#' @return list with lfqdata and protein annotation
#' @export
preprocess_FP_multi_site <- function(
    quant_data,
    fasta,
    annotation,
    pattern_contaminants = "^zz|^CON|Cont_",
    pattern_decoys = "^REV_|^rev_"){

  annot <- annotation$annot
  atable <- annotation$atable
  annot <- annot |> dplyr::mutate(
    raw.file = gsub("^x|.d.zip$|.raw$","",
                    (basename(annot[[atable$fileName]]))
    ))

  xx <- readr::read_tsv(quant_data)
  quant_idx_start <- grep(pattern = "ReferenceIntensity", x = colnames(xx))  + 1
  multiSite_long <- xx |>
    tidyr::pivot_longer(cols = all_of(quant_idx_start:ncol(xx)), values_to = "abundance", names_to = "channel")


  # join with anno again this should work now with Name # if not all samples are used in the dataset they would be removed here (to be tested)
  multiSite_long <- dplyr::inner_join(x = annotation$annot, y = multiSite_long)
  # add missing required parameters (qvalue)
  multiSite_long$qValue <- 1 - multiSite_long$MaxPepProb
  multiSite_long$nr_children  <- 1


  # Setup configuration manually for peptide analysis (phospho)
  atable$ident_Score = "MaxPepProb"
  atable$ident_qValue = "qValue"
  atable$nr_children = "nr_children"
  atable$hierarchy[["protein_Id"]] <- c("ProteinID")
  atable$hierarchy[["site"]] <- c("Index", "Peptide")
  atable$set_response("abundance")
  atable$hierarchyDepth <- 2

  # Preprocess data - aggregate proteins.
  config <- prolfqua::AnalysisConfiguration$new(atable)
  adata <- prolfqua::setup_analysis(multiSite_long, config)
  lfqdata <- prolfqua::LFQData$new(adata, config)
  lfqdata$remove_small_intensities(threshold = 1)


  # Create fasta annotation
  # Create Site Annotation
  site_annot <- multiSite_long |>
    dplyr::select(c("Index", "ProteinID", "Peptide", "SequenceWindow", "Start", "End", "MaxPepProb", "ReferenceIntensity")) |>
    dplyr::distinct()
  phosSite <- site_annot |> dplyr::rowwise() |> dplyr::mutate(siteinfo = gsub(ProteinID, "", Index))
  phosSite <- phosSite |>
    tidyr::separate_wider_delim(siteinfo, names = c(NA, "startModSite", "endModSite", "NumPhos", "LocalizedNumPhos", "PhosSites"), delim = "_",
                                too_few = "align_start")

  split_codes <- function(x) {
    if (is.na(x)) return(NA)
    return(gsub("([A-Z]\\d+)(?=[A-Z]\\d+)", "\\1;", x, perl = TRUE))
  }
  phosSite$PhosSites <- sapply(phosSite$PhosSites, split_codes)


  nrPep_exp <- multiSite_long |>
    dplyr::select(ProteinID, Peptide) |>
    dplyr::distinct() |>
    dplyr::group_by(ProteinID) |>
    dplyr::summarize(nrPeptides = dplyr::n()) |> dplyr::ungroup()

  fasta_annot <- prolfquapp::get_annot_from_fasta(fasta, pattern_decoys = pattern_decoys)
  fasta_annot <- dplyr::left_join(nrPep_exp, fasta_annot, by = c(ProteinID = "proteinname"), multiple = "all")
  fasta_annot <- fasta_annot |> dplyr::rename(description = fasta.header)
  fasta_annot2 <- dplyr::inner_join(fasta_annot, phosSite, by = "ProteinID")

  # Make names to match lfqdata
  fasta_annot2 <- fasta_annot2 |> dplyr::rename(!!lfqdata$config$table$hierarchy_keys_depth()[1] := !!rlang::sym("ProteinID"))
  fasta_annot2 <- fasta_annot2 |> dplyr::mutate(!!lfqdata$config$table$hierarchy_keys_depth()[2] := paste(!!rlang::sym("Index"),!!rlang::sym("Peptide"), sep = "~"))
  prot_annot <- prolfquapp::ProteinAnnotation$new(
    lfqdata ,
    fasta_annot2,
    description = "description",
    cleaned_ids = "protein_Id",
    full_id = "protein_Id",
    exp_nr_children = "nrPeptides",
    pattern_contaminants = pattern_contaminants,
    pattern_decoys = pattern_decoys
  )
  return(list(lfqdata = lfqdata , protein_annotation = prot_annot))
}



#' get report.tsv and fasta file location in folder
#' @return list with paths to data and fasta
#' @export
get_FP_combined_STY_files <- function(path){
  psm_file <- dir(path = path, pattern = "^combined_site_STY_.+\\.tsv", recursive = TRUE, full.names = TRUE)
  fasta.files <- grep("*.fasta$|*.fas$", dir(path = path, recursive = TRUE,full.names = TRUE), value = TRUE)
  fp.manifest <- grep("*.fp-manifest",  dir(path = path, recursive = TRUE,full.names = TRUE), value = TRUE)
  if (any(grepl("database[0-9]*.fasta$", fasta.files))) {
    fasta.files <- grep("database[0-9]*.fasta$", fasta.files, value = TRUE)
  }
  if (length(fasta.files) == 0) {
    logger::log_error("No fasta file found!")
    stop()
  }
  return(list(data = psm_file, fasta = fasta.files, fp.manifest = fp.manifest))
}


#' reads combined_site_STY file and converts to long format.
#' @return tidy data table
#' @export
#'
read_combined_STY_file <- function(file){
  xd <- readr::read_tsv(file, show_col_types = FALSE)
  colnames(xd) <- gsub("Localization Probability", "Localization_Probability", colnames(xd))
  colnames(xd) <- gsub("MaxLFQ Intensity", "MaxLFQ_Intensity", colnames(xd))
  xd <- xd |> dplyr::rename(BLP = "Best Localization_Probability")
  tidy_data <- xd |>
    tidyr::pivot_longer(
      cols = tidyselect::contains("Localization_Probability") | tidyselect::contains("Intensity") | tidyselect::contains("MaxLFQ_Intensity"),
      names_to = c("SampleName", ".value"),
      names_sep = " "
    )
  tidy_data <- tidy_data |> dplyr::rename(ProteinID = !!sym("Protein ID"))
  return(tidy_data)
}


#' preprocess FP multisite, filter by purity_threshold and PeptideProphetProb
#' @return list with lfqdata and protein annotation
#' @export
#' @param annotation_join_by column in annotation file
preprocess_FP_combined_STY <- function(
    quant_data,
    fasta,
    annotation,
    pattern_contaminants = "^zz|^CON|Cont_",
    pattern_decoys = "^REV_|^rev_",
    annotation_join_by = c("raw.file", "Name")
){

  annotation_join_by <- match.arg(annotation_join_by)
  pattern_contaminants = "^zz|^CON"
  pattern_decoys = "REV_"

  annot <- annotation$annot
  atable <- annotation$atable
  annot <- annot |> dplyr::mutate(
    raw.file = gsub("^x|.d.zip$|.raw$","",
                    (basename(normalize_path(annot[[atable$fileName]])))
    ))

  multiSite_long <- prolfquapp::read_combined_STY_file(quant_data)
  # join with anno again this should work now with Name # if not all samples are used in the dataset they would be removed here (to be tested)
  by = "SampleName"
  names(by) = annotation_join_by
  multiSite_long <- dplyr::inner_join(x = annotation$annot, y = multiSite_long, by = by)

  # add missing required parameters (qvalue)
  multiSite_long$qValue <- 0
  multiSite_long$nr_children  <- 1


  # Setup configuration manually for peptide analysis (phospho)
  atable$ident_Score = "Localization_Probability"
  atable$ident_qValue = "qValue"
  atable$nr_children = "nr_children"
  atable$hierarchy[["protein_Id"]] <- c("ProteinID")
  atable$hierarchy[["site"]] <- c("Index", "Peptide")
  atable$set_response("Intensity")
  atable$hierarchyDepth <- 2

  # Preprocess data - aggregate proteins.
  config <- prolfqua::AnalysisConfiguration$new(atable)
  adata <- prolfqua::setup_analysis(multiSite_long, config)

  lfqdata <- prolfqua::LFQData$new(adata, config)
  lfqdata$remove_small_intensities(threshold = 1)


  # Create fasta annotation
  # Create Site Annotation
  site_annot <- multiSite_long |>
    dplyr::select(c("Index", "Protein", "ProteinID", "Peptide", "BLP")) |>
    dplyr::distinct()

  phosSite <- site_annot |> dplyr::rowwise() |> dplyr::mutate(siteinfo = gsub(ProteinID, "", Index))
  phosSite <- phosSite |> dplyr::rowwise() |> dplyr::mutate(siteinfo = gsub("^_", "", siteinfo))

  nrPep_exp <- multiSite_long |>
    dplyr::select(ProteinID, Peptide) |>
    dplyr::distinct() |>
    dplyr::group_by(ProteinID) |>
    dplyr::summarize(nrPeptides = dplyr::n()) |> dplyr::ungroup()

  fasta_annot <- prolfquapp::get_annot_from_fasta(fasta, pattern_decoys = pattern_decoys)
  fasta_annot <- dplyr::left_join(nrPep_exp, fasta_annot, by = c(ProteinID = "proteinname"), multiple = "all")
  fasta_annot <- fasta_annot |> dplyr::rename(description = fasta.header)
  fasta_annot2 <- dplyr::inner_join(fasta_annot, phosSite, by = "ProteinID")

  # Make names to match lfqdata
  fasta_annot2 <- fasta_annot2 |> dplyr::rename(!!lfqdata$config$table$hierarchy_keys_depth()[1] := !!rlang::sym("ProteinID"))
  fasta_annot2 <- fasta_annot2 |> dplyr::mutate(!!lfqdata$config$table$hierarchy_keys_depth()[2] := paste(!!rlang::sym("Index"),!!rlang::sym("Peptide"), sep = "~"))
  prot_annot <- prolfquapp::ProteinAnnotation$new(
    lfqdata ,
    fasta_annot2,
    description = "description",
    cleaned_ids = "protein_Id",
    full_id = "protein_Id",
    exp_nr_children = "nrPeptides",
    pattern_contaminants = pattern_contaminants,
    pattern_decoys = pattern_decoys
  )
  return(list(lfqdata = lfqdata , protein_annotation = prot_annot))
}




