#' convert mzmine features to tidy table
#' @param data path to csv or data frame of mzMine features
#' @return A long-format tibble with one row per feature and data file.
#' @export
#' @examples
#' if(FALSE){
#' raw_df <- readr::read_csv("outputs-20250407T1707/mzmine/result_features.csv")
#' res <- tidy_mzMineFeatures(raw_df)
#' head(res)
#' }
tidy_mzMineFeatures <- function(data) {
  x <- if (is.character(data)) readr::read_csv(data) else data

  drop <- c("ion_identities", "compound_db_identity", "lipid_annotations", "molecular_networking")
  stopifnot(nrow(na.omit(dplyr::select(x, starts_with(drop)))) == 0)
  x <- dplyr::select(x, -starts_with(drop))

  colnames(x) <- gsub("^datafile:", "datafile_", colnames(x))
  feature_cols <- !grepl("^datafile_", colnames(x))
  colnames(x)[feature_cols] <- paste0("feature_", colnames(x)[feature_cols])

  xl <- x |>
    tidyr::pivot_longer(
      cols = starts_with("datafile_"),
      names_to = c("datafile", ".value"),
      names_pattern = "^(.+?):(.+)$"
    )
  colnames(xl) <- gsub(":", "_", colnames(xl))
  xl$datafile <- gsub("^datafile_", "", xl$datafile)
  xl <- xl |>
    tidyr::unite("metabolite_feature_Id", feature_id, feature_mz, feature_rt, feature_charge, remove = FALSE)
  return(xl)
}

#' get best feature annotation.
#' @param x data frame of feature annotations
#' @export
feature_annotation_get_best_score <- function(x) {
  x |>
    dplyr::group_by(id) |>
    dplyr::slice_max(score, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
}

#' get best feature annotation.
#' @param x data frame of feature annotations
#' @export
# nolint next: object_length_linter.
feature_annotation_collapse_to_single_row <- function(x) {
  x |>
    dplyr::group_by(id) |>
    dplyr::summarise(
      dplyr::across(.cols = setdiff(names(x), "id"), ~ stringr::str_c(as.character(.x), collapse = "; ")),
      .groups = "drop"
    )
}

#' get feature annotation.
#' @param x data frame of feature annotations
#' @param .loader function to select best annotations
#' @export
#' @examples
#' if (FALSE) {
#' x <- readr::read_csv("outputs-20250407T1707/mzmine/result_annotations.csv")
#' bestscore <- make_feature_annotation(x)
#' }
make_feature_annotation <- function(
  x,
  .loader = feature_annotation_get_best_score
) {
  .loader(x) |>
    dplyr::rename(annotation_rt = rt) |>
    tidyr::unite("description", compound_name, adduct, score, mol_formula, sep = ";", remove = FALSE)
}

#' get mzmine fliles
#' @param path path to data directory
#' @export
#' @examples
#' path <- "WU323671_mzMine_o35537_WpH9V2_neg_v2_result"
#' files <- get_mzMine_files(path)
get_mzMine_files <- function(path) {
  files <- dir(path = path, recursive = TRUE, full.names = TRUE)
  return(list(
    data = grep("*_features.csv$", files, value = TRUE),
    fasta = grep("*_annotations.csv$", files, value = TRUE)
  ))
}

#' preprocess mzMine input
#' @param quant_data path to mzMine features csv file
#' @param fasta_file path to annotations csv file
#' @param annotation annotation list from read_annotation
#' @param pattern_contaminants regex pattern for contaminants
#' @param pattern_decoys regex pattern for decoys
#' @param annotated if TRUE only keep annotated features
#' @param nr_peptides accepted for interface uniformity but ignored (mzMine features have one child per protein)
#' @return A list containing the prepared \code{LFQData} and
#'   \code{ProteinAnnotation} objects.
#' @export
#' @examples
#' if(FALSE){
#' annotation <- read_annotation(readr::read_tsv("outputs-20250407T1707/bfabric/input_dataset.tsv"), QC = TRUE)
#' files <- get_mzMine_files("outputs-20250407T1707/")
#' res <- preprocess_mzMine(files$data, files$fasta , annotation, annotated = TRUE)
#' }
preprocess_mzMine <- function(
  quant_data,
  fasta_file,
  annotation,
  pattern_contaminants = NULL,
  pattern_decoys = NULL,
  annotated = FALSE,
  # nr_peptides accepted for interface uniformity but ignored: mzMine features
  # have one child per protein, so a >= 2 cut would drop every feature.
  nr_peptides = 1
) {
  xdl <- tidy_mzMineFeatures(quant_data)
  m_annot <- make_feature_annotation(readr::read_csv(fasta_file))
  join <- if (annotated) dplyr::inner_join else dplyr::right_join
  xdl <- join(m_annot, xdl, by = c("id" = "feature_id"), relationship = "many-to-many")

  annot <- annotation$annot
  config <- annotation$atable$clone(deep = TRUE)
  annot$relative_path <- basename(annot$relative_path)
  nr <- sum(annot$relative_path %in% sort(unique(xdl$datafile)))
  logger::log_info("nr : ", nr, " files annotated out of ", length(unique(xdl$datafile)))
  stopifnot(nr > 0)
  config$hierarchy[["metabolite_feature_Id"]] <- "metabolite_feature_Id"
  config$set_response("area")
  byv <- c(stats::setNames("datafile", config$file_name), intersect(colnames(annot), colnames(xdl)))

  feature <- dplyr::inner_join(annot, xdl, by = byv, multiple = "all")
  lfqdata <- prolfqua::LFQData$new(prolfqua::setup_analysis(feature, config), config)
  lfqdata$remove_small_intensities()

  m_annot <- xdl |>
    dplyr::select("metabolite_feature_Id", "id", "feature_rt", "feature_mz", "feature_charge", "description") |>
    dplyr::distinct() |>
    dplyr::mutate(
      exp_children = 1,
      nrPeptides = 1,
      protein_length = 1,
      nr_tryptic_peptides = 1,
      IDcolumn = metabolite_feature_Id
    )
  prot_annot <- prolfquapp::ProteinAnnotation$new(
    lfqdata,
    m_annot,
    description = "description",
    cleaned_ids = "IDcolumn",
    full_id = "metabolite_feature_Id",
    exp_nr_children = "exp_children",
    pattern_contaminants = NULL,
    pattern_decoys = NULL
  )
  return(list(lfqdata = lfqdata, protein_annotation = prot_annot))
}
