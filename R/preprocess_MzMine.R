#' convert mzmine features to tidy table
#' @param data path to csv or data frame of mzMine features
#' @export
#' @examples
#' # example code
#' if(FALSE){
#' file_path = "outputs-20250407T1707/mzmine/result_features.csv"
#' raw_df <- readr::read_csv(file_path)
#' res <- tidy_mzMineFeatures(raw_df)
#'  file_path <- paste0(
#'    "/Users/witoldwolski/__checkout/prolfquapp/",
#'    "inst/application/mzMine/out_results_zip/",
#'    "mzmine/result_features.csv")
#' #file_path = "WU323671_mzMine_o35537_WpH9V2_neg_v2_result/mzmine/result_features.csv"
#' x <- readr::read_csv(file_path)
#' raw_df
#' res <- tidy_mzMineFeatures(raw_df)
#' head(res)
#'
#' }
tidy_mzMineFeatures <- function(data) {
  if("character" %in% class(data)){
    x <- readr::read_csv(data)
  } else if ("data.frame" %in% class(data)) {
    x <- data
  }

  xdrop <- x |> dplyr::select(
    #starts_with("alignment_scores"),
    starts_with("ion_identities"),
    starts_with("compound_db_identity"),
    starts_with("lipid_annotations"),
    starts_with("molecular_networking"))
  stopifnot(nrow(na.omit(xdrop)) == 0)


  x <- x |> dplyr::select(
    #-starts_with("alignment_scores"),
    -starts_with("ion_identities"),
    -starts_with("compound_db_identity"),
    -starts_with("lipid_annotations"),
    -starts_with("molecular_networking"))

  colnames(x) <- gsub("^datafile:", "datafile_", colnames(x))
  colnames(x)[!grepl("^datafile_",colnames(x))] <- paste0("feature_",colnames(x)[!grepl("^datafile_",colnames(x))])

  xl <- x |> tidyr::pivot_longer(cols = starts_with("datafile_"),
                                 names_to = c("datafile", ".value"),
                                 names_pattern = "^(.+?):(.+)$" )
  colnames(xl) <- gsub(":","_", colnames(xl))
  xl$datafile <- gsub("^datafile_","",xl$datafile)
  xl <- xl |> tidyr::unite("metabolite_feature_Id", feature_id, feature_mz, feature_rt, feature_charge, remove = FALSE )

  return(xl)
}

#' get best feature annotation.
#' @param x data frame of feature annotations
#' @export
feature_annotation_get_best_score <- function(x){
  best_per_id <- x |>
    dplyr::group_by(id) |>
    dplyr::slice_max(score, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  return(best_per_id)
}

#' get best feature annotation.
#' @param x data frame of feature annotations
#' @examples
#' if (FALSE) {
#'   x <- readr::read_csv("test/out_results_zip/mzmine/result_annotations.csv")
#'   x
#'   feature_annotation_collapse_to_single_row(x)
#' }
#' @export
feature_annotation_collapse_to_single_row <- function(x){
  collapsed_df <- x |>
    dplyr::group_by(id) |>
    dplyr::summarise(
      dplyr::across(
        .cols = setdiff(names(x), "id"),
        ~ stringr::str_c(as.character(.x), collapse = "; ")
      ),
      .groups = "drop"
    )
  return(collapsed_df)
}

#'
#' get feature annotation.
#' @param x data frame of feature annotations
#' @param .loader function to select best annotations
#' @export
#' @examples
#' if (FALSE) {
#' x <- readr::read_csv("outputs-20250407T1707/mzmine/result_annotations.csv")
#' bestscore <- make_feature_annotation(x)
#' dim(bestscore)
#' file_path = "outputs-20250407T1707/mzmine/result_features.csv"
#' raw_df <- read_csv(file_path)
#' res <- tidy_mzMineFeatures(raw_df)
#' intersect(colnames(res), colnames(bestscore))
#'
#' xx <- inner_join(bestscore, res, by = c("id" = "feature_id"),relationship = "many-to-many")
#'
#' dim(xx)
#' }
#' @export
#'
make_feature_annotation <- function(x, .loader = feature_annotation_get_best_score){
  res <- .loader(x)
  res <- res |> dplyr::rename(annotation_rt = rt)
  res <- res |> tidyr::unite("description", compound_name, adduct, score, mol_formula, sep = ";", remove = FALSE)
  return(res)
}



#' get mzmine fliles
#' @param path path to data directory
#' @export
#' @examples
#'
#' path <- "WU323671_mzMine_o35537_WpH9V2_neg_v2_result"
#' files <- get_mzMine_files(path)
get_mzMine_files <- function(path){
  feature <- grep("*_features.csv$", dir(path = path, recursive = TRUE, full.names = TRUE), value = TRUE)
  annot <- grep("*_annotations.csv$", dir(path = path, recursive = TRUE, full.names = TRUE), value = TRUE)
  return(list(data = feature, fasta = annot))
}


#' preprocess mzMine input
#' @param quant_data path to mzMine features csv file
#' @param fasta_file path to annotations csv file
#' @param annotation annotation list from read_annotation
#' @param pattern_contaminants regex pattern for contaminants
#' @param pattern_decoys regex pattern for decoys
#' @param annotated if TRUE only keep annotated features
#' @export
#'
#' @examples
#'
#' if(FALSE){
#' xd <- "outputs-20250407T1707/bfabric/input_dataset.tsv"
#' annot <- readr::read_tsv(xd)
#'
#' annotation <- read_annotation(annot, QC = TRUE)
#' xd <- "outputs-20250407T1707/"
#' files <- get_mzMine_files(path)
#' files
#' undebug(preprocess_mzMine)
#' res <- preprocess_mzMine(files$data, files$fasta , annotation)
#' dim(res$lfqdata$data)
#' res <- preprocess_mzMine(files$data, files$fasta , annotation, annotated = TRUE)
#' dim(res$lfqdata$data)
#' }
preprocess_mzMine <- function(
    quant_data,
    fasta_file,
    annotation,
    pattern_contaminants = NULL,
    pattern_decoys = NULL,
    annotated = FALSE
){
  xdl <- readr::read_csv(quant_data)
  xdl <- tidy_mzMineFeatures(xdl)
  annot <- readr::read_csv(fasta_file)
  annot <- make_feature_annotation(annot)

  if (annotated) {
    xdl <- dplyr::inner_join(annot,xdl, by = c("id" = "feature_id"),relationship = "many-to-many")
  } else{
    xdl <- dplyr::right_join(annot,xdl, by = c("id" = "feature_id"),relationship = "many-to-many")
  }

  annot <- annotation$annot
  atable <- annotation$atable
  annot$relative_path <- basename(annot$relative_path)
  nr <- sum(annot$relative_path %in% sort(unique(xdl$datafile)))
  logger::log_info("nr : ", nr, " files annotated out of ", length(unique(xdl$datafile)))
  stopifnot( nr > 0)
  # atable$sampleName = "file_id"
  atable$hierarchy[["metabolite_feature_Id"]] <- c("metabolite_feature_Id")
  atable$set_response("area")
  byv <- c("datafile")
  names(byv) <- atable$fileName
  byv <- c(byv, intersect(colnames(annot), colnames(xdl)))

  feature <- dplyr::inner_join(annot, xdl, by = byv, multiple = "all")
  config <- prolfqua::AnalysisConfiguration$new(atable)
  adata <- prolfqua::setup_analysis(feature, config)
  lfqdata <- prolfqua::LFQData$new(adata, config)
  lfqdata$remove_small_intensities()

  m_annot <- xdl |>
    dplyr::select("metabolite_feature_Id",
                  "id",
                  "feature_rt",
                  "feature_mz",
                  "feature_charge",
                  "description") |>
    dplyr::distinct()

  m_annot$exp_children <- 1
  m_annot$nrPeptides <- 1
  # c("protein_length", "nr_tryptic_peptides")
  m_annot$protein_length <- 1
  m_annot$nr_tryptic_peptides <- 1
  # handle not identified
  # m_annot$nr_compounds <- ifelse(m_annot$Checked, 2 ,1)

  m_annot <- m_annot |> dplyr::mutate(IDcolumn = metabolite_feature_Id)
  ProteinAnnotation$undebug("initialize")
  prot_annot <- prolfquapp::ProteinAnnotation$new(
    lfqdata ,
    m_annot,
    description = "description",
    cleaned_ids = "IDcolumn",
    full_id = "metabolite_feature_Id",
    exp_nr_children = "exp_children",
    pattern_contaminants = NULL,
    pattern_decoys = NULL
  )
  return(list(lfqdata = lfqdata , protein_annotation = prot_annot))
}
