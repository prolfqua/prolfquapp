#' get petpide.txt and fasta file location in folder
#' @param path path to data directory
#' @return list with paths to data and fasta
#' @export
get_MSstats_files <- function(path) {
  files <- dir(path = path, recursive = TRUE, full.names = TRUE)
  msstats.path <- grep("msstats.*\\.(csv|tsv)$", files, value = TRUE, ignore.case = TRUE)
  fasta.files <- grep("*\\.fasta$|*\\.fas$", files, ignore.case = TRUE, value = TRUE)
  if (any(grepl("database[0-9]*.fasta$", fasta.files))) {
    fasta.files <- grep("database[0-9]*.fasta$", fasta.files, value = TRUE)
  }
  if (length(fasta.files) == 0) {
    logger::log_error("No fasta file found!")
    stop()
  }
  if (length(msstats.path) > 1) {
    logger::log_warn("more then 1 msstats.tsv file found :", length(msstats.path), ". Returning first.")
  }
  list(data = msstats.path[1], fasta = fasta.files)
}


#' read MSstats.csv files and rollup to ProteinSequence level.
#' @param file path to MSstats csv file
#' @return A tibble with peptide-level intensities and child counts.
#' @export
read_msstats <- function(file) {
  readr::read_csv(file) |>
    dplyr::group_by(dplyr::across(c("ProteinName", "PeptideSequence", "IsotopeLabelType", "Run"))) |>
    dplyr::summarise(nr_children = dplyr::n(), Intensity = sum(Intensity, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(Intensity = ifelse(Intensity < 1e-10, NA, Intensity))
}

#' create dataset template from MSStats data.
#' @param files list with data and fasta file paths
#' @export
#'
dataset_template_MSSTATS <- function(files) {
  datasetannot <- prolfquapp::read_table_data(files$data) |>
    dplyr::select(raw.file = "Run", "Group" = "Condition", "Subject" = "BioReplicate") |>
    dplyr::distinct()
  datasetannot$Control <- ""
  tidyr::unite(datasetannot, "Name", "Group", "Subject", sep = "_", remove = FALSE)
}

# Shared MSstats reader. `fasta_key` is the FASTA column matched to `ProteinName`:
# "proteinname" (FragPipe DIA, cleaned accessions) or "fasta.id" (full FASTA ids).
.preprocess_msstats <- function(
  quant_data,
  fasta_file,
  annotation,
  pattern_contaminants,
  pattern_decoys,
  hierarchy_depth,
  nr_peptides,
  fasta_key
) {
  config <- annotation$atable$clone(deep = TRUE)
  annot <- annotation$annot |>
    dplyr::mutate(
      !!config$file_name := gsub("^x|\\.d\\.zip$|\\.raw$", "", basename(.data[[config$file_name]]))
    )

  peptide <- read_msstats(quant_data)
  peptide$nr_peptides <- 1
  nrPeptides_exp <- peptide |>
    dplyr::distinct(dplyr::across(c("ProteinName", "PeptideSequence"))) |>
    dplyr::count(dplyr::across("ProteinName"), name = "nrPeptides")

  # MSstats `PeptideSequence` may carry modifications depending on the upstream
  # converter; if so the count over-counts (under-filters, never over-drops).
  peptide <- prolfquapp::filter_by_peptide_count(peptide, "ProteinName", "PeptideSequence", nr_peptides)
  .stop_if_unannotated(annot[[config$file_name]], peptide$Run)

  peptide$qValue <- 0
  config$ident_q_value <- "qValue"
  config$hierarchy[["protein_Id"]] <- c("ProteinName")
  config$hierarchy[["peptide_Id"]] <- c("PeptideSequence")
  config$nr_children <- "nrPeptides"
  config$set_response("Intensity")
  config$hierarchy_depth <- hierarchy_depth

  apeptide <- dplyr::inner_join(annot, peptide, multiple = "all", by = stats::setNames("Run", config$file_name))
  .diagnose_sample_join(annot[[config$file_name]], peptide$Run, apeptide[[config$file_name]], "MSstats")

  adata <- prolfqua::setup_analysis(apeptide, config)
  lfqdata <- prolfqua::LFQData$new(adata, config)
  logger::log_info("Start reading fasta: ", fasta_file)
  fasta_annot <- get_annot_from_fasta(fasta_file, pattern_decoys = pattern_decoys)
  logger::log_info("Finished reading fasta: ", fasta_file)

  fasta_annot <- nrPeptides_exp |>
    dplyr::left_join(fasta_annot, by = c("ProteinName" = fasta_key)) |>
    dplyr::rename(!!lfqdata$relevant_hierarchy_keys()[1] := "ProteinName", description = "fasta.header")
  fpdia <- fasta_key == "proteinname"
  prot_annot <- prolfquapp::ProteinAnnotation$new(
    lfqdata,
    fasta_annot,
    description = "description",
    cleaned_ids = if (fpdia) "protein_Id" else "proteinname",
    full_id = if (fpdia) "fasta.id" else "protein_Id",
    exp_nr_children = "nrPeptides",
    pattern_contaminants = pattern_contaminants,
    pattern_decoys = pattern_decoys
  )
  lfqdata$remove_small_intensities()
  list(lfqdata = lfqdata, protein_annotation = prot_annot)
}

#' preprocess MSstats fragpipe
#' @inheritParams preprocess_MSstats
#' @return A list containing the prepared \code{LFQData} and
#'   \code{ProteinAnnotation} objects.
#' @export
#'
preprocess_MSstats_FPDIA <- function(
  quant_data,
  fasta_file,
  annotation,
  pattern_contaminants = "",
  pattern_decoys = "",
  hierarchy_depth = 1,
  nr_peptides = 1
) {
  .preprocess_msstats(
    quant_data,
    fasta_file,
    annotation,
    pattern_contaminants,
    pattern_decoys,
    hierarchy_depth,
    nr_peptides,
    fasta_key = "proteinname"
  )
}


#' preprocess MSstats file coming from FragPipe
#' @param quant_data path to MSstats csv file
#' @param fasta_file path to fasta file(s)
#' @param annotation annotation list from read_annotation
#' @param pattern_contaminants regex pattern for contaminants
#' @param pattern_decoys regex pattern for decoys
#' @param hierarchy_depth hierarchy depth for aggregation
#' @param nr_peptides minimum number of distinct (stripped) peptides per protein (>= 1, default 1)
#' @return A list containing the prepared \code{LFQData} and
#'   \code{ProteinAnnotation} objects.
#' @export
#'
preprocess_MSstats <- function(
  quant_data,
  fasta_file,
  annotation,
  pattern_contaminants = "^zz|^CON|Cont_",
  pattern_decoys = "^REV_|^rev_",
  hierarchy_depth = 1,
  nr_peptides = 1
) {
  .preprocess_msstats(
    quant_data,
    fasta_file,
    annotation,
    pattern_contaminants,
    pattern_decoys,
    hierarchy_depth,
    nr_peptides,
    fasta_key = "fasta.id"
  )
}
