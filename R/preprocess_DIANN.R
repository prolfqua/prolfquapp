.normalize_raw_file <- function(x) {
  gsub(
    "^x|\\.d\\.zip$|\\.d$|\\.raw$|\\.mzML$",
    "",
    basename(gsub("\\\\", "/", x))
  )
}


.diagnose_sample_join <- function(
  annotation_keys,
  quant_keys,
  matched_keys,
  context = "reader"
) {
  clean_keys <- function(x) sort(unique(as.character(x[!is.na(x)])))
  annotation_keys <- clean_keys(annotation_keys)
  annotation_missing <- setdiff(annotation_keys, clean_keys(matched_keys))
  if (length(annotation_missing) > 0) {
    logger::log_warn(
      "{context}: annotated files not found in quantification data: ",
      paste(annotation_missing, collapse = " ; ")
    )
  }
  invisible(list(
    annotation_missing = annotation_missing,
    quant_only = setdiff(clean_keys(quant_keys), annotation_keys)
  ))
}


# Stop unless at least one quantified file is annotated.
.stop_if_unannotated <- function(annotated, quantified) {
  quantified <- unique(quantified)
  nr <- sum(annotated %in% quantified)
  logger::log_info("nr : ", nr, " files annotated out of ", length(quantified))
  if (nr == 0) {
    stop(
      "No files are annotated. The annotation file is not compatible withe quant data."
    )
  }
}


# FASTA files in `path`, preferring database[0-9]*.fasta, without first-pass files.
.get_fasta_files <- function(path) {
  fasta.files <- grep(
    "*.fasta$|*.fas$",
    dir(path = path, recursive = TRUE, full.names = TRUE),
    value = TRUE
  )
  if (any(grepl("database[0-9]*.fasta$", fasta.files))) {
    fasta.files <- grep("database[0-9]*.fasta$", fasta.files, value = TRUE)
  }
  fasta.files <- fasta.files[!grepl("first-pass", fasta.files)]
  if (length(fasta.files) == 0) {
    logger::log_error("No fasta file found!")
    stop()
  }
  fasta.files
}


# ProteinAnnotation for protein groups led by a UniProt accession (DIA-NN,
# Spectronaut): peptides counted per group in `report`, joined to the FASTA.
.uniprot_protein_annotation <- function(
  lfqdata,
  report,
  protein_col,
  peptide_col,
  fasta_file,
  pattern_contaminants,
  pattern_decoys
) {
  nrPEP <- report |>
    dplyr::distinct(dplyr::across(dplyr::all_of(c(protein_col, peptide_col)))) |>
    dplyr::count(dplyr::across(dplyr::all_of(protein_col)), name = "nrPeptides")
  nrPEP$IDcolumn <- sub("[ ;].*", "", nrPEP[[protein_col]])

  logger::log_info("start reading fasta.")
  fasta_annot <- get_annot_from_fasta(
    fasta_file,
    pattern_decoys = pattern_decoys,
    isUniprot = TRUE
  )
  logger::log_info("reading fasta done, creating protein annotation.")
  prot_annot <- nrPEP |>
    dplyr::left_join(fasta_annot, by = c(IDcolumn = "proteinname")) |>
    dplyr::rename(
      description = "fasta.header",
      protein_Id = dplyr::all_of(protein_col)
    )
  protAnnot <- prolfquapp::ProteinAnnotation$new(
    lfqdata,
    prot_annot,
    description = "description",
    cleaned_ids = "IDcolumn",
    full_id = "fasta.id",
    exp_nr_children = "nrPeptides",
    pattern_contaminants = pattern_contaminants,
    pattern_decoys = pattern_decoys
  )
  logger::log_info("protein annotation done.")
  protAnnot
}


#' read DiaNN diann-output.tsv file
#'
#' filter for 2 peptides per protein, and for Q.Value < 0.01 (default)
#' @param data data frame of DIA-NN report
#' @param Lib.PG.Q.Value library protein group q-value threshold
#' @param PG.Q.Value protein group q-value threshold
#' @import data.table
#' @export
#' @examples
#' \dontrun{
#' xx <- readr::read_tsv("WU292720_report.tsv")
#' report2 <- prolfquapp::diann_read_output(xx)
#' nrow(report2)
#' }
#'
diann_read_output <- function(data, Lib.PG.Q.Value = 0.01, PG.Q.Value = 0.05) {
  report2 <- data |>
    dplyr::filter(
      .data$Lib.PG.Q.Value < !!Lib.PG.Q.Value,
      .data$PG.Q.Value < !!PG.Q.Value
    )
  # DIA-NN 2.x has a bare `Run` column, DIA-NN 1.x a full-path `File.Name`.
  run_col <- intersect(c("Run", "File.Name"), names(report2))[1]
  if (is.na(run_col)) {
    stop("DIA-NN report has neither 'Run' nor 'File.Name'")
  }
  report2$raw.file <- .normalize_raw_file(report2[[run_col]])
  report2$Protein.Group <- sub("zz\\|(.+)\\|.+", "\\1", report2$Protein.Group)
  report2
}


#' Create peptide level (stripped sequences) report by aggregating Precursor abundances.
#'
#' \code{\link{diann_read_output}}
#' @param report2 filtered DIA-NN report data frame
#' @return A peptide-level tibble with aggregated precursor abundances.
#' @export
#'
diann_output_to_peptide <- function(report2) {
  pg_quantity_col <- intersect(c("PG.Quantity", "PG.MaxLFQ"), names(report2))[1]
  if (is.na(pg_quantity_col)) {
    stop("No protein group quantity column found")
  }
  report2 |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(
      "raw.file",
      "Protein.Group",
      "Protein.Names",
      pg_quantity_col,
      "Stripped.Sequence"
    )))) |>
    dplyr::summarize(
      Peptide.Quantity = sum(.data$Precursor.Quantity, na.rm = TRUE),
      Peptide.Normalised = sum(.data$Precursor.Normalised, na.rm = TRUE),
      PEP = min(.data$PEP, na.rm = TRUE),
      nr_children = n(),
      .groups = "drop"
    )
}


#' get report.tsv and fasta file location in folder
#' @param path path to data directory
#' @return list with paths to data and fasta
#' @export
#' @examples
#' \dontrun{
#' x <- get_DIANN_files("inst/application/DIANN/2517219/")
#' }
get_DIANN_files <- function(path) {
  diann.path <- grep(
    "report\\.parquet$|report\\.tsv$|diann-output\\.tsv",
    dir(path = path, recursive = TRUE, full.names = TRUE),
    value = TRUE
  )
  # drop DIA-NN PTM site reports, which also end in "report.parquet"
  diann.path <- diann.path[!grepl("site_report\\.(parquet|tsv)$", diann.path)]
  # prefer the native DIA-NN 2.x parquet when both a parquet and a tsv are found
  parquet.path <- grep("\\.parquet$", diann.path, value = TRUE)
  if (length(parquet.path) > 0) {
    diann.path <- parquet.path
  }
  list(data = diann.path, fasta = .get_fasta_files(path))
}


# read a DIA-NN report from parquet (native 2.x) or tsv (legacy)
read_diann_report <- function(path) {
  if (grepl("\\.parquet$", path)) {
    arrow::read_parquet(path)
  } else {
    readr::read_tsv(path)
  }
}


#' preprocess DIANN ouput, filter by q_value and nr_peptides
#' @param quant_data path to quantification data file
#' @param fasta_file path to fasta file(s)
#' @param annotation annotation list from read_annotation
#' @param pattern_contaminants regex pattern for contaminants
#' @param pattern_decoys regex pattern for decoys
#' @param q_value q-value threshold for filtering
#' @param hierarchy_depth hierarchy depth for aggregation
#' @param nr_peptides minimum number of peptides per protein
#' @return list with lfqdata and protein annotation
#' @export
#' @examples
#' \dontrun{
#' x <- get_DIANN_files("inst/application/DIANN/2706527/")
#' annotation <- file.path("inst/application/DIANN/2706527/dataset.csv") |>
#'   readr::read_csv() |>
#'   prolfquapp::read_annotation(QC = TRUE)
#' xd <- preprocess_DIANN(x$data, x$fasta, annotation, nr_peptides = 2)
#' xd$lfqdata$hierarchy_counts()
#' }
preprocess_DIANN <- function(
  quant_data,
  fasta_file,
  annotation,
  pattern_contaminants = "^zz|^CON|Cont_",
  pattern_decoys = "^REV_|^rev",
  q_value = 0.01,
  hierarchy_depth = 1,
  nr_peptides = 1
) {
  config <- annotation$atable$clone(deep = TRUE)
  annot <- annotation$annot |>
    dplyr::mutate(raw.file = .normalize_raw_file(.data[[config$file_name]]))
  report2 <- prolfquapp::diann_read_output(
    read_diann_report(quant_data),
    Lib.PG.Q.Value = q_value,
    PG.Q.Value = q_value
  )
  if (nrow(report2) == 0) {
    stop(
      "DIA-NN report contains no rows after filtering at q_value = ",
      q_value,
      ". Check that the report has quantified precursors/protein groups ",
      "and that Lib.PG.Q.Value and PG.Q.Value pass the threshold.",
      call. = FALSE
    )
  }
  peptide <- prolfquapp::diann_output_to_peptide(report2)
  peptide$qValue <- peptide$PEP
  .stop_if_unannotated(annot$raw.file, peptide$raw.file)

  config$file_name <- "raw.file"
  config$nr_children <- "nr_children"
  config$ident_q_value <- "qValue"
  config$hierarchy[["protein_Id"]] <- c("Protein.Group")
  config$hierarchy[["peptide_Id"]] <- c("Stripped.Sequence")
  config$set_response("Peptide.Quantity")
  config$hierarchy_depth <- hierarchy_depth

  peptide <- prolfquapp::filter_by_peptide_count(
    peptide,
    "Protein.Group",
    "Stripped.Sequence",
    nr_peptides
  )
  quant_keys <- peptide$raw.file
  peptide <- dplyr::inner_join(annot, peptide, multiple = "all")
  .diagnose_sample_join(annot$raw.file, quant_keys, peptide$raw.file, "DIA-NN")
  adata <- prolfqua::setup_analysis(peptide, config)
  lfqdata <- prolfqua::LFQData$new(adata, config)
  lfqdata$remove_small_intensities()

  protAnnot <- .uniprot_protein_annotation(
    lfqdata,
    report2,
    "Protein.Group",
    "Stripped.Sequence",
    fasta_file,
    pattern_contaminants,
    pattern_decoys
  )
  list(lfqdata = lfqdata, protein_annotation = protAnnot)
}

#' create dataset template from DIANN outputs
#' @param files list with data and fasta file paths
#' @export
dataset_template_diann <- function(files) {
  data <- read_diann_report(files$data)
  logger::log_info("Files: ", files$data, " loaded. Starting filtering.")
  xx <- prolfquapp::diann_read_output(
    data,
    Lib.PG.Q.Value = 0.01,
    PG.Q.Value = 0.01
  )
  data.frame(
    raw.file = unique(xx$raw.file),
    Name = NA,
    Group = NA,
    Subject = NA,
    Control = NA
  )
}
