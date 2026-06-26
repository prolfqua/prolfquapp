get_nr_pep <- function(report) {
  nrPEP <- report |>
    dplyr::select(all_of(c("Protein.Group", "Stripped.Sequence"))) |>
    dplyr::distinct() |>
    dplyr::group_by(!!sym("Protein.Group")) |>
    dplyr::summarize(nrPeptides = dplyr::n())

  return(nrPEP)
}


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
  clean_keys <- function(x) {
    sort(unique(as.character(x[!is.na(x)])))
  }
  annotation_keys <- clean_keys(annotation_keys)
  quant_keys <- clean_keys(quant_keys)
  matched_keys <- clean_keys(matched_keys)

  annotation_missing <- setdiff(annotation_keys, matched_keys)
  quant_only <- setdiff(quant_keys, annotation_keys)

  if (length(annotation_missing) > 0) {
    logger::log_warn(
      "{context}: annotated files not found in quantification data: ",
      paste(annotation_missing, collapse = " ; ")
    )
  }

  invisible(list(
    annotation_missing = annotation_missing,
    quant_only = quant_only
  ))
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
  filter_PG <- function(PG, .Lib.PG.Q.Value = 0.01, .PG.Q.Value = 0.05) {
    PG <- PG |> dplyr::filter(.data$Lib.PG.Q.Value < .Lib.PG.Q.Value)
    PG <- PG |> dplyr::filter(.data$PG.Q.Value < .PG.Q.Value)
    return(PG)
  }

  report <- data
  report2 <- filter_PG(
    report,
    .Lib.PG.Q.Value = Lib.PG.Q.Value,
    .PG.Q.Value = PG.Q.Value
  )
  # DIA-NN 2.x main report carries a bare `Run` column (no `File.Name`);
  # DIA-NN 1.x carries `File.Name` (full path). Prefer `Run` and fall back to
  # `File.Name`. The `Run` value is already the bare basename, so basename(),
  # the backslash gsub, and the extension strip are no-ops on it.
  run_col <- if ("Run" %in% names(report2)) {
    "Run"
  } else if ("File.Name" %in% names(report2)) {
    "File.Name"
  } else {
    stop("DIA-NN report has neither 'Run' nor 'File.Name'")
  }
  report2$raw.file <- .normalize_raw_file(report2[[run_col]])
  report2$Protein.Group <- sub("zz\\|(.+)\\|.+", "\\1", report2$Protein.Group)
  return(report2)
}


#' Create peptide level (stripped sequences) report by aggregating Precursor abundances.
#'
#' \code{\link{diann_read_output}}
#' @param report2 filtered DIA-NN report data frame
#' @export
#'
diann_output_to_peptide <- function(report2) {
  pg_quantity_col <- if ("PG.Quantity" %in% names(report2)) {
    "PG.Quantity"
  } else if ("PG.MaxLFQ" %in% names(report2)) {
    "PG.MaxLFQ"
  } else {
    stop("No protein group quantity column found")
  }
  peptide <- report2 |>
    dplyr::group_by(
      !!!syms(c(
        "raw.file",
        "Protein.Group",
        "Protein.Names",
        pg_quantity_col,
        "Stripped.Sequence"
      ))
    ) |>
    dplyr::summarize(
      Peptide.Quantity = sum(.data$Precursor.Quantity, na.rm = TRUE),
      Peptide.Normalised = sum(.data$Precursor.Normalised, na.rm = TRUE),
      PEP = min(.data$PEP, na.rm = TRUE),
      nr_children = n(),
      .groups = "drop"
    )
  return(peptide)
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
  # DIA-NN PTM "site_report.parquet" files also end in "report.parquet" and
  # would otherwise be matched above; drop them so only the main report remains.
  diann.path <- diann.path[!grepl("site_report\\.(parquet|tsv)$", diann.path)]
  # prefer the native DIA-NN 2.x parquet when both a parquet and a tsv are found
  parquet.path <- grep("\\.parquet$", diann.path, value = TRUE)
  if (length(parquet.path) > 0) {
    diann.path <- parquet.path
  }
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
  return(list(data = diann.path, fasta = fasta.files))
}


#' read a DIA-NN report from parquet (native 2.x) or tsv (legacy)
#' @param path path to a DIA-NN report (.parquet or .tsv)
#' @return a data frame / tibble
#' @noRd
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
#'
#' annotation <- file.path("inst/application/DIANN/2706527/dataset.csv") |>
#'   readr::read_csv() |>
#'   prolfquapp::read_annotation(QC = TRUE)
#' x$fasta
#' undebug(preprocess_DIANN)
#' xd <- preprocess_DIANN(x$data, x$fasta, annotation)
#' xd$lfqdata$hierarchy_counts()
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
  annot <- annotation$annot
  config <- annotation$atable$clone(deep = TRUE)
  annot <- annot |>
    dplyr::mutate(
      raw.file = .normalize_raw_file(annot[[config$file_name]])
    )
  data <- read_diann_report(quant_data)
  report2 <- prolfquapp::diann_read_output(
    data,
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
  nrPEP <- get_nr_pep(report2)
  nrPEP$Protein.Group.2 <- sapply(nrPEP$Protein.Group, function(x) {
    unlist(strsplit(x, "[ ;]"))[1]
  })

  peptide <- prolfquapp::diann_output_to_peptide(report2)
  peptide$qValue <- peptide$PEP
  nr <- sum(annot$raw.file %in% sort(unique(peptide$raw.file)))
  logger::log_info(
    "nr : ",
    nr,
    " files annotated out of ",
    length(unique(peptide$raw.file))
  )
  if (nr == 0) {
    stop(
      "No files are annotated. The annotation file is not compatible withe quant data."
    )
  }

  config$file_name <- "raw.file"
  config$nr_children <- "nr_children"
  config$ident_q_value <- "qValue"
  config$hierarchy[["protein_Id"]] <- c("Protein.Group")
  config$hierarchy[["peptide_Id"]] <- c("Stripped.Sequence")
  config$set_response("Peptide.Quantity")
  config$hierarchy_depth <- hierarchy_depth

  if (nr_peptides > 1) {
    nrPEP <- nrPEP |> dplyr::filter(nrPeptides >= nr_peptides)
    peptide <- peptide[peptide$Protein.Group %in% nrPEP$Protein.Group, ]
  }

  annotation_keys <- unique(annot$raw.file)
  quant_keys <- unique(peptide$raw.file)
  peptide <- dplyr::inner_join(annot, peptide, multiple = "all")
  .diagnose_sample_join(
    annotation_keys = annotation_keys,
    quant_keys = quant_keys,
    matched_keys = peptide$raw.file,
    context = "DIA-NN"
  )
  adata <- prolfqua::setup_analysis(peptide, config)
  lfqdata <- prolfqua::LFQData$new(adata, config)
  lfqdata$remove_small_intensities()

  # build protein annotation
  logger::log_info("start reading fasta.")
  fasta_annot <- get_annot_from_fasta(
    fasta_file,
    pattern_decoys = pattern_decoys,
    isUniprot = TRUE
  )
  logger::log_info("reading fasta done, creating protein annotation.")
  prot_annot <- dplyr::left_join(
    nrPEP,
    fasta_annot,
    by = c(Protein.Group.2 = "proteinname")
  )
  prot_annot <- dplyr::rename(
    prot_annot,
    IDcolumn = "Protein.Group.2",
    description = "fasta.header",
    protein_Id = "Protein.Group"
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
  return(list(lfqdata = lfqdata, protein_annotation = protAnnot))
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
  datasetannot <- data.frame(
    raw.file = unique(xx$raw.file),
    Name = NA,
    Group = NA,
    Subject = NA,
    Control = NA
  )
  return(datasetannot)
}
