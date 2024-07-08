get_nr_pep <- function(report){
  nrPEP <- report |>
    dplyr::select(all_of(c("Protein.Group", "Stripped.Sequence"))) |>
    dplyr::distinct() |>
    dplyr::group_by(!!sym("Protein.Group")) |>
    dplyr::summarize(nrPeptides = dplyr::n())

  return(nrPEP)
}

#' DIANN helper functions
#' @param diann.path path to diann-output.tsv
#' @param fasta.file path to fasta file
#' @param nrPeptides peptide threshold
#' @param Q.Value q value threshold
#' @import data.table
#' @export
#' @examples
#'
#' x <- get_DIANN_files("inst/application/DIANN/2517219/")
#' xd <- read_DIANN_output(x$data, x$fasta)
#' debug(read_DIANN_output)
read_DIANN_output <- function(diann.path,
                              fasta.file,
                              Q.Value = 0.01,
                              isUniprot = TRUE,
                              rev = "REV_"
) {
  warning("read_DIANN_output is being DEPRECATED")
  data <- readr::read_tsv(diann.path)
  report2 <- prolfquapp::diann_read_output(data, Q.Value = Q.Value)
  if (nrow(report2) == 0) {
    return(NULL)
  }
  peptide <- prolfquapp::diann_output_to_peptide(report2)
  peptide$Protein.Group.2 <- sapply(peptide$Protein.Group, function(x){ unlist(strsplit(x, "[ ;]"))[1]} )
  # we need to add the fasta.header information.
  fasta_annot <- get_annot_from_fasta(fasta.file, rev = rev, isUniprot = isUniprot)
  message("Percent of Proteins with description:" ,mean(peptide$Protein.Group.2 %in% fasta_annot$proteinname) * 100)
  # add fasta headers.
  if (nrow(peptide) == 0) {
    return(NULL)
  }
  peptide <- dplyr::left_join(dtplyr::lazy_dt(peptide), dtplyr::lazy_dt(fasta_annot),
                              by = c( Protein.Group.2 = "proteinname")) |>
    dplyr::as_tibble()
  return(peptide)
}


#' read DiaNN diann-output.tsv file
#'
#' filter for 2 peptides per protein, and for Q.Value < 0.01 (default)
#' @import data.table
#' @export
#' @examples
#'
#' xx <- readr::read_tsv("inst/application/DIANN/2517219/out-2024-06-12/WU304486_report.tsv")
#' report2 <- prolfquapp::diann_read_output(xx)
#' report2$Protein.Group.2 <- sapply(report2$Protein.Group, function(x){ unlist(strsplit(x, "[ ;]"))[1]} )
#'
#'
diann_read_output <- function(data, Q.Value = 0.01){

  select_PG <- function(report){
    columns  <- c("File.Name",
                  "Protein.Group",
                  "Protein.Names",
                  "PG.Quantity",
                  "PG.MaxLFQ",
                  "PG.Q.Value",
                  "Global.PG.Q.Value",
                  "Lib.PG.Q.Value")
    columns %in% names(report)

    PG <- report |> dplyr::select( dplyr::all_of( columns )) |>
      dplyr::distinct() |>
      dplyr::ungroup()
    return(PG)
  }

  filter_PG <- function(PG,  Q.Value = 0.01){
    PG <- PG |> dplyr::filter(.data$Lib.PG.Q.Value < Q.Value)
    PG <- PG |> dplyr::filter(.data$PG.Q.Value < Q.Value)
  }

  report <- data
  PG <- select_PG(report)
  PG2 <- filter_PG(PG, Q.Value = Q.Value)
  PG2 <- PG2 |> dplyr::select(c("File.Name", "Protein.Group", "Protein.Names"))

  report2 <- dplyr::inner_join(dtplyr::lazy_dt(PG2), dtplyr::lazy_dt(report),
                               by = c("File.Name", "Protein.Group", "Protein.Names")) |>
    dplyr::as_tibble()

  report2$raw.file <- gsub("^x|.d.zip$|.d$|.raw$|.mzML$","",basename(gsub("\\\\","/",report2$File.Name)))
  report2$Protein.Group <- sub("zz\\|(.+)\\|.+", "\\1", report2$Protein.Group )
  return(report2)
}

#' Create peptide level (stripped sequences) report by aggregating Precursor abundances.
#'
#' \code{\link{diann_read_output}}
#' @export
#'
diann_output_to_peptide <- function(report2){
  peptide <- report2 |>
    dplyr::group_by(!!!syms(c("raw.file",
                              "Protein.Group",
                              "Protein.Names",
                              "PG.Quantity",
                              "Stripped.Sequence" ) )) |>
    dplyr::summarize(Peptide.Quantity = sum(.data$Precursor.Quantity, na.rm = TRUE),
                     Peptide.Normalised = sum(.data$Precursor.Normalised, na.rm = TRUE),
                     Peptide.Translated = sum(.data$Precursor.Translated, na.rm = TRUE),
                     Peptide.Ms1.Translated = sum(.data$Ms1.Translated, na.rm = TRUE),
                     PEP = min(.data$PEP, na.rm = TRUE),
                     nr_children = n()
                     ,.groups = "drop")
  return(peptide)
}


#' get report.tsv and fasta file location in folder
#' @return list with paths to data and fasta
#' @export
#' @examples
#' x <- get_DIANN_files("inst/application/DIANN/2517219/")
#'
get_DIANN_files <- function(path){
  diann.path <- grep("report\\.tsv$|diann-output\\.tsv", dir(path = path, recursive = TRUE, full.names = TRUE), value = TRUE)
  fasta.files <- grep("*.fasta$", dir(path = path, recursive = TRUE,full.names = TRUE), value = TRUE)
  if (any(grepl("database[0-9]*.fasta$", fasta.files))) {
    fasta.files <- grep("database[0-9]*.fasta$", fasta.files, value = TRUE)
  }
  if (length(fasta.files) == 0) {
    logger::log_error("No fasta file found!")
    stop()
  }
  return(list(data = diann.path, fasta = fasta.files))
}

#' preprocess DIANN ouput, filter by q_value and nr_peptides
#' @return list with lfqdata and protein annotation
#' @export
#' @examples
#'
#' x <- get_DIANN_files("inst/application/DIANN/2517219/")
#'
#' annotation <- file.path("inst/application/DIANN/2517219/dataset.csv") |>
#'  readr::read_csv() |> prolfquapp::read_annotation(QC = TRUE)
#'  x$fasta
#' #debug(preprocess_DIANN)
#' xd <- preprocess_DIANN(x$data, x$fasta, annotation)
#'
preprocess_DIANN <- function(quant_data,
                             fasta_file,
                             annotation,
                             q_value = 0.01,
                             pattern_contaminants = "^zz|^CON",
                             pattern_decoys = "REV_"){

  annot <- annotation$annot
  atable <- annotation$atable$clone(deep = FALSE)
  annot <- annot |> dplyr::mutate(
    raw.file = gsub("^x|.d.zip$|.raw$","",
                    (basename(annot[[atable$fileName]]))
    ))
  data <- readr::read_tsv(quant_data)
  report2 <- prolfquapp::diann_read_output(data, Q.Value = q_value)
  nrPEP <- get_nr_pep(report2)
  nrPEP$Protein.Group.2 <- sapply(nrPEP$Protein.Group, function(x){ unlist(strsplit(x, "[ ;]"))[1]} )

  peptide <- prolfquapp::diann_output_to_peptide(report2)
  peptide$qValue <- 1 - peptide$PEP
  nr <- sum(annot$raw.file %in% sort(unique(peptide$raw.file)))
  logger::log_info("nr : ", nr, " files annotated out of ", length(unique(peptide$raw.file)))
  if (nr == 0) { stop("No files are annotated. The annotation file is not compatible withe quant data.") }

  atable$fileName = "raw.file"
  atable$nr_children = "nr_children"
  atable$ident_qValue = "qValue"
  atable$hierarchy[["protein_Id"]] <- c("Protein.Group")
  atable$hierarchy[["peptide_Id"]] <- c("Stripped.Sequence")
  atable$set_response("Peptide.Quantity")
  atable$hierarchyDepth <- 1

  peptide <- dplyr::inner_join(annot, peptide, multiple = "all")
  config <- prolfqua::AnalysisConfiguration$new(atable)
  adata <- prolfqua::setup_analysis(peptide, config)
  lfqdata <- prolfqua::LFQData$new(adata, config)
  lfqdata$remove_small_intensities()

  # build protein annotation
  fasta_annot <- get_annot_from_fasta(fasta_file, rev = pattern_contaminants, isUniprot = TRUE)
  prot_annot <- dplyr::left_join(nrPEP, fasta_annot, by = c(Protein.Group.2 = "proteinname"))
  prot_annot <- dplyr::rename(prot_annot, IDcolumn = "Protein.Group.2",description = "fasta.header",protein_Id = "Protein.Group" )

  protAnnot <- prolfquapp::ProteinAnnotation$new(
    lfqdata , prot_annot, description = "description",
    cleaned_ids = "IDcolumn",
    full_id = "fasta.id",
    exp_nr_children = "nrPeptides",
    pattern_contaminants = pattern_contaminants,
    pattern_decoys = pattern_decoys
  )
  return(list(lfqdata = lfqdata , protein_annotation = protAnnot))
}


