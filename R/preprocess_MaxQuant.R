#' Methods for reading  MaxQuant outputs
#'
#' Convert MaxQuant outputs into tidy tables. For more details see functions listed in the see also section.
#'
#' @family MaxQuant
#' @name MaxQuant
NULL

# Read a MaxQuant txt file (or `filename` from a zip archive) unless `x` is already a data.frame;
# column names are lower-cased.
.read_MQ_txt <- function(x, filename) {
  if (is.character(x)) {
    x <- read.csv(
      if (grepl("\\.zip$", tolower(x))) unz(x, filename) else x,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE
    )
  }
  colnames(x) <- tolower(colnames(x))
  x
}

# Gather the per raw file columns starting with `prefix` into long format.
.gather_MQ <- function(data, id, prefix, value) {
  data |>
    dplyr::select(!!id := "id", dplyr::starts_with(prefix)) |>
    tidyr::gather(key = "raw.file", value = !!value, dplyr::starts_with(prefix)) |>
    dplyr::mutate(raw.file = gsub(prefix, "", .data$raw.file))
}

#' extract intensities and annotations from MQ proteinGroups.txt
#' @export
#' @keywords internal
#' @param MQProteinGroups data.frame generated with read.csv("peptide.txt",sep="\\t", stringsAsFactors=FALSE)
#' @family MaxQuant
#' @examples
#' protein_txt <- prolfqua::find_package_file("prolfquapp","samples/maxquant_txt/tiny2.zip")
#' protein_txt <- read.csv(
#'   unz(protein_txt, "proteinGroups.txt"),
#'   header = TRUE, stringsAsFactors = FALSE, sep = "\t")
#' mq_proteins <-tidyMQ_ProteinGroups(protein_txt)
tidyMQ_ProteinGroups <- function(MQProteinGroups) {
  MQProteinGroups <- .read_MQ_txt(MQProteinGroups, "proteinGroups.txt")
  meta <- dplyr::select(
    MQProteinGroups,
    "protein.ids" = "protein.ids",
    "majority.protein.ids" = "majority.protein.ids",
    "nr.peptides" = "peptides",
    "fasta.headers",
    "protein.group.id" = "id",
    "protein.score" = one_of("score")
  )
  meta <- meta |> dplyr::mutate(proteinID = gsub(";.*$", "", .data$majority.protein.ids))

  pint <- .gather_MQ(MQProteinGroups, "protein.group.id", "intensity.", "mq.protein.intensity")
  plfq <- .gather_MQ(MQProteinGroups, "protein.group.id", "lfq.intensity.", "mq.protein.lfq.intensity")
  pcount <- .gather_MQ(MQProteinGroups, "protein.group.id", "ms.ms.count.", "mq.protein.ms.ms.count")
  pint <- dplyr::inner_join(pint, plfq, by = c("protein.group.id", "raw.file")) |>
    dplyr::inner_join(pcount, by = c("protein.group.id", "raw.file"))
  return(dplyr::inner_join(meta, pint, by = "protein.group.id"))
}

#' read evidence file
#' @param Evidence MQ evidence file or zip archive with evidence file
#' @export
#' @keywords internal
#' @family MaxQuant
#' @examples
#' evidence_txt <- prolfqua::find_package_file("prolfquapp", "samples/maxquant_txt/tiny2.zip")
#' evidence_txt <- read.csv(
#'   unz(evidence_txt, "evidence.txt"),
#'   header = TRUE, stringsAsFactors = FALSE, sep = "\t")
#' mq_evidence <- tidyMQ_Evidence(evidence_txt)
tidyMQ_Evidence <- function(Evidence) {
  Evidence <- .read_MQ_txt(Evidence, "evidence.txt")
  res <- dplyr::select(
    Evidence,
    "evidence.id" = "id",
    "peptide.id",
    "raw.file",
    "protein.group.id" = "protein.group.ids",
    "mod.peptide.id" = "mod..peptide.id",
    "leading.razor.protein",
    "evidence.score" = "score",
    "delta.score",
    "pep",
    "calibrated.retention.time",
    "retention.time",
    "retention.length",
    "charge",
    "mass",
    "ms.ms.count",
    "ms.ms.scan.number",
    "evidence.intensity" = "intensity",
    "modifications",
    "modified.sequence",
    "missed.cleavages",
    "reverse"
  )
  res <- res |>
    dplyr::mutate(
      reverse = .data$reverse %in% "+",
      raw.file = tolower(.data$raw.file),
      proteotypic = !grepl(";", .data$protein.group.id)
    )
  return(tidyr::separate_rows(res, .data$protein.group.id, sep = ";", convert = TRUE))
}

#' parse MQ peptides.txt
#' @param MQPeptides data.frame generated with read.csv("peptide.txt",sep = "\\t", stringsAsFactors = FALSE)
#' @family MaxQuant
#' @export
#' @keywords internal
#' @examples
#' peptide_txt <- prolfqua::find_package_file("prolfquapp", "samples/maxquant_txt/tiny2.zip")
#' peptides_txt <- read.csv(
#'   unz(peptide_txt, "peptides.txt"),
#'   header = TRUE, stringsAsFactors = FALSE, sep = "\t")
#' mq_peptides <- tidyMQ_Peptides(peptides_txt)
tidyMQ_Peptides <- function(MQPeptides, proteotypic_only = TRUE) {
  MQPeptides <- .read_MQ_txt(MQPeptides, "peptides.txt")
  meta <- dplyr::select(
    MQPeptides,
    "peptide.id" = "id",
    "sequence",
    "proteins",
    "leading.razor.protein",
    "protein.group.id" = "protein.group.ids",
    "peptide.score" = "score",
    "pep",
    "ms.ms.count" = "ms.ms.count",
    dplyr::one_of("missed.cleavages"),
    "unique.groups" = "unique..groups.",
    "potential.contaminant" = ends_with("contaminant"),
    "reverse" = "reverse"
  ) |>
    dplyr::mutate(
      potential.contaminant = dplyr::case_when(
        .data$potential.contaminant == "" ~ FALSE,
        .data$potential.contaminant == "+" ~ TRUE
      ),
      unique.groups = dplyr::case_when(.data$unique.groups == "yes" ~ TRUE, .data$unique.groups == "no" ~ FALSE),
      reverse = dplyr::case_when(.data$reverse == "+" ~ TRUE, .data$reverse == "" ~ FALSE)
    )

  PepIntensities <- .gather_MQ(MQPeptides, "peptide.id", "intensity.", "peptide.intensity")
  if (any(startsWith(colnames(MQPeptides), "identification.type."))) {
    # if only one file no id type is provided
    PepIDType <- .gather_MQ(MQPeptides, "peptide.id", "identification.type.", "id.type")
    PepIntensities <- dplyr::inner_join(PepIntensities, PepIDType, by = c("peptide.id", "raw.file"))
  } else {
    PepIntensities$id.type <- "By MS/MS"
  }

  xx <- dplyr::inner_join(meta, PepIntensities, by = "peptide.id")
  xx$proteotypic <- !grepl(";", xx$protein.group.id)
  xx <- xx |>
    tidyr::separate_rows(.data$protein.group.id, sep = ";", convert = TRUE) |>
    dplyr::mutate(proteins = dplyr::if_else(.data$proteins == "", .data$leading.razor.protein, .data$proteins))
  xx$isotope <- "light"
  if (proteotypic_only) {
    xx <- xx |> dplyr::filter(.data$proteotypic == TRUE)
  }
  return(xx)
}

#' get petpide.txt and fasta file location in folder
#' @param path path to data directory
#' @return list with paths to data and fasta
#' @export
get_MQ_peptide_files <- function(path) {
  files <- dir(path = path, recursive = TRUE, full.names = TRUE)
  data <- grep("peptides.txt", files, value = TRUE)
  fasta.files <- grep("*.fasta$", files, ignore.case = TRUE, value = TRUE)
  if (any(grepl("database[0-9]*.fasta$", fasta.files))) {
    fasta.files <- grep("database[0-9]*.fasta$", fasta.files, value = TRUE)
  }
  if (length(fasta.files) == 0) {
    logger::log_error("No fasta file found!")
    stop()
  }
  if (length(data) == 0) {
    logger::log_error(
      "No peptides.txt file found in '",
      path,
      "'. MaxQuant preprocessor requires peptides.txt (not proteinGroups.txt)."
    )
    stop()
  }
  return(list(data = data, fasta = fasta.files))
}

#' create template dataset for MAXQUANT data
#' @param files list with data and fasta file paths
#' @family MaxQuant
#' @export
dataset_template_MAXQUANT <- function(files) {
  peptide <- prolfquapp::tidyMQ_Peptides(files$data, proteotypic_only = TRUE)
  return(data.frame(raw.file = unique(peptide$raw.file), name = NA, group = NA, subject = NA, CONTROL = NA))
}

#' preprocess MQ peptide
#' @param quant_data path to peptides.txt file
#' @param fasta_file path to fasta file(s)
#' @param annotation annotation list from read_annotation
#' @param pattern_contaminants regex pattern for contaminants
#' @param pattern_decoys regex pattern for decoys
#' @param hierarchy_depth hierarchy depth for aggregation
#' @param nr_peptides minimum number of distinct (stripped) peptides per protein (>= 1, default 1)
#' @return A list containing the prepared \code{LFQData} and
#'   \code{ProteinAnnotation} objects.
#' @export
preprocess_MQ_peptide <- function(
  quant_data,
  fasta_file,
  annotation,
  pattern_contaminants = "^zz|^CON|Cont_",
  pattern_decoys = "^REV_|^rev_",
  hierarchy_depth = 1,
  nr_peptides = 1
) {
  annot <- annotation$annot
  config <- annotation$atable$clone(deep = TRUE)
  annot[[config$file_name]] <- tolower(gsub("^x|\\.d\\.zip$|\\.raw$", "", basename(annot[[config$file_name]])))

  peptide <- prolfquapp::tidyMQ_Peptides(quant_data, proteotypic_only = TRUE)
  nrPeptides_exp <- peptide |>
    dplyr::distinct(leading.razor.protein, sequence) |>
    dplyr::count(leading.razor.protein, name = "nrPeptides")
  peptide <- prolfquapp::filter_by_peptide_count(peptide, "leading.razor.protein", "sequence", nr_peptides)
  annotated <- annot[[config$file_name]]
  found <- sort(unique(peptide$raw.file))
  nr <- sum(annotated %in% found)
  logger::log_info("nr : ", nr, " files annotated out of ", length(found))
  stopifnot(nr > 0)
  missing <- paste(setdiff(annotated, found), collapse = " ; ")
  logger::log_info("channels in annotation which are not in peptide.txt file : ", missing)
  extra <- paste(setdiff(found, annotated), collapse = " ; ")
  logger::log_info("channels in peptide.txt which are not in annotation file : ", extra)

  peptide$qValue <- 1 - peptide$pep
  config$ident_score <- "pep"
  config$ident_q_value <- "qValue"
  config$hierarchy[["protein_Id"]] <- "leading.razor.protein"
  config$hierarchy[["peptide_Id"]] <- "sequence"
  config$set_response("peptide.intensity")
  config$hierarchy_depth <- hierarchy_depth

  apeptide <- dplyr::inner_join(annot, peptide, multiple = "all", by = stats::setNames("raw.file", config$file_name))
  lfqdata <- prolfqua::LFQData$new(prolfqua::setup_analysis(apeptide, config), config)

  fasta_annot <- get_annot_from_fasta(fasta_file)
  fasta_annot <- dplyr::left_join(nrPeptides_exp, fasta_annot, by = c(leading.razor.protein = "fasta.id")) |>
    dplyr::rename(!!lfqdata$relevant_hierarchy_keys()[1] := "leading.razor.protein", description = "fasta.header")
  prot_annot <- prolfquapp::ProteinAnnotation$new(
    lfqdata,
    fasta_annot,
    description = "description",
    cleaned_ids = "proteinname",
    full_id = "protein_Id",
    exp_nr_children = "nrPeptides",
    pattern_contaminants = pattern_contaminants,
    pattern_decoys = pattern_decoys
  )
  lfqdata$remove_small_intensities()
  return(list(lfqdata = lfqdata, protein_annotation = prot_annot))
}
