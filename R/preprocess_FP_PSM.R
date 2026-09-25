#' Methods for reading Fragpipe outputs
#'
#' Convert FragPipe outputs into tidy tables. For more details see functions listed in the see also section.
#'
#' @family FragPipe
#' @name FragPipe
NULL

.read_FP_combined_protein <- function(combprot, ...) {
  if (is.character(combprot) && file.exists(combprot)) {
    tibble::as_tibble(read.csv(combprot, header = TRUE, sep = "\t", stringsAsFactors = FALSE, ...))
  } else if ("tbl_df" %in% class(combprot)) {
    combprot
  } else {
    stop(class(combprot), " not supported.")
  }
}

#' FragPipe read FragPipe combined protein files up to Version 15
#'
#' @export
#' @param combprot path to combined_protein.tsv file
#' @param intnames intensity column prefix
#' @param protIDcol default protein.group
#' @param subgroup default subgroup
#' @keywords internal
#' @family FragPipe
#' @examples
#' prottsv <- prolfqua::find_package_file("prolfquapp", "samples/FragPipe/combined_protein_small.tsv")
#' prot <- tidy_FragPipe_combined_protein_deprec(prottsv)
#' stopifnot( dim(prot) ==c(19980,27))
tidy_FragPipe_combined_protein_deprec <- function(
  combprot,
  intnames = c(
    "total.intensity",
    "unique.intensity",
    "razor.intensity",

    "total.ion.count",
    "unique.ion.count",
    "razor.ion.count",

    "total.spectral.count",
    "unique.spectral.count",
    "razor.spectral.count"
  ),
  protIDcol = "protein.group",
  subgroup = "subgroup",
  as_list = FALSE
) {
  Cprotein <- .read_FP_combined_protein(combprot)
  colnames(Cprotein) <- tolower(colnames(Cprotein))
  cnam <- colnames(Cprotein)
  cnam <- cnam[1:which(cnam == "summarized.razor.spectral.count")]
  message("annotation columns : ", paste(cnam, collapse = "\n"))
  annot <- Cprotein |> dplyr::select(all_of(cnam))

  res <- lapply(stats::setNames(nm = intnames), function(what) {
    Cprotein |>
      dplyr::select(protIDcol, subgroup, dplyr::ends_with(what)) |>
      tidyr::pivot_longer(cols = dplyr::ends_with(what), names_to = "raw.file", values_to = what) |>
      dplyr::mutate(raw.file = gsub(paste0("\\.", what, "$"), "", .data$raw.file))
  })
  if (as_list) {
    return(res)
  }
  dplyr::inner_join(annot, Reduce(dplyr::inner_join, res))
}

#' read combined_protein.tsv file for FragPipe Version 16 or newer
#' @export
#' @param combprot path to combined_protein.tsv file
#' @param as_list return as list
#' @return tidy dataframe or list with df (e.g. total.spectral.count or total.intensity etc).
#' @keywords internal
#' @family FragPipe
tidy_FragPipe_combined_protein <- function(
  combprot,
  as_list = FALSE,
  spcnames = c(
    "Total Spectral Count",
    "Unique Spectral Count",
    "Razor Spectral Count"
  ),
  intnames = c("Total Intensity", "Unique Intensity", "Razor Intensity"),
  maxlfqnames = c(
    "MaxLFQ Total Intensity",
    "MaxLFQ Unique Intensity",
    "MaxLFQ Razor Intensity"
  )
) {
  Cprotein <- .read_FP_combined_protein(combprot, check.names = FALSE)
  cnam <- gsub(" Spectral Count$", " Razor Spectral Count", colnames(Cprotein))
  cnam <- gsub(" Intensity$", " Razor Intensity", cnam)
  cnam <- gsub("Unique Razor ", "Unique ", cnam)
  cnam <- gsub("Total Razor ", "Total ", cnam)
  colnames(Cprotein) <- cnam
  cnam <- cnam[1:which(cnam == "Combined Total Spectral Count")]
  message("annotation columns : ", paste(cnam, collapse = "\n"))
  annot <- Cprotein |> dplyr::select(all_of(cnam))

  extract_long <- function(what, butNot = NULL) {
    message("DD: ", what)
    cols <- grep(paste0(what, "$"), colnames(Cprotein), value = TRUE)
    cols <- setdiff(cols, if (!is.null(butNot)) grep(butNot, colnames(Cprotein), value = TRUE))
    Cprotein |>
      dplyr::select(dplyr::all_of(c("Protein", cols))) |>
      tidyr::pivot_longer(cols = dplyr::ends_with(what), names_to = "raw.file", values_to = what) |>
      dplyr::mutate(raw.file = gsub(paste0("\\.", what, "$"), "", .data$raw.file))
  }
  res <- lapply(stats::setNames(nm = c(intnames, spcnames)), extract_long, butNot = "maxlfq")
  if (any(grepl(".MaxLFQ.", colnames(Cprotein)))) {
    res <- c(res, lapply(stats::setNames(nm = maxlfqnames), extract_long))
  }
  if (as_list) {
    return(res)
  }

  merged <- Reduce(function(x, y) dplyr::inner_join(x, y, multiple = "all"), res)
  merged <- dplyr::inner_join(annot, merged, multiple = "all")
  colnames(merged) <- tolower(make.names(colnames(merged)))
  return(merged)
}

#' read psm.tsv produced by FragPipe and convert into long format
#' @export
#' @param psm_files path(s) to psm.tsv file(s)
#' @param purity_threshold purity threshold default = 0.5
#' @param PeptideProphetProb default 0.9
#' @param abundance_threshold minimum abundance threshold
#' @param column_before_quants describes the last column before the
#'   quantitative values (not consistent across FP versions),
#'   default "Quan Usage"
#' @param aggregate aggregate spectra to psm level
#' @return A list with the long-format PSM data and expected peptide counts.
tidy_FragPipe_psm <- function(
  psm_files,
  purity_threshold = 0.5,
  PeptideProphetProb = 0.9,
  abundance_threshold = 0,
  column_before_quants = c("Quan Usage", "Mapped Proteins"),
  aggregate = TRUE
) {
  psm_long <- list()
  for (psm_file in psm_files) {
    psm <- readr::read_tsv(psm_file)
    column_before_quants <- tail(intersect(colnames(psm), column_before_quants), n = 1)
    if (!"Purity" %in% colnames(psm)) {
      warning("no Purity column in psm file!")
      psm <- psm |> dplyr::mutate(Purity = 1, .before = column_before_quants)
    }
    x <- which(colnames(psm) == column_before_quants)
    colnamesQuan <- colnames(psm)[(x + 1):ncol(psm)]
    probability_column <- intersect(c("PeptideProphet Probability", "Probability"), colnames(psm))
    psm_long[[psm_file]] <- psm |>
      dplyr::select(dplyr::all_of(c(
        "Spectrum",
        "Spectrum File",
        "Peptide",
        "Modified Peptide",
        "Charge",
        "Intensity",
        "Purity",
        "Protein",
        "Protein Description",
        Probability = probability_column,
        "Retention",
        "Calibrated Observed Mass",
        "Assigned Modifications",
        colnamesQuan
      ))) |>
      tidyr::pivot_longer(tidyselect::all_of(colnamesQuan), values_to = "abundance", names_to = "channel")
  }
  psm_long <- dplyr::bind_rows(psm_long)
  if (!is.null(abundance_threshold)) {
    psm_long <- dplyr::filter(psm_long, abundance > abundance_threshold)
  }
  nrPeptides_exp <- psm_long |> dplyr::distinct(Protein, Peptide) |> dplyr::count(Protein, name = "nrPeptides")

  colnames(psm_long) <- make.names(colnames(psm_long))
  psm_long <- dplyr::filter(psm_long, Purity > purity_threshold & Probability > PeptideProphetProb)
  if (aggregate) {
    drop <- c("Spectrum.File", "Spectrum", "Intensity", "Purity", "Retention", "Calibrated.Observed.Mass", "Charge")
    psm_long <- psm_long |>
      dplyr::select(-all_of(drop)) |>
      dplyr::group_by(dplyr::across(-c(abundance, Probability))) |>
      dplyr::summarize(
        nr_psm = n(),
        abundance = sum(abundance, na.rm = TRUE),
        Probability = max(Probability, na.rm = TRUE)
      )
  }
  return(list(data = psm_long, nrPeptides_exp = nrPeptides_exp))
}

.get_FP_files <- function(path, pattern) {
  data <- dir(path = path, pattern = pattern, recursive = TRUE, full.names = TRUE)
  fasta.files <- grep("*.fasta$", dir(path = path, recursive = TRUE, full.names = TRUE), value = TRUE)
  if (any(grepl("database[0-9]*.fasta$", fasta.files))) {
    fasta.files <- grep("database[0-9]*.fasta$", fasta.files, value = TRUE)
  }
  if (length(fasta.files) == 0) {
    logger::log_error("No fasta file found!")
    stop()
  }
  return(list(data = data, fasta = fasta.files))
}

#' get psm.tsv and fasta file location in folder
#' @param path path to data directory
#' @return list with paths to data and fasta
#' @export
get_FP_PSM_files <- function(path) {
  .get_FP_files(path, "psm.tsv")
}

#' preprocess FP psm, filter by purity_threshold and PeptideProphetProb
#' @param quant_data path to quantification data file(s)
#' @param fasta_file path to fasta file(s)
#' @param annotation annotation list from read_annotation
#' @param pattern_contaminants regex pattern for contaminants
#' @param pattern_decoys regex pattern for decoys
#' @param purity_threshold purity threshold for filtering
#' @param PeptideProphetProb PeptideProphet probability threshold
#' @param hierarchy_depth hierarchy depth for aggregation
#' @param nr_peptides minimum number of distinct (stripped) peptides per protein (>= 1, default 1)
#' @param parse_fun function for parsing PSM files
#' @return list with lfqdata and protein annotation
#' @export
preprocess_FP_PSM <- function(
  quant_data,
  fasta_file,
  annotation,

  pattern_contaminants = "^zz|^CON|Cont_",
  pattern_decoys = "^REV_|^rev_",
  purity_threshold = 0.5,
  PeptideProphetProb = 0.9,
  hierarchy_depth = 1,
  nr_peptides = 1,
  parse_fun = tidy_FragPipe_psm
) {
  annot <- annotation$annot
  config <- annotation$atable$clone(deep = TRUE)
  annot$raw.file <- gsub("^x|\\.d\\.zip$|\\.raw$", "", basename(annot[[config$file_name]]))

  psm <- parse_fun(quant_data)
  nrPeptides_exp <- psm$nrPeptides
  psm <- psm$data
  psm$qValue <- 1 - psm$Probability
  # FragPipe `Peptide` is the stripped sequence (`Modified.Peptide` carries the mods).
  psm <- prolfquapp::filter_by_peptide_count(psm, "Protein", "Peptide", nr_peptides)
  annotated <- annot[[config$file_name]]
  found <- sort(unique(psm$channel))
  nr <- sum(annotated %in% found)
  logger::log_info("nr : ", nr, " files annotated out of ", length(found))
  stopifnot(nr > 0)
  missing <- paste(setdiff(annotated, found), collapse = " ; ")
  logger::log_info("channels in annotation which are not in psm.tsv file : ", missing)
  extra <- paste(setdiff(found, annotated), collapse = " ; ")
  logger::log_info("channels in psm.tsv which are not in annotation file : ", extra)

  config$ident_score <- "Probability"
  config$ident_q_value <- "qValue"
  config$hierarchy[["protein_Id"]] <- "Protein"
  config$hierarchy[["peptide_Id"]] <- "Peptide"
  config$hierarchy[["mod_peptide_Id"]] <- c("Modified.Peptide", "Assigned.Modifications")
  config$set_response("abundance")
  if ("nr_psm" %in% colnames(psm)) {
    config$nr_children <- "nr_psm"
  }
  config$hierarchy_depth <- hierarchy_depth

  psma <- dplyr::inner_join(annot, psm, multiple = "all", by = stats::setNames("channel", config$file_name))
  lfqdata <- prolfqua::LFQData$new(prolfqua::setup_analysis(psma, config), config)

  fasta_annot <- dplyr::left_join(nrPeptides_exp, get_annot_from_fasta(fasta_file), by = c(Protein = "fasta.id")) |>
    dplyr::rename(!!lfqdata$relevant_hierarchy_keys()[1] := "Protein", description = "fasta.header")
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

#' get dataset annotation template
#' @param files list with data and fasta file paths
#' @return data.frame
#' @export
dataset_template_FP_TMT <- function(files) {
  channel <- unique(prolfquapp::tidy_FragPipe_psm(files$data)$data$channel)
  return(data.frame(channel = channel, Name = channel, group = NA, subject = NA, CONTROL = NA))
}
