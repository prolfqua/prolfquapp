#' Methods for reading Fragpipe outputs
#'
#' Convert FragPipe outputs into tidy tables. For more details see functions listed in the see also section.
#'
#' @family FragPipe
#' @name FragPipe
NULL

#'
#' read FragPipe generated MSstats formatted csv files.
#'
#' sanitize entries in the Bioreplicate and Condition columns
#'
#' @family FragPipe
#' @export
#' @param file MSstats formatted file
#' @keywords internal
tidy_FragPipe_MSstats_csv <- function(file){
  inputFile <- readr::read_csv(unz(file, filename = "MSstats.csv"))
  inputFile$BioReplicate <- paste("br", inputFile$BioReplicate, sep = "")
  inputFile$Condition <- make.names(inputFile$Condition)
  inputFile$pep <- 0
  return(inputFile)
}



#' FragPipe read FragPipe combined protein files up to Version 15
#'
#' @export
#' @param combprot path to combined_protein.tsv file
#' @param intnames intensity column prefix
#' @param protIDcol default protein.group
#' @param subgroup default subgroup
#'
#' @keywords internal
#'
#' @family FragPipe
#' @examples
#'
#' prottsv <- prolfqua::find_package_file("prolfquapp", "samples/FragPipe/combined_protein_small.tsv")
#'
#' prot <- tidy_FragPipe_combined_protein_deprec(prottsv)
#' stopifnot( dim(prot) ==c(19980,27))
tidy_FragPipe_combined_protein_deprec <- function(
    combprot, intnames = c("total.intensity",
                           "unique.intensity",
                           "razor.intensity",

                           "total.ion.count",
                           "unique.ion.count",
                           "razor.ion.count",

                           "total.spectral.count",
                           "unique.spectral.count",
                           "razor.spectral.count"),
    protIDcol = "protein.group",
    subgroup = "subgroup",
    as_list = FALSE) {
  if (is.character(combprot) && file.exists(combprot)) {
    Cprotein <- tibble::as_tibble(read.csv(combprot,
                                           header = TRUE, sep = "\t", stringsAsFactors = FALSE))

  } else if ("tbl_df" %in% class(combprot)) {
    Cprotein <- combprot
  } else {
    stop(class(combprot), " not supported.")
  }

  ### start processing
  colnames(Cprotein) <- tolower(colnames(Cprotein))
  cnam <- tolower(colnames(Cprotein))
  cnam <- cnam[1:which(cnam == "summarized.razor.spectral.count")]
  message("annotation columns : ", paste(cnam, collapse = "\n"))

  annot <- Cprotein |> dplyr::select(all_of(cnam))

  extractDataLong <- function(Cprotein, what = "total.intensity"){
    gg <- Cprotein |> dplyr::select( protIDcol, subgroup, dplyr::ends_with(what))
    gg <- gg |> tidyr::pivot_longer(cols = dplyr::ends_with(what), names_to = "raw.file",values_to = what)
    gg <- gg |> dplyr::mutate(raw.file = gsub(paste0(".",what,"$"),"", .data$raw.file))
    gg
  }

  res <- vector( mode = "list", length = length(intnames))
  names(res)  <- intnames

  for (i in seq_along(intnames)) {
    res[[intnames[i]]] <- extractDataLong(Cprotein, what = intnames[i] )
  }
  if (as_list) {
    return( res )
  }

  merged <- Reduce(dplyr::inner_join, res)
  merged <- dplyr::inner_join(annot, merged)


  return(merged)
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
    spcnames = c("Total Spectral Count",
                 "Unique Spectral Count",
                 "Razor Spectral Count"),
    intnames = c("Total Intensity",
                 "Unique Intensity",
                 "Razor Intensity"),
    maxlfqnames = c("MaxLFQ Total Intensity",
                    "MaxLFQ Unique Intensity",
                    "MaxLFQ Razor Intensity")
) {
  protIDcol = "Protein"
  if (is.character(combprot) && file.exists(combprot)) {
    Cprotein <- tibble::as_tibble(
      read.csv(combprot,
               header = TRUE,
               sep = "\t",
               stringsAsFactors = FALSE,
               check.names = FALSE))

  } else if ("tbl_df" %in% class(combprot)) {
    Cprotein <- combprot
  } else {
    stop(class(combprot), " not supported.")
  }

  cnam <- gsub(
    "Total Razor ", "Total ",
    gsub(
      "Unique Razor ","Unique ",
      gsub(
        " Intensity$"," Razor Intensity",
        gsub(
          " Spectral Count$"," Razor Spectral Count",colnames(Cprotein))
      )
    )
  )
  ### start processing
  cnam <- cnam
  colnames(Cprotein) <- cnam
  cnam <- cnam[1:which(cnam == "Combined Total Spectral Count")]

  message("annotation columns : ", paste(cnam, collapse = "\n"))
  annot <- Cprotein |> dplyr::select(all_of(cnam))
  colnames(Cprotein)

  extractDataLong <- function(Cprotein, what = "Total Intensity", butNot = NULL){
    cols <- colnames(Cprotein)
    cols <- setdiff( grep(paste0(what,"$"), cols, value = TRUE) , if (is.null(butNot)) {NULL} else { grep(butNot, cols, value = TRUE) })
    gg <- Cprotein |> dplyr::select( dplyr::all_of(protIDcol), dplyr::all_of(cols) )

    gg <- gg |> tidyr::pivot_longer(cols = dplyr::ends_with(what), names_to = "raw.file",values_to = what)
    gg <- gg |> dplyr::mutate(raw.file = gsub(paste0(".",what,"$"),"", .data$raw.file))
    gg
  }

  res <- vector( mode = "list", length = length(c(intnames, spcnames)))
  names(res)  <- c(intnames, spcnames)

  for (i in seq_along(c(intnames, spcnames))) {
    message("DD: ", c(intnames, spcnames)[i] )
    res[[c(intnames, spcnames)[i]]] <- extractDataLong(Cprotein, what = c(intnames, spcnames)[i], butNot = "maxlfq" )
  }

  if (sum(grepl(".MaxLFQ.", colnames(Cprotein))) > 0) {
    res_maxlfq <- vector( mode = "list", length(maxlfqnames))
    names(res_maxlfq)  <- maxlfqnames
    for (i in seq_along(maxlfqnames)) {
      message("DD: ", maxlfqnames[i] )
      res_maxlfq[[maxlfqnames[i] ]] <-  extractDataLong(Cprotein, what = maxlfqnames[i], butNot = NULL )
    }
    res <- c(res, res_maxlfq)
  }

  if (as_list) {
    return( res )
  }

  sql_inner_join <- function(x, y){
    dplyr::inner_join(x,y, multiple = "all")
  }
  merged <- Reduce( sql_inner_join , res )
  merged <- dplyr::inner_join( annot, merged , multiple = "all")
  colnames(merged) <- tolower(make.names(colnames(merged)))
  return( merged )
}


#' read psm.tsv produced by FragPipe and convert into long format
#' @export
#' @param psm_file path to psm.tsv file
#' @param purity_threshold purity threshold default = 0.5
#' @param PeptideProphetProb default 0.9
#' @param column_before_quants describes the last column before the quantitative values (this is not consistent with in different versions of FP, default "Quan Usage"
#' @param aggregate aggregate spectra to psm level
tidy_FragPipe_psm <- function(psm_file,
                              purity_threshold = 0.5,
                              PeptideProphetProb = 0.9,
                              abundance_threshold = 0,
                              column_before_quants =  c("Quan Usage" , "Mapped Proteins"),
                              aggregate = TRUE){


  psm <- readr::read_tsv(psm_file)
  column_before_quants <- intersect(colnames(psm), c("Quan Usage" , "Mapped Proteins"))
  column_before_quants <- tail(column_before_quants, n=1)
  if (!"Purity" %in% colnames(psm) ) {
    warning("no Purity column in psm file!")
    psm <- psm |> dplyr::mutate(Purity = 1, .before = column_before_quants)
  }
  x <- which(colnames(psm) == column_before_quants)
  colnamesQuan <- colnames(psm)[(x + 1):ncol(psm)]
  probability_column <- intersect(c("PeptideProphet Probability", "Probability"), colnames(psm))

  psm_relevant <- psm |> dplyr::select(
    dplyr::all_of(
      c(c("Spectrum",
          "Spectrum File",
          "Peptide",
          "Modified Peptide",
          "Charge",
          "Intensity",
          "Purity",
          "Protein",
          "Protein Description",
          Probability = probability_column,
          "Protein Description",
          "Retention",
          "Calibrated Observed Mass",
          "Assigned Modifications",
          "Charge"),
        colnamesQuan) ))

  psm_long <- psm_relevant |> tidyr::pivot_longer( tidyselect::all_of(colnamesQuan), values_to = "abundance", names_to = "channel")
  if (!is.null(abundance_threshold)) {
    psm_long <- dplyr::filter(psm_long, abundance > abundance_threshold)
  }

  nrPeptides_exp <- psm_long |>
    dplyr::select(Protein, Peptide) |>
    dplyr::distinct() |>
    dplyr::group_by(Protein) |>
    dplyr::summarize(nrPeptides = dplyr::n())

  colnames(psm_long) <- make.names(colnames(psm_long))
  psm_long <- dplyr::filter(psm_long, Purity > purity_threshold & Probability > PeptideProphetProb)

  if (aggregate) {
    psm_long <- psm_long |>
      dplyr::select(-all_of(c("Spectrum.File","Spectrum","Intensity","Purity","Retention","Calibrated.Observed.Mass","Charge"))) |>
      dplyr::group_by(dplyr::across(-c(abundance, Probability))) |>
      dplyr::summarize(nr_psm = n(), abundance = sum(abundance, na.rm = TRUE), Probability = max(Probability, na.rm = TRUE))
  }

  return(list(data = psm_long, nrPeptides_exp = nrPeptides_exp))
}



#' get psm.tsv and fasta file location in folder
#' @return list with paths to data and fasta
#' @export
get_FP_PSM_files <- function(path){

  psm_file <- dir(path = path, pattern = "psm.tsv", recursive = TRUE, full.names = TRUE)
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



#' preprocess FP psm, filter by purity_threshold and PeptideProphetProb
#' @return list with lfqdata and protein annotation
#' @export
preprocess_FP_PSM <- function(quant_data,
                              fasta_file,
                              annotation,
                              purity_threshold = 0.5,
                              PeptideProphetProb = 0.9,
                              column_before_quants = c("Quan Usage" , "Mapped Proteins"),
                              pattern_contaminants = "^zz|^CON|Cont_",
                              pattern_decoys = "^REV_|^rev_"){
  annot <- annotation$annot
  atable <- annotation$atable
  annot <- annot |> dplyr::mutate(
    raw.file = gsub("^x|.d.zip$|.raw$","",
                    (basename(annot[[atable$fileName]]))
    ))

  psm <- prolfquapp::tidy_FragPipe_psm(quant_data, column_before_quants = column_before_quants)
  nrPeptides_exp <- psm$nrPeptides # this is a data.frame
  psm <- psm$data
  psm$qValue <- 1 - psm$Probability

  nr <- sum(annot[[annotation$atable$fileName]] %in% sort(unique(psm$channel)))
  logger::log_info("nr : ", nr, " files annotated out of ", length(unique(psm$channel)))
  stopifnot(nr > 0)
  logger::log_info("channels in annotation which are not in psm.tsv file : ", paste(setdiff(annot[[annotation$atable$fileName]],sort(unique(psm$channel))), collapse = " ; ") )
  logger::log_info("channels in psm.tsv which are not in annotation file : ", paste(setdiff(sort(unique(psm$channel)),annot[[annotation$atable$fileName]]), collapse = " ; ") )

  atable$ident_Score = "Probability"
  atable$ident_qValue = "qValue"
  atable$hierarchy[["protein_Id"]] <- c("Protein")
  atable$hierarchy[["peptide_Id"]] <- c("Peptide")
  atable$hierarchy[["mod_peptide_Id"]] <- c("Modified.Peptide","Assigned.Modifications")
  atable$set_response("abundance")

  bycol <- c("channel")
  names(bycol) <- atable$fileName
  psma <- dplyr::inner_join(annot, psm, multiple = "all", by = bycol)
  config <- prolfqua::AnalysisConfiguration$new(atable)
  adata <- prolfqua::setup_analysis(psma, config)
  lfqdata <- prolfqua::LFQData$new(adata, config)

  # build rowAnnotation.
  fasta_annot <- get_annot_from_fasta(fasta_file, pattern_decoys = pattern_decoys)
  fasta_annot <- dplyr::left_join(nrPeptides_exp, fasta_annot, by = c("Protein" = "fasta.id"))

  fasta_annot <- fasta_annot |> dplyr::rename(!!lfqdata$config$table$hierarchy_keys_depth()[1] := !!sym("Protein"))
  fasta_annot <- fasta_annot |> dplyr::rename(description = fasta.header)

  prot_annot <- prolfquapp::ProteinAnnotation$new(
    lfqdata ,
    fasta_annot,
    description = "description",
    cleaned_ids = "proteinname",
    full_id = "protein_Id",
    exp_nr_children = "nrPeptides",
    pattern_contaminants = pattern_contaminants,
    pattern_decoys = pattern_decoys
  )

  lfqdata$remove_small_intensities()
  return(list(lfqdata = lfqdata , protein_annotation = prot_annot))
}


#' get report.tsv and fasta file location in folder
#' @return list with paths to data and fasta
#' @export
get_FP_multiSite_files <- function(path){
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


