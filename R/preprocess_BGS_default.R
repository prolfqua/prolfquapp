#' get BGS and fasta file location in folder
#' @param file path to BGS report file
#' @return list with paths to data and fasta
#' @export
read_BGS <- function(
  file = "Experiment1_Report_BGS Factory Report (Normal).tsv"
) {
  bgs <- readr::read_tsv(file)
  colnames(bgs) <- colnames(bgs) |>
    stringr::str_replace_all("[[:space:]\\(\\)\\-]", "_") |>
    stringr::str_replace_all("_+", "_") |>
    stringr::str_replace("_$", "")

  ctoselect <- c(
    "R.FileName",
    "PG.ProteinGroups",
    "PG.ProteinAccessions",
    "PG.Qvalue",
    "PG.QValue_Run_Wise",
    "PG.Quantity",
    "PEP.GroupingKey",
    "PEP.IsProteotypic",
    "PEP.RunEvidenceCount",
    "PEP.NrOfMissedCleavages",
    "EG.ModifiedSequence",
    "EG.Qvalue",
    "FG.Qvalue",
    "FG.Charge",
    "FG.Quantity"
  )
  bgs[, ctoselect]
}


#' get BGS and fasta file location in folder
#' @param path path to data directory
#' @param bgs_pattern glob pattern for BGS report files
#' @return list with paths to data and fasta
#' @export
#' @examples
#' \dontrun{
#' x <- get_DIANN_files("inst/application/DIANN/2517219/")
#' }
get_BGS_files <- function(
  path,
  bgs_pattern = "*BGS Factory Report \\(Normal\\).tsv|_Report.tsv"
) {
  diann.path <- grep(
    bgs_pattern,
    dir(path = path, recursive = TRUE, full.names = TRUE),
    value = TRUE
  )
  list(data = diann.path, fasta = .get_fasta_files(path))
}

#' create templte dataset for BGS data
#' @param files list with data and fasta file paths
#' @return data.frame
#' @export
dataset_template_BGS <- function(files) {
  xt <- readr::read_tsv(files$data)
  ds <- xt |>
    dplyr::select(
      raw.file = "R.FileName",
      Group = "R.Condition",
      name = "R.Replicate"
    )
  ds <- ds |>
    tidyr::unite("Name", c("Group", "name"), remove = FALSE) |>
    dplyr::distinct()
  ds$name <- NULL
  return(ds)
}

#' preprocess DIANN ouput, filter by q_value and nr_peptides
#' @param quant_data path to quantification data file
#' @param fasta_file path to fasta file(s)
#' @param annotation annotation list from read_annotation
#' @param pattern_contaminants regex pattern for contaminants
#' @param pattern_decoys regex pattern for decoys
#' @param q_value q-value threshold for filtering
#' @param hierarchy_depth hierarchy depth for aggregation
#' @param nr_peptides minimum number of distinct (stripped) peptides per protein (>= 1, default 1)
#' @return list with lfqdata and protein annotation
#' @export
#' @examples
#' \dontrun{
#' x <- get_BGS_files("DefaultParsing")
#' bgs <- read_BGS(x$data)
#' annot <- data.frame(raw.file = bgs$R.FileName |> unique(),
#'  Name = paste(c(rep("A",3),rep("B",3)),1:6, sep="_"),
#' group = c(rep("A",3),rep("B",3)))
#' annotation <- annot |> prolfquapp::read_annotation(QC = TRUE)
#' #debug(preprocess_BGS)
#' xd <- preprocess_BGS(x$data, x$fasta, annotation)
#' }
preprocess_BGS <- function(
  quant_data,
  fasta_file,
  annotation,
  pattern_contaminants = "^zz|^CON|Cont_",
  pattern_decoys = "^REV_|^rev",
  q_value = 0.01,
  hierarchy_depth = 2,
  nr_peptides = 1
) {
  config <- annotation$atable$clone(deep = TRUE)
  annot <- annotation$annot |>
    dplyr::mutate(
      raw.file = gsub("^x|\\.d\\.zip$|\\.raw$", "", basename(.data[[config$file_name]]))
    )
  bgs <- read_BGS(quant_data)
  # `PEP.GroupingKey` is the Spectronaut peptide-level grouping key
  # (`EG.ModifiedSequence` carries the mods).
  report2 <- prolfquapp::filter_by_peptide_count(
    bgs,
    "PG.ProteinGroups",
    "PEP.GroupingKey",
    nr_peptides
  )
  .stop_if_unannotated(annot$raw.file, report2$R.FileName)

  config$file_name <- "raw.file"
  config$ident_q_value <- "FG.Qvalue"
  config$hierarchy[["protein_Id"]] <- c("PG.ProteinGroups")
  config$hierarchy[["peptide_Id"]] <- c("PEP.GroupingKey")
  config$hierarchy[["elution_group"]] <- c("EG.ModifiedSequence", "FG.Charge")
  config$set_response("FG.Quantity")
  config$hierarchy_depth <- hierarchy_depth

  report2 <- dplyr::inner_join(
    annot,
    report2,
    multiple = "all",
    by = c("raw.file" = "R.FileName")
  )
  adata <- prolfqua::setup_analysis(report2, config)
  lfqdata <- prolfqua::LFQData$new(adata, config)
  lfqdata$remove_small_intensities()

  protAnnot <- .uniprot_protein_annotation(
    lfqdata,
    bgs,
    "PG.ProteinGroups",
    "PEP.GroupingKey",
    fasta_file,
    pattern_contaminants,
    pattern_decoys
  )
  list(lfqdata = lfqdata, protein_annotation = protAnnot)
}
