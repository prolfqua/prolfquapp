#' Find the apb-export file in a folder
#'
#' apb-export's \code{prolfqua} target writes one \code{.h5ad} file from the
#' output of any vendor APB supports. The file carries its protein annotation,
#' so no FASTA file is needed.
#' @param path path to data directory
#' @return list with \code{data}, the \code{.h5ad} file, and an empty \code{fasta}
#' @export
#' @examples
#' \dontrun{
#' x <- get_APB_files("data_dir/")
#' }
get_APB_files <- function(path) {
  h5ad <- grep("\\.h5ad$", dir(path = path, recursive = TRUE, full.names = TRUE), value = TRUE)
  if (length(h5ad) != 1) {
    stop("Expected one .h5ad file in '", path, "', found ", length(h5ad), ": ", paste(h5ad, collapse = ", "))
  }
  list(data = h5ad, fasta = character())
}


#' Preprocess an apb-export \code{.h5ad} file, filter by q-values and nr_peptides
#'
#' Reads the wide AnnData written by apb-export's \code{prolfqua} target, drops
#' every observation whose q-value exceeds the threshold of a layer the file
#' carries, and sums the precursors of each peptide, counting them in
#' \code{nr_children}. A file of peptides (MaxQuant) has no precursors to sum.
#' The analysis configuration (hierarchy, response, q-value and children
#' columns) comes from the file; samples and factors from \code{annotation}.
#' Protein annotation, with \code{protein_length} and \code{nr_tryptic_peptides}
#' when the file was exported with a FASTA, comes from the file's \code{var}; the
#' number of peptides per protein is counted after the q-value filter.
#' @param quant_data path to the \code{.h5ad} file
#' @param fasta_file unused; the file carries its protein annotation
#' @param annotation annotation list from read_annotation
#' @param pattern_contaminants regex pattern for contaminants
#' @param pattern_decoys regex pattern for decoys
#' @param q_values named list of thresholds, one per q-value layer; an
#'   observation is kept when its value is below the threshold. A threshold for
#'   a layer the file does not carry is not applied.
#' @param hierarchy_depth hierarchy depth for aggregation
#' @param nr_peptides minimum number of peptides per protein
#' @return list with lfqdata and protein annotation
#' @export
#' @examples
#' \dontrun{
#' x <- get_APB_files("data_dir/")
#' annotation <- readr::read_csv("data_dir/dataset.csv") |>
#'   prolfquapp::read_annotation(QC = TRUE)
#' xd <- preprocess_APB(x$data, x$fasta, annotation)
#' xd$lfqdata$hierarchy_counts()
#' }
preprocess_APB <- function(
  quant_data,
  fasta_file,
  annotation,
  pattern_contaminants = "^zz|^CON|Cont_",
  pattern_decoys = "^REV_|^rev",
  q_values = list(pg_qValue = 0.01, pg_qValue_experiment = 0.01),
  hierarchy_depth = 1,
  nr_peptides = 1
) {
  adata <- anndataR::read_h5ad(quant_data)
  validate_prolfquapp_anndata(adata)
  pmeta <- adata$uns[["prolfquapp"]]
  stored <- .anndata_analysis_configuration(pmeta$analysis_configuration)
  response <- stored$get_response()
  logger::log_info(
    "APB: ",
    pmeta$source_software,
    " export; layers present: ",
    paste(c(response, intersect(pmeta$layer_names, adata$layers_keys())), collapse = ", ")
  )

  long <- anndata_to_long(adata, stored)
  long <- .apb_filter_q_values(long[!is.na(long[[response]]), , drop = FALSE], q_values)
  if (nrow(long) == 0) {
    stop("APB file contains no observations below the q-value thresholds.", call. = FALSE)
  }
  # the file's hierarchy down to the peptide; its precursors are summed
  keys <- utils::head(stored$hierarchy_keys(), 2)
  file_name <- stored$file_name
  peptide <- long |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(file_name, keys)))) |>
    dplyr::summarize(
      dplyr::across(dplyr::all_of(response), sum),
      dplyr::across(dplyr::any_of(stored$ident_q_value), min),
      dplyr::across(dplyr::all_of(stored$nr_children), sum),
      .groups = "drop"
    )
  peptide[[file_name]] <- .normalize_raw_file(peptide[[file_name]])

  # the file's configuration; samples and factors come from the annotation
  config <- stored$clone(deep = TRUE)
  config$hierarchy <- stored$hierarchy[keys]
  config$hierarchy_depth <- hierarchy_depth
  atable <- annotation$atable
  config$sample_name <- atable$sample_name
  config$factors <- atable$factors
  config$factor_depth <- atable$factor_depth
  config$norm_value <- atable$norm_value
  annot <- annotation$annot
  annot[[file_name]] <- .normalize_raw_file(annot[[atable$file_name]])
  .stop_if_unannotated(annot[[file_name]], peptide[[file_name]])

  pa <- pmeta$protein_annotation
  nrPEP <- peptide |>
    dplyr::distinct(dplyr::across(dplyr::all_of(keys))) |>
    dplyr::count(dplyr::across(dplyr::all_of(keys[[1]])), name = pa$exp_nr_children)

  peptide <- prolfquapp::filter_by_peptide_count(peptide, keys[[1]], keys[[2]], nr_peptides)
  quant_keys <- peptide[[file_name]]
  peptide <- dplyr::inner_join(annot, peptide, by = file_name, multiple = "all")
  .diagnose_sample_join(annot[[file_name]], quant_keys, peptide[[file_name]], "APB")
  lfqdata <- prolfqua::LFQData$new(prolfqua::setup_analysis(peptide, config), config)
  lfqdata$remove_small_intensities()

  # var holds one row per feature; keep the columns describing its protein,
  # among them protein_length and nr_tryptic_peptides when exported with a FASTA
  var <- as.data.frame(adata$var)
  per_feature <- c(stored$hierarchy_keys()[-1], stored$nr_children, stored$ident_q_value)
  prot_annot <- dplyr::distinct(var[, setdiff(colnames(var), per_feature), drop = FALSE])
  rownames(prot_annot) <- NULL
  protAnnot <- prolfquapp::ProteinAnnotation$new(
    lfqdata,
    dplyr::left_join(nrPEP, prot_annot, by = keys[[1]]),
    description = pa$description,
    cleaned_ids = pa$cleaned_ids,
    full_id = pa$full_id,
    exp_nr_children = pa$exp_nr_children,
    pattern_contaminants = pattern_contaminants,
    pattern_decoys = pattern_decoys
  )
  list(lfqdata = lfqdata, protein_annotation = protAnnot)
}

# Keeps the observations below the threshold of each q-value layer present.
.apb_filter_q_values <- function(long, q_values) {
  for (layer in names(q_values)) {
    threshold <- q_values[[layer]]
    if (!layer %in% colnames(long)) {
      logger::log_info("APB: no ", layer, " layer; its threshold ", threshold, " is not applied.")
      next
    }
    before <- nrow(long)
    long <- long[!is.na(long[[layer]]) & long[[layer]] < threshold, , drop = FALSE]
    logger::log_info("APB: ", layer, " < ", threshold, " keeps ", nrow(long), " of ", before, " observations.")
  }
  long
}

#' create dataset template from an apb-export file
#' @param files list with the \code{.h5ad} file in \code{data}
#' @export
dataset_template_APB <- function(files) {
  adata <- anndataR::read_h5ad(files$data)
  file_name <- adata$uns[["prolfquapp"]]$analysis_configuration$file_name
  data.frame(
    raw.file = as.data.frame(adata$obs)[[file_name]],
    Name = NA,
    Group = NA,
    Subject = NA,
    Control = NA
  )
}
