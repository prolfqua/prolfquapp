# Orchestration helpers for CMD scripts
#
# Each function encapsulates the core logic of a CMD_*.R script,
# keeping optparse and file I/O in the CMD script itself.

# --- CMD_CONTRASTS helpers ----------------------------------------------------

#' Generate contrasts for a single-factor design
#'
#' Reads an annotation file, validates the control level exists,
#' and adds a CONTROL column (C for control, T for treatment).
#'
#' @param annotation_file path to annotation CSV/TSV/XLSX
#' @param control reference level name (e.g. "WT")
#' @param group group column name, or NULL to auto-detect
#' @return data.frame with CONTROL column added
#' @export
#' @examples
#' csv <- system.file("application/contrasts/scenario1_single_factor.csv",
#'   package = "prolfquapp")
#' result <- run_contrasts_single(csv, control = "WT")
#' table(result$group, result$CONTROL)
#'
run_contrasts_single <- function(annotation_file, control, group = NULL) {
  stopifnot(file.exists(annotation_file))

  res <- prolfquapp::read_annotation(annotation_file, QC = TRUE)
  annot <- res$annot
  group_col <- if (!is.null(group)) group else res$atable$factors[["G_"]]

  if (!control %in% annot[[group_col]]) {
    stop(
      "control '", control, "' not found in column '", group_col, "'.",
      call. = FALSE
    )
  }

  annot$CONTROL <- ifelse(annot[[group_col]] == control, "C", "T")
  annot
}

#' Generate contrasts for a two-factor design
#'
#' Reads an annotation file and adds ContrastName/Contrast columns
#' using \code{\link[prolfqua]{annotation_add_contrasts}}.
#'
#' @param annotation_file path to annotation CSV/TSV/XLSX
#' @param f1 primary factor column name
#' @param f2 secondary factor column name
#' @param interactions logical; include interaction contrasts? Default TRUE.
#' @return data.frame with ContrastName and Contrast columns added
#' @export
#' @examples
#' csv <- system.file("application/contrasts/scenario2_two_factor.csv",
#'   package = "prolfquapp")
#' result <- run_contrasts_twofactor(csv, f1 = "treatment", f2 = "time")
#' unique(result[!is.na(result$ContrastName), c("ContrastName", "Contrast")])
#'
run_contrasts_twofactor <- function(
  annotation_file,
  f1,
  f2,
  interactions = TRUE
) {
  stopifnot(file.exists(annotation_file))

  df <- prolfquapp::read_table_data(annotation_file)
  missing_cols <- setdiff(c(f1, f2), colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      "Column(s) not found: ", paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  res <- prolfqua::annotation_add_contrasts(
    df,
    primary_col = f1,
    secondary_col = f2,
    interactions = interactions
  )
  res$annot
}

# --- CMD_MAKE_YAML helper -----------------------------------------------------

#' Generate a DEA configuration list
#'
#' Creates a \code{\link{ProlfquAppConfig}} object, converts it to a
#' plain list, and reorders fields so that verbose/internal sections
#' appear at the bottom of the YAML output.
#'
#' @param project project ID
#' @param order order ID
#' @param workunit workunit ID
#' @param norm normalization method (e.g. "vsn", "none", "robscale")
#' @param model contrast facade method (see
#'   \code{names(prolfqua::FACADE_REGISTRY)})
#' @param outdir optional output directory; if it exists, stored in the
#'   config so downstream scripts know where to write results
#' @return named list suitable for \code{yaml::write_yaml()}
#' @export
#' @examples
#' cfg <- run_make_yaml(project = "p100", workunit = "WU123")
#' cfg$project_spec$workunit_Id
#'
run_make_yaml <- function(
  project = "",
  order = "",
  workunit = "",
  norm = "vsn",
  model = "lm_missing",
  outdir = NULL
) {
  GRP2 <- prolfquapp::make_DEA_config_R6(
    PROJECTID = project,
    ORDERID = order,
    WORKUNITID = workunit,
    Normalization = norm,
    model = model
  )
  GRP2$set_zipdir_name()
  if (!is.null(outdir) && dir.exists(outdir)) {
    GRP2$path <- outdir
  }
  cfg <- GRP2$as_list()

  # Move verbose/internal fields to bottom for readability
  fields_to_move <- c("ext_reader", "group")
  main <- cfg[!names(cfg) %in% fields_to_move]
  bottom <- cfg[names(cfg) %in% fields_to_move]
  c(main, bottom)
}

# --- CMD_QUANT_QC helper ------------------------------------------------------

#' Preprocess quantification data for QC reporting
#'
#' Resolves configuration (from YAML or defaults), reads the annotation,
#' and runs software-specific preprocessing. Returns everything needed
#' to construct a \code{\link{QC_generator}}.
#'
#' @param indir directory containing quantification output files
#' @param dataset path to annotation CSV/TSV/XLSX
#' @param software software key (e.g. "DIANN", "SIM"); looked up in
#'   \code{\link{prolfqua_preprocess_functions}}
#' @param yaml_file optional path to config YAML; if it exists, used
#'   instead of generating a default config
#' @param outdir output directory (used only when no yaml_file)
#' @param project project ID (used only when no yaml_file)
#' @param order order ID (used only when no yaml_file)
#' @param workunit workunit ID (used only when no yaml_file)
#' @return list with \code{xd} (preprocessed data), \code{files}
#'   (discovered file paths), and \code{config} (ProlfquAppConfig)
#' @export
run_qc_preprocess <- function(
  indir,
  dataset,
  software,
  yaml_file = NULL,
  outdir = "qc_dir",
  project = "",
  order = "",
  workunit = ""
) {
  if (!is.null(yaml_file) && file.exists(yaml_file)) {
    GRP2 <- prolfquapp::get_config(yaml_file)
  } else {
    GRP2 <- prolfquapp::make_DEA_config_R6(
      PATH = outdir,
      ORDERID = order,
      PROJECTID = project,
      WORKUNITID = workunit,
      application = software,
      prefix = "QC"
    )
  }

  if (!file.exists(dataset)) {
    stop("No annotation file found: ", dataset, call. = FALSE)
  }

  annotation <- prolfquapp::read_table_data(dataset) |>
    prolfquapp::read_annotation(QC = TRUE, repeated = FALSE)

  procsoft <- prolfquapp::preprocess_software(
    indir,
    annotation,
    preprocess_functions =
      prolfquapp::prolfqua_preprocess_functions[[software]],
    pattern_contaminants = GRP2$processing_options$pattern_contaminants,
    pattern_decoys = GRP2$processing_options$pattern_decoys
  )

  xd <- procsoft$xd
  xd$lfqdata$set_config_value("hierarchy_depth", 1)

  list(xd = xd, files = procsoft$files, config = GRP2)
}

# --- CMD_DEA_V2 helper -------------------------------------------------------

#' Run differential expression analysis pipeline
#'
#' Reads annotation, preprocesses quantification data, then runs the
#' full DEA pipeline: aggregation, transformation, model fitting, and
#' contrast computation. Returns the \code{\link{DEAnalyse}} object and
#' supporting data needed for report generation.
#'
#' @param indir directory containing quantification output files
#' @param dataset path to annotation CSV/TSV/XLSX
#' @param software software key as returned by
#'   \code{\link{get_procfuncs}} (e.g. "prolfquapp.DIANN")
#' @param config a \code{\link{ProlfquAppConfig}} object (already
#'   merged with CLI overrides via \code{\link{sync_opt_config}})
#' @return list with \code{deanalyse} (DEAnalyse), \code{xd}
#'   (preprocessed data), \code{annotation} (parsed annotation
#'   including contrasts), and \code{files} (discovered file paths)
#' @export
run_dea <- function(indir, dataset, software, config) {
  if (!file.exists(dataset)) {
    stop("Annotation file not found: ", dataset, call. = FALSE)
  }

  annotation <- prolfquapp::read_table_data(dataset) |>
    prolfquapp::read_annotation(prefix = config$group)

  pfuncs <- prolfquapp::get_procfuncs()
  if (!software %in% names(pfuncs)) {
    stop(
      "Software '", software, "' not found. Available: ",
      paste(names(pfuncs), collapse = ", "),
      call. = FALSE
    )
  }

  procsoft <- prolfquapp::preprocess_software(
    indir,
    annotation,
    preprocess_functions = pfuncs[[software]],
    pattern_contaminants =
      config$processing_options$pattern_contaminants,
    pattern_decoys = config$processing_options$pattern_decoys
  )

  xd <- procsoft$xd

  data_prep <- prolfquapp::ProteinDataPrep$new(
    xd$lfqdata, xd$protein_annotation, config
  )
  data_prep$cont_decoy_summary()
  data_prep$remove_cont_decoy()
  data_prep$aggregate()
  data_prep$transform_data()

  deanalyse <- data_prep$build_deanalyse(annotation$contrasts)
  deanalyse$build_default()
  deanalyse$get_annotated_contrasts()

  list(
    deanalyse = deanalyse,
    xd = xd,
    annotation = annotation,
    files = procsoft$files
  )
}
