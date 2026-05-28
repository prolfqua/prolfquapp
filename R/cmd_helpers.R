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
      "control '",
      control,
      "' not found in column '",
      group_col,
      "'.",
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
      "Column(s) not found: ",
      paste(missing_cols, collapse = ", "),
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
#'   \code{names(prolfqua::FACADE_REGISTRY)}) or "saint"
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
    preprocess_functions = prolfquapp::prolfqua_preprocess_functions[[
      software
    ]],
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

  default_model <- .resolve_facade_model(
    config$processing_options$model,
    config$processing_options$model_missing
  )
  model_entry <- prolfqua::lookup_facade(default_model)
  if (is.null(model_entry)) {
    stop("Unknown facade: ", default_model, call. = FALSE)
  }
  is_nested <- identical(model_entry$needs, "nested")
  saint_annot <- isTRUE(model_entry$needs_saint_annotation)

  is_peptide_reader <- grepl("_PEPTIDE$", software)
  if (is_nested && !is_peptide_reader) {
    stop(
      "Facade '",
      default_model,
      "' (needs=\"nested\") requires a peptide-level reader ",
      "(e.g. DIANN_PEPTIDE). Got '",
      software,
      "'.",
      call. = FALSE
    )
  }

  annotation <- prolfquapp::read_table_data(dataset) |>
    prolfquapp::read_annotation(prefix = config$group, SAINT = saint_annot)

  pfuncs <- prolfquapp::get_procfuncs()
  if (!software %in% names(pfuncs)) {
    stop(
      "Software '",
      software,
      "' not found. Available: ",
      paste(names(pfuncs), collapse = ", "),
      call. = FALSE
    )
  }

  procsoft <- prolfquapp::preprocess_software(
    indir,
    annotation,
    preprocess_functions = pfuncs[[software]],
    pattern_contaminants = config$processing_options$pattern_contaminants,
    pattern_decoys = config$processing_options$pattern_decoys
  )

  xd <- procsoft$xd

  data_prep <- prolfquapp::ProteinDataPrep$new(
    xd$lfqdata,
    xd$protein_annotation,
    config
  )
  data_prep$cont_decoy_summary()
  data_prep$remove_cont_decoy()

  if (is_nested) {
    lfq_peptide <- data_prep$lfq_data_peptide$get_copy()
    lfq_peptide$set_config_value("hierarchy_depth", 1)
    lfq_raw <- lfq_peptide
    lfq_model <- prolfquapp::transform_lfqdata(
      lfq_peptide,
      method = config$processing_options$transform
    )
    report_prep <- prolfquapp::ProteinDataPrep$new(
      lfq_peptide$get_copy(),
      data_prep$rowAnnot,
      config
    )
    report_prep$aggregate()
    report_prep$transform_data()
    lfq_raw$rename_response("abundance")
    lfq_model$rename_response("normalized_abundance")

    deanalyse <- DEAnalysePeptideToProtein$new(
      lfq_data = lfq_model,
      rowAnnot = data_prep$rowAnnot,
      prolfq_app_config = config,
      contrasts = annotation$contrasts,
      default_model = default_model,
      lfq_data_raw = lfq_raw,
      summary = data_prep$summary
    )
    deanalyse$build_default()
    deanalyse$get_annotated_contrasts()
    deanalyse$lfq_data <- report_prep$lfq_data_transformed
    deanalyse$lfq_data_raw <- report_prep$lfq_data
  } else {
    data_prep$aggregate()
    data_prep$transform_data()
    deanalyse <- data_prep$build_deanalyse(annotation$contrasts)
    deanalyse$build_default()
    deanalyse$get_annotated_contrasts()
  }

  list(
    deanalyse = deanalyse,
    xd = xd,
    annotation = annotation,
    files = procsoft$files
  )
}

write_dea_run_outputs <- function(result, config, opt, ymlfile) {
  deanalyse <- result$deanalyse
  xd <- result$xd
  annotation <- result$annotation
  files <- result$files

  logger::log_info("Processing done: ", opt$software)
  logger::log_info(paste(
    c(
      "Protein Annotation :\n",
      capture.output(print(xd$protein_annotation$get_summary()))
    ),
    collapse = "\n"
  ))
  logger::log_info(
    "ContrastNames: \n",
    paste(names(annotation$contrasts), collapse = "\n")
  )
  logger::log_info("END OF ANALYSIS")

  logger::log_info("CREATING DEAReportGenerator")
  reporter <- DEAReportGenerator$new(deanalyse, config, name = "")
  logger::log_info("Writing results to: ", config$get_zipdir())

  outdir <- reporter$write_DEA_all(
    boxplot = FALSE,
    markdown = "Grp2Analysis_V2_R6.Rmd",
    markdown_qc = "DiffExpQC_R6.Rmd"
  )

  arrow::write_parquet(
    deanalyse$lfq_data$data_long(),
    sink = file.path(config$get_result_dir(), "lfqdata_normalized.parquet")
  )
  cfg <- prolfqua::R6_extract_values(deanalyse$lfq_data$get_config())
  yaml::write_yaml(cfg, file.path(config$get_result_dir(), "lfqdata.yaml"))

  lfqdataIB <- xd$lfqdata$get_subset(xd$protein_annotation$clean(
    contaminants = config$processing_options$remove_cont,
    decoys = config$processing_options$remove_decoys
  ))

  ibaq_file <- file.path(
    reporter$resultdir,
    paste0("IBAQ_", opt$workunit, ".xlsx")
  )
  if (length(xd$lfqdata$relevant_hierarchy_keys()) == 1) {
    ibaq <- compute_IBAQ_values(lfqdataIB, xd$protein_annotation)
    writexl::write_xlsx(
      ibaq$data_wide()$data,
      path = ibaq_file
    )
  }
  outdir$data_files$ibaq_file <- ibaq_file

  tryCatch(
    {
      logger::log_info("Writing summarized experiment.")
      se_file <- file.path(reporter$resultdir, "SummarizedExperiment.rds")
      SE <- reporter$make_SummarizedExperiment()
      saveRDS(SE, file = se_file)

      logger::log_info("Rendering SE Quarto DEA report.")
      outdir$quarto_file <- render_quarto_se_report(
        se_file = se_file,
        output_dir = reporter$resultdir,
        output_file = paste0(reporter$fname, "_quarto.html")
      )
    },
    error = function(e) {
      logger::log_warn(
        "Skipping SummarizedExperiment/Quarto export: ",
        conditionMessage(e)
      )
    }
  )

  prolfquapp::write_index_html(outdir, result_dir = reporter$ZIPDIR)

  logger::log_info(
    "Creating directory with input files :",
    config$get_input_dir()
  )
  dir.create(config$get_input_dir())

  prolfquapp::copy_DEA_R6_Files(workdir = config$get_input_dir())
  prolfquapp::copy_shell_script(workdir = config$get_input_dir())

  file.copy(
    c(files$data, files$fasta, ymlfile, opt$dataset),
    config$get_input_dir()
  )

  logger::log_info(
    "Write yaml with parameters: ",
    file.path(config$get_input_dir(), "minimal.yaml")
  )

  yaml::write_yaml(
    prolfqua::R6_extract_values(config),
    file = file.path(config$get_input_dir(), "minimal.yaml")
  )

  invisible(outdir)
}

#' Run differential expression analysis for CompoundDiscoverer ZIP exports
#'
#' Reads the embedded prolfqua sample annotation and long feature table from a
#' CompoundDiscoverer ZIP export, then runs the same DEA preparation and model
#' pipeline as \code{\link{run_dea}}.
#'
#' @param input ZIP file path or directory containing one ZIP export
#' @param config a \code{\link{ProlfquAppConfig}} object
#' @param files optional file list from \code{\link{get_CD_export_files}}
#' @param subset_column optional long-table subset column to filter features
#' @return list with \code{deanalyse}, \code{xd}, \code{annotation}, and
#'   \code{files}
#' @export
run_dea_cd <- function(
  input = NULL,
  config,
  files = NULL,
  subset_column = NULL
) {
  if (is.null(files)) {
    if (is.null(input)) {
      stop("Either input or files must be supplied.", call. = FALSE)
    }
    files <- prolfquapp::get_CD_export_files(input)
  }

  cd <- prolfquapp::preprocess_CD_export(
    long_file = files$data,
    sample_file = files$samples,
    config = config,
    subset_column = subset_column
  )
  xd <- list(
    lfqdata = cd$lfqdata,
    protein_annotation = cd$protein_annotation
  )
  annotation <- cd$annotation

  data_prep <- prolfquapp::ProteinDataPrep$new(
    xd$lfqdata,
    xd$protein_annotation,
    config
  )
  data_prep$cont_decoy_summary()
  data_prep$remove_cont_decoy()
  if (
    length(xd$lfqdata$hierarchy_keys()) ==
      xd$lfqdata$get_config()$hierarchy_depth
  ) {
    data_prep$lfq_data <- data_prep$lfq_data_peptide
  } else {
    data_prep$aggregate()
  }
  data_prep$transform_data()

  deanalyse <- data_prep$build_deanalyse(annotation$contrasts)
  deanalyse$build_default()
  deanalyse$get_annotated_contrasts()

  list(
    deanalyse = deanalyse,
    xd = xd,
    annotation = annotation,
    files = files
  )
}
