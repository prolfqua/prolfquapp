# Orchestration helpers for CMD scripts: the core logic of each CMD_*.R script,
# while optparse and file I/O stay in the CMD script itself.

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
  group_col <- group %||% res$atable$factors[["G_"]]
  if (!control %in% annot[[group_col]]) {
    stop("control '", control, "' not found in column '", group_col, "'.", call. = FALSE)
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
run_contrasts_twofactor <- function(annotation_file, f1, f2, interactions = TRUE) {
  stopifnot(file.exists(annotation_file))
  df <- prolfquapp::read_table_data(annotation_file)
  missing_cols <- setdiff(c(f1, f2), colnames(df))
  if (length(missing_cols) > 0) {
    stop("Column(s) not found: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  prolfqua::annotation_add_contrasts(df, primary_col = f1, secondary_col = f2, interactions = interactions)$annot
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
#' @param nr_peptides minimum distinct peptides per protein (>= 1, default 1)
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
  model = "lm_impute",
  nr_peptides = 1,
  outdir = NULL
) {
  GRP2 <- prolfquapp::make_DEA_config_R6(
    PROJECTID = project,
    ORDERID = order,
    WORKUNITID = workunit,
    Normalization = norm,
    model = model,
    nr_peptides = nr_peptides
  )
  GRP2$set_zipdir_name()
  if (!is.null(outdir) && dir.exists(outdir)) {
    GRP2$path <- outdir
  }
  cfg <- GRP2$as_list()
  # Move verbose/internal fields to bottom for readability
  bottom <- names(cfg) %in% c("ext_reader", "group")
  c(cfg[!bottom], cfg[bottom])
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
#' @param flat_outdir when TRUE, write QC outputs directly into \code{outdir}
#'   without a dated subdir (default FALSE)
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
  workunit = "",
  flat_outdir = FALSE
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
  GRP2$flat_outdir <- isTRUE(flat_outdir)
  if (!file.exists(dataset)) {
    stop("No annotation file found: ", dataset, call. = FALSE)
  }
  annotation <- prolfquapp::read_table_data(dataset) |>
    prolfquapp::read_annotation(QC = TRUE, repeated = FALSE)
  procsoft <- prolfquapp::preprocess_software(
    indir,
    annotation,
    preprocess_functions = prolfquapp::prolfqua_preprocess_functions[[software]],
    pattern_contaminants = GRP2$processing_options$pattern_contaminants,
    pattern_decoys = GRP2$processing_options$pattern_decoys
  )
  procsoft$xd$lfqdata$set_config_value("hierarchy_depth", 1)
  list(xd = procsoft$xd, files = procsoft$files, config = GRP2)
}

# --- CMD_DEA_V2 helper -------------------------------------------------------

# Nested facades (e.g. firth_nested, lmer_nested) fit peptide-level data. Paired
# with a protein-level reader, switch to its registered "<reader>_PEPTIDE"
# counterpart (same reader, other hierarchy_depth) instead of failing.
.resolve_nested_reader <- function(software, is_nested, available, facade = "") {
  if (!is_nested || grepl("_PEPTIDE$", software)) {
    return(software)
  }
  peptide_software <- paste0(software, "_PEPTIDE")
  if (!peptide_software %in% available) {
    stop(
      "Facade '",
      facade,
      "' (needs=\"nested\") requires a peptide-level reader ",
      "(e.g. DIANN_PEPTIDE), but no peptide-level counterpart '",
      peptide_software,
      "' is registered for '",
      software,
      "'.",
      call. = FALSE
    )
  }
  logger::log_info(
    "Facade '",
    facade,
    "' (needs=\"nested\") requires a peptide-level reader; ",
    "switching software from '",
    software,
    "' to '",
    peptide_software,
    "'."
  )
  peptide_software
}

# Fit a nested facade on peptide-level data; the report data are the aggregated
# protein-level data.
.build_nested_deanalyse <- function(data_prep, config, contrasts, default_model) {
  lfq_peptide <- data_prep$lfq_data_peptide$get_copy()
  lfq_peptide$set_config_value("hierarchy_depth", 1)
  lfq_model <- prolfquapp::transform_lfqdata(lfq_peptide, method = config$processing_options$transform)
  report_prep <- prolfquapp::ProteinDataPrep$new(lfq_peptide$get_copy(), data_prep$rowAnnot, config)
  report_prep$aggregate()
  report_prep$transform_data()
  lfq_peptide$rename_response("abundance")
  lfq_model$rename_response("normalized_abundance")

  deanalyse <- DEAnalysePeptideToProtein$new(
    lfq_data = lfq_model,
    rowAnnot = data_prep$rowAnnot,
    prolfq_app_config = config,
    contrasts = contrasts,
    default_model = default_model,
    lfq_data_raw = lfq_peptide,
    summary = data_prep$summary
  )
  deanalyse$build_default()
  deanalyse$get_annotated_contrasts()
  deanalyse$lfq_data <- report_prep$lfq_data_transformed
  deanalyse$lfq_data_raw <- report_prep$lfq_data
  deanalyse
}

# Transform the prepared data and fit the default model.
.build_protein_deanalyse <- function(data_prep, contrasts) {
  data_prep$transform_data()
  deanalyse <- data_prep$build_deanalyse(contrasts)
  deanalyse$build_default()
  deanalyse$get_annotated_contrasts()
  deanalyse
}

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
#'   including contrasts), \code{files} (discovered file paths),
#'   \code{software} (resolved software key), and
#'   \code{requested_software} (software key requested by the caller)
#' @export
#' @examples
#' dataset <- system.file(
#'   "application/sim_test/dataset_sim.csv",
#'   package = "prolfquapp"
#' )
#' config <- make_DEA_config_R6(
#'   PATH = tempdir(),
#'   application = "SIM",
#'   model = "lm"
#' )
#' result <- run_dea(tempdir(), dataset, "prolfquapp.SIM", config)
#' names(result)
run_dea <- function(indir, dataset, software, config) {
  if (!file.exists(dataset)) {
    stop("Annotation file not found: ", dataset, call. = FALSE)
  }
  po <- config$processing_options
  default_model <- .resolve_facade_model(po$model, po$model_missing)
  model_entry <- prolfqua::lookup_facade(default_model)
  if (is.null(model_entry)) {
    stop("Unknown facade: ", default_model, call. = FALSE)
  }
  preprocess_functions <- prolfquapp::get_procfuncs()
  is_nested <- identical(model_entry$needs, "nested")
  config$software <- .resolve_nested_reader(software, is_nested, names(preprocess_functions), default_model)

  annotation <- prolfquapp::read_table_data(dataset) |>
    prolfquapp::read_annotation(prefix = config$group, SAINT = isTRUE(model_entry$needs_saint_annotation))
  if (!config$software %in% names(preprocess_functions)) {
    stop(
      "Software '",
      config$software,
      "' not found. Available: ",
      paste(names(preprocess_functions), collapse = ", "),
      call. = FALSE
    )
  }
  procsoft <- prolfquapp::preprocess_software(
    indir,
    annotation,
    preprocess_functions = preprocess_functions[[config$software]],
    pattern_contaminants = po$pattern_contaminants,
    pattern_decoys = po$pattern_decoys,
    nr_peptides = po$nr_peptides
  )
  xd <- procsoft$xd

  data_prep <- prolfquapp::ProteinDataPrep$new(xd$lfqdata, xd$protein_annotation, config)
  # Single filtering path: no quant filtering here. Contaminants are kept and
  # labelled; decoys are dropped only at the model fit (DEAnalyse$build_facade).
  # This just records the contaminant / decoy QC proportions.
  data_prep$cont_decoy_summary()

  # Route prolfqua's progress into the log: its progress bar is disabled on the
  # non-tty stderr of docker/slurm runs, freezing the log during long fits.
  # Don't clobber a reporter the caller already set.
  if (is.null(getOption("prolfqua.progress"))) {
    options(prolfqua.progress = function(i, total, label) {
      if (is.null(label) || identical(label, "")) {
        label <- "fit"
      }
      logger::log_info("{label}: {i}/{total} ({round(100 * i / total)}%)")
    })
  }
  fit_start <- Sys.time()
  logger::log_info("start fitting / contrasts for model: {default_model}")
  deanalyse <- if (is_nested) {
    .build_nested_deanalyse(data_prep, config, annotation$contrasts, default_model)
  } else {
    data_prep$aggregate()
    .build_protein_deanalyse(data_prep, annotation$contrasts)
  }
  fit_min <- round(as.numeric(difftime(Sys.time(), fit_start, units = "mins")), 2)
  logger::log_info(paste0("done fitting / contrasts for ", default_model, " in ", fit_min, " min"))

  list(
    deanalyse = deanalyse,
    xd = xd,
    annotation = annotation,
    files = procsoft$files,
    software = config$software,
    requested_software = software
  )
}

write_dea_run_outputs <- function(result, config, opt, ymlfile) {
  deanalyse <- result$deanalyse
  xd <- result$xd
  config$software <- result$software %||% opt$software

  logger::log_info("Processing done: ", config$software)
  logger::log_info(paste(
    c("Protein Annotation :\n", capture.output(print(xd$protein_annotation$get_summary()))),
    collapse = "\n"
  ))
  logger::log_info("ContrastNames: \n", paste(names(result$annotation$contrasts), collapse = "\n"))
  logger::log_info("END OF ANALYSIS")

  logger::log_info("CREATING DEAReportGenerator")
  reporter <- DEAReportGenerator$new(deanalyse, config, name = "")
  logger::log_info("Writing results to: ", config$get_zipdir())
  outdir <- reporter$write_DEA_all(boxplot = FALSE)

  arrow::write_parquet(
    deanalyse$lfq_data$data_long(),
    sink = file.path(config$get_result_dir(), "lfqdata_normalized.parquet")
  )
  cfg <- prolfqua::R6_extract_values(deanalyse$lfq_data$get_config())
  yaml::write_yaml(cfg, file.path(config$get_result_dir(), "lfqdata.yaml"))

  # compute_IBAQ_values joins the protein annotation itself, so contaminants
  # stay labelled and decoys need no explicit drop. It mutates its LFQData,
  # hence the copy.
  ibaq_file <- file.path(reporter$resultdir, paste0("IBAQ_", opt$workunit, ".xlsx"))
  ibaq <- NULL
  if (length(xd$lfqdata$relevant_hierarchy_keys()) == 1) {
    ibaq <- compute_IBAQ_values(xd$lfqdata$get_copy(), xd$protein_annotation)
    writexl::write_xlsx(ibaq$data_wide()$data, path = ibaq_file)
  }
  outdir$data_files$ibaq_file <- ibaq_file

  summarized_experiment <- reporter$make_SummarizedExperiment(ibaq = ibaq)
  # Writes AnnData.h5ad + SummarizedExperiment.rds + DEAnalyse.rds and renders
  # the Quarto reports; a failing report warns without aborting the run.
  logger::log_info("Writing AnnData and summarized experiment, then rendering Quarto reports.")
  reports <- render_dea_reports(reporter, summarized_experiment)
  outdir$data_files$anndata_file <- reports$anndata_file
  outdir$dea_file <- reports$dea_file
  outdir$qc_file <- reports$qc_file
  outdir$quarto_file <- reports$tabset_file
  outdir$sse_file <- reports$sse_file
  prolfquapp::write_index_html(outdir, result_dir = reporter$ZIPDIR)

  input_dir <- config$get_input_dir()
  logger::log_info("Creating directory with input files :", input_dir)
  dir.create(input_dir)
  prolfquapp::copy_shell_script(workdir = input_dir)
  file.copy(c(result$files$data, result$files$fasta, ymlfile, opt$dataset), input_dir)
  logger::log_info("Write yaml with parameters: ", file.path(input_dir, "minimal.yaml"))
  yaml::write_yaml(prolfqua::R6_extract_values(config), file = file.path(input_dir, "minimal.yaml"))
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
run_dea_cd <- function(input = NULL, config, files = NULL, subset_column = NULL) {
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
  xd <- list(lfqdata = cd$lfqdata, protein_annotation = cd$protein_annotation)
  data_prep <- prolfquapp::ProteinDataPrep$new(xd$lfqdata, xd$protein_annotation, config)
  # Single filtering path (see run_dea): QC only, no quant filtering here.
  data_prep$cont_decoy_summary()
  if (length(xd$lfqdata$hierarchy_keys()) == xd$lfqdata$get_config()$hierarchy_depth) {
    data_prep$lfq_data <- data_prep$lfq_data_peptide
  } else {
    data_prep$aggregate()
  }
  deanalyse <- .build_protein_deanalyse(data_prep, cd$annotation$contrasts)
  list(deanalyse = deanalyse, xd = xd, annotation = cd$annotation, files = files)
}
