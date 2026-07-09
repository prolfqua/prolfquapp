if (!require("optparse", quietly = TRUE)) {
  install.packages("optparse", dependencies = TRUE)
}

option_list <- list(
  optparse::make_option(
    c("-i", "--input"),
    type = "character",
    default = ".",
    help = "CompoundDiscoverer ZIP export or folder containing one export",
    metavar = "path"
  ),
  optparse::make_option(
    c("-y", "--yaml"),
    type = "character",
    default = NULL,
    help = "optional YAML configuration file",
    metavar = "character"
  ),
  optparse::make_option(
    c("-w", "--workunit"),
    type = "character",
    default = NULL,
    help = "workunit ID override",
    metavar = "character"
  ),
  optparse::make_option(
    c("-o", "--outdir"),
    type = "character",
    default = ".",
    help = "output directory",
    metavar = "character"
  ),
  optparse::make_option(
    c("-m", "--model"),
    type = "character",
    default = NULL,
    help = paste0(
      "contrast facade method (overrides config). Options: ",
      paste(names(prolfqua::FACADE_REGISTRY), collapse = ", ")
    ),
    metavar = "character"
  ),
  optparse::make_option(
    c("-n", "--normalization"),
    type = "character",
    default = NULL,
    help = "normalization method override: vsn, robscale, or none. Default without YAML: vsn",
    metavar = "character"
  ),
  optparse::make_option(
    c("--subset-columns"),
    type = "character",
    default = "auto",
    help = paste(
      "comma-separated long-table subset columns to run separately;",
      "'auto' uses all columns after Group, 'none' disables subsets"
    ),
    metavar = "character"
  ),
  optparse::make_option(
    c("--include-full"),
    type = "logical",
    default = TRUE,
    help = paste(
      "run the full unfiltered CD table in addition to subset runs;",
      "set to FALSE to generate subset-only outputs"
    ),
    metavar = "logical"
  ),
  optparse::make_option(
    c("--libPath"),
    type = "character",
    default = NULL,
    help = "optional R library path",
    metavar = "string"
  )
)

parser <- optparse::OptionParser(
  usage = "%prog --input compound_discoverer_export.zip --outdir results",
  option_list = option_list
)

if (length(commandArgs(TRUE)) == 0) {
  optparse::print_help(parser)
}

arguments <- optparse::parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

logger::log_appender(logger::appender_console)
prolfquapp::route_messages_to_logger()
logger::log_info(
  "LIBRARY PATHS (.libPaths()):",
  paste(.libPaths(), collapse = "\n")
)

if (!is.null(opt$libPath) && nchar(opt$libPath) > 0) {
  logger::log_info("Setting libPath: ", opt$libPath)
  .libPaths(opt$libPath)
  logger::log_info(
    "LIBRARY PATHS (.libPaths()):",
    paste(.libPaths(), collapse = "\n")
  )
}

library(prolfquapp)
logger::log_info("using : ", system.file(package = "prolfqua"))
logger::log_info("using : ", system.file(package = "prolfquapp"))

infer_ids <- function(input) {
  stem <- tools::file_path_sans_ext(basename(input))
  project <- sub("^.*p([0-9]+).*$", "\\1", stem)
  order <- sub("^.*o([0-9]+).*$", "\\1", stem)
  list(
    project = if (identical(project, stem)) "" else project,
    order = if (identical(order, stem)) "" else order,
    workunit = stem
  )
}

ids <- infer_ids(opt$input)
normalization <- if (is.null(opt$normalization)) {
  NULL
} else {
  match.arg(opt$normalization, c("vsn", "robscale", "none"))
}
if (!is.null(opt$yaml)) {
  if (!file.exists(opt$yaml)) {
    stop("YAML file not found: ", opt$yaml, call. = FALSE)
  }
  logger::log_info("YAML file read: ", opt$yaml)
  GRP2 <- prolfquapp::get_config(opt$yaml)
} else {
  logger::log_info("No YAML supplied; using default CD DEA configuration.")
  GRP2 <- prolfquapp::make_DEA_config_R6(
    PATH = opt$outdir,
    PROJECTID = ids$project,
    ORDERID = ids$order,
    WORKUNITID = ids$workunit,
    Normalization = if (is.null(normalization)) "vsn" else normalization,
    application = "CompoundDiscoverer",
    model = "lm_missing"
  )
}

GRP2$software <- "CompoundDiscoverer"
GRP2$path <- opt$outdir
if (!is.null(opt$workunit)) {
  logger::log_info("Setting workunit to: ", opt$workunit)
  GRP2$project_spec$workunit_Id <- opt$workunit
}
if (!is.null(opt$model)) {
  GRP2$processing_options$model <- opt$model
}
if (
  !is.null(normalization) &&
    !identical(GRP2$processing_options$transform, normalization)
) {
  logger::log_info("Setting normalization to: ", normalization)
  GRP2$processing_options$transform <- normalization
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
files <- prolfquapp::get_CD_export_files(opt$input)
on.exit(unlink(files$tempdir, recursive = TRUE), add = TRUE)

subset_columns_opt <- opt$subset_columns
if (is.null(subset_columns_opt)) {
  subset_columns_opt <- opt[["subset-columns"]]
}
if (!is.character(subset_columns_opt) || length(subset_columns_opt) == 0) {
  subset_columns_opt <- "auto"
}
include_full <- opt$include_full
if (is.null(include_full)) {
  include_full <- opt[["include-full"]]
}
include_full <- isTRUE(include_full)

subset_columns <- prolfquapp::get_CD_subset_columns(files$data)
if (!identical(subset_columns_opt, "auto")) {
  subset_columns <- if (identical(tolower(subset_columns_opt), "none")) {
    character()
  } else {
    trimws(strsplit(subset_columns_opt, ",", fixed = TRUE)[[1]])
  }
}
subset_columns <- subset_columns[nchar(subset_columns) > 0]

runs <- list()
if (length(subset_columns) == 0 || include_full) {
  runs[[length(runs) + 1]] <- list(label = NULL, subset_column = NULL)
}
for (subset_column in subset_columns) {
  runs[[length(runs) + 1]] <- list(
    label = prolfquapp::sanitize_CD_subset_name(subset_column),
    subset_column = subset_column
  )
}
if (length(runs) == 0) {
  runs[[1]] <- list(label = NULL, subset_column = NULL)
}

logger::log_info("Software: CompoundDiscoverer")
logger::log_info(
  "CD runs: ",
  paste(
    vapply(
      runs,
      function(run) if (is.null(run$label)) "full" else run$label,
      character(1)
    ),
    collapse = ", "
  )
)

run_one <- function(run, base_config) {
  run_config <- base_config$clone(deep = TRUE)
  run_config$set_zipdir_name()
  if (!is.null(run$label)) {
    run_config$zipdir_name <- paste0(run_config$zipdir_name, "_", run$label)
  }

  dir.create(run_config$get_zipdir(), showWarnings = FALSE, recursive = TRUE)

  current_time <- Sys.time()
  formatted_time <- format(current_time, "%Y%m%d%H%M")
  logfile <- paste0("prolfqua_", formatted_time, ".log")
  appender_combined <- logger::appender_tee(file.path(
    run_config$get_zipdir(),
    logfile
  ))
  logger::log_appender(appender_combined)
  logger::log_info(prolfquapp::capture_output(quote(lobstr::tree(opt))))
  logger::log_info(
    "Writing to output directory : ",
    run_config$get_zipdir(),
    " and file :",
    logfile
  )
  if (!is.null(run$subset_column)) {
    logger::log_info(
      "Running CD subset column: ",
      run$subset_column,
      " -> folder suffix: ",
      run$label
    )
  }

  logger::log_info("prolfquapp parameters : ")
  logger::log_info(paste(
    capture.output(lobstr::tree(prolfqua::R6_extract_values(run_config))),
    collapse = "\n"
  ))

  result <- tryCatch(
    prolfquapp::run_dea_cd(
      config = run_config,
      files = files,
      subset_column = run$subset_column
    ),
    error = function(e) {
      stack_trace <- capture.output(traceback())
      logger::log_error(conditionMessage(e), "\n")
      logger::log_error("Stack trace:\n")
      logger::log_error(
        paste(stack_trace, collapse = "\n"),
        "\n"
      )
      if (interactive()) stop(e) else quit(save = "no", status = 1)
    }
  )

  deanalyse <- result$deanalyse
  xd <- result$xd
  annotation <- result$annotation

  logger::log_info("Processing done: CompoundDiscoverer")
  logger::log_info(paste(
    c(
      "Metabolite Annotation :\n",
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
  reporter <- prolfquapp::DEAReportGenerator$new(
    deanalyse,
    run_config,
    name = ""
  )

  logger::log_info("Writing results to: ", run_config$get_zipdir())
  outdir <- reporter$write_DEA_all(
    boxplot = FALSE
  )

  cfg <- prolfqua::R6_extract_values(deanalyse$lfq_data$get_config())
  yaml::write_yaml(cfg, file.path(run_config$get_result_dir(), "lfqdata.yaml"))

  # Writes SummarizedExperiment.rds + DEAnalyse.rds and renders the Quarto
  # reports (primary R6 DEA report, SE-tabset overview, differential-expression
  # QC, and sample-size estimation). Each renders independently; a failure warns
  # without aborting the run.
  logger::log_info("Writing summarized experiment and rendering Quarto reports.")
  reports <- prolfquapp:::render_dea_reports(reporter)
  outdir$dea_file <- reports$dea_file
  outdir$qc_file <- reports$qc_file
  outdir$quarto_file <- reports$tabset_file
  outdir$sse_file <- reports$sse_file

  prolfquapp::write_index_html(outdir, result_dir = reporter$ZIPDIR)

  logger::log_info("Writing normalized LFQData Parquet export.")
  parquet <- prolfquapp:::write_parquet_isolated(
    deanalyse$lfq_data$data_long(),
    sink = file.path(run_config$get_result_dir(), "lfqdata_normalized.parquet")
  )
  if (isTRUE(parquet$success)) {
    logger::log_info("Parquet export written: ", parquet$sink)
  } else {
    logger::log_warn(
      "Parquet export failed with exit status ",
      parquet$exit_status,
      ". Continuing because the report outputs were already written."
    )
    if (length(parquet$output) > 0) {
      logger::log_warn(
        "Parquet export output:\n",
        paste(parquet$output, collapse = "\n")
      )
    }
  }

  logger::log_info(
    "Creating directory with input files :",
    run_config$get_input_dir()
  )
  dir.create(run_config$get_input_dir(), showWarnings = FALSE, recursive = TRUE)

  prolfquapp::copy_shell_script(workdir = run_config$get_input_dir())

  file.copy(
    c(files$zip, files$data, files$samples, opt$yaml),
    run_config$get_input_dir()
  )

  logger::log_info(
    "Write yaml with parameters: ",
    file.path(run_config$get_input_dir(), "minimal.yaml")
  )

  yaml::write_yaml(
    prolfqua::R6_extract_values(run_config),
    file = file.path(run_config$get_input_dir(), "minimal.yaml")
  )
}

for (run in runs) {
  run_one(run, GRP2)
}
