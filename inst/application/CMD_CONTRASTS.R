# CMD_CONTRASTS.R
#
# Generate contrast definitions from an annotation file.
#
# Scenario 1 - single factor:
#   Rscript CMD_CONTRASTS.R annotation.csv --control WT -o out.csv
#   Rscript CMD_CONTRASTS.R annotation.csv --control WT --group treatment -o out.csv
#
# Scenario 2 - two factor:
#   Rscript CMD_CONTRASTS.R annotation.csv --f1 treatment --f2 time -o out.csv
#   Rscript CMD_CONTRASTS.R annotation.csv --f1 treatment --f2 time --interactions FALSE

if (!require("optparse", quietly = TRUE)) {
  install.packages("optparse", dependencies = TRUE)
}

option_list <- list(
  optparse::make_option(
    c("--control"),
    type = "character",
    default = NULL,
    help = "[Scenario 1] Reference level of the group factor",
    metavar = "character"
  ),
  optparse::make_option(
    c("--group"),
    type = "character",
    default = NULL,
    help = "[Scenario 1] Group column name [auto-detected if omitted]",
    metavar = "character"
  ),
  optparse::make_option(
    c("--f1"),
    type = "character",
    default = NULL,
    help = "[Scenario 2] Primary factor column name",
    metavar = "character"
  ),
  optparse::make_option(
    c("--f2"),
    type = "character",
    default = NULL,
    help = "[Scenario 2] Secondary factor column name",
    metavar = "character"
  ),
  optparse::make_option(
    c("--interactions"),
    type = "logical",
    default = TRUE,
    help = "[Scenario 2] Include interaction contrasts [default: TRUE]",
    metavar = "logical"
  ),
  optparse::make_option(
    c("-o", "--output"),
    type = "character",
    default = NULL,
    help = "Output file (csv/tsv/xlsx); if omitted, only prints contrasts",
    metavar = "character"
  ),
  optparse::make_option(
    c("--libPath"),
    type = "character",
    default = NULL,
    help = "(optional) R library path",
    metavar = "string"
  )
)

parser <- optparse::OptionParser(
  usage = paste(
    "\n  Scenario 1 (single factor):  %prog annotation.csv --control WT [-o out.csv]",
    "\n  Scenario 2 (two factor):     %prog annotation.csv --f1 treatment --f2 time [-o out.csv]"
  ),
  option_list = option_list
)

if (length(commandArgs(TRUE)) == 0) {
  optparse::print_help(parser)
  quit(status = 0)
}

arguments <- optparse::parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
args <- arguments$args

logger::log_appender(logger::appender_console)

if (!is.null(opt$libPath) && nchar(opt$libPath) > 0) {
  .libPaths(c(opt$libPath, .libPaths()))
}

library(prolfquapp)
logger::log_info("prolfquapp : ", system.file(package = "prolfquapp"))
logger::log_info("prolfqua   : ", system.file(package = "prolfqua"))

# --- Validate ----------------------------------------------------------------

if (length(args) == 0) {
  logger::log_error("No annotation file provided.")
  optparse::print_help(parser)
  quit(status = 1)
}

annotation_file <- args[1]
if (!file.exists(annotation_file)) {
  logger::log_error("File not found: ", annotation_file)
  quit(status = 1)
}

two_factor <- !is.null(opt$f1) && !is.null(opt$f2)
one_factor <- !is.null(opt$control)

if (!one_factor && !two_factor) {
  logger::log_error(
    "Specify --control (scenario 1) or --f1 + --f2 (scenario 2)."
  )
  optparse::print_help(parser)
  quit(status = 1)
}

# =============================================================================
# Scenario 1: single factor — add CONTROL column (C/T)
# =============================================================================
if (one_factor) {
  res <- prolfquapp::read_annotation(annotation_file, QC = TRUE)
  annot <- res$annot
  group_col <- if (!is.null(opt$group)) {
    opt$group
  } else {
    res$atable$factors[["G_"]]
  }

  logger::log_info("Group column : ", group_col)
  logger::log_info(
    "Levels found : ",
    paste(unique(annot[[group_col]]), collapse = ", ")
  )

  if (!opt$control %in% annot[[group_col]]) {
    logger::log_error(
      "--control '",
      opt$control,
      "' not found in '",
      group_col,
      "'."
    )
    quit(status = 1)
  }

  annot$CONTROL <- ifelse(annot[[group_col]] == opt$control, "C", "T")
  annot_out <- annot

  logger::log_info(
    "CONTROL column added (C = '",
    opt$control,
    "', T = everything else):\n",
    paste(
      capture.output(print(table(annot_out[[group_col]], annot_out$CONTROL))),
      collapse = "\n"
    )
  )
}

# =============================================================================
# Scenario 2: two factor
# =============================================================================
if (two_factor) {
  df <- prolfquapp::read_table_data(annotation_file)

  missing_cols <- setdiff(c(opt$f1, opt$f2), colnames(df))
  if (length(missing_cols) > 0) {
    logger::log_error(
      "Column(s) not found: ",
      paste(missing_cols, collapse = ", ")
    )
    quit(status = 1)
  }

  logger::log_info(
    "f1 = '",
    opt$f1,
    "' levels: ",
    paste(unique(df[[opt$f1]]), collapse = ", ")
  )
  logger::log_info(
    "f2 = '",
    opt$f2,
    "' levels: ",
    paste(unique(df[[opt$f2]]), collapse = ", ")
  )

  res <- prolfqua::annotation_add_contrasts(
    df,
    primary_col = opt$f1,
    secondary_col = opt$f2,
    interactions = opt$interactions
  )
  annot_out <- res$annot
}

# --- Print (scenario 2 only) -------------------------------------------------

if (two_factor) {
  ct <- dplyr::distinct(annot_out[
    !is.na(annot_out$ContrastName),
    c("ContrastName", "Contrast")
  ])
  logger::log_info(
    "Suggested contrasts:\n",
    paste(ct$ContrastName, ct$Contrast, sep = " = ", collapse = "\n")
  )
}

# --- Write (optional) --------------------------------------------------------

if (!is.null(opt$output)) {
  prolfquapp::write_annotation_file(annot_out, opt$output)
  logger::log_info("Written to: ", opt$output)
}
