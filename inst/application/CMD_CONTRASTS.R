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

# --- Dispatch -----------------------------------------------------------------

if (one_factor) {
  annot_out <- prolfquapp::run_contrasts_single(
    annotation_file, opt$control, opt$group
  )
  logger::log_info(
    "CONTROL column added (C = '", opt$control, "'):\n",
    paste(
      capture.output(print(
        table(annot_out[[grep("group", colnames(annot_out), value = TRUE)[1]]],
              annot_out$CONTROL)
      )),
      collapse = "\n"
    )
  )
} else {
  annot_out <- prolfquapp::run_contrasts_twofactor(
    annotation_file, opt$f1, opt$f2, opt$interactions
  )
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
