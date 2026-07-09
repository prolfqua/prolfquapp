# generates template yaml file
if (!require("optparse", quietly = TRUE)) {
  install.packages("optparse", dependencies = TRUE)
}

prolfquapp::route_messages_to_logger()

option_list <- list(
  make_option(
    c("-n", "--norm"),
    type = "character",
    default = "vsn",
    help = "normalization method to use, vsn, none, or robscale",
    metavar = "character"
  ),
  make_option(
    c("-o", "--outdir"),
    type = "character",
    default = ".",
    help = "folder write yaml file to",
    metavar = "string"
  ),
  make_option(
    c("-y", "--yaml"),
    type = "character",
    default = "config.yaml",
    help = "yaml configuration file name",
    metavar = "character"
  ),
  make_option(
    c("-w", "--workunit"),
    type = "character",
    default = "",
    help = "workunit ID",
    metavar = "character"
  ),
  make_option(
    c("-O", "--order"),
    type = "character",
    default = "",
    help = "order ID",
    metavar = "character"
  ),
  make_option(
    c("-p", "--project"),
    type = "character",
    default = "",
    help = "project ID",
    metavar = "character"
  ),
  make_option(
    c("-s", "--software"),
    type = "character",
    default = "DIANN",
    help = "either DIANN, FP_TMT, MAXQUANT",
    metavar = "character"
  ),
  make_option(
    c("-m", "--model"),
    type = "character",
    default = "lm_missing",
    help = paste0(
      "contrast facade method. Options: ",
      paste(c(names(prolfqua::FACADE_REGISTRY), "saint"), collapse = ", ")
    ),
    metavar = "character"
  ),
  make_option(
    c("--nr_peptides"),
    type = "integer",
    default = 1,
    help = "minimum distinct peptides per protein (>= 1)",
    metavar = "N"
  )
)

parser <- optparse::OptionParser(
  usage = "%prog --norm vsn --workunit WUID332211",
  option_list = option_list
)
arguments <- optparse::parse_args(parser, positional_arguments = TRUE)


lobstr::tree(arguments)


opt <- arguments$options
ymlfile <- arguments$args
ymlfile <- if (length(ymlfile) == 0) {
  opt$yaml
} else {
  ymlfile
}

logger::log_info("writing yaml file : ", ymlfile)
cfg <- prolfquapp::run_make_yaml(
  project = opt$project,
  order = opt$order,
  workunit = opt$workunit,
  norm = opt$norm,
  model = opt$model,
  nr_peptides = opt$nr_peptides,
  outdir = opt$outdir
)

if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir)
}
yaml::write_yaml(cfg, file = file.path(opt$outdir, ymlfile))
