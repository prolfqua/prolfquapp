# generates template yaml file
if (!require("optparse", quietly = TRUE)) {
  install.packages("optparse", dependencies = TRUE)
}


option_list <- list(
  make_option(c("-n","--norm"), type = "character", default = "vsn",
                help = "normalization method to use, vsn, none, or robscale",
                metavar = "character"),
  make_option(c("-o", "--outdir"), type = "character", default = ".",
              help = "folder write yaml file to",
              metavar = "string"),
  make_option(c("-y", "--yaml"), type = "character", default = "config.yaml",
              help = "yaml configuration file name",
              metavar = "character"),
  make_option(c("-w", "--workunit"), type = "character", default = "",
              help = "workunit ID",
              metavar = "character"),
  make_option(c("-O", "--order"), type = "character", default = "",
              help = "order ID",
              metavar = "character"),
  make_option(c("-p", "--project"), type = "character", default = "",
              help = "project ID",
              metavar = "character"),
  make_option(c("-s", "--software"), type = "character", default = "DIANN",
              help = "either DIANN, FP_TMT, MAXQUANT",
              metavar = "character")

)

parser <- optparse::OptionParser(usage = "%prog --norm vsn --workunit WUID332211", option_list = option_list)
arguments <- optparse::parse_args(parser, positional_arguments = TRUE)


#parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
#arguments <- parse_args(parser, args = "test.yml", positional_arguments = TRUE)
lobstr::tree(arguments)


opt <- arguments$options
ymlfile <- arguments$args
ymlfile <- if ( length(ymlfile) == 0 ) { opt$yaml } else {ymlfile}

logger::log_info("writing yaml file : ", ymlfile)
GRP2 <- prolfquapp::make_DEA_config_R6(
  PROJECTID = opt$project,
  ORDERID = opt$order,
  WORKUNITID = opt$workunit,
  Normalization = opt$norm
  )
GRP2$set_zipdir_name()
if (!is.null(opt$outdir) && dir.exists(opt$outdir)) {
  GRP2$path <- opt$outdir
}
cfg <- prolfquapp::R6_extract_values(GRP2)
cfg <- GRP2$as_list()

# Define the fields that should be moved to the bottom
fields_to_move <- c("ext_reader", "group", "RES", "pop")
# Separate the fields into 'main' and 'bottom'
main_fields <- cfg[!names(cfg) %in% fields_to_move]
bottom_fields <- cfg[names(cfg) %in% fields_to_move]
# Combine the main fields and bottom fields in the desired order
cfg <- c(main_fields, bottom_fields)

if (!dir.exists(opt$outdir)) {dir.create(opt$outdir)}
yaml::write_yaml(cfg, file = file.path(opt$outdir , ymlfile))

