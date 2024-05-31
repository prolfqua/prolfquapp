# generates template yaml file
if (!require("optparse", quietly = TRUE)) {
  install.packages("optparse", dependencies = TRUE)
}


option_list <- list(
  make_option(c("-t","--trans"), type = "character", default = "vsn",
                help = "normalization method to use, vsn, none, or robscale",
                metavar = "character"),
  make_option(c("-d", "--outdir"), type = "character", default = ".",
              help = "folder write yaml file to",
              metavar = "string"),
  make_option(c("-y", "--yaml"), type = "character", default = "config.yaml",
              help = "yaml configuration file name",
              metavar = "character"),
  make_option(c("-w", "--workunit"), type = "character", default = "1234",
              help = "workunit ID",
              metavar = "character"),
  make_option(c("-o", "--order"), type = "character", default = "1234",
              help = "order ID",
              metavar = "character"),
  make_option(c("-p", "--project"), type = "character", default = "1234",
              help = "project ID",
              metavar = "character")

)

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list, add_help_option = TRUE)
arguments <- parse_args(parser, args = "test.yml", positional_arguments = TRUE)

opt <- arguments$options

ymlfile <- arguments$args
ymlfile


ymlfile <- if ( length(ymlfile) == 0 ) { opt$yaml } else {ymlfile}
logger::log_info("writing yaml file : ", ymlfile)
GRP2 <- prolfquapp::make_DEA_config_R6(
  ZIPDIR = "DEA",
  PROJECTID = opt$project,
  ORDERID = opt$order,
  WORKUNITID = opt$workunit,
  Normalization = opt$trans
  )
yaml::write_yaml(prolfquapp::R6_extract_values(GRP2), file = file.path(opt$outdir , ymlfile))
