if (!require("prolfqua", quietly = TRUE))
  remotes::install_github("prolfqua/prolfquapp", dependencies = TRUE)
if (!require("prolfqua", quietly = TRUE))
  remotes::install_github("prolfqua/prolfquapp", dependencies = TRUE)
if (!require("optparse", quietly = TRUE))
  install.packages("optparse", dependencies = TRUE)


library("optparse")
parser <- OptionParser()
parser <- add_option(parser, c("-d", "--dir"), type = "character", default = ".",
                     help = "folder containing fasta, diann-output.tsv and dataset.tsv file",
                     metavar = "string")
parser <- add_option(parser, c("-p", "--projectid"), type = "character", default = "1234",
                     help = "your project identifier",
                     metavar = "string")
parser <- add_option(parser, c("-o", "--output"), type = "character", default = "DEA",
                     help = "folder to write the results to.",
                     metavar = "string")
parser <- add_option(parser, c("--libPath"), type = "character", default = "DEA",
                     help = " (optional) R library path",
                     metavar = "string")

opt <- parse_args(parser)

