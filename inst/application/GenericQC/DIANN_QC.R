#R
# 20211210 WEW/CP make it work for WU272669
# 20230131 make it work for Application 312
# R_LIBS_SITE="/scratch/FRAGPIPEDIA_312/R_LIBS_V1/"; R --vanilla  < ~/slurmworker/R/fgcz_fragpipeDIA_DIANN_prolfqua_qc.R
##### QCs

library(dplyr)
library(tidyverse)
library(prolfquapp)
library(logger)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  path <- args[1]
  project_Id <- args[2]
  output_dir <- args[3]
  libPath <- args[3]

  logger::log_info(paste("libPath :" , .libPaths(), collapse = " ;"))

  if (is.na(libPath) & dir.exists(libPath) ) {
    logger::log_info(paste("Setting libPath:", libPath, collapse = " ;"))
    .libPaths(libPath)
  }
} else {
  message("please provide :\n",
          "path       : folder containing fasta, diann-output.tsv and dataset.tsv file \n",
          "project_id : b-fabric project id\n",
          "output_dir : folder to write the results to \n",
          "libPath.   : (optional) R library path\n"
          )
  path = "WU287241"
  project_Id = "o32211"
  output_dir = "drummmmm"
}

if (!dir.exists(output_dir)) { dir.create(output_dir) }

mdir <- function(path, pattern){
  dir(path, pattern, full.names = TRUE, recursive = TRUE)
}

fasta.file <- mdir(path,
                   pattern = "*.fasta$")
diann.output <- mdir(path,
                     pattern = "diann-output.tsv")
dataset.csv <- mdir(path,
                    pattern = "dataset.csv")


xx <- readr::read_tsv(diann.output)
peptide <- read_DIANN_output(
  diann.path = diann.output,
  fasta.file = fasta.file,
  nrPeptides = 1,
  Q.Value = 0.1)

prot_annot <- prolfquapp::dataset_protein_annot(
  peptide,
  "Protein.Group.2",
  protein_annot = "fasta.header")

# dataset.csv must either contain columns:
# `Relative Path` Name `Grouping Var` FileType
# else the report will be generated but without this information.
if (length(dataset.csv) > 0) {
  annotation <- read.csv(dataset.csv)
  annotation$inputresource.name <- tools::file_path_sans_ext(basename(annotation$Relative.Path))
  annotation$sample.name <- make.unique(annotation$Name)
} else{
  annotation <- data.frame(Relative.Path = unique(peptide$raw.file))
  annotation$inputresource.name <- tools::file_path_sans_ext(basename(annotation$Relative.Path))
  annotation$sample.name <- gsub("^[0-9]{8,8}_","", annotation$Relative.Path)
  annotation$Grouping.Var <- "None_Specified"
}

peptide <- dplyr::inner_join(
  annotation,
  peptide,
  by = c("inputresource.name" = "raw.file"),
  multiple = "all")

# MSFragger specific (moving target)
atable <- prolfqua::AnalysisTableAnnotation$new()
atable$fileName = "inputresource.name"
atable$sampleName = "sample.name"
atable$ident_qValue = "PEP"
atable$hierarchy[["protein_Id"]] <- c("Protein.Group.2")
atable$hierarchy[["peptide_Id"]] <- c("Stripped.Sequence")
atable$hierarchyDepth <- 1
atable$set_response("Peptide.Quantity")
atable$factors[["group"]] = "Grouping.Var"
atable$factorDepth <- 1


config <- prolfqua::AnalysisConfiguration$new(atable)
adata <- prolfqua::setup_analysis(peptide, config)
lfqdata <- prolfqua::LFQData$new(adata, config)
lfqdata$remove_small_intensities(threshold = 1)

TABLES2WRITE <- list()
TABLES2WRITE$peptide_wide <- left_join(prot_annot,
                                       lfqdata$to_wide()$data,
                                       multiple = "all")

TABLES2WRITE$annotation <- lfqdata$factors()
lfqdataProt <- prolfquapp::aggregate_data(lfqdata, agg_method = "medpolish")
TABLES2WRITE$proteins_wide <- left_join(prot_annot,
                                        lfqdataProt$to_wide()$data,
                                        multiple = "all")

summarizer <- lfqdata$get_Summariser()
precabund <- summarizer$percentage_abundance()
precabund <- inner_join(
  prot_annot,
  precabund,
  multiple = "all",
  by = lfqdata$config$table$hierarchy_keys_depth())

TABLES2WRITE$proteinAbundances <- precabund
writexl::write_xlsx(TABLES2WRITE, path = file.path(output_dir, "proteinAbundances.xlsx"))

file.copy(system.file("application/genericQC/test_NewQc.Rmd", package = "prolfquapp"),
          to = output_dir, overwrite = TRUE)

rmarkdown::render(file.path(output_dir,"test_NewQc.Rmd"), params = list(config = lfqdataProt$config$table,
                                                  precabund = precabund,
                                                  factors = TRUE))
ps = prolfqua::ProjectStructure$new(outpath = path,
                                    project_Id = "",
                                    workunit_Id = basename(getwd()),
                                    order_Id = "",
                                    inputAnnotation = NULL,
                                    inputData = NULL)
ps$create()
if (nrow(lfqdata$factors()) > 1) {
  file.copy(system.file("application/genericQC/QCandSSE.Rmd", package = "prolfquapp"),
            to = output_dir, overwrite = TRUE)

  rmarkdown::render(file.path(output_dir,"QCandSSE.Rmd"),
                    params = list(data = lfqdata$data, configuration = lfqdata$config, project_conf = ps, pep = FALSE),

                    )
} else{
  message("only a single sample: ", nrow(lfqdata$factors()))

}

