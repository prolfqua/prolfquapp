#R
# 20211210 WEW/CP make it work for WU272669
# 20230131 make it work for Application 312
# R_LIBS_SITE="/scratch/FRAGPIPEDIA_312/R_LIBS_V1/"; R --no-save  < ~/slurmworker/R/fgcz_fragpipeDIA_DIANN_prolfqua_qc.R
##### QCs

library(dplyr)
library(tidyverse)
library(prolfquapp)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  path <- args[1]
} else {
  path = "WU287241"
  path <- "."
}

mdir <- function(path, pattern){
  dir(path, pattern, full.names = TRUE, recursive = TRUE)
}

fasta.file <- mdir(path,
  pattern = "*.fasta$")
diann.output <- mdir(path,
  pattern = "diann-output.tsv")
dataset.csv <- mdir(path,
  pattern = "dataset1.csv")

xx <- readr::read_tsv(diann.output)
peptide <- read_DIANN_output(
  diann.path = diann.output,
  fasta.file = fasta.file,
  ,nrPeptides = 1,
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
  annotation$Grouping.Var <- "None"
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

tables2write <- list()
tables2write$peptide_wide <- left_join(prot_annot, lfqdata$to_wide()$data, multiple = "all")
tables2write$annotation <- lfqdata$factors()

lfqdataProt <- prolfquapp::aggregate_data(lfqdata, agg_method = "medpolish")
tables2write$proteins_wide <- left_join(prot_annot, lfqdataProt$to_wide()$data, multiple = "all")

logger::log_info("data aggregated: medpolish.")
ps <- prolfqua::ProjectStructure$new(outpath = path,
                                     project_Id = "",
                                     workunit_Id = basename(getwd()),
                                     order_Id = "",
                                     inputAnnotation = NULL,
                                     inputData = NULL)

ps$create()
prolfqua::LFQDataSummariser$undebug("percentage_abundance")
summarizer <- lfqdata$get_Summariser()
summarizer$plot_hierarchy_counts_sample()

precabund <- summarizer$percentage_abundance()


precabund <- inner_join(
  prot_annot,
  precabund,
  multiple = "all",
  by = lfqdata$config$table$hierarchy_keys_depth())

readr::write_tsv(precabund,file = "testData.tsv")
saveRDS(lfqdata, "lfqdata.Rds")
colnames(precabund)
precabund$IDcolumn <- NULL
precabund$isotopeLabel <- NULL
precabund$id <- NULL
precabund$medianArea <- NULL
precabund$nrNAs <- NULL

# pp <- prolfquapp::plot_abundance_vs_percent(
#   precabund,
#   cfg_table = lfqdata$config$table,
#   top_N = NULL,factors = FALSE,
#   cumulative = TRUE)

precabund_table <- precabund |> mutate(
  abundance_percent = signif(abundance_percent, n ),
  abundance_percent_cumulative = signif(abundance_percent_cumulative, n),
  percent_prot = signif(percent_prot, 3))

protID = lfqdata$config$table$hierarchy_keys_depth()
datax <- crosstalk::SharedData$new(precabund_table , key = as.formula(paste(" ~ ", protID )), group = "Choose Protein")
DT::datatable(datax)


params <- list()
params$plot <- pp
params$data <- datax

rmarkdown::render("test_NewQc.Rmd")

if (FALSE) {
  if (nrow(lfqdata$factors()) > 1) {
    rmarkdown::render("QCandSSE.Rmd", params = list(data = lfqdata$data, configuration = lfqdata$config, project_conf = ps, pep = FALSE))
  } else{
    message("only a single sample: ", nrow(lfqdata$factors()))

  }

}
