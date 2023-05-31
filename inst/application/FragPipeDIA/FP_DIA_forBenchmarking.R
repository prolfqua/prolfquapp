library(tidyverse)
library(prolfquapp)

#### inputdir <- "/Users/witoldwolski/__checkout/protriple/inst/FP_diann_DataBench/Triple_Proteome_Exp2_Grad30_staggered/"
inputdir <- "/Users/witoldwolski/__checkout/protriple/inst/FP_diann_DataBench/Triple_Proteome_Exp2_Grad90/"
#inputdir <- "/Users/witoldwolski/__checkout/protriple/inst/FP_diann_DataBench/Triple_Proteome_Exp2_Grad90_staggered/"


#inputdir <- "/Users/witoldwolski/__checkout/protriple/inst/diann_DataBench/Triple_Proteome_Exp2_Grad90_V2/"
inputdir <- "/Users/witoldwolski/__checkout/protriple/inst/diann_DataBench/Triple_Proteome_Exp2_Grad90_staggered_nondecoy_db_V2/"


zipdir <- basename(inputdir)
if (grepl("FP_diann", inputdir)) {
  print("FP")
  GRP2 <- prolfquapp::make_DEA_config(ZIPDIR = paste0("FP_diann_", zipdir))
} else {
  print("DIANN")
  GRP2 <- prolfquapp::make_DEA_config(ZIPDIR = paste0("diann_",zipdir))
}

GRP2$pop$transform = "none"

###
dir.create(GRP2$zipdir)
###
# reading foreign data
REPEATED <- TRUE

diann.path <- grep("diann-output.tsv", dir(inputdir,recursive = TRUE,full.names = TRUE), value = TRUE)
fasta.file <- grep("*.fasta", dir(inputdir,recursive = TRUE, full.names = TRUE), value = TRUE)
ds_file <-  grep("dataset.csv", dir(inputdir,recursive = TRUE, full.names = TRUE), value = TRUE)

peptide <- read_DIANN_output(
  diann.path = diann.path[1],
  fasta.file = fasta.file[1],
  ,nrPeptides = 1,
  Q.Value = 0.1)

prot_annot <- prolfquapp::dataset_protein_annot(
  peptide,
  "Protein.Group.2",
  protein_annot = "fasta.header")


annot <- read.csv(ds_file)
annot <- data.frame(lapply(annot, as.character))
annot <- annot |> dplyr::mutate(
  raw.file = gsub("^x|.d.zip$|.raw$","",
                  (basename(annot$Relative.Path))
  )
)
annot$CONTROL <- ifelse(annot$GroupingVar == "A", "C", "T")
GRP2 <- prolfquapp::dataset_extract_contrasts(annot, GRP2)
annot$raw.file[ !annot$raw.file %in% sort(unique(peptide$raw.file)) ]
nr <- sum(annot$raw.file %in% sort(unique(peptide$raw.file)))
logger::log_info("nr : ", nr, " files annotated")
annot$Relative.Path <- NULL

annot$Name <- paste(annot$Name,"_" , 1:nrow(annot), sep = "")
peptide <- dplyr::inner_join(annot, peptide, multiple = "all")
annotPep <- peptide

atable <- prolfqua::AnalysisTableAnnotation$new()
atable$fileName = "raw.file"
atable$hierarchy[["protein_Id"]] <- c("Protein.Group.2")
atable$hierarchy[["peptide_Id"]] <- c("Stripped.Sequence")
atable$set_response("Peptide.Quantity")
atable$hierarchyDepth <- 1

res <- prolfquapp::dataset_set_factors(atable, peptide)
atable <- res$atable
peptide <- res$msdata

# Preprocess data - aggregate proteins.
config <- prolfqua::AnalysisConfiguration$new(atable)

adata <- prolfqua::setup_analysis(peptide, config)


lfqdata <- prolfqua::LFQData$new(adata, config)
lfqdata$remove_small_intensities()
GRP2$pop$nrPeptides <- 1
logger::log_info("AGGREGATING PEPTIDE DATA!")

lfqdata <- prolfquapp::aggregate_data(lfqdata, agg_method = GRP2$pop$aggregate)
logger::log_info("data aggregated: {GRP2$pop$aggregate}.")
logger::log_info("END OF DATA TRANSFORMATION.")


prolfquapp::copy_DEA_DIANN()
grp <- prolfquapp::generate_DEA_reports(lfqdata, GRP2, prot_annot)
prolfquapp::render_DEA(grp[[1]], outpath = GRP2$zipdir, htmlname = "TestTheBest")

for (i in seq_along(grp)) {
  prolfquapp::write_DEA_all(grp[[i]], names(grp)[i], GRP2$zipdir , boxplot = FALSE)
}

library(prozor)
fasta <- prozor::readPeptideFasta(fasta.file[1])
View(annotPep)
pepinfo <- data.frame(peptideSeq = unique(annotPep$Stripped.Sequence))
dim(pepinfo)

annotPep <- prozor::annotatePeptides(pepinfo = pepinfo, fasta = fasta)
annotPep$proteinID |> unique() |> length()
dd <- prepareMatrix(res, peptideID = "peptideSeq" )

prozor:::.greedy2
undebug(greedy)
res <- greedy(dd)
x <- tibble(protein = unlist(res), peptide = names(res))
length(unique(x$protein))

xn <- x |> dplyr::group_by(protein) |> dplyr::summarize(n = n(),.groups = "drop")
xn <- xn |> filter(n > 1)
xn$protein |> unique() |> length()


