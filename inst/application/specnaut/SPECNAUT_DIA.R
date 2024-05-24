library(tidyverse)
library(prolfquapp)


inputdir <- "/Users/witoldwolski/__checkout/protriple/inst/specnaut_DataBench/forWeWat33/20230117_141504_p3000_DIA_tripleProteome_90minGradient_NOimputing_meansBackgroundSignal_Report_DIAservice_v2.xls"
inputdir <- "/Users/witoldwolski/__checkout/protriple/inst/specnaut_DataBench/forWeWat33/20230117_141554_p3000_DIA_tripleProteome_90MinGradien_defaultSettings_Report_DIAservice_v2.xls"
#inputdir <- "/Users/witoldwolski/__checkout/protriple/inst/specnaut_DataBench/nonstaggered/20230426_074926_p3000_DIA_tripleProteome_90MinGradien_defaultSettings_NONstaggered_defaultSettings_Report.xls"

fasta <- "/Users/witoldwolski/__checkout/protriple/inst/diann_DataBench/Triple_Proteome_Exp2_Grad90_V2/misc/fasta/fgcz_tripleProteome_MSV000090837_20221214.fasta"


fastadb <- prozor::readPeptideFasta(fasta)
GRP2 <- prolfquapp::make_DEA_config(ZIPDIR = paste0("specnaut_",zipdir))


GRP2$pop$transform = "none"
###
dir.create(GRP2$zipdir)
###
# reading foreign data
REPEATED <- TRUE

prec <- read_tsv(inputdir)
prec$R.FileName |> unique()
prec <- prec |> select(R.FileName,
                       PG.ProteinAccessions,
                       PG.ProteinDescriptions,
                       EG.StrippedSequence,
                       EG.ModifiedSequence,
                       EG.PrecursorId,
                       EG.Qvalue,
                       FG.Quantity)

prec <- prec |> tidyr::separate(PG.ProteinAccessions, "protID", sep=";", remove=FALSE)



hist(log(prec$FG.Quantity))
hist(prec$FG.Quantity[prec$FG.Quantity < 2000], breaks = 1000 , ylim = c(0,100))
abline(v = 50, col=2)
prec <- prec |> filter(FG.Quantity > 50)
hist(prec$FG.Quantity[prec$FG.Quantity < 2000], breaks = 1000 , ylim = c(0,100))

peptide <- prec |>
  group_by(R.FileName, PG.ProteinAccessions, protID, PG.ProteinDescriptions, EG.StrippedSequence) |>
  summarize(n = n(), Quantity = sum(FG.Quantity, na.rm = TRUE), Qvalue = min(EG.Qvalue, na.rm = TRUE)) |>
  ungroup()

is_grouped_df(peptide)


nrpepprot <- peptide |>
  select(PG.ProteinAccessions, EG.StrippedSequence) |>
  distinct() |>
  group_by(PG.ProteinAccessions) |>
  summarize(nrPeptides = n())

plot(table(nrpepprot$nrPeptides), xlim = c(1,100))


peptide <- dplyr::inner_join(nrpepprot, peptide, multiple = "all")
peptide$raw.file <- peptide$R.FileName

prot_annot <- prolfquapp::dataset_protein_annot(
  peptide,
  "protID",
  protein_annot = "PG.ProteinDescriptions",
  more_columns = c("PG.ProteinAccessions","nrPeptides"))

annot <- data.frame(raw.file = peptide$R.FileName |> unique())
annot <- annot |> dplyr::mutate(Name = gsub(".+_LFQ", "LFQ", raw.file))
annot <- annot |> dplyr::mutate(GroupingVar = dplyr::case_when(grepl("_A_", Name) ~ "A", TRUE ~ "B"))
annot$CONTROL <- ifelse(annot$GroupingVar == "A", "C", "T")


annot$raw.file[ !annot$raw.file %in% sort(unique(peptide$raw.file)) ]
nr <- sum(annot$raw.file %in% sort(unique(peptide$raw.file)))
logger::log_info("nr : ", nr, " files annotated")

peptide <- dplyr::inner_join(annot, peptide, multiple = "all")
peptide$identScore <- -log10(peptide$Qvalue)
GRP2 <- prolfquapp::dataset_set_factors_deprecated(annot, GRP2)

atable <- prolfqua::AnalysisTableAnnotation$new()
atable$fileName = "raw.file"
atable$hierarchy[["protein_Id"]] <- c("protID")
atable$hierarchy[["peptide_Id"]] <- c("EG.StrippedSequence")
atable$set_response("Quantity")
atable$ident_qValue <- "Qvalue"
atable$ident_Score <- "identScore"
atable$hierarchyDepth <- 1

res <- prolfquapp::dataset_set_factors(atable, peptide)
atable <- res$atable
atable$factor_keys()

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

prot_annot <- prolfquapp::ProteinAnnotation$new(lfqdata = lfqdata, prot_annot)

grp <- prolfquapp::generate_DEA_reports(lfqdata, GRP2, prot_annot)

prolfquapp::render_DEA(grp[[1]], outpath = GRP2$zipdir, htmlname = "TestTheBest")

for (i in seq_along(grp)) {
  prolfquapp::write_DEA_all(grp[[i]], names(grp)[i], GRP2$zipdir , boxplot = FALSE)
}

