library(prolfquapp)
library(tidyverse)

pathtoZip <- "qcdata/FragPipeTMT_o33075_QC_mini.zip"
files <- unzip(pathtoZip, list = TRUE)


ffasta <- grep(".fasta", files$Name , value = TRUE)
fpsm <- grep("psm.tsv", files$Name , value = TRUE)
fdataset <- grep("dataset.csv", files$Name, value = TRUE)

dataset <- read_csv(unz(pathtoZip, filename = fdataset))
psmFuLL$`Modified Peptide`
psmFuLL$
psmFuLL <- read_tsv(unz(pathtoZip,filename =  fpsm))
psm <- prolfqua::tidy_FragPipe_psm(unz(pathtoZip,filename =  fpsm))
nrowPSM <- nrow(psm)

fasta_annot <- get_annot_from_fasta(unz(pathtoZip, filename = ffasta))
psm <- dplyr::inner_join(psm, fasta_annot, by = c(Protein = "fasta.id"), multiple = "all")
stopifnot(nrow(psm) == nrowPSM)

prot_annot <- prolfquapp::dataset_protein_annot(
  psm,
  c("protein_Id" = "Protein"),
  protein_annot = "fasta.header",
  more_columns = "nrPeptides")
prot_annot |> filter(grepl("^zz",protein_Id)) |> head()

psm$qValue <- 1 - psm$PeptideProphet.Probability


atable <- prolfqua::AnalysisConfiguration$new()
atable$ident_Score = "PeptideProphet.Probability"
atable$ident_qValue = "qValue"
atable$fileName = "channel"
atable$opt_rt = "Retention"
atable$opt_mz = "psmcharge"

atable$hierarchy[["protein_Id"]] <- c("Protein")
atable$hierarchy[["peptide_Id"]] <- c("Peptide")
atable$hierarchy[["mod_peptide_Id"]] <- c("Modified.Peptide","Assigned.Modifications")
atable$hierarchy[["precursor"]] <- c("Modified.Peptide","Assigned.Modifications", "Charge")
atable$hierarchy[["Spectrum"]] <- c("Spectrum")

atable$factors[["experiment"]] <- "experiment"
atable$set_response("abundance")

psm$experiment <- "QC"
psm$psmcharge <- psm$Charge
config <- prolfqua::AnalysisConfiguration$new(atable)
adata <- prolfqua::setup_analysis(psm, config)

lfqdata <- prolfqua::LFQData$new(adata, config)
lfqdata$hierarchy_counts()

allmods <- data.frame(allmods = unique(lfqdata$data$mod_peptide_Id ))
tmp <- allmods |> separate(allmods, into = c("modSeq", "Assigned.Modifications"), sep = "~")

x <- tmp |> tidyr::separate_longer_delim("Assigned.Modifications", delim=",")

x <- x |> dplyr::mutate(Assigned.Modifications = trimws(Assigned.Modifications))
x <- x |> dplyr::mutate(Modification = gsub("^[0-9]+","", Assigned.Modifications))
x |> head()
x |> dplyr::filter(grepl("N-term\\(42", x$Modification))
tx <- x$Modification |> table()
tx
plot(tx)

xN <- x |> filter(Modification == "N-term(229.1629)")
nrNtermMod <- xN$modSeq |> unique() |> length()
nrPeptides <- x$modSeq |> unique() |> length()


xL <- x |> filter(Modification == "K(229.1629)")
nrLysMod <- xL$modSeq |> unique() |> length()
nrLysPeptides <-  grepl("L",unique(x$modSeq)) |> sum()



