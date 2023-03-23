library(rlang)
ymlfile <- "config.yaml"
GRP2 <- prolfquapp::read_yaml(ymlfile)

###
dir.create(GRP2$zipdir)
path <- "diann-output.tsv"

report2 <- prolfqua::diann_read_output(path,nrPeptides = 2, Q.Value = 0.01)


report2$raw.file <- gsub("^x|.d.zip$|.d$|.raw$|.mzML$","",basename(gsub("\\\\","/",report2$File.Name)))
peptide <- prolfqua::diann_output_to_peptide(report2)
peptide$Protein.Group.2 <- sapply(peptide$Protein.Group, function(x){ unlist(strsplit(x, " "))[1]} )
# we need to add the fasta.header information.

library(seqinr)
#### code specific to Uniprot

fasta.file <- dir(".", pattern = "*.fasta")

fasta <- prozor::readPeptideFasta(fasta.file)
fasta_annot <- prolfqua::matrix_to_tibble(
  data.frame(annot = sapply(fasta, seqinr::getAnnot)), preserve_row_names = NULL
)
fasta_annot <- fasta_annot |> tidyr::separate(.data$annot,
                                              c("proteinname","fasta.header"),
                                              sep = "\\s", extra = "merge")
fasta_annot <- fasta_annot |> dplyr::mutate(proteinname = gsub(">","", .data$proteinname) )
if (TRUE) {
  fasta_annot <- fasta_annot |> dplyr::mutate(proteinname = gsub(".+\\|(.+)\\|.*","\\1", .data$proteinname) )
}

# remove duplicated id's
fasta_annot <- fasta_annot[!duplicated(fasta_annot$proteinname),]

stopifnot(mean(peptide$Protein.Group.2 %in% fasta_annot$proteinname) > 0.9)
# add fasta headers.
peptide <- dplyr::left_join(peptide, fasta_annot, by = c( Protein.Group.2 = "proteinname"))
mean(is.na(peptide$fasta.header))


dsf <- "dataset.csv"
annot <- read.csv(dsf)
annot <- data.frame(lapply(annot, as.character))
annot <- annot |> dplyr::mutate(
  raw.file = gsub("^x|.d.zip$|.raw$","",
                  (basename(annot$Relative.Path))
  )
)
GRP2 <- prolfquapp::dataset_extract_contrasts(annot, GRP2)

annot$raw.file[ !annot$raw.file %in% sort(unique(peptide$raw.file)) ]
nr <- sum(annot$raw.file %in% sort(unique(peptide$raw.file)))

logger::log_info("nr : ", nr, " files annotated")
annot$Relative.Path <- NULL
peptide <- dplyr::inner_join(annot, peptide, multiple = "all")


atable <- prolfqua::AnalysisTableAnnotation$new()
atable$fileName = "raw.file"
atable$hierarchy[["protein_Id"]] <- c("Protein.Group.2")
atable$hierarchy[["peptide_Id"]] <- c("Stripped.Sequence")
atable$set_response("Peptide.Quantity")
atable$hierarchyDepth <- 1

res <- prolfquapp::dataset_set_factors(atable, peptide)
atable <- res$atable
peptide <- res$msdata

prot_annot <- prolfquapp::dataset_protein_annot(
  peptide,
  "Protein.Group.2",
  protein_annot = "fasta.header")

# Preprocess data - aggregate proteins.
config <- prolfqua::AnalysisConfiguration$new(atable)
adata <- prolfqua::setup_analysis(peptide, config)
lfqdata <- prolfqua::LFQData$new(adata, config)
lfqdata$remove_small_intensities()
GRP2$pop$nrPeptides <- 2
logger::log_info("AGGREGATING PEPTIDE DATA!")

lfqdata <- prolfquapp::aggregate_data(lfqdata, agg_method = GRP2$pop$aggregate)
logger::log_info("data aggregated: {GRP2$pop$aggregate}.")
lfqdata$factors()
logger::log_info("END OF DATA TRANSFORMATION.")


prolfquapp::copy_DEA_FragPipe_DIA()
undebug(prolfquapp::generate_DEA_reports)
grp <- prolfquapp::generate_DEA_reports(lfqdata, GRP2, prot_annot)
for (i in seq_along(grp)) {
  prolfquapp::write_DEA_all(grp[[i]], names(grp)[i], GRP2$zipdir )
}
