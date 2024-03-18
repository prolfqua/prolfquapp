prolfquapp::copy_DEA_DIANN()



### Annotation
path <- "input_ex/"
annot <- readxl::read_xlsx(file.path(path,"o33926_input_files.xlsx"))
colnames(annot) <- make.names(colnames(annot))
annot <- annot |> dplyr::select(StudyFile = Study.File.ID, File.Name, Material = Meterial,Gender,Genotype)
annot <- annot |> tidyr::unite(Group, Material, Gender, Genotype, remove = FALSE)

annot$File.Name <- gsub("\\", "/",annot$File.Name , fixed = TRUE)
annot <- annot |> dplyr::mutate(File.Name = basename(File.Name))
annot <- annot |> dplyr::filter(Material != "none")

contrast <- c(
  CFvsMR_plasma = "(G_plasma_Female_CF + G_plasma_Male_CF)/2 - (G_plasma_Female_MR + G_plasma_Male_MR)/2",
  CVvsMR_male_plasma = "G_plasma_Male_CF - G_plasma_Male_MR",
  CVvsMR_female_plasma = "G_plasma_Female_CF - G_plasma_Female_MR",

  CFvsMR_urine = "(G_urine_Female_CF + G_urine_Male_CF)/2 - (G_urine_Female_MR + G_urine_Male_MR)/2",
  CVvsMR_male_urine = "G_urine_Male_CF - G_urine_Male_MR",
  CVvsMR_female_urine = "G_urine_Female_CF - G_urine_Female_MR"

)

annot$ContrastName <- c(names(contrast), rep(NA, nrow(annot) - length(contrast)))
annot$Contrast <- c(contrast,rep(NA, nrow(annot) - length(contrast)))

annot <- annot |> dplyr::rename(Name = StudyFile)
writexl::write_xlsx(annot, path = "annotation.xlsx")

annotation <- prolfquapp::read_annotation(annot)
#######

GRP2 <- prolfquapp::make_DEA_config_R6(ZIPDIR = "o33926_Metabo")
dir.create(GRP2$zipdir)

in_file <- file.path(path,"o33926_Areas_IDs_MS1andRT_or_MS2.xlsx")
xd <- prolfquapp::preprocess_CD(in_file, annotation = annotation)

logger::log_info("run analysis")
grp <- prolfquapp::generate_DEA_reports2(xd$lfqdata, GRP2, xd$protein_annotation, annotation$contrasts)


logger::log_info("write results and html reports")
prolfquapp::write_DEA_all(grp[[1]], names(grp)[1], GRP2$zipdir , boxplot = FALSE)

logger::log_info("write results and summarized experiment")
undebug( prolfquapp::make_SummarizedExperiment)
SE <- prolfquapp::make_SummarizedExperiment(grp[[1]])
saveRDS(SE, file = file.path(GRP2$zipdir, paste0("DE_", names(grp)[1]) , paste0("SummarizedExperiment",".rds") ))





