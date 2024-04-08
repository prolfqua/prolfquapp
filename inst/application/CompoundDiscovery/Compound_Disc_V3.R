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

compound <- readxl::read_xlsx("NewOutput/o33926_Areas_IDs_MS1andRT_or_MS2_Compounds.xlsx")
nrow(compound)
mzCloud <- readxl::read_xlsx("NewOutput/o33926_Areas_IDs_MS1andRT_or_MS2_mzCloud.xlsx")
View(mzCloud)

mzVault <- readxl::read_xlsx("NewOutput/o33926_Areas_IDs_MS1andRT_or_MS2_mzVault.xlsx")
head(mzVault)
dim(mzVault)


