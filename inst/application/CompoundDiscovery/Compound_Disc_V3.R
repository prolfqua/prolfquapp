library(tidyverse)
library(prolfquapp)
prolfquapp::copy_DEA_DIANN()


files <- dir()
annot <- readxl::read_xlsx(grep("_InputFiles.xlsx$",files, value = TRUE))
colnames(annot) <- make.names(colnames(annot))
str(annot)

annot <- dplyr::select(annot, any_of(c("Study.File.ID", "File.Name","Group","SampleGroup")))
annot$Group |> table() |> t() |> t()
annot <- annot |> dplyr::mutate(CONTROL = case_when(Group == "a" ~ "C" , TRUE ~ "T"))



annot$File.Name <- gsub("\\", "/",annot$File.Name , fixed = TRUE)
annot <- annot |> dplyr::mutate(File.Name = basename(File.Name))

annot <- annot |> dplyr::rename(Name = Study.File.ID)
writexl::write_xlsx(annot, path = "annotation.xlsx")

annotation <- prolfquapp::read_annotation(annot)
annotation$annot
#######

GRP2 <- prolfquapp::make_DEA_config_R6()
dir.create(GRP2$zipdir)


in_file <- file.path(path,"test_Compounds.xlsx")
xd <- prolfquapp::preprocess_CD(in_file, annotation = annotation,.func_massage = )



logger::log_info("run analysis")
grp <- prolfquapp::generate_DEA_reports2(xd$lfqdata, GRP2, xd$protein_annotation, annotation$contrasts)


logger::log_info("write results and html reports")
prolfquapp::write_DEA_all(
  grp,
  "Groups_vs_Controls",
  GRP2$zipdir,
  boxplot = FALSE,
  markdown = "_Grp2Analysis_V2.Rmd")

logger::log_info("write results and summarized experiment")
SE <- prolfquapp::make_SummarizedExperiment(grp)
saveRDS(SE, file = file.path(GRP2$zipdir,
                             paste0("Results_DEA_WU", grp$project_spec$workunit_Id) ,
                             paste0("SummarizedExperiment",".rds") ))

### put all inputs into indir

inputs <- file.path(GRP2$zipdir,
                    paste0("Inputs_DEA_WU", GRP2$project_spec$workunit_Id))
dir.create(inputs)
#prolfquapp::copy_DEA_DIANN(workdir = inputs, run_script = TRUE)
#file.copy(files$data, inputs)
#file.copy(files$fasta, inputs)





