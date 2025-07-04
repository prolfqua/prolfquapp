# example code
#'
pep <- prolfqua::sim_lfq_data_peptide_config(Nprot = 10)
#'
pep <- prolfqua::LFQData$new(pep$data, pep$config)
pA <- data.frame(protein_Id = unique(pep$data$protein_Id))
pA <- pA |> dplyr::mutate(fasta.annot = paste0(pA$protein_Id, "_description"))
pA <- prolfquapp::ProteinAnnotation$new(pep,row_annot = pA ,description = "fasta.annot")
GRP2 <- prolfquapp::make_DEA_config_R6()
GRP2$processing_options$diff_threshold = 0.2
#'
GRP2$processing_options$transform <- "robscale"
pep$factors()
contrasts = c("AVsC" = "group_A - group_Ctrl", BVsC = "group_B - group_Ctrl")
# prolfquapp::DEAnalyse$debug("build_model_glm_peptide")
deanalyse <- prolfquapp::DEAnalyse$new(pep, pA, GRP2, contrasts)

xx <- deanalyse$build_model_glm_peptide()



deanalyse$lfq_data_peptide$hierarchy_counts()
deanalyse$cont_decoy_summary()
deanalyse$prolfq_app_config$processing_options$remove_cont = TRUE
deanalyse$remove_cont_decoy()
deanalyse$aggregate()
pl <- deanalyse$get_aggregation_plots(exp_nr_children = 2)
print(pl$plots[[3]])

deanalyse$transform_data()
mod <- deanalyse$build_model_linear_protein()
contlm <- deanalyse$get_contrasts_linear_protein()
#'

merged <- deanalyse$get_contrasts_merged_protein()
stopifnot(nrow(merged$get_contrasts()) == 20)

deanalyse$lfq_data$complete_cases()
str <- deanalyse$get_strategy_glm_prot()
x <- deanalyse$build_model_glm_protein()

x$modelDF$linear_model[[1]]
x$modelDF$linear_model[[2]]
x$modelDF$linear_model[[3]]


xprot <- deanalyse$get_contrasts_glm_protein()
xprot$get_contrasts()

#'
stopifnot(nrow(merged$get_contrasts()) == 20)
#'
xpep <- deanalyse$get_contrasts_glm_peptide()

xprot$get_Plotter()$volcano()
xpep$get_Plotter()$volcano()
sr <- deanalyse$lfq_data_peptide$get_Summariser()
#'
#'
#'
deanalyse$filter_contrasts()
#'
xd <- deanalyse$filter_data()
xd <- deanalyse$contrasts_to_Grob()
bb <- deanalyse$get_boxplots()
bx <- deanalyse$get_boxplots_contrasts()
dev.off()
grid::grid.draw(bx$bxpl_grobs[[1]])
# deanalyse$write_boxplots_contrasts("test.pdf")
#'
