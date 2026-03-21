# Example: facade-based DEA using ProteinDataPrep + DEAnalyse
#'
pep <- prolfqua::sim_lfq_data_peptide_config(Nprot = 100)
pep <- prolfqua::LFQData$new(pep$data, pep$config)
pA <- data.frame(protein_Id = unique(pep$data$protein_Id))
pA <- pA |> dplyr::mutate(fasta.annot = paste0(pA$protein_Id, "_description"))
pA <- prolfquapp::ProteinAnnotation$new(
  pep,
  row_annot = pA,
  description = "fasta.annot"
)
GRP2 <- prolfquapp::make_DEA_config_R6()
GRP2$processing_options$diff_threshold <- 0.2
GRP2$processing_options$transform <- "robscale"
#'
contrasts <- c("AVsC" = "group_A - group_Ctrl", BVsC = "group_B - group_Ctrl")
#'
# ---- Data preparation ----
data_prep <- prolfquapp::ProteinDataPrep$new(pep, pA, GRP2)
data_prep$cont_decoy_summary()
data_prep$remove_cont_decoy()
data_prep$aggregate()
data_prep$transform_data()
#'
# ---- Build DEAnalyse with default facade (lm_missing) ----
deanalyse <- data_prep$build_deanalyse(contrasts)
deanalyse$build_default()
stopifnot(nrow(deanalyse$contrast_results[[deanalyse$default_model]]$get_contrasts()) == 200)
#'
# ---- Annotate and filter contrasts ----
deanalyse$get_annotated_contrasts()
deanalyse$filter_contrasts()
#'
# ---- Build additional facades ----
deanalyse$build_facade("lm")
deanalyse$build_facade("limma")
#'
# ---- Access facade plotter ----
cpl <- deanalyse$contrast_results[["lm_missing"]]$get_Plotter()
cpl$volcano()
#'
