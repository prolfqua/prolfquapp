# find files
get_files <- function(path){
  diann.path <- grep("report\\.tsv$|diann-output\\.tsv", dir(path = path, recursive = TRUE, full.names = TRUE), value = TRUE)
  fasta.files <- grep("*.fasta$", dir(path = path, recursive = TRUE,full.names = TRUE), value = TRUE)
  if (any(grepl("database[0-9]*.fasta$", fasta.files))) {
    fasta.files <- grep("database[0-9]*.fasta$", fasta.files, value = TRUE)
  }
  if (length(fasta.files) == 0) {
    logger::log_error("No fasta file found!")
    stop()
  }
  return(list(data = diann.path, fasta = fasta.files))
}





preprocess_DIANN <- function(quant_data,
                             fasta_file,
                             annotation,
                             q_value = 0.01,
                             nrPeptides = 1){

  annot <- annotation$annot
  atable <- annotation$atable
  annot <- annot |> dplyr::mutate(
    raw.file = gsub("^x|.d.zip$|.raw$","",
                    (basename(annot[[atable$fileName]]))
    ))


  peptide <- prolfquapp::read_DIANN_output(
    diann.path = quant_data,
    fasta.file = fasta_file,
    nrPeptides = 1,
    Q.Value = q_value)

  prot_annot <- prolfquapp::dataset_protein_annot(
    peptide,
    c("protein_Id" = "Protein.Group"),
    protein_annot = "fasta.header",
    more_columns = c("nrPeptides", "fasta.id", "Protein.Group.2")
  )

  nr <- sum(annot$raw.file %in% sort(unique(peptide$raw.file)))
  logger::log_info("nr : ", nr, " files annotated out of ", length(unique(peptide$raw.file)))

  atable$fileName = "raw.file"
  atable$hierarchy[["protein_Id"]] <- c("Protein.Group")
  atable$hierarchy[["peptide_Id"]] <- c("Stripped.Sequence")
  atable$set_response("Peptide.Quantity")
  atable$hierarchyDepth <- 1

  peptide <- dplyr::inner_join(annot, peptide, multiple = "all")
  config <- prolfqua::AnalysisConfiguration$new(atable)
  adata <- prolfqua::setup_analysis(peptide, config)
  lfqdata <- prolfqua::LFQData$new(adata, config)
  lfqdata$remove_small_intensities()
  protAnnot <- prolfqua::ProteinAnnotation$new(lfqdata , prot_annot)
  return(list(lfqdata = lfqdata , protein_annotation = protAnnot))
}






make_DEA_report2 <- function(lfqdata,
                             protAnnot,
                             GRP2
) {
  ### Do some type of data normalization (or do not)

  transformed <- prolfquapp::transform_lfqdata(
    lfqdata,
    method = GRP2$pop$transform,
    internal = GRP2$pop$internal
  )


  allProt <- nrow( protAnnot$row_annot )
  GRP2$RES <- list()
  GRP2$RES$Summary <- data.frame(
    totalNrOfProteins = allProt,
    percentOfContaminants = round(protAnnot$annotate_contaminants(GRP2$processing_options$pattern_contaminants)/allProt * 100 , digits = 2),
    percentOfFalsePositives  = round(protAnnot$annotate_decoys(GRP2$processing_options$pattern_decoys)/allProt * 100 , digits = 2),
    NrOfProteinsNoDecoys = protAnnot$nr_clean()
  )
  GRP2$RES$rowAnnot <- protAnnot

  if (GRP2$processing_options$remove_cont || GRP2$processing_options$remove_cont) {
    lfqdata <- lfqdata$get_subset(protAnnot$clean(contaminants = GRP2$processing_options$remove_cont, decoys = GRP2$processing_options$remove_decoys))
    transformed <- transformed$get_subset(protAnnot$clean())
    logger::log_info(
      paste0("removing contaminants and reverse sequences with patterns:",
             GRP2$processing_options$pattern_contaminants,
             GRP2$processing_options$pattern_decoys ))
  }

  lfqdata$rename_response("protein_abundance")
  transformed$rename_response("normalized_protein_abundance")

  GRP2$RES$lfqData <- lfqdata
  GRP2$RES$transformedlfqData <- transformed

  ################## Run Modelling ###############
  factors <- transformed$config$table$factor_keys()[!grepl("^control", transformed$config$table$factor_keys() , ignore.case = TRUE)]

  if (is.null(GRP2$processing_options$interaction) || !GRP2$processing_options$interaction ) {
    formula <- paste0(transformed$config$table$get_response(), " ~ ",
                      paste(factors, collapse = " + "))
  } else {
    formula <- paste0(transformed$config$table$get_response(), " ~ ",
                      paste(factors, collapse = " * "))
  }

  message("FORMULA :",  formula)
  GRP2$RES$formula <- formula
  formula_Condition <-  prolfqua::strategy_lm(formula)
  # specify model definition
  modelName  <- "Model"

  mod <- prolfqua::build_model(
    transformed,
    formula_Condition,
    subject_Id = transformed$config$table$hierarchy_keys() )

  logger::log_info("fitted model with formula : {formula}")
  GRP2$RES$models <- mod

  contr <- prolfqua::Contrasts$new(mod, GRP2$pop$Contrasts,
                                   modelName = "Linear_Model")
  conrM <- prolfqua::ContrastsModerated$new(
    contr)

  if (is.null(GRP2$processing_options$missing) || GRP2$processing_options$missing ) {
    mC <- prolfqua::ContrastsMissing$new(
      lfqdata = transformed,
      contrasts = GRP2$pop$Contrasts,
      modelName = "Imputed_Mean")

    conMI <- prolfqua::ContrastsModerated$new(
      mC)

    res <- prolfqua::merge_contrasts_results(conrM, conMI)
    GRP2$RES$contrMerged <- res$merged
    GRP2$RES$contrMore <- res$more

  } else {
    GRP2$RES$contrMerged <- conrM
    GRP2$RES$contrMore <- NULL
  }


  datax <- GRP2$RES$contrMerged$get_contrasts()
  datax <- dplyr::inner_join(GRP2$RES$rowAnnot$row_annot, datax, multiple = "all")
  GRP2$RES$contrastsData  <- datax

  datax <- datax |>
    dplyr::filter(.data$FDR < GRP2$processing_options$FDR_threshold &
                    abs(.data$diff) > GRP2$processing_options$diff_threshold )
  GRP2$RES$contrastsData_signif <- datax
  return(GRP2)
}


generate_DEA_reports2 <- function(lfqdata, GRP2, prot_annot, Contrasts) {
  # Compute all possible 2 Grps to avoid specifying reference.
  res <- list()
  GRP2$pop$Contrasts <- Contrasts
  logger::log_info("CONTRAST : ", paste( GRP2$pop$Contrasts, collapse = " "))
  lfqwork <- lfqdata$get_copy()
  # lfqwork$data <- lfqdata$data |> dplyr::filter(.data$Group_ %in% levels$Group_)
  grp2 <- make_DEA_report2(
    lfqwork,
    prot_annot,
    GRP2)
  return(list("Groups_vs_Controls" = grp2) )
}

