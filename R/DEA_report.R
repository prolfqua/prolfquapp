# Load necessary libraries

DEAresult <- R6::R6Class(
  "DEAresult",
  public = list(
  )
)

#' will replace make_DEA_report
#' @export
#'
#' @examples
#' # example code
#'
#' pep <- prolfqua::sim_lfq_data_protein_config(Nprot = 100)
#' pep <- prolfqua::LFQData$new(pep$data, pep$config)
#' pA <- data.frame(protein_Id = unique(pep$data$protein_Id))
#' pA <- pA |> dplyr::mutate(fasta.annot = paste0(pA$protein_Id, "_description"))
#' pA <- prolfquapp::ProteinAnnotation$new(pep,row_annot = pA ,description = "fasta.annot")
#' GRP2 <- prolfquapp::make_DEA_config_R6()
#' GRP2$processing_options$transform <- "robscale"
#' pep$factors()
#' GRP2$pop <- list(Contrasts = c("AVsC" = "group_A - group_Ctrl", BVsC = "group_B - group_Ctrl"))
#' DEAnalyse$debug("transform_data")
#' deanalyse <- DEAnalyse$new(pep, pA, GRP2)
#' deanalyse$cont_decoy_summary()
#' deanalyse$GRP2$processing_options$remove_cont = TRUE
#' deanalyse$remove_cont_decoy()
#' deanalyse$transform_data()
#' deanalyse$create_model_formula()
#' deanalyse$build_model()
#' deanalyse$compute_contrasts()
#' deanalyse$filter_contrasts()
#' deanalyse$RES
DEAnalyse <- R6::R6Class(
  "DEAnalyse",
  public = list(
    #' @field lfqdata LFQData
    lfqdata = NULL,
    #' @field rowAnnot ProteinAnnotation
    rowAnnot = NULL,
    #'@field GRP2 ProlfquAppConfig
    GRP2 = NULL,
    #' @field reference_proteins reference proteins to use for internal normalization
    reference_proteins = NULL,
    decoy_summary = NULL,
    lfqData = NULL,
    transformedlfqData = NULL,
    formula = character(),
    models = data.frame(),
    contrasts = list(),

    #' initialize
    #' @param lfqdata lfqdata
    #' @param rowAnnot ProteinAnnotation
    #' @param GRP2 ProlfquAppConfig
    initialize = function(lfqdata, rowAnnot, GRP2) {
      self$lfqdata <- lfqdata
      self$rowAnnot <- rowAnnot
      self$GRP2 <- GRP2
      self$RES <- list()
    },

    cont_decoy_summary = function() {
      self$rowAnnot$get_summary()
    },
    remove_cont_decoy = function() {
        self$lfqdata <- self$lfqdata$get_subset(self$rowAnnot$clean(
          contaminants = self$GRP2$processing_options$remove_cont,
          decoys = self$GRP2$processing_options$remove_decoys
        ))
        logger::log_info("removing contaminants and reverse sequences with patterns: ",
                         self$GRP2$processing_options$pattern_contaminants,
                         self$GRP2$processing_options$pattern_decoys)
    },
    transform_data = function() {
      transformed <- prolfquapp::transform_lfqdata(
        self$lfqdata,
        method = self$GRP2$processing_options$transform,
        internal = self$reference_proteins
      )
      self$lfqdata$rename_response("protein_abundance")
      transformed$rename_response("normalized_protein_abundance")
      self$transformedlfqData <- transformed
      invisible(transformed)
    },
    # static
    create_model_formula = function(prlconfig, interaction = FALSE) {
      factors <- prlconfig$table$factor_keys_depth()[
        !grepl("^control", prlconfig$table$factor_keys_depth() , ignore.case = TRUE)
      ]
      # model with or without interactions
      if (interaction ) {
        formula <- paste0(prlconfig$table$get_response(), " ~ ",
                          paste(factors, collapse = " + "))
      } else {
        formula <- paste0(prlconfig$table$get_response(), " ~ ",
                          paste(factors, collapse = " * "))
      }
      logger::log_info("fitted model with formula : {formula}")
      return(formula)
    },
    #static
    build_model = function(transformedlfqData, interaction = FALSE) {
      formula <- create_model_formula(transformedlfqData$config, interaction = interaction)
      formula_Condition <-  prolfqua::strategy_lm(formula)
      models <- prolfqua::build_model(
        transformed,
        formula_Condition,
        subject_Id = transformed$config$table$hierarchy_keys() )
      return(models)
    },
    #static
    compute_contrasts = function(models, contrasts, modelName = "Linear_Model") {
      contr <- prolfqua::Contrasts$new(models,
                                       contrasts,
                                       modelName = "Linear_Model")
      conrM <- prolfqua::ContrastsModerated$new(contr)
      return(conrM)
    },
    #static
    compute_contrasts_missing = function(transformedlfqData, contrasts){
      mC <- prolfqua::ContrastsMissing$new(
        lfqdata = transformedlfqData,
        contrasts = contrasts,
        modelName = "Imputed_Mean")
      conMI <- prolfqua::ContrastsModerated$new(mC)
      return(conMI)
    },
    #static
    filter_contrasts = function(contr, FDR_threshold = 0.1, diff_threshol = 1){
      datax <- contr$get_contrasts()
      datax <- datax |>
        dplyr::filter(.data$FDR < FDR_threshold &
                        abs(.data$diff) > diff_threshold )
      invisible(datax)
    })
)

