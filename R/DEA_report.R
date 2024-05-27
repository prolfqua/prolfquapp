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
#' # DEAnalyse$debug("transform_data")
#' deanalyse <- DEAnalyse$new(pep, pA, GRP2)
#' deanalyse$cont_decoy_summary()
#' deanalyse$GRP2$processing_options$remove_cont = TRUE
#' deanalyse$remove_cont_decoy()
#' deanalyse$transform_data()
#' deanalyse$create_model_formula()
#' deanalyse$build_model_linear()
#' deanalyse$compute_contrasts_linear()
#' #deanalyse$filter_contrasts()
#' deanalyse$RES
DEAnalyse <- R6::R6Class(
  "DEAnalyse",
  public = list(
    #' @field lfqdata LFQData
    lfqdata = NULL,
    #' @field rowAnnot ProteinAnnotation
    rowAnnot = NULL,
    contrasts = character(),
    #'@field GRP2 ProlfquAppConfig
    GRP2 = NULL,


    #' @field reference_proteins reference proteins to use for internal normalization
    reference_proteins = NULL, # to use for internal calibration
    decoy_summary = NULL,
    lfqData = NULL,
    transformedlfqData = NULL,
    formula = character(),
    models = list(),
    contrastRes = list(),

    #' initialize
    #' @param lfqdata lfqdata
    #' @param rowAnnot ProteinAnnotation
    #' @param GRP2 ProlfquAppConfig
    initialize = function(lfqdata, rowAnnot, GRP2, contrasts) {
      self$lfqdata <- lfqdata
      self$rowAnnot <- rowAnnot
      self$GRP2 <- GRP2
      self$contrasts
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
    create_model_formula = function() {
      prlconfig <- self$transformedlfqData$config
      interaction <- self$GRP2$processing_options$interaction
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
      self$formula <- formula
      return(formula)
    },
    build_model_linear = function(modelName = "modelLinear") {
      formula <- self$create_model_formula()
      formula_Condition <-  prolfqua::strategy_lm(formula)
      models <- prolfqua::build_model(
        self$transformedlfqData,
        formula_Condition)
      self$models[[modelName]] <- models
      return(models)
    },
    build_model_glm = function(){
      self$models[["glm"]] <- models
    },
    compute_contrasts_linear = function(modelName = "modelLinear") {
      contr <- prolfqua::Contrasts$new(self$models[["modelLinear"]],
                                       self$contrasts,
                                       modelName = modelName)
      conrM <- prolfqua::ContrastsModerated$new(contr)
      self$contrastRes[[modelName]] <- conrM
      invisible(conrM)
    },
    #static
    compute_contrasts_missing = function(modelName = "Imputed_Mean"){
      mC <- prolfqua::ContrastsMissing$new(
        lfqdata = transformedlfqData,
        contrasts = contrasts,
        modelName = modelName)
      conMI <- prolfqua::ContrastsModerated$new(mC)
      self$contrastRes[[modelName]] <- conMI
      return(conMI)
    },
    #static
    filter_contrasts = function(FDR_threshold = 0.1, diff_threshol = 1,
                                modelName = "modelLinear"){
      datax <- self$contrastRes[[modelName]]$get_contrasts()
      datax <- datax |>
        dplyr::filter(.data$FDR < FDR_threshold &
                        abs(.data$diff) > diff_threshold )
      invisible(datax)
    })
)

