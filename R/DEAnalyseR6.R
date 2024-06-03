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
#' contrasts = c("AVsC" = "group_A - group_Ctrl", BVsC = "group_B - group_Ctrl")
#' #DEAnalyse$debug("filter_contrasts")
#' deanalyse <- DEAnalyse$new(pep, pA, GRP2, contrasts)
#' deanalyse$cont_decoy_summary()
#' deanalyse$GRP2$processing_options$remove_cont = TRUE
#' deanalyse$remove_cont_decoy()
#' deanalyse$transform_data()
#' deanalyse$create_model_formula()
#' mod <- deanalyse$build_model_linear()
#' class(mod)
#' contlm <- deanalyse$compute_contrasts_linear()
#' class(contlm)
#' contlm$modelName == "modelLinear_moderated"
#'
#' xd <-deanalyse$filter_contrasts(FDR_threshold = 0.3, diff_threshold = 0.1)
#' xd <- deanalyse$filter_data(FDR_threshold = 0.3, diff_threshold = 0.1)
#' dim(xd)
#' bb <- deanalyse$get_boxplots()
#' bb$boxplot[[1]]
DEAnalyse <- R6::R6Class(
  "DEAnalyse",
  public = list(
    #' @field lfqdata LFQData
    lfqdata = NULL,
    #' @field rowAnnot ProteinAnnotation
    rowAnnot = NULL,
    #' @field contrasts vector with contrasts
    contrasts = character(),
    #'@field GRP2 ProlfquAppConfig
    GRP2 = NULL,


    #' @field reference_proteins reference proteins to use for internal normalization
    reference_proteins = NULL, # to use for internal calibration
    #' @field decoy_summary decoy_summary
    decoy_summary = NULL,
    #' @field lfqData lfqData
    lfqData = NULL,
    #' @field lfqData_transformed todo
    lfqData_transformed = NULL,
    #' @field lfqData_subset todo
    lfqData_subset = NULL,

    #' @field formula todo
    formula = character(),
    #' @field models todo
    models = list(),
    #' @field contrastRes todo
    contrastRes = list(),

    #' @description
    #' initialize
    #' @param lfqdata lfqdata
    #' @param rowAnnot ProteinAnnotation
    #' @param GRP2 ProlfquAppConfig
    #' @param contrasts vector with contrasts
    initialize = function(lfqdata, rowAnnot, GRP2, contrasts) {
      stopifnot("LFQData" %in% class(lfqdata))
      stopifnot("ProteinAnnotation" %in% class(rowAnnot))
      stopifnot("ProlfquAppConfig" %in% class(GRP2))
      stopifnot(length(contrasts) >= 1)
      self$lfqdata <- lfqdata
      self$rowAnnot <- rowAnnot
      self$GRP2 <- GRP2
      self$contrasts <- contrasts
    },
    #' @description
    #' count number of decoys
    cont_decoy_summary = function() {
      self$rowAnnot$get_summary()
    },
    #' @description
    #' remove contaminants and decoys
    remove_cont_decoy = function() {
      self$lfqdata <- self$lfqdata$get_subset(self$rowAnnot$clean(
        contaminants = self$GRP2$processing_options$remove_cont,
        decoys = self$GRP2$processing_options$remove_decoys
      ))
      logger::log_info("removing contaminants and reverse sequences with patterns: ",
                       self$GRP2$processing_options$pattern_contaminants,
                       self$GRP2$processing_options$pattern_decoys)
    },
    #' @description
    #' transform data
    transform_data = function() {
      transformed <- prolfquapp::transform_lfqdata(
        self$lfqdata,
        method = self$GRP2$processing_options$transform,
        internal = self$reference_proteins
      )
      self$lfqdata$rename_response("protein_abundance")
      transformed$rename_response("normalized_protein_abundance")
      self$lfqData_transformed <- transformed
      invisible(transformed)
    },
    #' @description
    #' static create model formula
    create_model_formula = function() {
      prlconfig <- self$lfqData_transformed$config
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
    #' @description
    #' fit linear model
    #' @param modelName modelLinear
    build_model_linear = function(modelName = "modelLinear") {
      formula <- self$create_model_formula()
      formula_Condition <-  prolfqua::strategy_lm(formula)
      models <- prolfqua::build_model(
        self$lfqData_transformed,
        formula_Condition)
      self$models[[modelName]] <- models
      return(models)
    },
    #' @description
    #' fit generalized linear model
    build_model_glm = function(){
      self$models[["glm"]] <- models
    },
    #' @description
    #' compute contrasts linear
    #' @param modelName of model default modelLinear
    compute_contrasts_linear = function(modelName = "modelLinear") {
      contr <- prolfqua::Contrasts$new(self$models[[modelName]],
                                       self$contrasts,
                                       modelName = modelName)
      conrM <- prolfqua::ContrastsModerated$new(contr)
      self$contrastRes[[modelName]] <- conrM
      invisible(conrM)
    },
    #' @description
    #' compute missing contrasts
    #' @param modelName Imputed_Mean
    compute_contrasts_missing = function(modelName = "Imputed_Mean"){
      mC <- prolfqua::ContrastsMissing$new(
        lfqdata = self$lfqData_transformed,
        contrasts = self$contrasts,
        modelName = modelName)
      conMI <- prolfqua::ContrastsModerated$new(mC)
      self$contrastRes[[modelName]] <- conMI
      return(conMI)
    },
    #' @description
    #' filter contrasts for therehold
    #' @param FDR_threshold FDR threshold
    #' @param diff_threshol diff threshold
    #' @param modelName modelLinear
    filter_contrasts = function(FDR_threshold = 0.1, diff_threshold = 1,
                                modelName = "modelLinear"){
      datax <- self$contrastRes[[modelName]]$get_contrasts()
      datax <- datax |>
        dplyr::filter(.data$FDR < FDR_threshold &
                        abs(.data$diff) > diff_threshold )
      invisible(datax)
    },
    #' @description
    #' filter transformed lfq data for significant proteins.
    #' @param FDR_threshold FDR threshold
    #' @param diff_threshol diff threshold
    #' @param modelName modelLinear
    filter_data = function(FDR_threshold = 0.1,
                           diff_threshold = 1,
                           modelName = "modelLinear")
    {
      dx <- self$filter_contrasts(FDR_threshold, diff_threshold, modelName = modelName)
      self$lfqData_subset <- self$lfqData_transformed$get_subset(dx)
      return(self$lfqData_subset)
    },
    #' @description
    #' create boxplots
    #'
    get_boxplots = function(){
      return(self$lfqData_subset$get_Plotter()$boxplots())
    }

  )


)

