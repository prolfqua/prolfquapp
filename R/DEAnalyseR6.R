custom_round <- function(arr) {
  cr <- function(x){
    if (x == 0) {
      return(0)
    } else if (abs(x) >= 1) {
      return(round(x, 2))
    } else {
      return(signif(x, 2))
    }}
  return(vapply(arr, cr, numeric(1)))
}

#' will replace make_DEA_report
#' @export
#'
#' @examples
#' # example code
#'
#' pep <- prolfqua::sim_lfq_data_protein_config(Nprot = 100)
#'
#' pep <- prolfqua::LFQData$new(pep$data, pep$config)
#' pA <- data.frame(protein_Id = unique(pep$data$protein_Id))
#' pA <- pA |> dplyr::mutate(fasta.annot = paste0(pA$protein_Id, "_description"))
#' pA <- prolfquapp::ProteinAnnotation$new(pep,row_annot = pA ,description = "fasta.annot")
#' GRP2 <- prolfquapp::make_DEA_config_R6()
#' GRP2$processing_options$diff_threshold = 0.2
#'
#' GRP2$processing_options$transform <- "robscale"
#' pep$factors()
#' contrasts = c("AVsC" = "group_A - group_Ctrl", BVsC = "group_B - group_Ctrl")
#' # DEAnalyse$debug("filter_contrasts")
#' deanalyse <- DEAnalyse$new(pep, pA, GRP2, contrasts)
#' deanalyse$cont_decoy_summary()
#' deanalyse$prolfq_app_config$processing_options$remove_cont = TRUE
#' deanalyse$remove_cont_decoy()
#' deanalyse$transform_data()
#' deanalyse$create_model_formula()
#' mod <- deanalyse$build_model_linear()
#' contlm <- deanalyse$get_contrasts_linear()
#'
#' deanalyse$filter_contrasts()
#'
#' xd <- deanalyse$filter_data()
#' xd <- deanalyse$contrasts_to_Grob()
#' bb <- deanalyse$get_boxplots()
#' bx <- deanalyse$get_boxplots_contrasts()
#' grid::grid.draw(bx$bxpl_grobs[[1]])
#' # deanalyse$write_boxplots_contrasts("test.pdf")
#' deanalyse$build_model_glm_protein()
#'
DEAnalyse <- R6::R6Class(
  "DEAnalyse",
  public = list(
    #' @field lfqdata LFQData
    lfqdata = NULL,
    #' @field rowAnnot ProteinAnnotation
    rowAnnot = NULL,
    #' @field contrasts vector with contrasts
    contrasts = character(),
    #'@field prolfq_app_config ProlfquAppConfig
    prolfq_app_config = NULL,


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
    #' @field FDR_threshold fdr threshold
    FDR_threshold = 0.1,
    #' @field diff_threshold diff_threshold
    diff_threshold = 1,
    #' @field m1_linear linearModel
    m1_linear = "linearModel",
    #' @field m2_missing imputedModel
    m2_missing = "imputedModel",
    #' @field m3_merged mergedModel
    m3_merged = "mergedModel",
    #' @field m4_glm_protein m4_glm_protein
    m4_glm_protein = "glmModel",
    #' @field m4_glm_peptide m4_glm_peptide
    m4_glm_peptide = "glmModel",
    #' @field default_model default_model
    default_model = character(),
    #' @description
    #' initialize
    #' @param lfqdata lfqdata
    #' @param rowAnnot ProteinAnnotation
    #' @param prolfq_app_config ProlfquAppConfig
    #' @param contrasts vector with contrasts
    #' @param FDR_threshold FDR_threshold
    #' @param diff_threshold diff_threshold
    initialize = function(lfqdata,
                          rowAnnot,
                          prolfq_app_config,
                          contrasts,
                          default_model = "mergedModel") {
      stopifnot(default_model %in% c(self$m3_merged, self$m2_missing, self$m1_linear, self$m4_glm))
      self$default_model = default_model
      stopifnot("LFQData" %in% class(lfqdata))
      stopifnot("ProteinAnnotation" %in% class(rowAnnot))
      stopifnot("ProlfquAppConfig" %in% class(prolfq_app_config))
      stopifnot(length(contrasts) >= 1)
      self$lfqdata <- lfqdata
      self$rowAnnot <- rowAnnot
      self$prolfq_app_config <- prolfq_app_config
      self$contrasts <- contrasts
      self$FDR_threshold = prolfq_app_config$processing_options$FDR_threshold
      self$diff_threshold = prolfq_app_config$processing_options$diff_threshold
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
        contaminants = self$prolfq_app_config$processing_options$remove_cont,
        decoys = self$prolfq_app_config$processing_options$remove_decoys
      ))
      logger::log_info("removing contaminants and reverse sequences with patterns: ",
                       self$prolfq_app_config$processing_options$pattern_contaminants,
                       self$prolfq_app_config$processing_options$pattern_decoys)
    },
    #' @description
    #' transform data
    transform_data = function() {
      transformed <- prolfquapp::transform_lfqdata(
        self$lfqdata,
        method = self$prolfq_app_config$processing_options$transform,
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
      interaction <- self$prolfq_app_config$processing_options$interaction
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
    build_model_linear = function() {
      if (is.null(self$models[[self$m1_linear]]) ) {
        formula <- self$create_model_formula()
        formula_Condition <-  prolfqua::strategy_lm(formula)
        models <- prolfqua::build_model(
          self$lfqData_transformed,
          formula_Condition)
        self$models[[self$m1_linear]] <- models
      }
      return(self$models[[self$m1_linear]])
    },
    #' @description
    #' fit generalized linear model
    build_model_glm_protein = function(){
      if (is.null(self$models[[self$m4_glm_protein]])) {
        formula <- self$create_model_formula()
        modelFunction <-  prolfqua::strategy_glm(formula,
                                                 family = stats::binomial,
                                                 multiplier = 1.5,
                                                 offset = 1)
        models <- prolfqua::build_model(self$lfqData_transformed , modelFunction)
        self$models[[self$m4_glm_protein]] <- models
      }
      return(self$models[[self$m4_glm_protein]])
    },
    #' @description
    #' fit generalized linear model
    build_model_glm_peptide = function(){
      if (is.null(self$models[[self$m4_glm_peptide]])) {
        formula <- self$create_model_formula_peptide()
        modelFunction <-  strategy_glm(formula,
                                       family = stats::quasibinomial,
                                       multiplier = 1.2,
                                       offset = 1)
        models <- prolfqua::build_model(self$lfqData_peptide , modelFunction)
        self$models[[self$m4_glm_protein]] <- models
      }
      return(self$models[[self$m4_glm_protein]])
    },

    #' @description
    #' compute contrasts linear
    get_contrasts_linear = function() {
      private$get_contrasts(self$m1_linear)
    },
    #' @description
    #' get contrasts from glm model
    get_contrasts_glm = function(){
      private$get_contrasts(self$m4_glm_protein)
    },
    #' @description
    #' compute missing contrasts
    get_contrasts_missing = function(){
      if (is.null(self$contrastRes[[self$m2_missing]])) {
        mC <- prolfqua::ContrastsMissing$new(
          lfqdata = self$lfqData_transformed,
          contrasts = self$contrasts,
          modelName = self$m2_missing)
        conMI <- prolfqua::ContrastsModerated$new(mC)
        self$contrastRes[[self$m2_missing]] <- conMI
      }
      return(self$contrastRes[[self$m2_missing]])
    },

    #' @description
    #' merge contrasts
    get_contrasts_merged = function(){
      if (is.null(self$contrastRes[[self$m3_merged]])) {
        self$get_contrasts_linear()
        self$get_contrasts_missing()

        self$contrastRes[[self$m3_merged]] <- prolfqua::merge_contrasts_results(
          self$contrastRes[[self$m1_linear]],
          self$contrastRes[[self$m2_missing]])$merged

      }
      invisible(self$contrastRes[[self$m3_merged]])
    },
    #' @description
    #' filter contrasts for threshold
    filter_contrasts = function(){
      if (is.null(self$contrastRes[[self$default_model]])) {
        self$get_contrasts_merged()
      }
      datax <- self$contrastRes[[self$default_model]]$get_contrasts()
      datax <- datax |>
        dplyr::filter(.data$FDR < self$FDR_threshold &
                        abs(.data$diff) > self$diff_threshold )
      invisible(datax)
    },
    #' @description
    #' filter transformed lfq data for significant proteins.
    filter_data = function(){
      dx <- self$filter_contrasts()
      self$lfqData_subset <- self$lfqData_transformed$get_subset(dx)
      return(self$lfqData_subset)
    },
    #' @description
    #' create boxplots
    #'
    get_boxplots = function(){
      return(self$lfqData_subset$get_Plotter()$boxplots())
    },
    #' @description
    #' create boxplots
    #'
    contrasts_to_Grob = function(){
      datax <- self$filter_contrasts()
      xdn <- datax |> dplyr::nest_by(protein_Id)

      stats2grob <- function(data) {
        res <- data |>
          dplyr::select(contrast, diff, statistic, FDR) |>
          dplyr::mutate(
            diff = custom_round(diff),
            statistic = custom_round(statistic),
            FDR = custom_round(FDR)
          )
        res <- res |> gridExtra::tableGrob()
        res
      }
      grobs <- vector(mode = "list", length = length(nrow(xdn)))
      for (i in seq_len(nrow(xdn))) {
        grobs[[i]] <- stats2grob(xdn$data[[i]])
      }
      xdn$grobs <- grobs
      return(xdn)
    },
    #' @description
    #' get box with contrast information
    get_boxplots_contrasts = function(){
      ctrG <- self$contrasts_to_Grob()
      bp <- self$get_boxplots()
      stopifnot(nrow(ctrG) == nrow(bp))
      res <- vector(mode = "list", length = nrow(ctrG))
      for (i in seq_len(nrow(ctrG))) {
        res[[i]] <- gridExtra::arrangeGrob(
          bp$boxplot[[i]],
          ctrG$grobs[[i]], nrow = 2, heights = c(2/3, 1/3))
      }
      ctrG$bxpl_grobs = res
      return(ctrG)
    },
    #' @description
    #' write boxplots contrasts to file
    #' @param filename filename to write to
    write_boxplots_contrasts = function(filename){
      ctrG <- self$get_boxplots_contrasts()
      pdf(file = filename)
      for (i in seq_along(ctrG$bxpl_grobs)) {
        grid::grid.newpage()
        grid::grid.draw(ctrG$bxpl_grobs[[i]])
      }
      dev.off()
    }
  ),
  private = list(
    get_contrasts = function(contrastName){
      if (is.null(self$contrastRes[[contrastName]])) {
        self$build_model_linear()
        contr <- prolfqua::Contrasts$new(self$models[[contrastName]],
                                         self$contrasts,
                                         modelName = contrastName)
        conrM <- prolfqua::ContrastsModerated$new(contr)
        self$contrastRes[[contrastName]] <- conrM
      }
      invisible(self$contrastRes[[contrastName]])
    }

  )
)


