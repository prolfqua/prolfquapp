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


# DEAnalyse ----
#' will replace make_DEA_report
#' @export
#'
#' @examples
#' # example code
#'
#' pep <- prolfqua::sim_lfq_data_peptide_config(Nprot = 100)
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
#' #DEAnalyse$debug("get_contrasts_glm_peptide")
#' #DEAnalyse$debug("build_model_glm_protein")
#' deanalyse <- prolfquapp::DEAnalyse$new(pep, pA, GRP2, contrasts)
#' deanalyse$lfq_data_peptide$hierarchy_counts()
#' deanalyse$cont_decoy_summary()
#' deanalyse$prolfq_app_config$processing_options$remove_cont = TRUE
#' deanalyse$remove_cont_decoy()
#' deanalyse$aggregate()
#' pl <- deanalyse$get_aggregation_plots(exp_nr_children = 10)
#' print(pl$plots[[3]])
#' deanalyse$transform_data()
#' mod <- deanalyse$build_model_linear_protein()
#' contlm <- deanalyse$get_contrasts_linear_protein()
#'
#'
#' merged <- deanalyse$get_contrasts_merged_protein()
#' stopifnot(nrow(merged$get_contrasts()) == 200)
#' #deanalyse$create_model_formula()
#' #deanalyse$build_model_glm_protein()
#' #deanalyse$build_model_glm_peptide()
#' xprot <- deanalyse$get_contrasts_glm_protein()
#' stopifnot(nrow(merged$get_contrasts()) == 200)
#' xpep <- deanalyse$get_contrasts_glm_peptide()
#' xprot$get_contrasts()
#' xprot$get_Plotter()$volcano()
#' xpep$get_Plotter()$volcano()
#' sr <- deanalyse$lfq_data_peptide$get_Summariser()
#'
#'
#'
#' deanalyse$filter_contrasts()
#'
#' xd <- deanalyse$filter_data()
#' xd <- deanalyse$contrasts_to_Grob()
#' bb <- deanalyse$get_boxplots()
#' bx <- deanalyse$get_boxplots_contrasts()
#' grid::grid.draw(bx$bxpl_grobs[[1]])
#' # deanalyse$write_boxplots_contrasts("test.pdf")

DEAnalyse <- R6::R6Class(
  "DEAnalyse",
  public = list(

    #'@field prolfq_app_config ProlfquAppConfig
    prolfq_app_config = NULL,

    #' @field lfq_data_peptide LFQData peptide level
    lfq_data_peptide = NULL,
    #' @field lfq_data LFQData
    lfq_data = NULL,
    #' @field lfq_data_transformed todo
    lfq_data_transformed = NULL,
    #' @field lfq_data_subset todo
    lfq_data_subset = NULL,
    #' @field aggregator aggregator
    aggregator = NULL,

    #' @field rowAnnot ProteinAnnotation
    rowAnnot = NULL,
    #' @field contrasts vector with contrasts
    contrasts = character(),
    #' @field FDR_threshold fdr threshold
    FDR_threshold = 0.1,
    #' @field diff_threshold diff_threshold
    diff_threshold = 1,


    #' @field reference_proteins reference proteins to use for internal normalization
    reference_proteins = NULL, # to use for internal calibration

    #' @field formula todo
    formula = character(),
    #' @field formula glm peptide todo
    formula_glm_peptide = character(),

    #' @field models todo
    models = list(),
    #' @field contrast_results todo
    contrast_results = list(),

    #
    #' @field m1_linear linearModel
    m1_linear = "linearModel",
    #' @field m2_missing imputedModel
    m2_missing = "imputedModel",
    #' @field m3_merged mergedModel
    m3_merged = "mergedModel",
    #' @field m4_glm_protein m4_glm_protein
    m4_glm_protein = "glmModel",
    #' @field m4_glm_peptide m4_glm_peptide
    m4_glm_peptide = "glmModelPeptide",
    #' @field default_model default_model
    default_model = character(),
    #' @description
    #' initialize
    #' @param lfq_data lfq_data
    #' @param rowAnnot ProteinAnnotation
    #' @param prolfq_app_config ProlfquAppConfig
    #' @param contrasts vector with contrasts
    #' @param FDR_threshold FDR_threshold
    #' @param diff_threshold diff_threshold
    initialize = function(lfq_data_peptide,
                          rowAnnot,
                          prolfq_app_config,
                          contrasts,
                          default_model = "mergedModel") {
      stopifnot(default_model %in% c(self$m3_merged, self$m2_missing, self$m1_linear, self$m4_glm))
      self$default_model = default_model
      stopifnot("LFQData" %in% class(lfq_data_peptide))
      stopifnot("ProteinAnnotation" %in% class(rowAnnot))
      stopifnot("ProlfquAppConfig" %in% class(prolfq_app_config))
      stopifnot(length(contrasts) >= 1)
      self$lfq_data_peptide <- lfq_data_peptide
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
      self$lfq_data_peptide <- self$lfq_data_peptide$get_subset(self$rowAnnot$clean(
        contaminants = self$prolfq_app_config$processing_options$remove_cont,
        decoys = self$prolfq_app_config$processing_options$remove_decoys
      ))
      logger::log_info("removing contaminants and reverse sequences with patterns: ",
                       self$prolfq_app_config$processing_options$pattern_contaminants,
                       self$prolfq_app_config$processing_options$pattern_decoys)
    },
    #' @description
    #' aggregate peptide data
    aggregate = function(){
      agg_method = self$prolfq_app_config$processing_options$aggregate
      logger::log_info("AGGREGATING PEPTIDE DATA: {agg_method}.")

      lfqdata_peptide <- self$lfq_data_peptide

      if (length(lfqdata_peptide$config$table$hierarchy_keys()) == lfqdata_peptide$config$table$hierarchyDepth) {
        warning('nothing to aggregate from, returning unchanged data.')
        self$lfq_data <- lfqdata_peptide
        invisible(self$lfq_data)
      }

      if (agg_method == "topN") {
        self$aggregator <- lfqdata_peptide$get_Aggregator()
        self$aggregator$sum_topN(N = N)
        self$lfq_data <- self$aggregator$lfq_agg

      } else if (agg_method == "lmrob" || agg_method == "medpolish") {

        transformed_peptide <- lfqdata_peptide$get_Transformer()$intensity_array(log)$lfq
        self$aggregator <- transformed_peptide$get_Aggregator()

        if (agg_method == "lmrob" ) {
          self$aggregator$lmrob()
        } else if (agg_method == "medpolish") {
          self$aggregator$medpolish()
        }

        lfq_data <- self$aggregator$lfq_agg
        tr <- lfq_data$get_Transformer()
        tr <- tr$intensity_array(exp, force = TRUE)
        lfq_data <- tr$lfq
        lfq_data$is_transformed(FALSE)
        self$lfq_data <- lfq_data

      } else {
        logger::log_warn("no such aggregator {agg_method}.")
      }
      logger::log_info("END OF PROTEIN AGGREGATION")
      invisible(self$lfq_data)

    },
    #' @description
    #' get aggregation plots
    #'
    #' @param exp_nr_children nr children to filter; default >=2
    #'

    get_aggregation_plots = function(exp_nr_children = 2){
      subset <- self$rowAnnot$filter_by_nr_children(exp_nr_children = exp_nr_children)
      res <- self$aggregator$plot(subset)
      return(res)
    },
    #' @description
    #' write aggregation plots
    #' @param exp_nr_children nr children to filter; default >=2
    #'
    write_aggregation_plots = function(exp_nr_children = 2){
      subset <- self$rowAnnot$filter_by_nr_children(exp_nr_children = exp_nr_children)
      self$aggregator$write_plots(self$prolfq_app_config$zipdir, subset)
    },

    #' @description
    #' transform data
    transform_data = function() {
      transformed <- prolfquapp::transform_lfqdata(
        self$lfq_data,
        method = self$prolfq_app_config$processing_options$transform,
        internal = self$reference_proteins
      )
      self$lfq_data$rename_response("protein_abundance")
      transformed$rename_response("normalized_protein_abundance")
      self$lfq_data_transformed <- transformed
      invisible(transformed)
    },
    #' @description
    #' create model formula
    #'
    create_model_formula = function(){
      prlconfig <- self$lfq_data_transformed$config
      return(private$create_formula(prlconfig))
    },


    #' @description
    #' fit linear model
    build_model_linear_protein = function() {
      if (is.null(self$models[[self$m1_linear]]) ) {
        formula <- self$create_model_formula()
        formula_Condition <-  prolfqua::strategy_lm(formula)
        models <- prolfqua::build_model(
          self$lfq_data_transformed,
          formula_Condition)
        self$models[[self$m1_linear]] <- models
      }
      return(self$models[[self$m1_linear]])
    },
    #' @description
    #' fit generalized linear model
    build_model_glm_protein = function(){
      if (is.null(self$models[[self$m4_glm_protein]])) {
        lfq <- self$lfq_data_transformed

        lfq$complete_cases()
        lfq$data <- lfq$data |>
          dplyr::mutate(binresp =
                          factor(ifelse(is.na(!!sym(lfq$response())), 0, 1)))
        formula <- private$create_formula(lfq$config, response = "binresp")

        modelFunction <-  prolfqua::strategy_glm(formula,
                                                 family = stats::binomial,
                                                 multiplier = 1.2,
                                                 offset = 1)

        models <- prolfqua::build_model(lfq , modelFunction)
        self$models[[self$m4_glm_protein]] <- models
      }
      return(self$models[[self$m4_glm_protein]])
    },
    #' @description
    #' fit generalized linear model
    build_model_glm_peptide = function(){
      if (is.null(self$models[[self$m4_glm_peptide]])) {
        lfq <- self$lfq_data_peptide

        lfq$complete_cases()
        lfq$data <- lfq$data |>
          dplyr::mutate(binresp =
                          factor(ifelse(is.na(!!sym(lfq$response())), 0, 1)))

        formula <- private$create_formula(lfq$config, response = "binresp")
        # block for peptides
        hkey <- tail(lfq$config$table$hierarchy_keys(), n = 1)
        formula <- paste0(formula, "+", hkey)

        modelFunction <-  prolfqua::strategy_glm(
          formula,
          family = stats::quasibinomial,
          multiplier = 1.2,
          offset = 1)
        models <- prolfqua::build_model(lfq , modelFunction)
        self$models[[self$m4_glm_peptide]] <- models
      }
      return(self$models[[self$m4_glm_peptide]])
    },

    #' @description
    #' compute contrasts linear
    get_contrasts_linear_protein = function() {
      self$build_model_linear_protein()
      private$get_contrasts(self$m1_linear)
    },
    #' @description
    #' get contrasts from glm model
    get_contrasts_glm_peptide = function(){
      self$build_model_glm_peptide()
      private$get_contrasts(self$m4_glm_peptide)
    },
    #' @description
    #' get contrasts from glm model for peptides
    get_contrasts_glm_protein = function(){
      self$build_model_glm_protein()
      private$get_contrasts(self$m4_glm_protein)
    },
    #' @description
    #' compute missing contrasts
    get_contrasts_missing_protein = function(){
      if (is.null(self$contrast_results[[self$m2_missing]])) {
        mC <- prolfqua::ContrastsMissing$new(
          lfqdata = self$lfq_data_transformed,
          contrasts = self$contrasts,
          modelName = self$m2_missing)
        conMI <- prolfqua::ContrastsModerated$new(mC)
        self$contrast_results[[self$m2_missing]] <- conMI
      }
      return(self$contrast_results[[self$m2_missing]])
    },

    #' @description
    #' merge contrasts
    get_contrasts_merged_protein = function(){
      if (is.null(self$contrast_results[[self$m3_merged]])) {
        self$get_contrasts_linear_protein()
        self$get_contrasts_missing_protein()

        self$contrast_results[[self$m3_merged]] <- prolfqua::merge_contrasts_results(
          self$contrast_results[[self$m1_linear]],
          self$contrast_results[[self$m2_missing]])$merged

      }
      invisible(self$contrast_results[[self$m3_merged]])
    },
    #' @description
    #' filter contrasts for threshold
    filter_contrasts = function(){
      if (is.null(self$contrast_results[[self$default_model]])) {
        stop("no default model contrasts yet:", self$default_model)
      }
      datax <- self$contrast_results[[self$default_model]]$get_contrasts()
      datax <- datax |>
        dplyr::filter(.data$FDR < self$FDR_threshold &
                        abs(.data$diff) > self$diff_threshold )
      invisible(datax)
    },
    #' @description
    #' filter transformed lfq data for significant proteins.
    filter_data = function(){
      dx <- self$filter_contrasts()
      self$lfq_data_subset <- self$lfq_data_transformed$get_subset(dx)
      return(self$lfq_data_subset)
    },
    #' @description
    #' create boxplots
    #'
    get_boxplots = function(){
      if (is.null(self$lfq_data_subset)) {
        self$filter_data()
      }
      return(self$lfq_data_subset$get_Plotter()$boxplots())
    },
    #' @description
    #' create boxplots
    #'
    contrasts_to_Grob = function(){
      datax <- self$filter_contrasts()
      hkeys <- self$lfq_data_transformed$config$table$hierarchy_keys_depth()
      xdn <- datax |> dplyr::nest_by(!!!syms(hkeys))

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
      pb <- progress::progress_bar$new(total = nrow(xdn))
      for (i in seq_len(nrow(xdn))) {
        pb$tick()
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
      pb <- progress::progress_bar$new(total = nrow(ctrG))
      for (i in seq_len(nrow(ctrG))) {
        res[[i]] <- gridExtra::arrangeGrob(
          bp$boxplot[[i]],
          ctrG$grobs[[i]], nrow = 2, heights = c(2/3, 1/3))
        pb$tick()
      }
      ctrG$bxpl_grobs = res
      return(ctrG)
    },
    #' @description
    #' write boxplots contrasts to file
    #' @param filename filename to write to
    write_boxplots_contrasts = function(filename = "boxplots"){
      ctrG <- self$get_boxplots_contrasts()
      filename <- paste0(filename, "_FDR_", self$FDR_threshold, "_diff_", self$diff_threshold, ".pdf")
      logger::log_info("start writing boxplots into file : ", filename)
      pdf(file = file.path(self$prolfq_app_config$zipdir, filename))
      pb <- progress::progress_bar$new(total = length(ctrG$bxpl_grobs))
      for (i in seq_along(ctrG$bxpl_grobs)) {
        pb$tick()
        grid::grid.newpage()
        grid::grid.draw(ctrG$bxpl_grobs[[i]])
      }
      dev.off()
    }
  ),
  private = list(
    get_contrasts = function(contrastName){
      if (is.null(self$contrast_results[[contrastName]])) {
        if (!is.null(self$models[[contrastName]])) {
          contr <- prolfqua::Contrasts$new(self$models[[contrastName]],
                                           self$contrasts,
                                           modelName = contrastName)
          conrM <- prolfqua::ContrastsModerated$new(contr)
          self$contrast_results[[contrastName]] <- conrM
        } else {
          stop("No model for :", contrastName)
        }
      }
      invisible(self$contrast_results[[contrastName]])
    },
    create_formula = function(prlconfig, response = prlconfig$table$get_response()) {

      interaction <- self$prolfq_app_config$processing_options$interaction
      factors <- prlconfig$table$factor_keys_depth()[
        !grepl("^control", prlconfig$table$factor_keys_depth() , ignore.case = TRUE)
      ]
      # model with or without interactions
      if (interaction ) {
        formula <- paste0(response, " ~ ",
                          paste(factors, collapse = " + "))
      } else {
        formula <- paste0(response, " ~ ",
                          paste(factors, collapse = " * "))
      }
      logger::log_info("fitted model with formula : {formula}")
      self$formula <- formula
      return(formula)
    }


  )
)


