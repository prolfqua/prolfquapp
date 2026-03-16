# DEAnalyse ----
#' Differential expression analysis engine
#'
#' Takes a prepared \code{ProteinDataPrep} object and runs statistical modelling:
#' fits linear and GLM models, computes moderated contrasts, merges results.
#'
#' @export
#'
#' @examples
#' pep <- prolfqua::sim_lfq_data_peptide_config(Nprot = 100)
#' pep <- prolfqua::LFQData$new(pep$data, pep$config)
#' pA <- data.frame(protein_Id = unique(pep$data$protein_Id))
#' pA <- pA |> dplyr::mutate(fasta.annot = paste0(pA$protein_Id, "_description"))
#' pA <- prolfquapp::ProteinAnnotation$new(pep, row_annot = pA, description = "fasta.annot")
#' GRP2 <- prolfquapp::make_DEA_config_R6()
#' GRP2$processing_options$diff_threshold <- 0.2
#' GRP2$processing_options$transform <- "robscale"
#'
#' contrasts <- c("AVsC" = "group_A - group_Ctrl", BVsC = "group_B - group_Ctrl")
#'
#' data_prep <- prolfquapp::ProteinDataPrep$new(pep, pA, GRP2)
#' data_prep$cont_decoy_summary()
#' data_prep$remove_cont_decoy()
#' data_prep$aggregate()
#' data_prep$transform_data()
#'
#' deanalyse <- prolfquapp::DEAnalyse$new(data_prep, contrasts)
#' mod <- deanalyse$build_model_linear_protein()
#' contlm <- deanalyse$get_contrasts_linear_protein()
#' merged <- deanalyse$get_contrasts_merged_protein()
#' stopifnot(nrow(merged$get_contrasts()) == 200)
#'
DEAnalyse <- R6::R6Class(
  "DEAnalyse",
  public = list(
    #' @field prolfq_app_config ProlfquAppConfig
    prolfq_app_config = NULL,

    #' @field lfq_data_peptide LFQData peptide level
    lfq_data_peptide = NULL,
    #' @field lfq_data LFQData protein level
    lfq_data = NULL,
    #' @field lfq_data_transformed normalized LFQData
    lfq_data_transformed = NULL,

    #' @field rowAnnot ProteinAnnotation
    rowAnnot = NULL,
    #' @field contrasts vector with contrasts
    contrasts = character(),
    #' @field FDR_threshold FDR threshold
    FDR_threshold = 0.1,
    #' @field diff_threshold difference threshold
    diff_threshold = 1,

    #' @field summary data.frame with contaminant/decoy summary
    summary = NULL,
    #' @field annotated_contrasts contrasts joined with row annotations
    annotated_contrasts = NULL,
    #' @field annotated_contrasts_signif significant annotated contrasts
    annotated_contrasts_signif = NULL,

    #' @field formula model formula
    formula = character(),

    #' @field models list of fitted models
    models = list(),
    #' @field contrast_results list of contrast results
    contrast_results = list(),

    #' @field m1_linear Linear_Model
    m1_linear = "Linear_Model",
    #' @field m2_missing Imputed_Mean
    m2_missing = "Imputed_Mean",
    #' @field m3_merged mergedModel
    m3_merged = "mergedModel",
    #' @field m4_glm_protein glmModel
    m4_glm_protein = "glmModel",
    #' @field m4_glm_peptide glmModelPeptide
    m4_glm_peptide = "glmModelPeptide",
    #' @field default_model default model key
    default_model = character(),

    #' @description
    #' Initialize DEAnalyse from a prepared ProteinDataPrep object
    #' @param data_prep ProteinDataPrep object with aggregated and normalized data
    #' @param contrasts named vector of contrast definitions
    #' @param default_model which model to use for final results
    initialize = function(data_prep, contrasts, default_model = "mergedModel") {
      stopifnot("ProteinDataPrep" %in% class(data_prep))
      stopifnot(
        default_model %in%
          c(
            self$m3_merged,
            self$m2_missing,
            self$m1_linear,
            self$m4_glm_protein
          )
      )
      stopifnot(length(contrasts) >= 1)
      self$default_model <- default_model
      self$lfq_data_peptide <- data_prep$lfq_data_peptide
      self$lfq_data <- data_prep$lfq_data
      self$lfq_data_transformed <- data_prep$lfq_data_transformed
      self$rowAnnot <- data_prep$rowAnnot
      self$prolfq_app_config <- data_prep$prolfq_app_config
      self$summary <- data_prep$summary
      self$contrasts <- contrasts
      self$FDR_threshold <- data_prep$prolfq_app_config$processing_options$FDR_threshold
      self$diff_threshold <- data_prep$prolfq_app_config$processing_options$diff_threshold
    },

    #' @description
    #' Create model formula from transformed data config
    create_model_formula = function() {
      prlconfig <- self$lfq_data_transformed$config
      return(private$create_formula(prlconfig))
    },

    #' @description
    #' Fit linear model at protein level
    build_model_linear_protein = function() {
      if (is.null(self$models[[self$m1_linear]])) {
        formula <- self$create_model_formula()
        formula_Condition <- prolfqua::strategy_lm(formula)
        models <- prolfqua::build_model(
          self$lfq_data_transformed,
          formula_Condition
        )
        self$models[[self$m1_linear]] <- models
      }
      return(self$models[[self$m1_linear]])
    },

    #' @description
    #' Get GLM strategy for protein-level missingness model
    get_strategy_glm_prot = function() {
      lfq <- self$lfq_data_transformed
      formula <- private$create_formula(lfq$config, response = "binresp")
      modelFunction <- prolfqua::strategy_glm(
        formula,
        family = stats::binomial,
        multiplier = 1.2,
        offset = 1
      )
      modelFunction <- prolfqua::strategy_logistf(formula)
      return(modelFunction)
    },

    #' @description
    #' Fit generalized linear model at protein level
    build_model_glm_protein = function() {
      if (is.null(self$models[[self$m4_glm_protein]])) {
        lfq <- self$lfq_data_transformed
        lfq$complete_cases()
        lfq$data <- lfq$data |>
          dplyr::mutate(
            binresp = ifelse(is.na(!!sym(lfq$response())), 0, 1)
          )
        modelFunction <- self$get_strategy_glm_prot()
        models <- prolfqua::build_model(lfq, modelFunction)
        self$models[[self$m4_glm_protein]] <- models
      }
      return(self$models[[self$m4_glm_protein]])
    },

    #' @description
    #' Fit generalized linear model at peptide level (not yet implemented)
    build_model_glm_peptide = function() {
      stop(
        "Not yet implemented: build_model_glm_peptide references undefined variables (istar, models2)"
      )
    },

    #' @description
    #' Compute moderated contrasts from linear model
    get_contrasts_linear_protein = function() {
      self$build_model_linear_protein()
      private$get_contrasts(self$m1_linear)
    },

    #' @description
    #' Compute moderated contrasts from GLM peptide model
    get_contrasts_glm_peptide = function() {
      self$build_model_glm_peptide()
      private$get_contrasts(self$m4_glm_peptide)
    },

    #' @description
    #' Compute moderated contrasts from GLM protein model
    get_contrasts_glm_protein = function() {
      self$build_model_glm_protein()
      private$get_contrasts(self$m4_glm_protein)
    },

    #' @description
    #' Compute moderated contrasts from missing-value imputation model
    get_contrasts_missing_protein = function() {
      if (is.null(self$contrast_results[[self$m2_missing]])) {
        mC <- prolfqua::ContrastsMissing$new(
          lfqdata = self$lfq_data_transformed,
          contrasts = self$contrasts,
          modelName = self$m2_missing
        )
        conMI <- prolfqua::ContrastsModerated$new(mC)
        self$contrast_results[[self$m2_missing]] <- conMI
      }
      return(self$contrast_results[[self$m2_missing]])
    },

    #' @description
    #' Merge linear and missing-value contrasts (or use linear only if model_missing = FALSE)
    get_contrasts_merged_protein = function() {
      if (is.null(self$contrast_results[[self$m3_merged]])) {
        model_missing <- self$prolfq_app_config$processing_options$model_missing
        if (is.null(model_missing) || model_missing) {
          self$get_contrasts_linear_protein()
          self$get_contrasts_missing_protein()
          self$contrast_results[[
            self$m3_merged
          ]] <- prolfqua::merge_contrasts_results(
            self$contrast_results[[self$m1_linear]],
            self$contrast_results[[self$m2_missing]]
          )$merged
        } else {
          self$get_contrasts_linear_protein()
          self$contrast_results[[self$m3_merged]] <- self$contrast_results[[
            self$m1_linear
          ]]
        }
      }
      invisible(self$contrast_results[[self$m3_merged]])
    },

    #' @description
    #' Join default model contrasts with protein row annotations
    get_annotated_contrasts = function() {
      if (is.null(self$contrast_results[[self$default_model]])) {
        stop("no default model contrasts yet: ", self$default_model)
      }
      datax <- self$contrast_results[[self$default_model]]$get_contrasts()
      datax <- dplyr::inner_join(
        self$rowAnnot$row_annot,
        datax,
        multiple = "all"
      )
      self$annotated_contrasts <- datax

      datax_signif <- datax |>
        dplyr::filter(
          .data$FDR < self$FDR_threshold &
            abs(.data$diff) > self$diff_threshold
        )
      self$annotated_contrasts_signif <- datax_signif
      invisible(self$annotated_contrasts)
    },

    #' @description
    #' Return contrast rows passing FDR and difference thresholds
    filter_contrasts = function() {
      if (is.null(self$contrast_results[[self$default_model]])) {
        stop("no default model contrasts yet:", self$default_model)
      }
      datax <- self$contrast_results[[self$default_model]]$get_contrasts()
      datax <- datax |>
        dplyr::filter(
          .data$FDR < self$FDR_threshold &
            abs(.data$diff) > self$diff_threshold
        )
      invisible(datax)
    }
  ),
  private = list(
    get_contrasts = function(contrastName) {
      if (is.null(self$contrast_results[[contrastName]])) {
        if (!is.null(self$models[[contrastName]])) {
          contr <- prolfqua::Contrasts$new(
            self$models[[contrastName]],
            self$contrasts,
            modelName = contrastName
          )
          conrM <- prolfqua::ContrastsModerated$new(contr)
          self$contrast_results[[contrastName]] <- conrM
        } else {
          stop("No model for :", contrastName)
        }
      }
      invisible(self$contrast_results[[contrastName]])
    },
    create_formula = function(prlconfig, response = prlconfig$get_response()) {
      interaction <- self$prolfq_app_config$processing_options$interaction
      factors <- prlconfig$factor_keys_depth()[
        !grepl("^control", prlconfig$factor_keys_depth(), ignore.case = TRUE)
      ]
      if (is.null(interaction) || !interaction) {
        formula <- paste0(response, " ~ ", paste(factors, collapse = " + "))
      } else {
        formula <- paste0(response, " ~ ", paste(factors, collapse = " * "))
      }
      logger::log_info("fitted model with formula : {formula}")
      self$formula <- formula
      return(formula)
    }
  )
)
