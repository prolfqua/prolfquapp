# DEAnalyse ----
#' Differential expression analysis engine using prolfqua facade classes
#'
#' Takes prepared LFQData (at the correct hierarchy level for the chosen facade)
#' and runs statistical modelling via prolfqua's ContrastsFacade classes
#' or the prolfquasaint SAINTexpress adapter.
#'
#' The caller (e.g. \code{ProteinDataPrep$build_deanalyse()}) is responsible for
#' providing data at the right level: aggregated protein-level for most facades,
#' or nested peptide-level for \code{lmer}/\code{ropeca}.
#'
#' @export
#'
#' @examples
#' pep <- prolfqua::sim_lfq_data_peptide_config(Nprot = 100)
#' pep <- prolfqua::LFQData$new(pep$data, pep$config)
#' pA <- data.frame(protein_Id = unique(pep$data_long()$protein_Id))
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
#' deanalyse <- data_prep$build_deanalyse(contrasts)
#' deanalyse$build_default()
#' stopifnot(nrow(deanalyse$contrast_results[[deanalyse$default_model]]$get_contrasts()) == 200)
#'
DEAnalyse <- R6::R6Class(
  "DEAnalyse",
  public = list(
    #' @field prolfq_app_config ProlfquAppConfig
    prolfq_app_config = NULL,

    #' @field lfq_data LFQData to model (transformed, at correct hierarchy level)
    lfq_data = NULL,
    #' @field lfq_data_raw raw (untransformed) LFQData for reporting
    lfq_data_raw = NULL,

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

    #' @field contrast_results named list of facade objects
    contrast_results = list(),
    #' @field saint_input SAINTexpress input tables for model = "saint"
    saint_input = NULL,
    #' @field saint_result SAINTexpress raw result list for model = "saint"
    saint_result = NULL,
    #' @field default_model facade registry key for the default model
    default_model = "lm_missing",

    #' @description
    #' Initialize DEAnalyse
    #' @param lfq_data LFQData to model (transformed, at correct hierarchy level)
    #' @param rowAnnot ProteinAnnotation object
    #' @param prolfq_app_config ProlfquAppConfig object
    #' @param contrasts named vector of contrast definitions
    #' @param default_model facade registry key (default "lm_missing")
    #' @param lfq_data_raw raw (untransformed) LFQData for reporting (optional)
    #' @param summary data.frame with contaminant/decoy summary (optional)
    initialize = function(
      lfq_data,
      rowAnnot,
      prolfq_app_config,
      contrasts,
      default_model = "lm_missing",
      lfq_data_raw = NULL,
      summary = NULL
    ) {
      default_model <- .resolve_facade_model(
        default_model,
        prolfq_app_config$processing_options$model_missing
      )
      stopifnot(default_model %in% .valid_facade_models())
      if (!.is_saint_model(default_model)) {
        stopifnot(length(contrasts) >= 1)
      }
      self$lfq_data <- lfq_data
      self$lfq_data_raw <- lfq_data_raw
      self$rowAnnot <- rowAnnot
      self$prolfq_app_config <- prolfq_app_config
      self$contrasts <- contrasts
      self$default_model <- default_model
      self$summary <- summary
      self$FDR_threshold <- prolfq_app_config$processing_options$FDR_threshold
      self$diff_threshold <- prolfq_app_config$processing_options$diff_threshold
    },

    #' @description
    #' Build a facade by registry key
    #' @param name facade registry key (e.g. "lm", "lm_missing", "limma")
    #' @param modelstr model formula string; auto-generated if NULL
    #' @return the facade object (invisibly)
    build_facade = function(name, modelstr = NULL) {
      if (!is.null(self$contrast_results[[name]])) {
        return(invisible(self$contrast_results[[name]]))
      }
      if (.is_saint_model(name)) {
        saint <- .build_saint_contrast_result(
          self$lfq_data,
          self$rowAnnot,
          spc = FALSE,
          engine = "r"
        )
        self$contrast_results[[name]] <- saint$contrast
        self$saint_input <- saint$input
        self$saint_result <- saint$result
        self$formula <- "SAINTexpress intensity model"
        return(invisible(saint$contrast))
      }
      entry <- prolfqua::FACADE_REGISTRY[[name]]
      if (is.null(entry)) {
        stop("Unknown facade: ", name)
      }
      if (is.null(modelstr)) {
        modelstr <- private$create_modelstr()
      }

      facade_class <- utils::getFromNamespace(entry$class, "prolfqua")
      facade <- facade_class$new(self$lfq_data, modelstr, self$contrasts)

      self$contrast_results[[name]] <- facade
      self$formula <- paste(self$lfq_data$response(), modelstr)
      invisible(facade)
    },

    #' @description
    #' Build the default facade (as set in default_model)
    build_default = function() {
      self$build_facade(self$default_model)
    },

    #' @description
    #' Join default model contrasts with protein row annotations
    get_annotated_contrasts = function() {
      if (is.null(self$contrast_results[[self$default_model]])) {
        stop("no default model contrasts yet: ", self$default_model)
      }
      contrast_obj <- self$contrast_results[[self$default_model]]
      datax <- contrast_obj$get_contrasts()
      datax_signif <- if (.is_saint_model(self$default_model)) {
        contrast_obj$get_ora(
          up = TRUE,
          FDR_threshold = self$FDR_threshold,
          diff_threshold = self$diff_threshold
        )
      } else {
        private$filter_significant_contrasts(datax)
      }
      if (.is_saint_model(self$default_model)) {
        datax_signif <- dplyr::inner_join(
          self$rowAnnot$row_annot,
          datax_signif,
          multiple = "all"
        )
      }
      datax <- dplyr::inner_join(
        self$rowAnnot$row_annot,
        datax,
        multiple = "all"
      )
      self$annotated_contrasts <- datax
      if (!.is_saint_model(self$default_model)) {
        datax_signif <- dplyr::inner_join(
          self$rowAnnot$row_annot,
          datax_signif,
          multiple = "all"
        )
      }
      self$annotated_contrasts_signif <- datax_signif
      invisible(self$annotated_contrasts)
    },

    #' @description
    #' Return contrast rows passing FDR and difference thresholds
    filter_contrasts = function() {
      if (is.null(self$contrast_results[[self$default_model]])) {
        stop("no default model contrasts yet:", self$default_model)
      }
      contrast_obj <- self$contrast_results[[self$default_model]]
      if (.is_saint_model(self$default_model)) {
        datax <- contrast_obj$get_ora(
          up = TRUE,
          FDR_threshold = self$FDR_threshold,
          diff_threshold = self$diff_threshold
        )
      } else {
        datax <- contrast_obj$get_contrasts()
        datax <- private$filter_significant_contrasts(datax)
      }
      invisible(datax)
    }
  ),
  private = list(
    filter_significant_contrasts = function(datax) {
      if (.is_saint_model(self$default_model)) {
        return(
          datax |>
            dplyr::filter(
              .data$BFDR < self$FDR_threshold &
                abs(.data$log2_EFCs) > self$diff_threshold
            )
        )
      }
      datax |>
        dplyr::filter(
          .data$FDR < self$FDR_threshold &
            abs(.data$diff) > self$diff_threshold
        )
    },

    create_modelstr = function() {
      interaction <- self$prolfq_app_config$processing_options$interaction
      factors <- self$lfq_data$relevant_factor_keys()[
        !grepl(
          "^control",
          self$lfq_data$relevant_factor_keys(),
          ignore.case = TRUE
        )
      ]
      sep <- if (is.null(interaction) || !interaction) " + " else " * "
      modelstr <- paste0("~ ", paste(factors, collapse = sep))
      logger::log_info("model formula: {self$lfq_data$response()} {modelstr}")
      modelstr
    }
  )
)
