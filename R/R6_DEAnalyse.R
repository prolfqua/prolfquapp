# DEAnalyse ----
#' Differential expression analysis engine using prolfqua facade classes
#'
#' Takes prepared LFQData (at the correct hierarchy level for the chosen facade)
#' and runs statistical modelling via prolfqua's ContrastsFacade classes.
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
    initialize = function(lfq_data, rowAnnot, prolfq_app_config,
                          contrasts, default_model = "lm_missing",
                          lfq_data_raw = NULL, summary = NULL) {
      stopifnot(default_model %in% names(prolfqua::FACADE_REGISTRY))
      stopifnot(length(contrasts) >= 1)
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
      entry <- prolfqua::FACADE_REGISTRY[[name]]
      if (is.null(entry)) stop("Unknown facade: ", name)
      if (is.null(modelstr)) modelstr <- private$create_modelstr()

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
    create_modelstr = function() {
      interaction <- self$prolfq_app_config$processing_options$interaction
      factors <- self$lfq_data$relevant_factor_keys()[
        !grepl("^control", self$lfq_data$relevant_factor_keys(), ignore.case = TRUE)
      ]
      sep <- if (is.null(interaction) || !interaction) " + " else " * "
      modelstr <- paste0("~ ", paste(factors, collapse = sep))
      logger::log_info("model formula: {self$lfq_data$response()} {modelstr}")
      modelstr
    }
  )
)
