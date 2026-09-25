# Legacy YAML alias: "prolfqua" maps to the lm facade, or to lm_impute (LOD
# refit of failed fits, flagged "lod_imputed" in estimate_type) when
# processing_options$model_missing is set.
.resolve_facade_model <- function(model, model_missing = FALSE) {
  if (identical(model, "prolfqua")) {
    if (isTRUE(model_missing)) "lm_impute" else "lm"
  } else {
    model
  }
}

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
    default_model = "lm_impute",

    #' @description
    #' Initialize DEAnalyse
    #' @param lfq_data LFQData to model (transformed, at correct hierarchy level)
    #' @param rowAnnot ProteinAnnotation object
    #' @param prolfq_app_config ProlfquAppConfig object
    #' @param contrasts named vector of contrast definitions
    #' @param default_model facade registry key (default "lm_impute")
    #' @param lfq_data_raw raw (untransformed) LFQData for reporting (optional)
    #' @param summary data.frame with contaminant/decoy summary (optional)
    initialize = function(
      lfq_data,
      rowAnnot,
      prolfq_app_config,
      contrasts,
      default_model = "lm_impute",
      lfq_data_raw = NULL,
      summary = NULL
    ) {
      po <- prolfq_app_config$processing_options
      default_model <- .resolve_facade_model(default_model, po$model_missing)
      entry <- prolfqua::lookup_facade(default_model)
      if (is.null(entry)) {
        stop("Unknown facade: ", default_model)
      }
      if (!isTRUE(entry$needs_saint_annotation)) {
        stopifnot(length(contrasts) >= 1)
      }
      self$lfq_data <- lfq_data
      self$lfq_data_raw <- lfq_data_raw
      # Stamp the app's decoy / contaminant patterns onto the modelling config
      # here, at the modelling boundary, so build_facade's targets-only gate
      # fires even when aggregation/transform handed back a fresh config.
      self$lfq_data$set_config_value("pattern_decoys", po$pattern_decoys)
      self$lfq_data$set_config_value("pattern_contaminants", po$pattern_contaminants)
      self$rowAnnot <- rowAnnot
      self$prolfq_app_config <- prolfq_app_config
      self$contrasts <- contrasts
      self$default_model <- default_model
      self$summary <- summary
      self$FDR_threshold <- po$FDR_threshold
      self$diff_threshold <- po$diff_threshold
    },

    #' @description
    #' Build a facade by registry key via \code{prolfqua::lookup_facade()}, so
    #' facades registered by downstream packages (e.g. \code{"saint"}) work
    #' like the built-in ones. Backends with \code{needs_saint_annotation = TRUE}
    #' receive \code{row_annot} from \code{self$rowAnnot}.
    #' @param name facade registry key (e.g. "lm", "lm_missing", "limma", "saint")
    #' @param modelstr model formula string; auto-generated if NULL, ignored by SAINT-style backends
    #' @return the facade object (invisibly)
    build_facade = function(name, modelstr = NULL) {
      if (!is.null(self$contrast_results[[name]])) {
        return(invisible(self$contrast_results[[name]]))
      }
      entry <- prolfqua::lookup_facade(name)
      if (is.null(entry)) {
        stop("Unknown facade: ", name)
      }
      facade_class <- utils::getFromNamespace(entry$class, entry$package %||% "prolfqua")

      # Targets-only fit: drop decoys before the fit so they never enter the
      # shared variance pool (limma prior, DEqMS trend). A non-NULL
      # pattern_decoys (even "") opts in, as in prolfqua::build_contrast_analysis().
      # self$lfq_data keeps decoys so the abundance export still carries them.
      model_lfq <- self$lfq_data
      if (!is.null(model_lfq$get_config()$pattern_decoys)) {
        top <- model_lfq$hierarchy_keys()[1]
        n_before <- length(unique(model_lfq$data_long()[[top]]))
        model_lfq <- model_lfq$remove_decoys()
        n_dropped <- n_before - length(unique(model_lfq$data_long()[[top]]))
        if (n_dropped > 0) {
          logger::log_info(
            "targets-only fit: dropped {n_dropped} decoy {top} before modelling ",
            "(kept in raw data for export)."
          )
        }
      }

      if (isTRUE(entry$needs_saint_annotation)) {
        facade <- facade_class$new(model_lfq, modelstr = NULL, contrasts = NULL, row_annot = self$rowAnnot$row_annot)
        self$formula <- "SAINTexpress intensity model"
      } else {
        modelstr <- modelstr %||% private$create_modelstr()
        private$validate_group_coverage()
        facade <- facade_class$new(model_lfq, modelstr, self$contrasts)
        self$formula <- paste(model_lfq$response(), modelstr)
      }

      self$contrast_results[[name]] <- facade
      invisible(facade)
    },

    #' @description
    #' Build the default facade (as set in default_model)
    build_default = function() {
      self$build_facade(self$default_model)
    },

    #' @description
    #' Join default-model contrasts, and those passing \code{filter_contrasts()},
    #' with protein row annotations.
    get_annotated_contrasts = function() {
      datax_signif <- self$filter_contrasts()
      datax <- self$contrast_results[[self$default_model]]$get_contrasts()
      self$annotated_contrasts <- .join_annotation(self$rowAnnot$row_annot, datax, self$rowAnnot$pID)
      self$annotated_contrasts_signif <- .join_annotation(self$rowAnnot$row_annot, datax_signif, self$rowAnnot$pID)
      invisible(self$annotated_contrasts)
    },

    #' @description
    #' Return contrast rows passing FDR and difference thresholds, as given by
    #' the contrast object's \code{filter_significant()} (one-sided for
    #' directional backends such as SAINT).
    filter_contrasts = function() {
      contrast_obj <- self$contrast_results[[self$default_model]]
      if (is.null(contrast_obj)) {
        stop("no default model contrasts yet: ", self$default_model)
      }
      invisible(contrast_obj$filter_significant(
        FDR_threshold = self$FDR_threshold,
        diff_threshold = self$diff_threshold
      ))
    }
  ),
  private = list(
    # Model factors: the relevant factor keys except control columns.
    model_factors = function() {
      grep("^control", self$lfq_data$relevant_factor_keys(), ignore.case = TRUE, invert = TRUE, value = TRUE)
    },

    create_modelstr = function() {
      interaction <- self$prolfq_app_config$processing_options$interaction
      sep <- if (is.null(interaction) || !interaction) " + " else " * "
      modelstr <- paste0("~ ", paste(private$model_factors(), collapse = sep))
      logger::log_info("model formula: {self$lfq_data$response()} {modelstr}")
      modelstr
    },

    # Abort early when a contrast references a group without samples after
    # matching the annotation to the quantification data; otherwise prolfqua
    # fails later with a cryptic "subscript out of bounds" in linfct
    # construction. Only the populated LFQData, not the annotation, shows this.
    validate_group_coverage = function() {
      factor_cols <- private$model_factors()
      if (length(factor_cols) == 0 || length(self$contrasts) == 0) {
        return(invisible(NULL))
      }
      data <- self$lfq_data$data_long(na.omit = TRUE)
      # Level tokens with data, named like contrast terms: paste0(factor_key, level).
      level_counts <- unlist(lapply(intersect(factor_cols, colnames(data)), function(col) {
        present <- unique(stats::na.omit(data[[col]]))
        counts <- vapply(present, function(level) sum(data[[col]] == level, na.rm = TRUE), integer(1))
        stats::setNames(counts, paste0(col, present))
      }))
      for (i in seq_along(self$contrasts)) {
        contrast <- self$contrasts[[i]]
        referenced <- tryCatch(all.vars(str2lang(contrast)), error = function(e) character())
        for (token in unique(unlist(strsplit(referenced, ":", fixed = TRUE)))) {
          if (any(startsWith(token, factor_cols)) && !token %in% names(level_counts)) {
            stop(
              "Contrast '",
              names(self$contrasts)[i] %||% paste0("contrast_", i),
              "' (",
              contrast,
              ") cannot be computed: group '",
              token,
              "' has 0 samples after matching the annotation to the ",
              "quantification data. Groups with data: ",
              paste(sprintf("%s [n=%d]", names(level_counts), level_counts), collapse = ", "),
              ". Check that the quantification report contains the raw files ",
              "for every group in the annotation (a group whose runs are ",
              "missing from the report is dropped during annotation).",
              call. = FALSE
            )
          }
        }
      }
      invisible(NULL)
    }
  )
)

# DEAnalysePeptideToProtein ----
#' Differential expression analysis from peptide input to protein output
#'
#' Runs facades that consume peptide-level measurements but emit protein-level
#' contrasts. The input \code{LFQData} keeps peptide hierarchy columns, while its
#' active hierarchy depth points to the protein level.
#'
#' @export
DEAnalysePeptideToProtein <- R6::R6Class(
  "DEAnalysePeptideToProtein",
  inherit = DEAnalyse,
  public = list(
    #' @description
    #' Build the default peptide-to-protein facade.
    build_default = function() {
      # Nested facades take a fixed-effects-only modelstr; lmer adds random effects itself.
      self$build_facade(self$default_model, modelstr = private$create_modelstr())
    }
  )
)
