# Legacy YAML alias: the user-facing config keyword "prolfqua" maps to
# the lm facade, or lm_impute when processing_options$model_missing is
# set. The lm_impute facade refits failed/singular per-protein fits at
# LOD with borrowed variance; modelName is the facade key ("lm_impute")
# and rescued rows are flagged in the estimate_type column
# ("lod_imputed") so the rescue is visible in the contrast output. The
# older alias target lm_missing relied on ContrastsMissing (group-mean
# substitution) which is now deprecated.
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
      default_model <- .resolve_facade_model(
        default_model,
        prolfq_app_config$processing_options$model_missing
      )
      entry <- prolfqua::lookup_facade(default_model)
      if (is.null(entry)) {
        stop("Unknown facade: ", default_model)
      }
      if (!isTRUE(entry$needs_saint_annotation)) {
        stopifnot(length(contrasts) >= 1)
      }
      self$lfq_data <- lfq_data
      self$lfq_data_raw <- lfq_data_raw
      # Enable the quant-layer targets-only fit: stamp the app's decoy /
      # contaminant patterns onto the modelling LFQData config so build_facade's
      # gate fires. prolfquapp bypasses prolfqua::build_contrast_analysis, and
      # aggregation/transform may hand back a fresh config, so we stamp here at
      # the modelling boundary rather than relying on upstream propagation.
      po <- prolfq_app_config$processing_options
      self$lfq_data$set_config_value("pattern_decoys", po$pattern_decoys)
      self$lfq_data$set_config_value("pattern_contaminants", po$pattern_contaminants)
      self$rowAnnot <- rowAnnot
      self$prolfq_app_config <- prolfq_app_config
      self$contrasts <- contrasts
      self$default_model <- default_model
      self$summary <- summary
      self$FDR_threshold <- prolfq_app_config$processing_options$FDR_threshold
      self$diff_threshold <- prolfq_app_config$processing_options$diff_threshold
    },

    #' @description
    #' Build a facade by registry key. Dispatches through
    #' \code{prolfqua::lookup_facade()} so any facade registered by a
    #' downstream package (e.g. \code{prolfquasaint::ContrastsSAINTFacade}
    #' registered as \code{"saint"}) is reachable the same way as the
    #' built-in prolfqua facades. SAINT-style backends that need the
    #' protein annotation (registry attribute
    #' \code{needs_saint_annotation = TRUE}) receive \code{row_annot}
    #' from \code{self$rowAnnot}.
    #' @param name facade registry key (e.g. "lm", "lm_missing", "limma",
    #'   "saint")
    #' @param modelstr model formula string; auto-generated if NULL.
    #'   Ignored by facades whose backend derives contrasts from
    #'   annotation (e.g. SAINT).
    #' @return the facade object (invisibly)
    build_facade = function(name, modelstr = NULL) {
      if (!is.null(self$contrast_results[[name]])) {
        return(invisible(self$contrast_results[[name]]))
      }
      entry <- prolfqua::lookup_facade(name)
      if (is.null(entry)) {
        stop("Unknown facade: ", name)
      }
      facade_class <- utils::getFromNamespace(
        entry$class,
        entry$package %||% "prolfqua"
      )

      # Targets-only fit (F1): drop decoys immediately before constructing the
      # facade, so they never enter the fit or the shared variance pool (limma
      # prior / DEqMS variance-count trend) where they would perturb *target*
      # q-values. Gated on pattern_decoys being non-NULL -- the SAME opt-in
      # semantics as prolfqua::build_contrast_analysis(): a non-NULL pattern
      # (including "", which the shared detector reads as "defaults only") means
      # opt in; NULL means the decoy machinery is off. The app maps an empty
      # REVpattern to NULL upstream, so "cleared pattern" == off. self$lfq_data
      # keeps decoys so the abundance export still carries them (NA stats).
      model_lfq <- self$lfq_data
      rev_pat <- model_lfq$get_config()$pattern_decoys
      if (!is.null(rev_pat)) {
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
        facade <- facade_class$new(
          model_lfq,
          modelstr = NULL,
          contrasts = NULL,
          row_annot = self$rowAnnot$row_annot
        )
        self$formula <- "SAINTexpress intensity model"
      } else {
        if (is.null(modelstr)) {
          modelstr <- private$create_modelstr()
        }
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
    #' Join default-model contrasts with protein row annotations.
    #' Significance filtering is delegated to
    #' \code{contrast_obj$filter_significant()}; backends with
    #' \code{ContrastConfiguration$significance_directional = TRUE}
    #' (e.g. SAINT) get one-sided filtering automatically.
    get_annotated_contrasts = function() {
      if (is.null(self$contrast_results[[self$default_model]])) {
        stop("no default model contrasts yet: ", self$default_model)
      }
      contrast_obj <- self$contrast_results[[self$default_model]]
      datax <- contrast_obj$get_contrasts()
      datax_signif <- contrast_obj$filter_significant(
        FDR_threshold = self$FDR_threshold,
        diff_threshold = self$diff_threshold
      )
      hkeys <- self$lfq_data$hierarchy_keys()
      datax <- .join_annotation(self$rowAnnot$row_annot, datax, hkeys)
      datax_signif <- .join_annotation(self$rowAnnot$row_annot, datax_signif, hkeys)
      self$annotated_contrasts <- datax
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
      invisible(
        contrast_obj$filter_significant(
          FDR_threshold = self$FDR_threshold,
          diff_threshold = self$diff_threshold
        )
      )
    }
  ),
  private = list(
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
    },

    # Abort early with an informative error when a contrast references a group
    # that has no samples after matching the annotation to the quantification
    # data. Without this guard the run fails much later inside prolfqua's
    # linfct construction with a cryptic "subscript out of bounds", because the
    # missing factor level has no model coefficient. The annotation file alone
    # cannot reveal this: a group only disappears once its raw files are found
    # to be absent from the quant report, so the check reads the populated
    # LFQData rather than the annotation.
    validate_group_coverage = function() {
      factor_cols <- self$lfq_data$relevant_factor_keys()
      factor_cols <- factor_cols[
        !grepl("^control", factor_cols, ignore.case = TRUE)
      ]
      if (length(factor_cols) == 0 || length(self$contrasts) == 0) {
        return(invisible(NULL))
      }
      data <- self$lfq_data$data_long(na.omit = TRUE)
      # Reconstruct the level tokens that have data, mirroring how contrasts are
      # built: paste0(factor_key, level), e.g. "G_OPENTRON".
      present_tokens <- character(0)
      level_counts <- integer(0)
      for (f in factor_cols) {
        if (!f %in% colnames(data)) {
          next
        }
        levels_present <- unique(stats::na.omit(data[[f]]))
        tokens <- paste0(f, levels_present)
        counts <- vapply(
          levels_present,
          function(l) sum(data[[f]] == l, na.rm = TRUE),
          integer(1)
        )
        present_tokens <- c(present_tokens, tokens)
        names(counts) <- tokens
        level_counts <- c(level_counts, counts)
      }
      for (i in seq_along(self$contrasts)) {
        contrast_str <- self$contrasts[[i]]
        contrast_name <- names(self$contrasts)[i] %||% paste0("contrast_", i)
        referenced <- tryCatch(
          all.vars(str2lang(contrast_str)),
          error = function(e) character(0)
        )
        # Interaction terms appear as a single backticked symbol "A:B".
        referenced <- unique(unlist(strsplit(referenced, ":", fixed = TRUE)))
        for (token in referenced) {
          # Only classify tokens that look like <factor_key><level>.
          if (!any(startsWith(token, factor_cols))) {
            next
          }
          if (!(token %in% present_tokens)) {
            present_summary <- paste(
              sprintf("%s [n=%d]", names(level_counts), level_counts),
              collapse = ", "
            )
            stop(
              "Contrast '", contrast_name, "' (", contrast_str,
              ") cannot be computed: group '", token,
              "' has 0 samples after matching the annotation to the ",
              "quantification data. Groups with data: ",
              if (length(present_summary)) present_summary else "<none>",
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
      # All nested facades (lmer_nested, ropeca_nested, firth_nested,
      # limpa_nested) accept a fixed-effects-only modelstr. The lmer facade
      # augments it internally with random effects derived from the LFQData.
      self$build_facade(
        self$default_model,
        modelstr = private$create_modelstr()
      )
    }
  )
)
