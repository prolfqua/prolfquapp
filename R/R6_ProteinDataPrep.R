#' ProteinDataPrep
#'
#' Handles data preparation for differential expression analysis:
#' contaminant/decoy filtering, peptide-to-protein aggregation, and normalization.
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
#' GRP2$processing_options$transform <- "robscale"
#'
#' data_prep <- prolfquapp::ProteinDataPrep$new(pep, pA, GRP2)
#' data_prep$cont_decoy_summary()
#' data_prep$remove_cont_decoy()
#' data_prep$aggregate()
#' data_prep$transform_data()
#'
ProteinDataPrep <- R6::R6Class(
  "ProteinDataPrep",
  public = list(
    #' @field prolfq_app_config ProlfquAppConfig
    prolfq_app_config = NULL,
    #' @field lfq_data_peptide LFQData peptide level
    lfq_data_peptide = NULL,
    #' @field lfq_data LFQData protein level (after aggregation)
    lfq_data = NULL,
    #' @field lfq_data_transformed normalized LFQData
    lfq_data_transformed = NULL,
    #' @field aggregator aggregator object
    aggregator = NULL,
    #' @field rowAnnot ProteinAnnotation
    rowAnnot = NULL,
    #' @field summary data.frame with contaminant/decoy summary
    summary = NULL,

    #' @description
    #' Initialize ProteinDataPrep
    #' @param lfq_data_peptide LFQData object at peptide level
    #' @param rowAnnot ProteinAnnotation object
    #' @param prolfq_app_config ProlfquAppConfig object
    initialize = function(lfq_data_peptide, rowAnnot, prolfq_app_config) {
      stopifnot("LFQData" %in% class(lfq_data_peptide))
      stopifnot("ProteinAnnotation" %in% class(rowAnnot))
      stopifnot("ProlfquAppConfig" %in% class(prolfq_app_config))
      self$lfq_data_peptide <- lfq_data_peptide
      self$rowAnnot <- rowAnnot
      self$prolfq_app_config <- prolfq_app_config
    },

    #' @description
    #' Count number of contaminants and decoys
    cont_decoy_summary = function() {
      allProt <- nrow(self$rowAnnot$row_annot)
      self$summary <- data.frame(
        totalNrOfProteins = allProt,
        percentOfContaminants = round(self$rowAnnot$annotate_contaminants() / allProt * 100, digits = 2),
        percentOfFalsePositives = round(self$rowAnnot$annotate_decoys() / allProt * 100, digits = 2),
        NrOfProteinsNoDecoys = self$rowAnnot$nr_clean()
      )
      self$summary
    },

    #' @description
    #' Remove contaminants and decoys from peptide data
    remove_cont_decoy = function() {
      self$lfq_data_peptide <- self$lfq_data_peptide$get_subset(self$rowAnnot$clean(
        contaminants = self$prolfq_app_config$processing_options$remove_cont,
        decoys = self$prolfq_app_config$processing_options$remove_decoys
      ))
      logger::log_info(
        "removing contaminants and reverse sequences with patterns: ",
        self$prolfq_app_config$processing_options$pattern_contaminants,
        self$prolfq_app_config$processing_options$pattern_decoys
      )
    },

    #' @description
    #' Aggregate peptide data to protein level
    aggregate = function() {
      agg_method <- self$prolfq_app_config$processing_options$aggregate
      logger::log_info("AGGREGATING PEPTIDE DATA: {agg_method}.")
      lfqdata_peptide <- self$lfq_data_peptide

      if (length(lfqdata_peptide$config$hierarchy_keys()) == lfqdata_peptide$config$hierarchyDepth) {
        warning("nothing to aggregate from, returning unchanged data.")
        self$lfq_data <- lfqdata_peptide
        return(invisible(self$lfq_data))
      }

      if (agg_method == "topN") {
        self$aggregator <- lfqdata_peptide$get_Aggregator()
        self$aggregator$sum_topN(N = N)
        self$lfq_data <- self$aggregator$lfq_agg
      } else if (agg_method == "lmrob" || agg_method == "medpolish") {
        transformed_peptide <- lfqdata_peptide$get_Transformer()$intensity_array(log)$lfq
        self$aggregator <- transformed_peptide$get_Aggregator()

        if (agg_method == "lmrob") {
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
    #' Get aggregation plots
    #' @param exp_nr_children minimum number of peptides per protein; default 2
    get_aggregation_plots = function(exp_nr_children = 2) {
      subset <- self$rowAnnot$filter_by_nr_children(exp_nr_children = exp_nr_children)
      res <- self$aggregator$plot(subset)
      return(res)
    },

    #' @description
    #' Write aggregation plots to file
    #' @param exp_nr_children minimum number of peptides per protein; default 2
    write_aggregation_plots = function(exp_nr_children = 2) {
      subset <- self$rowAnnot$filter_by_nr_children(exp_nr_children = exp_nr_children)
      self$aggregator$write_plots(self$prolfq_app_config$zipdir, subset)
    },

    #' @description
    #' Transform and normalize protein-level data
    transform_data = function() {
      transformed <- prolfquapp::transform_lfqdata(
        self$lfq_data,
        method = self$prolfq_app_config$processing_options$transform
      )
      self$lfq_data$rename_response("abundance")

      if (length(self$prolfq_app_config$processing_options$internal) > 0) {
        x <- transformed$hierarchy()
        mm <- colnames(x)[1]
        x <- x |> dplyr::filter(!!dplyr::sym(mm) %in% self$prolfq_app_config$processing_options$internal)
        if (nrow(x) == 0) {
          what <- paste(self$prolfq_app_config$processing_options$internal, collapse = ",")
          warning("not in list : ", what)
        } else {
          xs <- transformed$get_subset(x)
          tr <- transformed$get_Transformer()
          transformed <- tr$center_to_reference(xs)$lfq
        }
      }

      transformed$rename_response("normalized_abundance")
      self$lfq_data_transformed <- transformed
      invisible(transformed)
    }
  )
)
