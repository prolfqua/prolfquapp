# Load necessary libraries

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
#' pA <- prolfqua::ProteinAnnotation$new(pep,row_annot = pA ,description = "fasta.annot")
#' GRP2 <- prolfquapp::make_DEA_config_R6()
#' GRP2$processing_options$transform <- "robscale"
#' pep$factors()
#' GRP2$pop <- list(Contrasts = c("AVsC" = "group_A - group_Ctrl", BVsC = "group_B - group_Ctrl"))
#' deanalyse <- DEAnalyse$new(pep, pA, GRP2)
#' deanalyse$cont_decoy_summary()
#' deanalyse$GRP2$processing_options$remove_cont = TRUE
#' deanalyse$remove_cont_decoy()
#' deanalyse$transform_data()
#' deanalyse$create_model_formula()
#' deanalyse$build_model()
#' deanalyse$compute_contrasts()
#' deanalyse$filter_contrasts()
#' deanalyse$RES
DEAnalyse <- R6::R6Class(
  "DEAnalyse",
  public = list(
    #' @field lfqdata `LFQData`
    lfqdata = NULL,
    #' @field rowAnnot `ProteinAnnotation`
    rowAnnot = NULL,
    #'@field GRP2 `ProlfquAppConfig`
    GRP2 = NULL,
    RES = list(),
    #' initialize
    #' @param lfqdata lfqdata
    #' @param rowAnnot ProteinAnnotation
    #' @param GRP2 ProlfquAppConfig
    initialize = function(lfqdata, rowAnnot, GRP2) {
      self$lfqdata <- lfqdata
      self$rowAnnot <- rowAnnot
      self$GRP2 <- GRP2
      self$RES <- list()
    },
    cont_decoy_summary = function() {
      allProt <- nrow(self$rowAnnot$row_annot)
      contdecoySummary <- data.frame(
        totalNrOfProteins = allProt,
        percentOfContaminants = round(self$rowAnnot$annotate_contaminants(
          self$GRP2$processing_options$pattern_contaminants) / allProt * 100, digits = 2),
        percentOfFalsePositives = round(self$rowAnnot$annotate_decoys(
          self$GRP2$processing_options$pattern_decoys) / allProt * 100, digits = 2),
        NrOfProteinsNoDecoys = self$rowAnnot$nr_clean()
      )
      self$RES$Summary <- contdecoySummary
      return(contdecoySummary)
    },
    remove_cont_decoy = function() {
      if (self$GRP2$processing_options$remove_cont || self$GRP2$processing_options$remove_decoys) {
        self$lfqdata <- self$lfqdata$get_subset(self$rowAnnot$clean(
          contaminants = self$GRP2$processing_options$remove_cont,
          decoys = self$GRP2$processing_options$remove_decoys
        ))
        logger::log_info("removing contaminants and reverse sequences with patterns: ",
                         self$GRP2$processing_options$pattern_contaminants,
                         self$GRP2$processing_options$pattern_decoys)
      }
    },
    transform_data = function() {
      transformed <- prolfquapp::transform_lfqdata(
        self$lfqdata,
        method = self$GRP2$processing_options$transform,
        internal = self$GRP2$pop$internal
      )
      self$lfqdata$rename_response("protein_abundance")
      transformed$rename_response("normalized_protein_abundance")
      self$RES$lfqData <- self$lfqdata
      self$RES$transformedlfqData <- transformed
      invisible(transformed)
    },
    create_model_formula = function() {
      transformed <- self$RES$transformedlfqData
      factors <- transformed$config$table$factor_keys_depth()[
        !grepl("^control", transformed$config$table$factor_keys_depth() , ignore.case = TRUE)
      ]
      # model with or without interactions
      if (is.null(self$GRP2$processing_options$interaction) || !self$GRP2$processing_options$interaction ) {
        formula <- paste0(transformed$config$table$get_response(), " ~ ",
                          paste(factors, collapse = " + "))
      } else {
        formula <- paste0(transformed$config$table$get_response(), " ~ ",
                          paste(factors, collapse = " * "))
      }
      message("FORMULA :",  formula)
      self$RES$formula <- formula
      return(formula)
    },
    build_model = function() {
      transformed <- self$RES$transformedlfqData
      formula <- self$RES$formula
      formula_Condition <-  prolfqua::strategy_lm(formula)
      self$RES$models <- prolfqua::build_model(
        transformed,
        formula_Condition,
        subject_Id = transformed$config$table$hierarchy_keys() )
      logger::log_info("fitted model with formula : {formula}")
      invisible(self$RES$models)
    },
    compute_contrasts = function() {
      contr <- prolfqua::Contrasts$new(self$RES$models,
                                       self$GRP2$pop$Contrasts,
                                       modelName = "Linear_Model")
      conrM <- prolfqua::ContrastsModerated$new(
        contr)

      if (is.null(self$GRP2$processing_options$missing) || self$GRP2$processing_options$missing ) {
        transformed <- self$RES$transformedlfqData
        mC <- prolfqua::ContrastsMissing$new(
          lfqdata = transformed,
          contrasts = self$GRP2$pop$Contrasts,
          modelName = "Imputed_Mean")
        conMI <- prolfqua::ContrastsModerated$new(
          mC)
        res <- prolfqua::merge_contrasts_results(conrM, conMI)
        self$RES$contrMerged <- res$merged
        self$RES$contrMore <- res$more

      } else {
        self$RES$contrMerged <- conrM
        self$RES$contrMore <- NULL
      }
      invisible(self$RES$contrMerged)
    },
    filter_contrasts = function(){
      datax <- self$RES$contrMerged$get_contrasts()
      datax <- datax |>
        dplyr::filter(.data$FDR < self$GRP2$processing_options$FDR_threshold &
                        abs(.data$diff) > self$GRP2$processing_options$diff_threshold )
      GRP2$RES$contrastsData_signif <- datax
      invisible(datax)
    })
)

