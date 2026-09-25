# QC_generator------
#' QC_generator
#' @export
QC_generator <- R6::R6Class(
  "QC_generator",
  public = list(
    #' @field lfqdata lfqdata
    lfqdata = NULL,
    #' @field lfqdata_prot lfqdata_prot
    lfqdata_prot = NULL,
    #' @field lfqdata_prot_IBAQ lfqdata_prot_IBAQ
    lfqdata_prot_IBAQ = NULL,
    #' @field lfqdata_prot_transformed lfqdata_prot_transformed (VSN normalized)
    lfqdata_prot_transformed = NULL,
    #' @field protein_annotation protein_annotation
    protein_annotation = NULL,
    #' @field lfqdata_peptide lfqdata_peptide
    lfqdata_peptide = NULL,
    #' @field lfqdata_peptide_transformed lfqdata_peptide_transformed (VSN normalized)
    lfqdata_peptide_transformed = NULL,
    #' @field output_dir output_dir
    output_dir = NULL,
    #' @field GRP2 GRP2
    GRP2 = NULL,
    #' @field TABLES2WRITE TABLES2WRITE
    TABLES2WRITE = list(),

    #' @field links links
    links = list(),
    #' @description
    #' initialize
    #' @param lfqdata LFQData object
    #' @param protein_annotation ProteinAnnotation object
    #' @param prolfquapp_config ProlfquAppConfig object
    initialize = function(lfqdata, protein_annotation, prolfquapp_config) {
      self$GRP2 <- prolfquapp_config
      self$lfqdata <- lfqdata
      self$protein_annotation <- protein_annotation
      self$output_dir <- self$GRP2$get_zipdir()
      self$TABLES2WRITE <- list()
    },
    #' @description
    #' get peptides in wide format
    #' @return peptide data in wide format
    get_peptides_wide = function() {
      lfqdata <- self$lfqdata$get_copy()
      lfqdata$set_config_value("hierarchy_depth", min(2, length(self$lfqdata$hierarchy_keys())))
      self$lfqdata_peptide <- prolfquapp::aggregate_data(lfqdata, agg_method = "medpolish")
      invisible(private$annotate(self$lfqdata$data_wide()$data))
    },
    #' @description
    #' get VSN-transformed peptide data
    #' @return VSN-transformed peptide LFQData
    get_peptides_transformed = function() {
      if (is.null(self$lfqdata_peptide_transformed)) {
        self$lfqdata_peptide_transformed <- prolfquapp::transform_lfqdata(self$lfqdata, method = "vsn")
      }
      invisible(self$lfqdata_peptide_transformed)
    },
    #' @description
    #' get VSN-transformed peptide data in wide format
    #' @return VSN-transformed peptide data in wide format
    get_peptides_transformed_wide = function() {
      private$annotate(self$get_peptides_transformed()$data_wide()$data)
    },
    #' @description
    #' get annotation data
    #' @return annotation data.frame
    get_annotation = function() {
      invisible(self$lfqdata$factors())
    },
    #' @description
    #' get protein data
    #' @return protein LFQData
    get_prot_data = function() {
      if (is.null(self$lfqdata_prot)) {
        self$lfqdata_prot <- prolfquapp::aggregate_data(self$lfqdata, agg_method = "medpolish")
      }
      invisible(self$lfqdata_prot)
    },
    #' @description
    #' get protein data in wide format
    #' @return protein data in wide format
    get_prot_wide = function() {
      private$annotate_with_nr_children(self$get_prot_data())
    },
    #' @description
    #' get VSN-transformed protein data
    #' @return VSN-transformed protein LFQData
    get_prot_transformed = function() {
      if (is.null(self$lfqdata_prot_transformed)) {
        self$lfqdata_prot_transformed <- prolfquapp::transform_lfqdata(self$get_prot_data(), method = "vsn")
      }
      invisible(self$lfqdata_prot_transformed)
    },
    #' @description
    #' get VSN-transformed protein data in wide format
    #' @return VSN-transformed protein data in wide format
    get_prot_transformed_wide = function() {
      private$annotate(self$get_prot_transformed()$data_wide()$data)
    },
    #' @description
    #' get IBAQ protein data
    #' @return IBAQ protein LFQData
    get_prot_IBAQ = function() {
      relevant_columns <- c("protein_length", "nr_tryptic_peptides")
      if (!all(relevant_columns %in% colnames(self$protein_annotation$row_annot))) {
        warning("skipping IBAQ computation, no:", paste(relevant_columns, collapse = "; "))
      } else if (is.null(self$lfqdata_prot_IBAQ)) {
        self$lfqdata_prot_IBAQ <- prolfquapp::compute_IBAQ_values(self$lfqdata, self$protein_annotation)
      }
      invisible(self$lfqdata_prot_IBAQ)
    },
    #' @description
    #' get protein abundance per group
    #' @return protein abundance per group
    get_protein_per_group_abundance = function() {
      invisible(self$get_prot_IBAQ()$get_Summariser()$percentage_abundance())
    },
    #' @description
    #' get protein abundance per group with row annotation
    #' @return protein abundance per group with annotation
    get_protein_per_group_abundance_with_row_annot = function() {
      private$annotate_per_group(self$get_protein_per_group_abundance())
    },
    #' @description
    #' get protein abundance per group in wide format
    #' @return protein abundance per group in wide format
    get_protein_per_group_abundance_wide = function() {
      precabund_data_wide <- self$get_protein_per_group_abundance() |>
        tidyr::pivot_wider(
          id_cols = self$lfqdata$hierarchy_keys()[1],
          names_from = interaction,
          values_from = c(
            nrReplicates,
            nrMeasured,
            nrNAs,
            sd,
            var,
            meanAbundance,
            medianAbundance,
            CV,
            id,
            abundance_percent,
            abundance_percent_cumulative,
            percent_prot
          )
        )
      private$annotate_per_group(precabund_data_wide)
    },
    #' @description
    #' get IBAQ protein data in wide format
    #' @return IBAQ protein data in wide format
    get_prot_IBAQ_wide = function() {
      ibaq <- self$get_prot_IBAQ()
      if (is.null(ibaq)) NULL else private$annotate_with_nr_children(ibaq)
    },
    #' @description
    #' get list of all tables
    #' @return list of tables
    get_list = function() {
      TABLES2WRITE <- list()
      TABLES2WRITE$peptide_wide <- self$get_peptides_wide()
      TABLES2WRITE$peptide_VSN_normalized <- self$get_peptides_transformed_wide()
      TABLES2WRITE$annotation <- self$get_annotation()
      TABLES2WRITE$prot_medpolish_estimate <- self$get_prot_wide()
      TABLES2WRITE$prot_VSN_normalized <- self$get_prot_transformed_wide()
      TABLES2WRITE$prot_IBAQ_estimate <- self$get_prot_IBAQ_wide()
      TABLES2WRITE$prot_IBAQ_per_group_stats <- self$get_protein_per_group_abundance_wide()
      return(TABLES2WRITE)
    },
    #' @description
    #' write tables to xlsx file
    write_xlsx = function() {
      xlsxfile <- file.path(self$output_dir, paste0("proteinAbundances_", self$GRP2$project_spec$workunit_Id, ".xlsx"))
      writexl::write_xlsx(self$get_list(), path = xlsxfile)
      self$links[["QC_XLSX"]] <- xlsxfile
    },
    #' @description
    #' copy dataset/annotation file to output directory
    #' @param dataset_path path to the dataset file
    copy_dataset = function(dataset_path) {
      if (file.exists(dataset_path)) {
        dest_path <- file.path(self$output_dir, basename(dataset_path))
        file.copy(dataset_path, dest_path, overwrite = TRUE)
        self$links[["DATASET"]] <- dest_path
        logger::log_info("Copied dataset to: ", dest_path)
      } else {
        logger::log_warn("Dataset file not found: ", dataset_path)
      }
    },
    #' @description
    #' render QC protein abundances report
    render_QC_protein_abundances = function() {
      tryCatch(
        {
          pap_file <- file.path(self$output_dir, "proteinAbundances.rds")
          saveRDS(self, file = pap_file)
          render_quarto_protein_abundances_report(
            pap_file = pap_file,
            output_dir = self$output_dir,
            output_file = "QC_ProteinAbundances_tabset.html",
            project_info = private$project_info(),
            factors = TRUE
          )
        },
        error = function(e) logger::log_warn("Skipping QC protein-abundances Quarto report: ", conditionMessage(e))
      )
      self$links[["QC_ABUNDANCES"]] <- file.path(self$output_dir, "QC_ProteinAbundances_tabset.html")
    },
    #' @description
    #' render sample size QC report
    render_sample_size_QC = function() {
      if (nrow(self$get_prot_data()$factors()) > 1) {
        tryCatch(
          {
            qc_data_file <- file.path(self$output_dir, "QC_sampleSizeEstimation.rds")
            prot <- self$get_prot_data()
            saveRDS(list(data = prot$data_long(), configuration = prot$get_config()), file = qc_data_file)
            software <- self$GRP2$software
            target_type <- if (grepl("MZMINE", software)) {
              "metabolite"
            } else if (grepl("PEPTIDE", software)) {
              "peptide"
            } else {
              "protein"
            }
            render_quarto_qc_sse_report(
              qc_data_file = qc_data_file,
              output_dir = self$output_dir,
              output_file = "QCandSSE_tabset.html",
              project_conf = private$project_info(),
              target_type = target_type
            )
          },
          error = function(e) logger::log_warn("Skipping QC sample-size Quarto report: ", conditionMessage(e))
        )
      } else {
        message("only a single sample: ", nrow(self$get_prot_data()$factors()))
      }
      self$links[["QC_SAMPLE_SIZE"]] <- file.path(self$output_dir, "QCandSSE_tabset.html")
    },
    #' @description
    #' render index HTML file
    render_index_html = function() {
      title <- paste0("QC Results for WU : ", self$GRP2$project_spec$workunit_Id, " and input : ", self$GRP2$software)
      self$links <- self$links[order(names(self$links) == "QC_XLSX")]
      items <- sprintf("<li><a href='%s'>%s</a></li>", vapply(self$links, basename, ""), names(self$links))
      str <- c("<!DOCTYPE html>", "<html>", "<head>", paste0("<title>", title, "</title>"), "</head>", "<body>")
      str <- c(str, paste0("<h1>", title, "</h1>"), "<ul>", items, "</ul>", "</body>", "</html>")
      cat(str, file = file.path(self$output_dir, "index.html"), sep = "\n")
    },
    #' @description
    #' render index markdown file
    render_index_md = function() {
      self$links <- self$links[order(names(self$links) == "QC_XLSX")]
      str <- c(
        paste0(
          "# QC Results for WU : ",
          self$GRP2$project_spec$workunit_Id,
          ", and input : ",
          self$GRP2$software,
          "\n"
        ),
        "\n## Available Reports\n",
        sprintf("- [%s](%s)", names(self$links), vapply(self$links, basename, ""))
      )
      cat(str, file = file.path(self$output_dir, "index.md"), sep = "\n")
    },
    #' @description
    #' get protein per group small wide format
    #' @return protein per group data in small wide format
    get_protein_per_group_small_wide = function() {
      tableconfig <- self$get_prot_IBAQ()$get_config()
      protID <- tableconfig$hierarchy_keys_depth()
      factor_keys <- tableconfig$factor_keys_depth()
      value_cols <- c("nrMeasured", "meanAbundance", "abundance_percent")
      dplyr::inner_join(self$protein_annotation$row_annot, self$get_protein_per_group_abundance(), by = protID) |>
        dplyr::mutate(abundance_percent = signif(abundance_percent, 2)) |>
        dplyr::select(all_of(c(protID, "nrPeptides", factor_keys, value_cols, "description"))) |>
        tidyr::pivot_wider(names_from = all_of(factor_keys), values_from = all_of(value_cols))
    }
  ),
  private = list(
    annotate = function(data) {
      dplyr::left_join(self$protein_annotation$row_annot, data, multiple = "all")
    },
    annotate_per_group = function(data) {
      by <- self$get_prot_IBAQ()$relevant_hierarchy_keys()
      invisible(dplyr::inner_join(self$protein_annotation$row_annot, data, multiple = "all", by = by))
    },
    annotate_with_nr_children = function(lfqdata) {
      dplyr::inner_join(
        private$annotate(lfqdata$data_wide()$data),
        lfqdata$data_wide(value = lfqdata$nr_children_col())$data,
        by = c(lfqdata$relevant_hierarchy_keys(), "isotopeLabel")
      )
    },
    # Plain list of scalar project identifiers (empty fields become NULL) for the Quarto report headers.
    project_info = function() {
      as_id <- function(x) if (length(x) >= 1 && nzchar(as.character(x)[[1]])) as.character(x)[[1]]
      ps <- self$GRP2$project_spec
      list(
        project_Id = as_id(ps$project_Id),
        project_name = as_id(ps$project_name),
        order_Id = as_id(ps$order_Id),
        workunit_Id = as_id(ps$workunit_Id),
        input_URL = as_id(ps$input_URL),
        software = as_id(self$GRP2$software)
      )
    }
  )
)
