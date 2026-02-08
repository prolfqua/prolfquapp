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
      lfqdata$config$table$hierarchyDepth <- min(2, length(self$lfqdata$config$table$hierarchyKeys()))
      self$lfqdata_peptide <- prolfquapp::aggregate_data(lfqdata, agg_method = "medpolish")
      peptide_wide <- dplyr::left_join(self$protein_annotation$row_annot,
        self$lfqdata$to_wide()$data,
        multiple = "all"
      )
      invisible(peptide_wide)
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
      lfqdata_pep_tr <- self$get_peptides_transformed()
      if (is.null(lfqdata_pep_tr)) {
        return(NULL)
      }
      peptide_wide <- dplyr::left_join(
        self$protein_annotation$row_annot,
        lfqdata_pep_tr$to_wide()$data,
        multiple = "all"
      )
      return(peptide_wide)
    },
    #' @description
    #' get annotation data
    #' @return annotation data.frame
    get_annotation = function() {
      annotation <- self$lfqdata$factors()
      invisible(annotation)
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
      lfqdata_prot <- self$get_prot_data()
      proteins_wide <- dplyr::left_join(
        self$protein_annotation$row_annot,
        lfqdata_prot$to_wide()$data,
        multiple = "all"
      )

      # Get nr_children data using helper method
      nr_children_data <- private$get_nr_children_data(lfqdata_prot)

      proteins_wide <- dplyr::inner_join(
        proteins_wide,
        nr_children_data,
        by = c(lfqdata_prot$config$table$hierarchy_keys_depth(), "isotopeLabel")
      )
      return(proteins_wide)
    },
    #' @description
    #' get VSN-transformed protein data
    #' @return VSN-transformed protein LFQData
    get_prot_transformed = function() {
      if (is.null(self$lfqdata_prot_transformed)) {
        lfqdata_prot <- self$get_prot_data()
        self$lfqdata_prot_transformed <- prolfquapp::transform_lfqdata(lfqdata_prot, method = "vsn")
      }
      invisible(self$lfqdata_prot_transformed)
    },
    #' @description
    #' get VSN-transformed protein data in wide format
    #' @return VSN-transformed protein data in wide format
    get_prot_transformed_wide = function() {
      lfqdata_prot_tr <- self$get_prot_transformed()
      if (is.null(lfqdata_prot_tr)) {
        return(NULL)
      }
      proteins_wide <- dplyr::left_join(
        self$protein_annotation$row_annot,
        lfqdata_prot_tr$to_wide()$data,
        multiple = "all"
      )
      return(proteins_wide)
    },
    #' @description
    #' get IBAQ protein data
    #' @return IBAQ protein LFQData
    get_prot_IBAQ = function() {
      relevant_columns <- c("protein_length", "nr_tryptic_peptides")
      if (is.null(self$lfqdata_prot_IBAQ) && all(relevant_columns %in% colnames(self$protein_annotation$row_annot))) {
        self$lfqdata_prot_IBAQ <- prolfquapp::compute_IBAQ_values(self$lfqdata, self$protein_annotation)
      } else if (!all(relevant_columns %in% colnames(self$protein_annotation$row_annot))) {
        warning("skipping IBAQ computation, no:", paste(relevant_columns, collapse = "; "))
      }
      invisible(self$lfqdata_prot_IBAQ)
    },
    #' @description
    #' get protein abundance per group
    #' @return protein abundance per group
    get_protein_per_group_abundance = function() {
      summarizer <- self$get_prot_IBAQ()$get_Summariser()
      precabund <- summarizer$percentage_abundance()
      invisible(precabund)
    },
    #' @description
    #' get protein abundance per group with row annotation
    #' @return protein abundance per group with annotation
    get_protein_per_group_abundance_with_row_annot = function() {
      summarizer <- self$get_prot_IBAQ()$get_Summariser()
      precabund <- summarizer$percentage_abundance()
      precabund <- dplyr::inner_join(
        self$protein_annotation$row_annot,
        precabund,
        multiple = "all",
        by = self$get_prot_IBAQ()$config$table$hierarchy_keys_depth()
      )
      invisible(precabund)
    },
    #' @description
    #' get protein abundance per group in wide format
    #' @return protein abundance per group in wide format
    get_protein_per_group_abundance_wide = function() {
      precabund <- self$get_protein_per_group_abundance()
      precabund_data_wide <- precabund |>
        tidyr::pivot_wider(
          id_cols = self$lfqdata$config$table$hierarchy_keys()[1],
          names_from = interaction,
          values_from = c(nrReplicates, nrMeasured, nrNAs, sd, var, meanAbundance, medianAbundance, CV, id, abundance_percent, abundance_percent_cumulative, percent_prot)
        )
      precabund <- dplyr::inner_join(
        self$protein_annotation$row_annot,
        precabund_data_wide,
        multiple = "all",
        by = self$get_prot_IBAQ()$config$table$hierarchy_keys_depth()
      )

      invisible(precabund)
    },
    #' @description
    #' get IBAQ protein data in wide format
    #' @return IBAQ protein data in wide format
    get_prot_IBAQ_wide = function() {
      if (!is.null(self$get_prot_IBAQ())) {
        IBAQ_abundances <-
          dplyr::left_join(self$protein_annotation$row_annot,
            self$get_prot_IBAQ()$to_wide()$data,
            multiple = "all"
          )

        # Get nr_children data using helper method
        nr_children_data <- private$get_nr_children_data(self$get_prot_IBAQ())

        IBAQ_abundances <- dplyr::inner_join(
          IBAQ_abundances,
          nr_children_data,
          by = c(self$get_prot_IBAQ()$config$table$hierarchy_keys_depth(), "isotopeLabel")
        )
        return(IBAQ_abundances)
      } else {
        return(NULL)
      }
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
      xlsxfile <- file.path(
        self$output_dir,
        paste0("proteinAbundances_", self$GRP2$project_spec$workunit_Id, ".xlsx")
      )
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
      file.copy(system.file("application/GenericQC/QC_ProteinAbundances.Rmd", package = "prolfquapp"),
        to = self$output_dir, overwrite = TRUE
      )
      if (TRUE) {
        rmarkdown::render(file.path(self$output_dir, "QC_ProteinAbundances.Rmd"),
          params = list(
            pap = self,
            project_info = self$GRP2$project_spec,
            factors = TRUE
          ),
          output_file = "proteinAbundances.html"
        )
      } else {
        str <- c(
          "<!DOCTYPE html>",
          "<html>",
          "<head>", "<title>",
          "There is a problem",
          "</title>", "</head>",
          "<body>",
          "<h1>",
          paste0("the input file :", files$data, " is empty"),
          "</h1>",
          "</body>",
          "</html>"
        )
        cat(str, file = file.path(self$output_dir, "proteinAbundances.html"), sep = "\n")
      }
      self$links[["QC_ABUNDANCES"]] <- file.path(self$output_dir, "proteinAbundances.html")
    },
    #' @description
    #' render sample size QC report
    render_sample_size_QC = function() {
      if (nrow(self$get_prot_data()$factors()) > 1) {
        file.copy(system.file("application/GenericQC/QCandSSE.Rmd", package = "prolfquapp"),
          to = self$output_dir, overwrite = TRUE
        )
        rmarkdown::render(file.path(self$output_dir, "QCandSSE.Rmd"),
          params = list(
            data = self$get_prot_data()$data,
            configuration = self$get_prot_data()$config,
            project_conf = GRP2$project_spec,
            target_type = private$get_target_type()
          ),
          output_file = "QC_sampleSizeEstimation.html"
        )
      } else {
        message("only a single sample: ", nrow(self$get_prot_data()$factors()))
      }
      self$links[["QC_SAMPLE_SIZE"]] <- file.path(self$output_dir, "QC_sampleSizeEstimation.html")
    },
    #' @description
    #' render index HTML file
    render_index_html = function() {
      str <- c(
        "<!DOCTYPE html>",
        "<html>",
        "<head>",
        paste0(
          "<title>QC Results for WU : ",
          self$GRP2$project_spec$workunit_Id, " and input : ",
          self$GRP2$software, "</title>"
        ),
        "</head>",
        "<body>",
        paste0("<h1>QC Results for WU : ", self$GRP2$project_spec$workunit_Id, " and input : ", self$GRP2$software, "</h1>"),
        "<ul>"
      )
      # Sort links to ensure QC_XLSX is last
      sorted_links <- names(self$links)
      if ("QC_XLSX" %in% sorted_links) {
        sorted_links <- c(sorted_links[sorted_links != "QC_XLSX"], "QC_XLSX")
        self$links <- self$links[sorted_links]
      }
      # Add links
      for (name in names(self$links)) {
        link_path <- basename(self$links[[name]])
        str <- c(
          str,
          paste0("<li><a href='", link_path, "'>", name, "</a></li>")
        )
      }

      str <- c(
        str,
        "</ul>",
        "</body>",
        "</html>"
      )

      cat(str, file = file.path(self$output_dir, "index.html"), sep = "\n")
      # self$links[["INDEX"]] = file.path(self$output_dir, "index.html")
    },
    #' @description
    #' render index markdown file
    render_index_md = function() {
      str <- c(
        paste0("# QC Results for WU : ", self$GRP2$project_spec$workunit_Id, ", and input : ", self$GRP2$software, "\n"),
        "\n## Available Reports\n"
      )

      sorted_links <- names(self$links)
      if ("QC_XLSX" %in% sorted_links) {
        sorted_links <- c(sorted_links[sorted_links != "QC_XLSX"], "QC_XLSX")
        self$links <- self$links[sorted_links]
      }

      # Add links
      for (name in names(self$links)) {
        link_path <- basename(self$links[[name]])
        str <- c(
          str,
          paste0("- [", name, "](", link_path, ")")
        )
      }

      cat(str, file = file.path(self$output_dir, "index.md"), sep = "\n")
    },
    #' @description
    #' get protein per group small wide format
    #' @return protein per group data in small wide format
    get_protein_per_group_small_wide = function() {
      n <- 2
      precabund <- self$get_protein_per_group_abundance()
      tableconfig <- self$get_prot_IBAQ()$config$table
      protID <- tableconfig$hierarchy_keys_depth()
      precabund <- dplyr::inner_join(self$protein_annotation$row_annot,
        precabund,
        by = protID
      )

      precabund_table <- precabund |> dplyr::mutate(
        abundance_percent = signif(abundance_percent, n),
        abundance_percent_cumulative = signif(abundance_percent_cumulative, n),
        percent_prot = signif(percent_prot, 3)
      )
      precabund_table <- precabund_table |>
        dplyr::select(
          all_of(c(
            protID, "nrPeptides",
            tableconfig$factor_keys_depth(),
            "nrMeasured", "meanAbundance", "abundance_percent", "description"
          ))
        )
      factors <- TRUE
      if (factors) {
        precabund_table <- precabund_table |>
          tidyr::pivot_wider(
            names_from = tableconfig$factor_keys_depth(),
            values_from = c("nrMeasured", "meanAbundance", "abundance_percent")
          )
      } else {
        precabund_table <- dplyr::select(precabund_table, -all_of(tableconfig$factor_keys_depth()))
      }
      return(precabund_table)
    }
  ),
  private = list(
    get_target_type = function() {
      if (grepl("MZMINE", self$GRP2$software)) {
        return("metabolite")
      } else if (grepl("PEPTIDE", self$GRP2$software)) {
        return("peptide")
      } else {
        return("protein")
      }
    },

    # Helper method to get nr_children data
    # @param lfqdata LFQData object to get nr_children data from
    # @return data frame with nr_children data
    get_nr_children_data = function(lfqdata) {
      # Get nr_children data using the configured column name
      nr_children_col_name <- lfqdata$config$table$nr_children
      nr_children_data <- lfqdata$to_wide(value = nr_children_col_name)$data
      return(nr_children_data)
    }
  )
)
