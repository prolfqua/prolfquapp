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
    #' @field protein_annotation protein_annotation
    protein_annotation = NULL,
    #' @field lfqdata_peptide lfqdata_peptide
    lfqdata_peptide = NULL,
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
    #' @param lfq_data lfq_data
    #' @param protein_annotation protein_annotation
    #' @param prolfquapp_config prolfquapp_config
    initialize = function(lfqdata, protein_annotation, prolfquapp_config) {
      self$GRP2 = prolfquapp_config
      self$lfqdata = lfqdata
      self$protein_annotation = protein_annotation
      self$output_dir = self$GRP2$get_zipdir()
      self$TABLES2WRITE = list()

    },
    get_peptides_wide = function(){
      lfqdata <- self$lfqdata$get_copy()
      lfqdata$config$table$hierarchyDepth <- min(2, length(self$lfqdata$config$table$hierarchyKeys()))
      self$lfqdata_peptide <- prolfquapp::aggregate_data(lfqdata, agg_method = "medpolish")
      peptide_wide <- dplyr::left_join(self$protein_annotation$row_annot,
                                       self$lfqdata$to_wide()$data,
                                       multiple = "all")
      invisible(peptide_wide)
    },
    get_annotation = function(){
      annotation <- self$lfqdata$factors()
      invisible(annotation)
    },
    get_prot_data = function(){
      if (is.null(self$lfqdata_prot)) {
        self$lfqdata_prot <- prolfquapp::aggregate_data(self$lfqdata, agg_method = "medpolish")
      }
      invisible(self$lfqdata_prot)
    },
    get_prot_wide = function(){
      lfqdata_prot = self$get_prot_data()
      proteins_wide <- dplyr::left_join(
        self$protein_annotation$row_annot,
        lfqdata_prot$to_wide()$data,
        multiple = "all")

      proteins_wide <- dplyr::inner_join(
        proteins_wide,
        lfqdata_prot$to_wide(value = "nr_children")$data,
        by = c(lfqdata_prot$config$table$hierarchy_keys_depth(),"isotopeLabel")
        , suffix = c("_abundance","_nr_children"))
      return(proteins_wide)
    },
    get_prot_IBAQ = function(){
      relevant_columns <- c("protein_length", "nr_tryptic_peptides")
      if ( is.null(self$lfqdata_prot_IBAQ ) && all(relevant_columns %in% colnames(self$protein_annotation$row_annot)) ) {
        self$lfqdata_prot_IBAQ <- prolfquapp::compute_IBAQ_values(self$lfqdata, self$protein_annotation)
      } else if (!all(relevant_columns %in% colnames(self$protein_annotation$row_annot)) ) {
        warning("skipping IBAQ computation, no:", paste(relevant_columns, collapse = "; "))
      }
      invisible(self$lfqdata_prot_IBAQ )
    },
    get_protein_per_group_abundance = function(){
      summarizer <- self$get_prot_IBAQ()$get_Summariser()
      precabund <- summarizer$percentage_abundance()
      invisible(precabund)
    },
    get_protein_per_group_abundance_with_row_annot = function(){
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
    get_protein_per_group_abundance_wide = function(){
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
        by = self$get_prot_IBAQ()$config$table$hierarchy_keys_depth())

      invisible(precabund)
    },

    get_prot_IBAQ_wide = function(){
      if (!is.null(self$get_prot_IBAQ())) {
      IBAQ_abundances <-
        dplyr::left_join(self$protein_annotation$row_annot,
                         self$get_prot_IBAQ()$to_wide()$data,
                         multiple = "all")
      IBAQ_abundances <- dplyr::inner_join(
        IBAQ_abundances,
        self$get_prot_IBAQ()$to_wide(value = "nr_children")$data,
        by = c(self$get_prot_IBAQ()$config$table$hierarchy_keys_depth(),"isotopeLabel")
        , suffix = c("_iBAQ","_nr_children"))
      return(IBAQ_abundances)
      } else {return(NULL)}
    },

    get_list = function(){
      TABLES2WRITE <- list()
      TABLES2WRITE$peptide_wide <- self$get_peptides_wide()
      TABLES2WRITE$annotation <- self$get_annotation()
      TABLES2WRITE$prot_medpolish_estimate <- self$get_prot_wide()
      TABLES2WRITE$prot_IBAQ_estimate <- self$get_prot_IBAQ_wide()
      TABLES2WRITE$prot_IBAQ_per_group_stats <- self$get_protein_per_group_abundance_wide()
      return(TABLES2WRITE)
    },
    write_xlsx = function(){
      xlsxfile = file.path(self$output_dir,
                           paste0("proteinAbundances_",self$GRP2$project_spec$workunit_Id,".xlsx"))
      writexl::write_xlsx(self$get_list(), path = xlsxfile)
      self$links[["QC_XLSX"]] = xlsxfile
    },

    render_QC_protein_abundances = function(){
      file.copy(system.file("application/GenericQC/QC_ProteinAbundances.Rmd", package = "prolfquapp"),
                to = self$output_dir, overwrite = TRUE)
      if (TRUE) {
        rmarkdown::render(file.path(self$output_dir,"QC_ProteinAbundances.Rmd"),
                          params = list(pap = self,
                                        project_info = self$GRP2$project_spec,
                                        factors = TRUE),
                          output_file = "proteinAbundances.html")
      } else {
        str <- c("<!DOCTYPE html>",
                 "<html>",
                 "<head>","<title>",
                 "There is a problem",
                 "</title>","</head>",
                 "<body>",
                 "<h1>",
                 paste0("the input file :" , files$data , " is empty"),
                 "</h1>",
                 "</body>",
                 "</html>")
        cat(str, file = file.path(self$output_dir,"proteinAbundances.html"), sep = "\n")
      }
      self$links[["QC_ABUNDANCES"]] = file.path(self$output_dir,"proteinAbundances.html")
    },
    render_sample_size_QC = function(){
      if (nrow(self$get_prot_data()$factors()) > 1) {
        file.copy(system.file("application/GenericQC/QCandSSE.Rmd", package = "prolfquapp"),
                  to = self$output_dir, overwrite = TRUE)
        rmarkdown::render(file.path(self$output_dir,"QCandSSE.Rmd"),
                          params = list(data = self$get_prot_data()$data,
                                        configuration = self$get_prot_data()$config,
                                        project_conf = GRP2$project_spec,
                                        pep = FALSE),
                          output_file = "QC_sampleSizeEstimation.html"
        )
      } else{
        message("only a single sample: ", nrow(self$get_prot_data()$factors()))
      }
      self$links[["QC_SAMPLE_SIZE"]] = file.path(self$output_dir,"QC_sampleSizeEstimation.html")
    },

    render_index_html = function(){

      str <- c("<!DOCTYPE html>",
               "<html>",
               "<head>",
               paste0("<title>QC Results for WU : ", self$GRP2$project_spec$workunit_Id, " and input : ", self$GRP2$software, "</title>"),
               "</head>",
               "<body>",
               paste0("<h1>QC Results for WU : ", self$GRP2$project_spec$workunit_Id, " and input : ", self$GRP2$software, "</h1>"),
               "<ul>")
      # Sort links to ensure QC_XLSX is last
      sorted_links <- names(self$links)
      if ("QC_XLSX" %in% sorted_links) {
        sorted_links <- c(sorted_links[sorted_links != "QC_XLSX"], "QC_XLSX")
        self$links <- self$links[sorted_links]
      }
      # Add links
      for(name in names(self$links)) {
        link_path <- basename(self$links[[name]])
        str <- c(str,
                paste0("<li><a href='", link_path, "'>", name, "</a></li>"))
      }

      str <- c(str,
               "</ul>",
               "</body>",
               "</html>")

      cat(str, file = file.path(self$output_dir, "index.html"), sep = "\n")
      #self$links[["INDEX"]] = file.path(self$output_dir, "index.html")
    },


    render_index_md = function(){
      str <- c(paste0("# QC Results for WU : ", self$GRP2$project_spec$workunit_Id,  ", and input : ", self$GRP2$software, "\n"),
               "\n## Available Reports\n")

      sorted_links <- names(self$links)
      if ("QC_XLSX" %in% sorted_links) {
        sorted_links <- c(sorted_links[sorted_links != "QC_XLSX"], "QC_XLSX")
        self$links <- self$links[sorted_links]
      }

      # Add links
      for(name in names(self$links)) {
        link_path <- basename(self$links[[name]])
        str <- c(str,
                paste0("- [", name, "](", link_path, ")"))
      }

      cat(str, file = file.path(self$output_dir, "index.md"), sep = "\n")
    },

    get_protein_per_group_small_wide = function(){
      n = 2
      precabund = self$get_protein_per_group_abundance()
      tableconfig = self$get_prot_IBAQ()$config$table
      protID <- tableconfig$hierarchy_keys_depth()
      precabund <- dplyr::inner_join(self$protein_annotation$row_annot,
                                     precabund,
                                     by = protID)

      precabund_table <- precabund |> dplyr::mutate(
        abundance_percent = signif(abundance_percent, n ),
        abundance_percent_cumulative = signif(abundance_percent_cumulative, n),
        percent_prot = signif(percent_prot, 3))
      precabund_table <- precabund_table |>
        dplyr::select(
          all_of(c(protID, "nrPeptides",
                   tableconfig$factor_keys_depth(),
                   "nrMeasured", "meanAbundance", "abundance_percent" ,"description")))
      factors = TRUE
      if (factors) {
        precabund_table <- precabund_table |>
          tidyr::pivot_wider(
            names_from = tableconfig$factor_keys_depth(),
            values_from = c("nrMeasured", "meanAbundance", "abundance_percent"))
      } else {
        precabund_table <- dplyr::select(precabund_table, -all_of(tableconfig$factor_keys_depth()))
      }
      return(precabund_table)
    }
  )
)
