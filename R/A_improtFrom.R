#' @importFrom ggplot2 aes ggplot geom_point facet_wrap scale_y_continuous
#' @importFrom ggrepel geom_label_repel
#' @importFrom dtplyr lazy_dt
#' @importFrom rlang sym syms .data
#' @importFrom tidyselect starts_with all_of ends_with one_of
#' @importFrom dplyr distinct n
#' @importFrom lobstr tree
#' @importFrom UpSetR upset
#' @importFrom vsn justvsn
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom arrow write_parquet
#' @importFrom optparse make_option
#' @importFrom pander pander
#' @importFrom bookdown render_book
#' @importFrom crosstalk SharedData
#' @importFrom purrr map
#' @importFrom readr read_tsv
#' @importFrom readxl read_excel
#' @importFrom stringi stri_replace_all_fixed
#' @importFrom stringr str_replace
#' @importFrom tibble tibble
#' @importFrom tidyr pivot_longer
#' @importFrom writexl write_xlsx
#' @importFrom yaml read_yaml
#' @importFrom stats as.formula na.omit rlnorm rt
#' @importFrom utils getFromNamespace installed.packages packageVersion read.csv tail write.table

NULL

#' set library path with logging
#' @param lib_path path to R library directory
#' @export
set_lib_path <- function(lib_path){
  if (!is.null(lib_path) && dir.exists(lib_path) ) {
    logger::log_info(paste("Setting libPath:", lib_path, collapse = " ;"))
    .libPaths(lib_path)
    logger::log_info(.libPaths(), sep = "\n")
  }
}
