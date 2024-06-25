#' @importFrom ggplot2 aes ggplot geom_point facet_wrap  scale_y_continuous
#' @importFrom ggrepel geom_label_repel
#' @importFrom dtplyr lazy_dt
#' @importFrom rlang sym syms
#' @importFrom tidyselect starts_with all_of
#' @importFrom dplyr distinct n

NULL

#' set library path with logging
#' @export
set_lib_path <- function(lib_path){
  if (!is.null(lib_path) && dir.exists(lib_path) ) {
    logger::log_info(paste("Setting libPath:", lib_path, collapse = " ;"))
    .libPaths(lib_path)
    logger::log_info(.libPaths(), sep = "\n")
  }
}
