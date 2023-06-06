#' processing options R6 class
#' @export
#' @family ProlfquAppConfig
ProcessingOptions <- R6::R6Class(
  "ProcessingOptions",
  public = list(
    #' @field transform data transformation method
    transform = "vsn",
    #' @field aggregate protein abundance estimation method
    aggregate = "medpolish",
    #' @field diff_threshold difference threshold
    diff_threshold = 1,
    #' @field FDR_threshold FDR threshold
    FDR_threshold = 0.1,
    #' @field remove_cont should contaminants be removed
    remove_cont = TRUE,
    #' @field remove_decoys should decoys be removed
    remove_decoys = TRUE,
    #' @field pattern_decoys decoy patterns
    pattern_decoys = "^REV",
    #' @field pattern_contaminants pattern contaminants
    pattern_contaminants = "^CON|^zz"
  )
)


#' project specification R6 class
#' @export
#' @family ProlfquAppConfig
ProjectSpec <- R6::R6Class(
  "ProjectSpec",
  public = list(
    #' @field projectID project ID
    projectID =  as.integer(23689),
    #' @field projectName project name
    projectName = "",
    #' @field orderID order ID
    orderID = as.integer(24227),
    #' @field workunitID workunit ID
    workunitID =  as.integer(289177),
    #' @field inputID input id
    inputID =  as.integer(2294093),
    #' @field inputURL input URL
    inputURL = "https://fgcz-bfabric.uzh.ch/bfabric/userlab/show-resource.html?id=2294093"
  )
)


#' R6 class representing
#'
#' @export
#' @family ProlfquAppConfig
#' @examples
#'
#'
#' r6obj_config <- ProlfquAppConfig$new(ProcessingOptions$new(), ProjectSpec$new())
#' xx <- R6_extract_values(r6obj_config)
#' yaml::write_yaml(xx, file = "test.yaml")
#' config <- yaml::read_yaml(file = "test.yaml")
ProlfquAppConfig <- R6::R6Class(
  "ProlfquAppConfig",
  public = list(
    #' @field processing_options ProcessingOption R6 class
    processing_options = NULL,
    #' @field project_spec Project Spec R6 class
    project_spec = NULL,
    #' @field software name of input software
    software = character(),
    #' @field zipdir results should go to zipdir
    zipdir = character(),
    #' @description
    #' set procession options and project spec
    #' @param processing_options instance of ProjectOptions
    #' @param project_spec instance of ProjectSpec
    #' @param zipdir where to store results
    #' @param software name of input software
    initialize = function(processing_options, project_spec, zipdir = NULL, software = "DIANN"){
      self$software = software
      self$zipdir <- if (!is.null(zipdir)) {
        zipdir
      } else {
        paste0("C",project_spec$projectID ,"WU",project_spec$workunitID)
      }
      self$processing_options = processing_options
      self$project_spec = project_spec
    }
  )
)




#' set arguments in list config to r6obj
#' @export
#' @family ProlfquAppConfig
#'
set_list_to_R6 <- function(config_list, r6obj_config){
  for (n in seq_along(config)) {
    if (class(config[[n]]) == "list") {
      message(paste0("setting fields in :", names(config)[n], "\n"))
      r6component <- r6obj_config[[names(config)[n]]]
      set_config(config[[n]], r6component)
    } else {
      cat(n , ":" , names(config)[n], " = ", config[[n]], "\n")
      r6obj_config[[names(config)[n]]] <- config[[n]]
    }
  }
}

#' extract R6 class fields as list
#' @export
#' @family ProlfquAppConfig
R6_extract_values <- function(r6class){
  tmp <- sapply(r6class, class)
  slots <- tmp[ !tmp %in% c("environment", "function")]
  res <- list()
  for (i in names(slots)) {
    if ("R6" %in% class(r6class[[i]])) {
      res[[i]]  <- R6_extract_values(r6class[[i]])
    }else{
      res[[i]] <- r6class[[i]]
    }
  }
  return(res)
}










