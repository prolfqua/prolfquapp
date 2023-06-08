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
    projectID =  integer(),
    #' @field projectName project name
    projectName = "",
    #' @field orderID order ID
    orderID = integer(),
    #' @field workunitID workunit ID
    workunitID =  integer(),
    #' @field inputID input id
    inputID =  integer(),
    #' @field inputURL input URL
    inputURL = "https://fgcz-bfabric.uzh.ch/bfabric/"
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



#' create GRP2 configuration.
#' Use this function if there is no Yaml Input.
#' @param patternDecoys default "REV_"
#' @param patternContaminants default "zz_"
#' @export
#' @family ProlfquAppConfig
#' @examples
#' DEAconfig <- make_DEA_config_R6()
#' R6list <- R6_extract_values(DEAconfig)
make_DEA_config_R6 <- function(
    ZIPDIR  = ".",
    PROJECTID = "",
    ORDERID ="",
    WORKUNITID ="",
    Normalization = c("vsn", "quantile", "robscale"),
    aggregation = c("medpolish" , "top3", "lmrob"),
    Diffthreshold = 1,
    FDRthreshold = 0.1,
    removeContaminants = FALSE,
    removeDecoys = FALSE,
    patternDecoys = "REV_",
    patternContaminants = "zz",
    application = "FragPipeTMT" ){

  pop <- ProcessingOptions$new()
  pop$pattern_contaminants = patternContaminants
  pop$pattern_decoys = patternDecoys
  pop$remove_cont = removeContaminants
  pop$remove_decoys = removeDecoys
  pop$FDR_threshold = FDRthreshold
  pop$diff_threshold = Diffthreshold
  pop$aggregate = aggregation
  pop$transform = Normalization
  ps <- ProjectSpec$new()
  ps$orderID = ORDERID
  ps$projectID = PROJECTID
  ps$workunitID = WORKUNITID

  r6obj_config <- ProlfquAppConfig$new(pop, ps)
  r6obj_config$zipdir = ZIPDIR
  r6obj_config$software = application

  return(r6obj_config)
}







