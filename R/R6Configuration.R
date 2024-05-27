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
    pattern_contaminants = "^CON|^zz",
    #' @field nr_peptides number of peptides
    nr_peptides = 1,
    #' @field interaction model with interactions default FALSE
    interaction = FALSE,
    #' @field missing model missigness, default TRUE
    missing = TRUE,
    #' @field model name of the model to use "prolfqua", "SE", "ROPECA", default "prolfqua"
    model = "prolfqua",
    #' @field other list with additional options
    other = NULL
  )
)


#' project specification R6 class
#' @export
#' @family ProlfquAppConfig
ProjectSpec <- R6::R6Class(
  "ProjectSpec",
  public = list(
    #' @field projectID project ID
    project_Id =  integer(),
    #' @field project_name project name
    project_name = "",
    #' @field orderID order ID
    order_Id = integer(),
    #' @field workunitID workunit ID
    workunit_Id =  integer(),
    #' @field inputURL input URL
    input_URL = "https://fgcz-bfabric.uzh.ch/bfabric/"
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
#' yaml::write_yaml(xx, file = file.path(tempdir(),"test.yaml"))
#' config <- yaml::read_yaml(file = file.path(tempdir(),"test.yaml"))
#'
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
    #' @field pop optional processing options
    pop = list(),
    #' @field RES resuls
    RES = list(),
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
        paste0("C",project_spec$project_Id ,"WU",project_spec$workunit_Id)
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

#' read minimal yaml
#' @export
#' @examples
#'
#' DEAconfig <- make_DEA_config_R6(ZIPDIR = "DEA", WORKUNITID = "3333")
#' configList <- R6_extract_values(DEAconfig)
#' stopifnot(class(configList) == "list")
#' configList$zipdir
#' config <- list_to_R6_app_config(configList)
#' stopifnot("ProlfquAppConfig" %in% class(config))
#' stopifnot(config$zipdir == configList$zipdir)
list_to_R6_app_config <- function(dd){
  popR6 <- ProcessingOptions$new()
  pop <- dd$processing_options
  for (i in names(pop)) {
    popR6[[i]] <- pop[[i]]
  }
  psR6 <- ProjectSpec$new()
  ps <- dd$project_spec
  for (i in names(ps)) {
    psR6[[i]] <- ps[[i]]
  }
  r6obj_config <- ProlfquAppConfig$new(popR6, psR6)
  r6obj_config$zipdir = dd$zipdir
  r6obj_config$software = dd$software
  return(r6obj_config)
}


#' create GRP2 configuration.
#' Use this function if there is no Yaml Input.
#' @param patternDecoys default "^REV_"
#' @param patternContaminants default "^zz_"
#' @export
#' @family ProlfquAppConfig
#' @examples
#' DEAconfig <- make_DEA_config_R6()
#' R6list <- R6_extract_values(DEAconfig)
#'
make_DEA_config_R6 <- function(
    ZIPDIR  = ".",
    PROJECTID = "",
    ORDERID ="",
    WORKUNITID ="",
    Normalization = c("none", "vsn", "quantile", "robscale"),
    aggregation = c("medpolish" , "top3", "lmrob"),
    diff_threshold = 1,
    FDR_threshold = 0.1,
    nr_peptides = 1,
    removeContaminants = FALSE,
    removeDecoys = FALSE,
    patternDecoys = "^REV_",
    patternContaminants = "^zz",
    application = "FragPipeTMT" ){

  Normalization <- match.arg(Normalization)
  aggregation <- match.arg(aggregation)

  pop <- ProcessingOptions$new()
  pop$pattern_contaminants = patternContaminants
  pop$pattern_decoys = patternDecoys
  pop$remove_cont = removeContaminants
  pop$remove_decoys = removeDecoys
  pop$FDR_threshold = FDR_threshold
  pop$diff_threshold = diff_threshold
  pop$aggregate = aggregation
  pop$transform = Normalization
  pop$nr_peptides = nr_peptides
  ps <- ProjectSpec$new()
  ps$order_Id = ORDERID
  ps$project_Id = PROJECTID
  ps$workunit_Id = WORKUNITID

  r6obj_config <- ProlfquAppConfig$new(pop, ps)
  pi <- if (PROJECTID != "") { paste0("_PI_", PROJECTID)} else {NULL}
  r6obj_config$zipdir = paste0(ZIPDIR, pi , "_OI_",ORDERID,"_WU_",WORKUNITID, "_", Normalization)
  r6obj_config$software = application

  return(r6obj_config)
}


#' read yaml file
#' @export
#' @return list with applications parameters
#' @examples
#' if(FALSE){
#'   yfile <- prolfqua::find_package_file("prolfquapp", "/inst/application/DIANN/myYamls.zip")
#'   file.exists(yfile)
#'   yfiles <- dir(yfile,recursive = TRUE,full.names = TRUE)
#'   config <- read_BF_yamlR6(yfiles[1])
#' }
#'
read_BF_yamlR6 <- function(ymlfile, application = "FragPipeTMT" ) {
  yml = yaml::read_yaml(ymlfile)

  WORKUNITID = yml$job_configuration$workunit_id
  PROJECTID = yml$job_configuration$project_id
  ORDERID = yml$job_configuration$order_id
  ORDERID <- if (is.null(ORDERID)) { PROJECTID }else{ ORDERID }
  ZIPDIR = paste0("C",ORDERID,"WU",WORKUNITID)

  ps <- ProjectSpec$new()
  ps$order_Id = ORDERID
  ps$project_Id = PROJECTID
  ps$workunit_Id = WORKUNITID
  ps$project_name <- ""
  idxzip <- grep("[0-9]{7,7}.zip",yml$application$input[[1]])

  ps$input_URL <- yml$job_configuration$input[[1]][[idxzip]]$resource_url

  #at least 2 peptides per protein
  pop <- ProcessingOptions$new()
  pop$transform <- yml$application$parameters$`3|Normalization`
  pop$aggregate <- "medpolish"
  pop$diff_threshold <- as.numeric(yml$application$parameters$`4|Difference_threshold`)
  pop$FDR_threshold <- as.numeric(yml$application$parameters$`5|FDR_threshold`)

  pop$remove_cont <- if (yml$application$parameters$`6|remConDec` == "true") { TRUE } else { FALSE }
  pop$remove_decoys <- if (yml$application$parameters$`6|remConDec` == "true") { TRUE } else { FALSE }

  pop$pattern_decoys <- yml$application$parameters$`7|REVpattern`
  pop$pattern_contaminants <- yml$application$parameters$`8|CONpattern`

  r6obj_config <- ProlfquAppConfig$new(pop, ps)
  r6obj_config$zipdir = ZIPDIR
  r6obj_config$software = application

  return(r6obj_config)
}



if (FALSE ) {
  yfile <- file.path(find.package("prolfquapp") , "/application/DIANN/myYamls.zip")
  file.exists(yfile)
  xx <- unzip(yfile, list= TRUE)
  yfiles <- grep(".yml$", xx$Name, value = TRUE)

  res <- list()
  for (file in yfiles) {
    config <- read_BF_yamlR6(unz(yfile, file))
    x <- (R6_extract_values(config))
    df <- data.frame(unlist(x))
    names(df)[1] <- basename(file)
    res[[basename(file)]] <- as.data.frame(t(df))
  }

}





