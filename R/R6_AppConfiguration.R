# ProcessingOptions ----
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
    #' @field model_missing model missigness, default TRUE
    model_missing = TRUE,
    #' @field model name of the model to use "prolfqua", "SE", "ROPECA", default "prolfqua"
    model = "prolfqua",
    #' @field other list with additional options
    other = NULL
  )
)


# ProjectSpec -----
#' project specification R6 class
#' @export
#' @family ProlfquAppConfig
ProjectSpec <- R6::R6Class(
  "ProjectSpec",
  public = list(
    #' @field project_Id project ID
    project_Id =  character(),
    #' @field project_name project name
    project_name = "",
    #' @field order_Id order ID
    order_Id = character(),
    #' @field workunit_Id workunit ID
    workunit_Id =  character(),
    #' @field input_URL input URL
    input_URL = "https://fgcz-bfabric.uzh.ch/bfabric/"
  )
)

# ExternalReader -----
#' project specification R6 class
#' @export
#' @family ProlfquAppConfig
#' @examples
#' ExternalReader$new()
#'
ExternalReader <- R6::R6Class(
  "ExternalReader",
  public = list(
    #' @field get_files get_files
    get_files = character(),
    #' @field preprocess preprocess
    preprocess = character(),
    #' @field extra_args extra_args
    extra_args = "list()"
    )
)


# ProlfquAppConfig -----
#' R6 class representing
#'
#' @export
#' @family ProlfquAppConfig
#' @examples
#'
#'
#' r6obj_config <- ProlfquAppConfig$new(ProcessingOptions$new(), ProjectSpec$new(), ExternalReader$new())
#' xx <- prolfqua::R6_extract_values(r6obj_config)
#' yaml::write_yaml(xx, file = file.path(tempdir(),"test.yaml"))
#' config <- yaml::read_yaml(file = file.path(tempdir(),"test.yaml"))
#'
#' r6obj_config$set_zipdir_name()
#'
#' r6obj_config$get_zipdir()
#' r6obj_config$get_result_dir()
#' r6obj_config$get_input_dir()
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
    #' @field prefix either QC or DEA
    prefix = character(),
    #' @field zipdir_name results should go to zipdir_name
    zipdir_name = character(),
    #' @field name description
    path = character(),
    #' @field pop optional processing options
    pop = list(),
    #' @field RES resuls
    RES = list(),

    #' @field group name
    group = "G_",

    #' @field pp_name Name of input
    ext_reader = NULL,
    #' @description
    #' set procession options and project spec
    #' @param processing_options instance of ProjectOptions
    #' @param project_spec instance of ProjectSpec
    #' @param zipdir_name where to store results
    #' @param software name of input software
    initialize = function(processing_options,
                          project_spec,
                          ext_reader,
                          zipdir_name = ".",
                          path = ".",
                          software = "DIANN",
                          prefix = "DEA"){
      self$software = software
      self$processing_options = processing_options
      self$project_spec = project_spec
      self$ext_reader = ext_reader
      self$zipdir_name <- zipdir_name
      self$path = path
      self$prefix = prefix
    },
    set_zipdir_name = function(){
        pi <- if (length(self$project_spec$project_Id) == 0 || self$project_spec$project_Id == "") {
          NULL
        } else { paste0("_PI", self$project_spec$project_Id) }
        oi <- if (length(self$project_spec$order_Id) == 0 || self$project_spec$order_Id == "") {
          NULL
        } else { paste0("_O",self$project_spec$order_Id) }
        wu <- if (length(self$project_spec$workunit_Id) == 0 || self$project_spec$workunit_Id == "") {
          NULL
        } else {paste0("_WU",self$project_spec$workunit_Id)}
        res <- paste0(
          self$prefix,
          "_",format(Sys.Date(), "%Y%m%d"),
          pi,
          oi,
          wu,
          "_", self$processing_options$transform)
      self$zipdir_name = res
      return(res)
    },
    get_zipdir = function(){
      return(file.path(self$path, self$zipdir_name))
    },
    get_result_dir = function(){
      tmp <- file.path( self$get_zipdir() , paste0("Results_WU_", self$project_spec$workunit_Id))
      return(tmp)
    },
    get_input_dir = function(){
      tmp <- file.path(self$get_zipdir(), paste0("Inputs_WU_", self$project_spec$workunit_Id))
      return(tmp)
    },
    as_list = function(){
      res <- prolfqua::R6_extract_values(self)
      return(res)
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




#' read minimal yaml
#' @export
#' @examples
#'
#' DEAconfig <- make_DEA_config_R6( WORKUNITID = "3333")
#' configList <- prolfqua::R6_extract_values(DEAconfig)
#' stopifnot(class(configList) == "list")
#' old <- configList$zipdir_name
#' config <- list_to_R6_app_config(configList)
#' stopifnot(config$zipdir_name == old)
#' stopifnot("ProlfquAppConfig" %in% class(config))
#' stopifnot(config$zipdir_name == configList$zipdir_name)
#'
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
  extR6 <- ExternalReader$new()
  ext <- dd$ext_reader
  for (i in names(ext)) {
    extR6[[i]] <- ext[[i]]
  }
  r6obj_config <- ProlfquAppConfig$new(popR6, psR6, extR6)
  r6obj_config$zipdir_name = dd$zipdir_name
  r6obj_config$software = dd$software
  r6obj_config$group = dd$group
  r6obj_config$path = dd$path
  r6obj_config$prefix = dd$prefix

  if (is.null(r6obj_config$zipdir_name)) {
    r6obj_config$set_zipdir_name()
  }
  return(r6obj_config)
}


#' create GRP2 configuration.
#' Use this function if there is no Yaml Input.
#' @param patternDecoys default "^REV_"
#' @param patternContaminants default "^zz_"
#' @export
#' @family ProlfquAppConfig
#' @examples
#'
#' DEAconfig <- make_DEA_config_R6(ORDERID = "1234", WORKUNITID = "1234")
#' DEAconfig$set_zipdir_name()
#' DEAconfig$get_zipdir()
#' DEAconfig$get_result_dir()
#' DEAconfig$get_input_dir()
#' R6list <- prolfqua::R6_extract_values(DEAconfig)
#'
#'
make_DEA_config_R6 <- function(
    PATH = ".",
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
    patternDecoys = "^REV_|^rev_",
    patternContaminants = "^zz|^CON|Cont_",
    application = "DIANN",
    prefix = "DEA"
    ){

  Normalization <- match.arg(Normalization)
  aggregation <- match.arg(aggregation)

  ext <- ExternalReader$new()

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

  r6obj_config <- ProlfquAppConfig$new(pop, ps, ext, prefix = prefix)
  r6obj_config$set_zipdir_name()
  r6obj_config$software = application
  r6obj_config$path <- PATH
  return(r6obj_config)
}


#' read yaml file
#' @export
#' @return list with applications parameters
#' @examples
#' if (FALSE) {
#'   yfile <- prolfqua::find_package_file("prolfquapp", "application/DIANN/config.yaml")
#'   file.exists(yfile)
#'   config <- read_BF_yamlR6(yfile)
#' }
#'
read_BF_yamlR6 <- function(ymlfile, application = "DIANN" ) {
  yml = yaml::read_yaml(ymlfile)

  WORKUNITID = yml$job_configuration$workunit_id
  PROJECTID = yml$job_configuration$project_id
  ORDERID = yml$job_configuration$order_id
  ORDERID <- if (is.null(ORDERID)) { PROJECTID }else{ ORDERID }

  ps <- ProjectSpec$new()
  ps$order_Id = ORDERID
  ps$project_Id = PROJECTID
  ps$workunit_Id = WORKUNITID
  ps$project_name <- ""
  # idxzip <- grep("[0-9]{7,7}.zip|DIANN_Result_WU[0-9]{6,6}.zip",yml$application$input[[1]], ignore.case = TRUE)
  # ps$input_URL <- yml$job_configuration$input[[1]][[idxzip]]$resource_url

  #at least 2 peptides per protein
  pop <- ProcessingOptions$new()
  pop$transform <- yml$application$parameters$`3|Normalization`
  pop$aggregate <- "medpolish"
  pop$diff_threshold <- as.numeric(yml$application$parameters$`4|Difference_threshold`)
  pop$FDR_threshold <- as.numeric(yml$application$parameters$`5|FDR_threshold`)

  pop$remove_cont <- if (yml$application$parameters$`6|remConDec` == "true") { TRUE } else { FALSE }
  pop$remove_decoys <- if (yml$application$parameters$`6|remConDec` == "true") { TRUE } else { FALSE }

  pop$pattern_decoys <- yml$application$parameters$`7|REVpattern`
  pop$pattern_decoys <- if (pop$pattern_decoys == "") { NULL } else {pop$pattern_decoys}
  pop$pattern_contaminants <- yml$application$parameters$`8|CONpattern`
  pop$pattern_contaminants <- if (pop$pattern_contaminants == "") { NULL } else (pop$pattern_contaminants)
  ext <- ExternalReader$new()
  r6obj_config <- ProlfquAppConfig$new(pop, ps, ext)
  r6obj_config$set_zipdir_name()
  r6obj_config$software = application

  return(r6obj_config)
}


#' get configuration from yaml if exists
#' @export
#' @examples
#'
#'
#' get_config()
get_config <- function(yamlfile, WORKUNITID =  "HelloWorld", ORDERID = "123") {
  if (missing(yamlfile)) {
    GRP2 <- prolfquapp::make_DEA_config_R6(
      PROJECTID = as.character(ORDERID) ,ORDERID = as.character(ORDERID), WORKUNITID = WORKUNITID )
  } else if (file.exists(yamlfile)) {
    xx <- yaml::read_yaml(yamlfile)
    if (!is.null(xx$project_spec)) {
      logger::log_info("prolfquapp yaml")
      GRP2 <- list_to_R6_app_config(xx)
      GRP2$set_zipdir_name()
    } else {
      GRP2 <- yamlfile |> prolfquapp::read_BF_yamlR6(application = "DIANN")
      GRP2$set_zipdir_name()
      logger::log_info("bfabric yaml")
    }
  } else {
    stop("no such file :", yamlfile )
  }
  return(GRP2)
}



if (FALSE) {
  yfile <- file.path(find.package("prolfquapp") , "/application/DIANN/myYamls.zip")
  file.exists(yfile)
  xx <- unzip(yfile, list = TRUE)
  yfiles <- grep(".yml$", xx$Name, value = TRUE)

  res <- list()
  for (file in yfiles) {
    config <- read_BF_yamlR6(unz(yfile, file))
    x <- (prolfqua::R6_extract_values(config))
    df <- data.frame(unlist(x))
    names(df)[1] <- basename(file)
    res[[basename(file)]] <- as.data.frame(t(df))
  }

}





