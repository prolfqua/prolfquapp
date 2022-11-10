#' read yaml file
#' @export
#' @return list with applications parameters
read_yaml <- function(ymlfile, application = "FragPipeTMT" ) {
  yml = yaml::read_yaml(ymlfile)

  WORKUNITID = yml$job_configuration$workunit_id
  PROJECTID = yml$job_configuration$project_id
  ORDERID = yml$job_configuration$order_id
  ORDERID <- if (is.null(ORDERID)) { PROJECTID }else{ ORDERID }
  ZIPDIR = paste0("C",ORDERID,"WU",WORKUNITID)

  GRP2 <- list()
  GRP2$Bfabric <- list()
  GRP2$Bfabric$projectID <- PROJECTID


  GRP2$Bfabric$projectName <- "" # workunit name in the future.
  GRP2$Bfabric$orderID <- ORDERID

  GRP2$Bfabric$workunitID <- WORKUNITID

  idxzip <- grep("*.zip",yml$application$input[[1]])

  GRP2$Bfabric$inputID <- yml$job_configuration$input[[1]][[idxzip]]$resource_id
  GRP2$Bfabric$inputURL <- yml$job_configuration$input[[1]][[idxzip]]$resource_url

  #at least 2 peptides per protein
  GRP2$pop <- list()
  GRP2$pop$transform <- yml$application$parameters$`3|Normalization`

  GRP2$pop$aggregate <- "medpolish"
  GRP2$pop$Diffthreshold <- as.numeric(yml$application$parameters$`4|Difference_threshold`)
  GRP2$pop$FDRthreshold <- as.numeric(yml$application$parameters$`5|FDR_threshold`)

  GRP2$pop$remove <- if (yml$application$parameters$`6|remConDec` == "true") { TRUE } else { FALSE }
  GRP2$pop$revpattern <- yml$application$parameters$`7|REVpattern`
  GRP2$pop$contpattern <- yml$application$parameters$`8|CONpattern`

  GRP2$Software <- application
  GRP2$zipdir <- ZIPDIR
  return(GRP2)
}
