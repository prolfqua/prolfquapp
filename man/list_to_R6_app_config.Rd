% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R6_AppConfiguration.R
\name{list_to_R6_app_config}
\alias{list_to_R6_app_config}
\title{read minimal yaml and convert to R6 object}
\usage{
list_to_R6_app_config(dd)
}
\arguments{
\item{dd}{list containing configuration data}
}
\value{
ProlfquAppConfig R6 object
}
\description{
read minimal yaml and convert to R6 object
}
\examples{

DEAconfig <- make_DEA_config_R6(WORKUNITID = "3333")
configList <- prolfqua::R6_extract_values(DEAconfig)
stopifnot(class(configList) == "list")
old <- configList$zipdir_name
config <- list_to_R6_app_config(configList)
stopifnot(config$zipdir_name == old)
stopifnot("ProlfquAppConfig" \%in\% class(config))
stopifnot(config$zipdir_name == configList$zipdir_name)

}
