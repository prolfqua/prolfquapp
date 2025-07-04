% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_function_in_packages.r
\name{get_procfuncs}
\alias{get_procfuncs}
\title{Get all processing functions from all packages}
\usage{
get_procfuncs(
  function_name = "prolfqua_preprocess_functions",
  prefix = "prolfqua"
)
}
\arguments{
\item{function_name}{The name of the function to get}

\item{prefix}{The prefix of the package names}
}
\value{
A list of processing functions
}
\description{
Get all processing functions from all packages
}
\examples{
get_procfuncs()
get_procfuncs("prolfqua_preprocess_functions", "prolfqua")
get_procfuncs("prolfqua_preprocess_functions", "prolfquapp")
get_procfuncs("prolfqua_preprocess_functions", "prolfquapp")
get_procfuncs("prolfqua_preprocess_functions", "xdx")

}
