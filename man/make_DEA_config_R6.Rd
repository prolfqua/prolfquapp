% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R6_AppConfiguration.R
\name{make_DEA_config_R6}
\alias{make_DEA_config_R6}
\title{create GRP2 configuration.
Use this function if there is no Yaml Input.}
\usage{
make_DEA_config_R6(
  PATH = ".",
  PROJECTID = "",
  ORDERID = "",
  WORKUNITID = "",
  Normalization = c("none", "vsn", "quantile", "robscale"),
  aggregation = c("medpolish", "top3", "lmrob"),
  diff_threshold = 1,
  FDR_threshold = 0.1,
  nr_peptides = 1,
  removeContaminants = FALSE,
  removeDecoys = FALSE,
  patternDecoys = "^REV_|^rev_",
  patternContaminants = "^zz|^CON|Cont_",
  application = "DIANN",
  prefix = "DEA"
)
}
\arguments{
\item{patternDecoys}{default "^REV_"}

\item{patternContaminants}{default "^zz_"}
}
\description{
create GRP2 configuration.
Use this function if there is no Yaml Input.
}
\examples{

DEAconfig <- make_DEA_config_R6(ORDERID = "1234", WORKUNITID = "1234")
DEAconfig$set_zipdir_name()
DEAconfig$get_zipdir()
DEAconfig$get_result_dir()
DEAconfig$get_input_dir()
R6list <- prolfqua::R6_extract_values(DEAconfig)


}
\seealso{
Other ProlfquAppConfig: 
\code{\link{ExternalReader}},
\code{\link{ProcessingOptions}},
\code{\link{ProjectSpec}},
\code{\link{ProlfquAppConfig}},
\code{\link{set_list_to_R6}()}
}
\concept{ProlfquAppConfig}
