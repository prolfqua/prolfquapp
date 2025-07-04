% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess_MaxQuant.R, R/tidyMS_MaxQuant.R
\name{tidyMQ_Evidence}
\alias{tidyMQ_Evidence}
\title{read evidence file}
\usage{
tidyMQ_Evidence(Evidence)

tidyMQ_Evidence(Evidence)
}
\arguments{
\item{Evidence}{MQ evidence file or zip archive with evidence file}
}
\description{
read evidence file

read evidence file
}
\examples{

evidence_txt <- prolfqua::find_package_file("prolfquapp", "samples/maxquant_txt/tiny2.zip")
evidence_txt <- read.csv(unz(evidence_txt,"evidence.txt"), header=TRUE, stringsAsFactors = FALSE, sep="\t")
mq_evidence <- tidyMQ_Evidence(evidence_txt)

evidence_txt <- prolfqua::find_package_file("prolfquapp", "samples/maxquant_txt/tiny2.zip")
evidence_txt <- read.csv(unz(evidence_txt,"evidence.txt"), header=TRUE, stringsAsFactors = FALSE, sep="\t")
mq_evidence <- tidyMQ_Evidence(evidence_txt)
}
\seealso{
Other MaxQuant: 
\code{\link{MaxQuant}},
\code{\link{tidyMQ_Peptides}()},
\code{\link{tidyMQ_ProteinGroups}()}

Other MaxQuant: 
\code{\link{MaxQuant}},
\code{\link{tidyMQ_Peptides}()},
\code{\link{tidyMQ_ProteinGroups}()}
}
\concept{MaxQuant}
\keyword{internal}
