% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess_FP_PSM.R
\name{tidy_FragPipe_combined_protein}
\alias{tidy_FragPipe_combined_protein}
\title{read combined_protein.tsv file for FragPipe Version 16 or newer}
\usage{
tidy_FragPipe_combined_protein(
  combprot,
  as_list = FALSE,
  spcnames = c("Total Spectral Count", "Unique Spectral Count", "Razor Spectral Count"),
  intnames = c("Total Intensity", "Unique Intensity", "Razor Intensity"),
  maxlfqnames = c("MaxLFQ Total Intensity", "MaxLFQ Unique Intensity",
    "MaxLFQ Razor Intensity")
)
}
\arguments{
\item{combprot}{path to combined_protein.tsv file}

\item{as_list}{return as list}
}
\value{
tidy dataframe or list with df (e.g. total.spectral.count or total.intensity etc).
}
\description{
read combined_protein.tsv file for FragPipe Version 16 or newer
}
\seealso{
Other FragPipe: 
\code{\link{FragPipe}},
\code{\link{tidy_FragPipe_MSstats_csv}()},
\code{\link{tidy_FragPipe_combined_protein_deprec}()}
}
\concept{FragPipe}
\keyword{internal}
