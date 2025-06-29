% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess_DIANN.R
\name{diann_read_output}
\alias{diann_read_output}
\title{read DiaNN diann-output.tsv file}
\usage{
diann_read_output(data, Lib.PG.Q.Value = 0.01, PG.Q.Value = 0.05)
}
\description{
filter for 2 peptides per protein, and for Q.Value < 0.01 (default)
}
\examples{
\dontrun{
xx <- readr::read_tsv("WU292720_report.tsv")
report2 <- prolfquapp::diann_read_output(xx)
nrow(report2)
}

}
