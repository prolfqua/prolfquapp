% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_DEA_report2.R
\name{bfabric_url_builder}
\alias{bfabric_url_builder}
\title{build bfabric urls}
\usage{
bfabric_url_builder(project_spec)
}
\description{
build bfabric urls
}
\examples{


ps <- ProjectSpec$new()
ps$project_Id <- 32258
ps$order_Id <- 34628
ps$workunit_Id <- 302212
bfabric_url_builder(ps)

ps <- ProjectSpec$new()
ps$order_Id <- 34628
ps$workunit_Id <- 302212
bfabric_url_builder(ps)

}
