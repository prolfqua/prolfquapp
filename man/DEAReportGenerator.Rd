% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R6_DEAReportGenerator.R
\name{DEAReportGenerator}
\alias{DEAReportGenerator}
\title{DEAReportGenerator}
\description{
DEAReportGenerator

DEAReportGenerator
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{lfqdata}}{LFQData object containing the quantitative data}

\item{\code{GRP2}}{ProlfquAppConfig object containing analysis configuration}

\item{\code{prot_annot}}{ProteinAnnotation object}

\item{\code{Contrasts}}{list of contrasts for differential expression analysis}

\item{\code{fname}}{filename for DEA results}

\item{\code{qcname}}{filename for QC results}

\item{\code{resultdir}}{directory for storing results}

\item{\code{ZIPDIR}}{zip directory path}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-DEAReportGenerator-new}{\code{DEAReportGenerator$new()}}
\item \href{#method-DEAReportGenerator-write_DEA_all}{\code{DEAReportGenerator$write_DEA_all()}}
\item \href{#method-DEAReportGenerator-render_DEA}{\code{DEAReportGenerator$render_DEA()}}
\item \href{#method-DEAReportGenerator-make_boxplots}{\code{DEAReportGenerator$make_boxplots()}}
\item \href{#method-DEAReportGenerator-prep_result_list}{\code{DEAReportGenerator$prep_result_list()}}
\item \href{#method-DEAReportGenerator-make_SummarizedExperiment}{\code{DEAReportGenerator$make_SummarizedExperiment()}}
\item \href{#method-DEAReportGenerator-clone}{\code{DEAReportGenerator$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-DEAReportGenerator-new"></a>}}
\if{latex}{\out{\hypertarget{method-DEAReportGenerator-new}{}}}
\subsection{Method \code{new()}}{
Initialize DEAReportGenerator with data and configuration
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DEAReportGenerator$new(lfqdata, GRP2, prot_annot, Contrasts)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{lfqdata}}{LFQData object containing quantitative data}

\item{\code{GRP2}}{ProlfquAppConfig object with analysis configuration}

\item{\code{prot_annot}}{ProteinAnnotation object}

\item{\code{Contrasts}}{list of contrasts for differential expression analysis}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-DEAReportGenerator-write_DEA_all"></a>}}
\if{latex}{\out{\hypertarget{method-DEAReportGenerator-write_DEA_all}{}}}
\subsection{Method \code{write_DEA_all()}}{
Write all DEA results to files
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DEAReportGenerator$write_DEA_all()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-DEAReportGenerator-render_DEA"></a>}}
\if{latex}{\out{\hypertarget{method-DEAReportGenerator-render_DEA}{}}}
\subsection{Method \code{render_DEA()}}{
Render DEA report using R Markdown
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DEAReportGenerator$render_DEA(
  htmlname,
  markdown = "_Grp2Analysis.Rmd",
  word = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{htmlname}}{name for the output HTML file}

\item{\code{markdown}}{path to the R Markdown template file}

\item{\code{word}}{logical, if TRUE output Word document, otherwise HTML}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-DEAReportGenerator-make_boxplots"></a>}}
\if{latex}{\out{\hypertarget{method-DEAReportGenerator-make_boxplots}{}}}
\subsection{Method \code{make_boxplots()}}{
Generate boxplots for quality control
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DEAReportGenerator$make_boxplots()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-DEAReportGenerator-prep_result_list"></a>}}
\if{latex}{\out{\hypertarget{method-DEAReportGenerator-prep_result_list}{}}}
\subsection{Method \code{prep_result_list()}}{
Prepare result list with all analysis outputs
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DEAReportGenerator$prep_result_list()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
list containing all analysis results
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-DEAReportGenerator-make_SummarizedExperiment"></a>}}
\if{latex}{\out{\hypertarget{method-DEAReportGenerator-make_SummarizedExperiment}{}}}
\subsection{Method \code{make_SummarizedExperiment()}}{
Create SummarizedExperiment object from analysis results
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DEAReportGenerator$make_SummarizedExperiment(
  strip = "~lfq~light",
  .url_builder = bfabric_url_builder
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{strip}}{pattern to strip from rownames}

\item{\code{.url_builder}}{function to build URLs for bfabric}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
SummarizedExperiment object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-DEAReportGenerator-clone"></a>}}
\if{latex}{\out{\hypertarget{method-DEAReportGenerator-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DEAReportGenerator$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
