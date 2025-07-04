% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R6_AnnotationProcessor.R
\name{AnnotationProcessor}
\alias{AnnotationProcessor}
\title{AnnotationProcessor}
\description{
AnnotationProcessor

AnnotationProcessor
}
\examples{

# AnnotationProcessor$debug("read_annotation")
ap <- AnnotationProcessor$new(prefix = "G_")

annot <- data.frame(
file = c("a1.raw","a2.raw","a3.raw","a4.raw"),
group = c("a","a","b","b"),
CONTROL = c("C","C","T","T"),
Subject = c("X","Y","X","Y"))
ap$check_annotation(annot)
af <- annot
af$file <- NULL
testthat::expect_error(ap$check_annotation(af), "column starting with :")
af <- annot
af$group <- NULL
testthat::expect_error(ap$check_annotation(af),"column starting with :")
aa <- ap$read_annotation(annot)
stopifnot(length(aa$atable$factor_keys_depth()) == 2)
stopifnot(all(c("atable", "annot", "contrasts") \%in\% names(aa)))
stopifnot(aa$contrasts == "G_b - G_a")
af <- annot
af$CONTROL <- NULL
testthat::expect_error(ap$check_annotation(af),"you must specify a CONTROL column")
af <- annot
af$Subject <- NULL
testthat::expect_warning(ap$check_annotation(af),"column starting with")


# should not throw exception since QC does not require group or subject
ap <- AnnotationProcessor$new(QC = TRUE)
af <- annot
# af$group <- NULL
af$CONTROL <- NULL
af$Subject <- NULL
ap$check_annotation(af)
aa <- ap$read_annotation(af)

stopifnot(aa$atable$factor_keys() == "G_")
stopifnot(aa$atable$factors == "group")
aa <- ap$read_annotation(annot)
aa$atable$fileName
aa$atable$sampleName
as <- annot
as$sample <- c("s1","s2","s3","s4")
aa <- ap$read_annotation(annot)
aa$atable$sampleName
stopifnot(is.null(aa$annotation))

annot <- data.frame(
file = c("a1.raw","a2.raw","a3.raw","a4.raw"),
Name = c("a1.raw","a2.raw","a3.raw","a4.raw"),
"Grouping Var" = c("a","a","b","b"),
CONTROL = c("C","C","T","T"),
Subject = c("X","Y","X","Y"))
ax <- ap$read_annotation(annot)
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{QC}}{is it a QC run}

\item{\code{prefix}}{name for one factor designs}

\item{\code{repeated}}{is it a repeated measurement}

\item{\code{SAINT}}{is it a AP MS experiment, then use Bait_ as prefix}

\item{\code{file_pattern}}{colnames for file}

\item{\code{grouping_pattern}}{colnames grouping variable}

\item{\code{subject_pattern}}{colnames for pairing variable}

\item{\code{control_pattern}}{contrast specification columns}

\item{\code{control_col_pattern}}{columns which contains C or T.}

\item{\code{sample_name_pattern}}{sample name column}

\item{\code{strict}}{should name check be strict}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-AnnotationProcessor-new}{\code{AnnotationProcessor$new()}}
\item \href{#method-AnnotationProcessor-check_annotation}{\code{AnnotationProcessor$check_annotation()}}
\item \href{#method-AnnotationProcessor-read_annotation}{\code{AnnotationProcessor$read_annotation()}}
\item \href{#method-AnnotationProcessor-extract_contrasts}{\code{AnnotationProcessor$extract_contrasts()}}
\item \href{#method-AnnotationProcessor-add_contrasts_vec}{\code{AnnotationProcessor$add_contrasts_vec()}}
\item \href{#method-AnnotationProcessor-clone}{\code{AnnotationProcessor$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-AnnotationProcessor-new"></a>}}
\if{latex}{\out{\hypertarget{method-AnnotationProcessor-new}{}}}
\subsection{Method \code{new()}}{
initialize
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AnnotationProcessor$new(
  QC = FALSE,
  prefix = "G_",
  repeated = TRUE,
  SAINT = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{QC}}{default FALSE}

\item{\code{prefix}}{default "G_"}

\item{\code{repeated}}{default TRUE}

\item{\code{SAINT}}{default FALSE}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-AnnotationProcessor-check_annotation"></a>}}
\if{latex}{\out{\hypertarget{method-AnnotationProcessor-check_annotation}{}}}
\subsection{Method \code{check_annotation()}}{
check annotation
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AnnotationProcessor$check_annotation(annot)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{annot}}{annotation}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-AnnotationProcessor-read_annotation"></a>}}
\if{latex}{\out{\hypertarget{method-AnnotationProcessor-read_annotation}{}}}
\subsection{Method \code{read_annotation()}}{
read annotation
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AnnotationProcessor$read_annotation(dsf)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{dsf}}{either dataframe or file path.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-AnnotationProcessor-extract_contrasts"></a>}}
\if{latex}{\out{\hypertarget{method-AnnotationProcessor-extract_contrasts}{}}}
\subsection{Method \code{extract_contrasts()}}{
check annotation
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AnnotationProcessor$extract_contrasts(annot, group)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{annot}}{annotation}

\item{\code{group}}{group column e.g. group}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-AnnotationProcessor-add_contrasts_vec"></a>}}
\if{latex}{\out{\hypertarget{method-AnnotationProcessor-add_contrasts_vec}{}}}
\subsection{Method \code{add_contrasts_vec()}}{
add vector of contrasts to annot table
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AnnotationProcessor$add_contrasts_vec(annot, Contrasts)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{annot}}{annotation}

\item{\code{Contrasts}}{vector with contrasts}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-AnnotationProcessor-clone"></a>}}
\if{latex}{\out{\hypertarget{method-AnnotationProcessor-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AnnotationProcessor$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
