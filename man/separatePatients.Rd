% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/separatePatients.R
\name{separatePatients}
\alias{separatePatients}
\title{Separates the patients into two groups based on the score of interest}
\usage{
separatePatients(
  applyFileUsed = "classification",
  metricToUse = "CS",
  verbose = TRUE
)
}
\arguments{
\item{applyFileUsed}{A character value indicating which \verb{*.apply} file to use for patient separation. If set to \code{classification}, then the Classification Signature \verb{*.apply} file located in the \code{output > apply} folder is used. If set to \code{custom}, then each \verb{*.apply} file from the \verb{output > apply > customSignatures > fullTables} folder are iteratively used.}

\item{metricToUse}{A character value indicating which metric to use for patients separation. If set to \code{CS}, then the Clinical Scores for each patient will be used. If set to \code{QS}, then the Quantitative Scores for each patient will be used for patients separation.}

\item{verbose}{A boolean determining if the function should regularly output its progression in the console (if set to \code{TRUE}).}
}
\value{
Generated output files are exported to the \code{output > class} folder. Their precise amounts and organization varies according to the \code{applyFileUsed} argument (see the online tutorial for more details).
}
\description{
This function separates the patients of a dataset into two groups: one which lowly express the desired metric and the other one which highly express the desired metric.
}
