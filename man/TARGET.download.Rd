% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TARGET.download.R
\name{TARGET.download}
\alias{TARGET.download}
\title{Downloads the data (RNA-Seq + clinical) from the desired TARGET project}
\usage{
TARGET.download(projectID = NULL, verbose = TRUE)
}
\arguments{
\item{projectID}{A character value determining which project to use.}

\item{verbose}{A boolean determining if the function should regularly output its progression in the console (if set to \code{TRUE}).}
}
\value{
Generated output files are fist exported to the \code{temp} folder before being moved to the \code{output > data} folder.
}
\description{
This function downloads the data (RNA-Seq + clinical) from the desired GDC TARGET project.
}
