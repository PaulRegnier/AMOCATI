% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resetWorkspace.R
\name{resetWorkspace}
\alias{resetWorkspace}
\title{Reset the workspace}
\usage{
resetWorkspace(eraseEntireRMemory = FALSE, verbose = TRUE)
}
\arguments{
\item{eraseEntireRMemory}{A boolean value indicating if the entire R memory should be erased (if set to \code{TRUE}) or not (if set to \code{FALSE}).}

\item{verbose}{A boolean determining if the function should regularly output its progression in the console (if set to \code{TRUE}).}
}
\description{
This function resets the workspace that AMOCATI will use for its tasks. It can also erase the whole R memory.
}
