% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diffreg.R
\name{.CalcGseaStat}
\alias{.CalcGseaStat}
\title{Calculating random walking score}
\usage{
.CalcGseaStat(ranklist, gs, alpha = NULL)
}
\arguments{
\item{ranklist}{A named vector for all genes sorted by correlation}

\item{gs}{A character vector for the set of targets of interest}

\item{alpha}{A numeric variable for gsea parameter}
}
\value{
A walking score vector
}
\description{
Calculating random walking score
}
