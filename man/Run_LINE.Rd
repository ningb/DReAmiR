% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pri_target_graph.R
\name{Run_LINE}
\alias{Run_LINE}
\title{Given a list of correlation matrcies, return the SpQN normalized matrices}
\usage{
Run_LINE(
  le.gene,
  dim = 5,
  edge.sign = c("Negative", "Positive", "Both"),
  do.dim.reduction = TRUE
)
}
\arguments{
\item{le.gene}{A list of correlation matrices with miRNA on columns and mRNA on rows}

\item{dim}{An integer indicating number embedded features returned from LINE}

\item{edge.sign}{An character for edge list weight; if negative, the reverse will be use}

\item{do.dim.reduction}{A logical variable, whether to return dimension reduction results}
}
\value{
A matrix of graph embedding features for all the leading edge genes
}
\description{
Given a list of correlation matrcies, return the SpQN normalized matrices
}
