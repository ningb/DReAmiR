% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/miSpQN.R
\name{normalize_correlation}
\alias{normalize_correlation}
\title{Main function to remove mean-correlation relationship between miRNA and mRNA}
\usage{
normalize_correlation(
  cor_mat,
  ave_exp1,
  ave_exp2,
  ngrp = 20,
  size_grp = c(1000, 150),
  ref_grp = 18
)
}
\arguments{
\item{cor_mat}{A correlation matrix for a group with mRNA on the rows and miRNA on the columns}

\item{ave_exp1}{A numeric vecotr with same length as the nrow(cor_mat) for average mRNA expressions}

\item{ave_exp2}{A numeric vecotr with same length as the ncol(cor_mat) for average miRNA expressions}

\item{ngrp}{An integer for number of inner bins. Default is 20}

\item{size_grp}{A numeric vector of length 2 for number of mRNA and miRNA per outer group. Default c(1000, 150)}

\item{ref_grp}{An integer for the index of bin to use as reference. Default is 18}
}
\value{
A normalized correlation matrix with the same dimension as cor_mat
}
\description{
Main function to remove mean-correlation relationship between miRNA and mRNA
}
