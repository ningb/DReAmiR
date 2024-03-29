% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/miSpQN.R
\name{get_condition_exp}
\alias{get_condition_exp}
\title{Calculate IQR and Q1 per bin x bin; Calculate expression mean per bin}
\usage{
get_condition_exp(cor_mat, ave_exp1, ave_exp2)
}
\arguments{
\item{cor_mat}{A correlation matrix for a group with mRNA on the rows and miRNA on the columns}

\item{ave_exp1}{A numeric vecotr with same length as the nrow(cor_mat) for average mRNA expressions}

\item{ave_exp2}{A numeric vecotr with same length as the ncol(cor_mat) for average miRNA expressions}
}
\value{
a list for IQR/Q1/average expression for mRNA and miRNA
}
\description{
Calculate IQR and Q1 per bin x bin; Calculate expression mean per bin
}
