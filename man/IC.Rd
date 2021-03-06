% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Information_Criterion.R
\name{IC}
\alias{IC}
\title{Calculate the Information Criterion for a robust model}
\usage{
IC(vector, data, family, tcc)
}
\arguments{
\item{vector}{a boolean vector indicating which variables to include.}

\item{data}{a data frame containing the variables in the model.}

\item{family}{a description of the error distribution and link function to be used in the model.}

\item{tcc}{the tuning constant c in Huber's psi-function.}
}
\description{
Calculate the Information Criterion for a robust model
}
