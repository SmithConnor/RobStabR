% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s_values.R
\name{s_values}
\alias{s_values}
\title{Compute the subtractive lack-of-fit measures}
\usage{
s_values(weights, data, n, p, family, coef, wald, dev, tcc)
}
\arguments{
\item{weights}{a vector contains=ing a subset of observations.}

\item{data}{a data frame containing the variables in the model.}

\item{n}{Number of observations.}

\item{p}{Number of variables.}

\item{family}{a description of the error distribution and link function to be used in the model.}

\item{coef}{A TRUE/FALSE value to indicate whether to evaluate RobStab using regression coefficients.}

\item{wald}{A TRUE/FALSE value to indicate whether to evaluate RobStab using wald statisitcs.}

\item{dev}{A TRUE/FALSE value to indicate whether to evaluate RobStab using deviances.}

\item{tcc}{the tuning constant c in Huber's psi-function.}
}
\description{
Compute the subtractive lack-of-fit measures
}
