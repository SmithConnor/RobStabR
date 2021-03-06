% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/candidate_space.R
\name{model_space}
\alias{model_space}
\title{Find the candidate model space after B bootstraps}
\usage{
model_space(
  data,
  B,
  m,
  nStrata = 8,
  family,
  k = 1,
  resid = "pearson",
  coef = TRUE,
  wald = TRUE,
  dev = TRUE,
  bootstraps = NA,
  tcc
)
}
\arguments{
\item{data}{a data frame containing the variables in the model.}

\item{B}{The number of bootstrap re-samples to run.}

\item{m}{The m used for the m-out-of-n boostrap.}

\item{nStrata}{Number of strat used for the stratified re-sampling.}

\item{family}{a description of the error distribution and link function to be used in the model.}

\item{k}{The minimum number of times a mondel must be identified through resampling to be included in the rediced candidate space.}

\item{resid}{the type of residuals to be used for the stratified re-sampling.}

\item{coef}{A TRUE/FALSE value to indicate whether to evaluate RobStab using regression coefficients.}

\item{wald}{A TRUE/FALSE value to indicate whether to evaluate RobStab using wald statisitcs.}

\item{dev}{A TRUE/FALSE value to indicate whether to evaluate RobStab using deviances.}

\item{bootstraps}{Recreact results using previously used resamples.}

\item{tcc}{the tuning constant c in Huber's psi-function.}
}
\description{
Find the candidate model space after B bootstraps
}
