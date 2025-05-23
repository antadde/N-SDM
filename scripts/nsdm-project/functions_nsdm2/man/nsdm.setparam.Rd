% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsdm.setparam.R
\name{nsdm.setparam}
\alias{nsdm.setparam}
\title{nsdm.setparam}
\usage{
nsdm.setparam(
  model_name,
  param_grid,
  covariate_names,
  tmp_path,
  ncov.esm = NULL,
  comb.esm = NULL,
  weights = 1,
  nthreads
)
}
\arguments{
\item{model_name}{A character string specifying the modeling algorithm to use (e.g., "glm", "gam", "rf", "gbm", "max", "esm").}

\item{param_grid}{A character string indicating the complete file path where the \code{.xlsx} document containing the parameter settings is stored.}

\item{covariate_names}{A character vector specifying the names of covariates to be included in the models.}

\item{tmp_path}{A character string specifying the directory where temporary outputs should be stored.}

\item{ncov.esm}{(Optional) A numeric value specifying the total number of covariates evaluated when using the ESM (Ensemble of Small Models) framework.}

\item{comb.esm}{(Optional) A numeric value indicating the number of covariates combined in each individual small model under the ESM framework.}

\item{weights}{A numeric vector specifying the weights to be applied in the models. If a single value is provided, no weighting is applied.}

\item{nthreads}{A numeric value specifying the number of cores to be used for parallel computations.}
}
\value{
A list containing elements of class \code{'multi.input'}, each representing a model configuration to be fitted.
}
\description{
Define model parameter settings for various modeling algorithms.
}
\author{
Antoine Adde (\email{antoine.adde@eawag.ch})
}
