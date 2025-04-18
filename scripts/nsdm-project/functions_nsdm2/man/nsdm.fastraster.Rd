% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsdm.fastraster.R
\name{nsdm.fastraster}
\alias{nsdm.fastraster}
\title{nsdm.fastraster}
\usage{
nsdm.fastraster(files, nsplits = parallel::detectCores() - 1)
}
\arguments{
\item{files}{A character vector containing file paths to raster layers (\code{.tif}).}

\item{nsplits}{A numeric value specifying the number of CPU cores to use for parallel loading.
Defaults to \code{parallel::detectCores() - 1} if NULL.}
}
\value{
A \code{SpatRaster}.
}
\description{
Parallelized loading of raster datasets using multiple cores.
}
\author{
Antoine Adde (antoine.adde@eawag.ch)
}
