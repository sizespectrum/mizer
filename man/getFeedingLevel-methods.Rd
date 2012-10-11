\docType{methods}
\name{getFeedingLevel}
\alias{getFeedingLevel}
\alias{getFeedingLevel,MizerParams,matrix,numeric-method}
\alias{getFeedingLevel,MizerSim,missing,missing-method}
\title{getFeedingLevel method for the size based model}
\arguments{
  \item{object}{A MizerParams or MizerSim object}

  \item{n}{A matrix of species abundance (species x size).
  Only used if \code{object} argument is of type
  MizerParam.}

  \item{n_pp}{A vector of the background abundance by size.
  Only used if \code{object} argument is of type
  MizerParam.}

  \item{time_range}{The time range (either a vector of
  values, a vector of min and max time, or a single value)
  to average the abundances over. Default is the whole time
  range. Only used if \code{object} argument is of type
  MizerSim.}

  \item{.drop}{should extra dimensions of length 1 in the
  output be dropped, simplifying the output. Defaults to
  TRUE}
}
\value{
  A two dimensional array (predator species x predator
  size)
}
\description{
  Calculates the amount of food consumed by a predator by
  predator size based on food availability, search volume
  and maximum intake
}

