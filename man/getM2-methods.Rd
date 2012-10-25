\docType{methods}
\name{getM2}
\alias{getM2}
\alias{getM2,MizerParams,matrix,numeric-method}
\alias{getM2,MizerSim,missing,missing-method}
\title{getM2 method for the size based model}
\arguments{
  \item{object}{An \code{MizerParams} object}

  \item{n}{A matrix of species abundance (species x size)}

  \item{n_pp}{A vector of the background abundance by size}

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
  A two dimensional array (prey species x prey size)
}
\description{
  Calculates the predation mortality on each prey species
  by prey size
}

