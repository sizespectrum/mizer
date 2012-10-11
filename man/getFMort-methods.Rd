\docType{methods}
\name{getFMort}
\alias{getFMort}
\alias{getFMort,MizerParams,matrix-method}
\alias{getFMort,MizerParams,numeric-method}
\alias{getFMort,MizerSim,missing-method}
\title{Get the total fishing mortality from all fishing gears by time, species and size}
\arguments{
  \item{object}{A \code{MizerParams} object or an
  \code{MizerSim} object}

  \item{effort}{The effort of each fishing gear. Only
  needed if the object argument is of class
  \code{MizerParams}. Can be an two dimensional array (time
  x gear), a vector of length equal to the number of gears,
  or a single numeric value (each gear has the same
  effort). If the object argument is of class
  \code{MizerSim} then the effort slot of the
  \code{MizerSim} object is used and this argument is not
  used.}

  \item{time_range}{The time range (either a vector of
  values, a vector of min and max time, or a single value)
  to average the abundances over. Default is the whole time
  range. Only used if \code{object} argument is of type
  MizerSim.}
}
\value{
  An array. If effort argument has a time dimension, output
  array has three dimensions (time x species x size). If
  effort argument does not have a time dimension, output
  array has two dimensions (species x size).
}
\description{
  Calculates the fishing mortality from all grars by
  species and size at each time step in the \code{effort}
  argument. Fishing mortality from each gear = catchability
  x selectivity x effort Total fishing mortality is just
  the sum from each gear
}

