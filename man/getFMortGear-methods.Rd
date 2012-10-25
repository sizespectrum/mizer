\docType{methods}
\name{getFMortGear}
\alias{getFMortGear}
\alias{getFMortGear,MizerParams,matrix-method}
\alias{getFMortGear,MizerParams,numeric-method}
\alias{getFMortGear,MizerSim,missing-method}
\title{Get the fishing mortality by time, gear, species and size}
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
}
\value{
  An array. If effort argument has a time dimension, output
  array has four dimensions (time x gear x species x size).
  If effort argument does not have a time dimension, output
  array has three dimensions (gear x species x size).
}
\description{
  Calculates the fishing mortality by gear, species and
  size at each time step in the \code{effort} argument.
  Fishing mortality = catchability x selectivity x effort
}

