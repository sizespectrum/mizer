\docType{methods}
\name{getFeedingLevel}
\alias{getFeedingLevel}
\alias{getFeedingLevel,MizerParams,matrix,numeric-method}
\alias{getFeedingLevel,MizerSim,missing,missing-method}
\title{getFeedingLevel method for the size based model}
\arguments{
  \item{object}{A \code{MizerParams} or \code{MizerSim}
  object}

  \item{n}{A matrix of species abundance (species x size).
  Only used if \code{object} argument is of type
  \code{MizerParams}.}

  \item{n_pp}{A vector of the background abundance by size.
  Only used if \code{object} argument is of type
  \code{MizerParams}.}

  \item{time_range}{The time range (either a vector of
  values, a vector of min and max time, or a single value)
  to average the abundances over. Default is the whole time
  range. Only used if \code{object} argument is of type
  \code{MizerSim}.}

  \item{drop}{should extra dimensions of length 1 in the
  output be dropped, simplifying the output. Defaults to
  TRUE}
}
\description{
  Calculates the amount of food consumed by a predator by
  predator size based on food availability, search volume
  and maximum intake This method is used by the
  \code{\link{project}} method for performing simulations.
}
\note{
  If a \code{MizerParams} object is passed in, the method
  returns a two dimensional array (predator species x
  predator size) based on the abundances also passed in. If
  a \code{MizerSim} object is passed in, the method returns
  a three dimensional array (time step x predator species x
  predator size) with the feeding level calculated at every
  time step in the simulation.
}
\examples{
data(species_params_gears)
data(inter)
params <- MizerParams(species_params_gears, inter)
# With constant fishing effort for all gears for 20 time steps
sim <- project(params, t_max = 20, effort = 0.5)
# Get the feeding level at one time step
n <- sim@n[21,,]
n_pp <- sim@n_pp[21,]
getFeedingLevel(params,n,n_pp)
# Get the feeding level at all saved time steps
getFeedingLevel(sim)
}

