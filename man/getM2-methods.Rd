\docType{methods}
\name{getM2}
\alias{getM2}
\alias{getM2,MizerParams,matrix,numeric-method}
\alias{getM2,MizerSim,missing,missing-method}
\title{getM2 method for the size based model}
\arguments{
  \item{object}{A \code{MizerParams} or \code{MizerSim}
  object.}

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

  \item{drop}{Only used when object is of type
  \code{MizerSim}. Should dimensions of length 1 in the
  output be dropped, simplifying the output. Defaults to
  TRUE}
}
\description{
  Calculates the total predation mortality on each prey
  species by prey size. This method is used by the
  \code{\link{project}} method for performing simulations.
}
\note{
  If a \code{MizerParams} object is passed in, the method
  returns a two dimensional array (prey species x prey
  size) based on the abundances also passed in. If a
  \code{MizerSim} object is passed in, the method returns a
  three dimensional array (time step x prey species x prey
  size) with the predation mortality calculated at every
  time step in the simulation.
}
\examples{
data(species_params_gears)
data(inter)
params <- MizerParams(species_params_gears, inter)
# With constant fishing effort for all gears for 20 time steps
sim <- project(params, t_max = 20, effort = 0.5)
# Get M2 at one time step
n <- sim@n[21,,]
n_pp <- sim@n_pp[21,]
getM2(params,n,n_pp)
# Get M2 at all saved time steps
getM2(sim)
}
\seealso{
  \code{\link{getPredRate}} and \code{\link{project}}.
}

