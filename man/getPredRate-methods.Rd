\docType{methods}
\name{getPredRate}
\alias{getPredRate}
\alias{getPredRate,MizerParams,matrix,numeric-method}
\alias{getPredRate-method}
\title{getPredRate method for the size based model}
\arguments{
  \item{object}{A \code{MizerParams} object.}

  \item{n}{A matrix of species abundance (species x size).}

  \item{n_pp}{A vector of the background abundance by
  size.}

  \item{feeding_level}{The current feeding level
  (optional). A matrix of size no. species x no. size bins.
  If not supplied, is calculated internally using the
  \code{getFeedingLevel()} method.}
}
\value{
  A three dimensional array (predator species x predator
  size x prey size)
}
\description{
  Calculates the predation rate of each predator species at
  size on prey size. This method is used by the
  \code{\link{project}} method for performing simulations.
  In the simulations, it is combined with the interaction
  matrix (see \code{\link{MizerParams}}) to calculate the
  realised predation mortality (see \code{\link{getM2}}).
}
\examples{
data(NS_species_params_gears)
data(inter)
params <- MizerParams(NS_species_params_gears, inter)
# With constant fishing effort for all gears for 20 time steps
sim <- project(params, t_max = 20, effort = 0.5)
# Get the feeding level at one time step
n <- sim@n[21,,]
n_pp <- sim@n_pp[21,]
getPredRate(params,n,n_pp)
}
\seealso{
  \code{\link{project}}, \code{\link{getM2}},
  \code{\link{getFeedingLevel}} and
  \code{\link{MizerParams}}
}

