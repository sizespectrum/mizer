\docType{methods}
\name{getEGrowth}
\alias{getEGrowth}
\alias{getEGrowth,MizerParams,matrix,numeric-method}
\title{getEGrowth method for the size based model}
\arguments{
  \item{object}{A \code{MizerParams} object.}

  \item{n}{A matrix of species abundance (species x size).}

  \item{n_pp}{A vector of the background abundance by
  size.}
}
\value{
  A two dimensional array (prey species x prey size)
}
\description{
  Calculates the energy available by species and size for
  growth after metabolism and movement have been accounted
  for. Used by the \code{project} method for performing
  simulations.
}
\examples{
data(species_params_gears)
data(inter)
params <- MizerParams(species_params_gears, inter)
# Project with constant fishing effort for all gears for 20 time steps
sim <- project(params, t_max = 20, effort = 0.5)
# Get the energy at a particular time step
getEGrowth(params,sim@n[21,,],sim@n_pp[21,])
}
\seealso{
  \code{\link{project}}
}

