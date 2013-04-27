\docType{methods}
\name{getRDI}
\alias{getRDI}
\alias{getRDI,MizerParams,matrix,numeric-method}
\title{getRDI method for the size based model}
\arguments{
  \item{object}{A \code{MizerParams} object.}

  \item{n}{A matrix of species abundance (species x size).}

  \item{n_pp}{A vector of the background abundance by
  size.}

  \item{sex_ratio}{Proportion of the population that is
  female. Default value is 0.5.}
}
\value{
  A numeric vector the length of the number of species
}
\description{
  Calculates the density independent recruitment (flux
  entering the smallest size class of each species) before
  density dependence, by species. Used by the
  \code{project} method for performing simulations.
}
\examples{
data(species_params_gears)
data(inter)
params <- MizerParams(species_params_gears, inter)
# Project with constant fishing effort for all gears for 20 time steps
sim <- project(params, t_max = 20, effort = 0.5)
# Get the recruitment at a particular time step
getRDI(params,sim@n[21,,],sim@n_pp[21,])
}
\seealso{
  \code{\link{project}}
}

