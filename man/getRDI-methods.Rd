\docType{methods}
\name{getRDI}
\alias{getRDI}
\alias{getRDI,MizerParams,matrix,numeric,matrix-method}
\alias{getRDI,MizerParams,matrix,numeric,missing-method}
\title{getRDI method for the size based model}
\arguments{
  \item{object}{A \code{MizerParams} object.}

  \item{n}{A matrix of species abundance (species x size).}

  \item{n_pp}{A vector of the background abundance by
  size.}

  \item{e_spawning}{The energy available for spawning
  (optional). A matrix of size no. species x no. size bins.
  If not supplied, is calculated internally using the
  \code{getESpawning()} method.}

  \item{sex_ratio}{Proportion of the population that is
  female. Default value is 0.5.}
}
\value{
  A numeric vector the length of the number of species
}
\description{
  Calculates the density independent recruitment (total egg
  production) before density dependence, by species. Used
  by the \code{project} method for performing simulations.
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

