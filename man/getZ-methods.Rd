\docType{methods}
\name{getZ}
\alias{getZ}
\alias{getZ,MizerParams,matrix,numeric,numeric-method}
\title{getZ method for the size based model}
\arguments{
  \item{object}{A \code{MizerParams} object.}

  \item{n}{A matrix of species abundance (species x size).}

  \item{n_pp}{A vector of the background abundance by
  size.}

  \item{effort}{A numeric vector of the effort by gear or a
  single numeric effort value which is used for all gears.}
}
\value{
  A two dimensional array (prey species x prey size).
}
\description{
  Calculates the total mortality on each species by size
  from predation mortality (M2), background mortality (M)
  and fishing mortality for a single time step.
}
\examples{
data(species_params_gears)
data(inter)
params <- MizerParams(species_params_gears, inter)
# Project with constant fishing effort for all gears for 20 time steps
sim <- project(params, t_max = 20, effort = 0.5)
# Get the total mortality at a particular time step
getZ(params,sim@n[21,,],sim@n_pp[21,],effort=0.5)
}
\seealso{
  \code{\link{getM2}}, \code{\link{getFMort}}
}

