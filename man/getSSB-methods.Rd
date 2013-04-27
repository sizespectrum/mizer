\docType{methods}
\name{getSSB}
\alias{getSSB}
\alias{getSSB,MizerSim-method}
\title{Calculate the SSB of species}
\arguments{
  \item{object}{An object of class \code{MizerSim}.}
}
\value{
  An array containing the SSB (time x species)
}
\description{
  Calculates the spawning stock biomass (SSB) through time
  of the species in the \code{MizerSim} class. SSB is
  calculated as the total mass of all mature individuals.
}
\examples{
data(species_params_gears)
data(inter)
params <- MizerParams(species_params_gears, inter)
# With constant fishing effort for all gears for 20 time steps
sim <- project(params, t_max = 20, effort = 0.5)
getSSB(sim)
}

