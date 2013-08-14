\docType{methods}
\name{getYieldGear}
\alias{getYieldGear}
\alias{getYieldGear,MizerSim-method}
\alias{getYieldGear-method}
\title{Calculate the total yield per gear and species}
\arguments{
  \item{object}{An object of class \code{MizerSim}.}
}
\value{
  An array containing the total yield (time x gear x
  species)
}
\description{
  Calculates the total yield per gear and species at each
  simulation time step.
}
\examples{
data(species_params_gears)
data(inter)
params <- MizerParams(species_params_gears, inter)
# With constant fishing effort for all gears for 20 time steps
sim <- project(params, t_max = 20, effort = 0.5)
getYieldGear(sim)
}
\seealso{
  \code{\link{getYield}}
}

