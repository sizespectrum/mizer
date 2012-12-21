\docType{methods}
\name{getYieldGear}
\alias{getYieldGear}
\alias{getYieldGear,MizerSim-method}
\title{Calculate the total yield per gear and species}
\arguments{
  \item{object}{An object of class \code{MizerSim}}
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
\dontrun{
data(species_params_gears)
data(inter)
params <- MizerParams(species_params_gears, inter)
sim <- project(params, effort=1, t_max=10)
y <- getYieldGear(sim)
}
}

