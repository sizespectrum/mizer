\docType{methods}
\name{getYield}
\alias{getYield}
\alias{getYield,MizerSim-method}
\alias{getYield-method}
\title{Calculate the total yield of each species}
\arguments{
  \item{object}{An object of class \code{MizerSim}.}
}
\value{
  An array containing the total yield (time x species)
}
\description{
  Calculates the total yield of each species across all
  gears at each simulation time step.
}
\examples{
\dontrun{
data(NS_species_params_gears)
data(inter)
params <- MizerParams(NS_species_params_gears, inter)
sim <- project(params, effort=1, t_max=10)
y <- getYield(sim)
}
}
\seealso{
  \code{\link{getYieldGear}}
}

